using namespace std;

TH1 * Correct_DataEff(std::string histoname,TH1 * Histo, TH1 * LATcorr);
TH1 * DivideHisto(TH1 *Histo1, TH1 *Histo2);

class DatavsMC
{

private:
	//LAT corr.
	TH1 * LATcorr_R  = NULL; 	
	TH1 * LATcorr_TOF= NULL;
	TH1 * LATcorr_NaF= NULL;
	TH1 * LATcorr_Agl= NULL; 	
        
	TH1 * Correction_R  ;
        TH1 * Correction_TOF;
        TH1 * Correction_NaF;
	TH1 * Correction_Agl;
	
	std::string Basename;

	int latzones  = 0;
	int selections= 0;
	int mc_types  = 0;
public:
	Efficiency * MCEff;
	Efficiency * DataEff;

	Efficiency * DataEff_corr;

	//creation constructor
		
	DatavsMC(std::string basename, int n=1, int S=1, int mcs=1){
		if(S==1){	
			if(mcs == 1){ 
				MCEff    = new Efficiency((basename + "_MC"  ).c_str());
				DataEff  = new Efficiency((basename + "_Data").c_str(),n);
			}
			else{
				MCEff    = new Efficiency((basename + "_MC"  ).c_str(),mcs);  	
				DataEff  = new Efficiency((basename + "_Data").c_str(),n);
			}

		}
		else{	
			if(mcs == 1){ 
				MCEff    = new Efficiency((basename + "_MC"  ).c_str(),S);
				DataEff  = new Efficiency((basename + "_Data").c_str(),n,S);
			}
			else{
				MCEff    = new Efficiency((basename + "_MC"  ).c_str(),mcs,S);  	
				DataEff  = new Efficiency((basename + "_Data").c_str(),n,S);
			}

		}


		latzones   = n;
		selections = S;
		mc_types   = mcs;	
	
		Basename = basename;

	}
 
	//reading constructors
	DatavsMC(TFile *file, std::string basename, int mcs=1){
		MCEff   = new Efficiency(file,(basename + "_MC"  ).c_str());
		DataEff = new Efficiency(file,(basename + "_Data").c_str());

		DataEff_corr  = new Efficiency((basename + "_Data").c_str());

		latzones   = DataEff -> beforeR -> GetNbinsY();
		selections = DataEff -> beforeR -> GetNbinsZ();
		mc_types   = mcs;
		
		Correction_R    =  (TH1 *) MCEff -> beforeR   -> Clone();	
		Correction_TOF  =  (TH1 *) MCEff -> beforeTOF -> Clone();
		Correction_NaF  =  (TH1 *) MCEff -> beforeNaF -> Clone();
		Correction_Agl  =  (TH1 *) MCEff -> beforeAgl -> Clone();
		
		Basename = basename;

	}	

	void Write();
	void Assign_LatCorr(TH1 * LatcorrR, TH1 * LatcorrTOF, TH1 * LatcorrNaF, TH1 * LatcorrAgl)
	{
		LATcorr_R   = (TH1*)LatcorrR   ->Clone();
                LATcorr_TOF = (TH1*)LatcorrTOF ->Clone();
                LATcorr_NaF = (TH1*)LatcorrNaF ->Clone();
                LATcorr_Agl = (TH1*)LatcorrAgl ->Clone();

	};
	
	void Eval_Corrected_DataEff();
	void Eval_DandMC_Eff();
	void DivideHisto(TH1 *Histo1, TH1 *Histo2, TH1 * Correction);
	void Eval_Corrections();	
	
	TH1 * GetCorrection_R()  { return Correction_R  ;};
	TH1 * GetCorrection_TOF(){ return Correction_TOF;};
	TH1 * GetCorrection_NaF(){ return Correction_NaF;};
	TH1 * GetCorrection_Agl(){ return Correction_Agl;}; 
};


void DatavsMC::Write(){
	MCEff 	-> Write();
	DataEff -> Write();
	return;
}

void DatavsMC::Eval_Corrected_DataEff(){

	if(selections == 1) {
		DataEff_corr -> beforeR   = ProjectionXtoTH1F((TH2F*)DataEff -> beforeR  , (Basename + "1_R" ),0,latzones);
		DataEff_corr -> beforeTOF = ProjectionXtoTH1F((TH2F*)DataEff -> beforeTOF, (Basename + "1"   ),0,latzones);
		DataEff_corr -> beforeNaF = ProjectionXtoTH1F((TH2F*)DataEff -> beforeNaF, (Basename + "1NaF"),0,latzones);
		DataEff_corr -> beforeAgl = ProjectionXtoTH1F((TH2F*)DataEff -> beforeAgl, (Basename + "1Agl"),0,latzones);
	}

	else {
		DataEff_corr -> beforeR   = (TH2F*)((TH2F *)((TH3F*)DataEff -> beforeR  ) -> Project3D("zx")) -> Clone();
		DataEff_corr -> beforeTOF = (TH2F*)((TH2F *)((TH3F*)DataEff -> beforeTOF) -> Project3D("zx")) -> Clone();
		DataEff_corr -> beforeNaF = (TH2F*)((TH2F *)((TH3F*)DataEff -> beforeNaF) -> Project3D("zx")) -> Clone();
		DataEff_corr -> beforeAgl = (TH2F*)((TH2F *)((TH3F*)DataEff -> beforeAgl) -> Project3D("zx")) -> Clone();
	}
	
	if(!LATcorr_R   ) cout<<"ERROR: Lat. corr for R   histos not assigned"<<endl;
	if(!LATcorr_TOF ) cout<<"ERROR: Lat. corr for TOF histos not assigned"<<endl;
	if(!LATcorr_NaF ) cout<<"ERROR: Lat. corr for NaF histos not assigned"<<endl;
	if(!LATcorr_Agl ) cout<<"ERROR: Lat. corr for Agl histos not assigned"<<endl;
	
	DataEff_corr -> afterR   =  Correct_DataEff( (Basename + "2_R" ),	DataEff -> afterR  ,	LATcorr_R   		);
        DataEff_corr -> afterTOF =  Correct_DataEff( (Basename + "2"   ),	DataEff -> afterTOF,	LATcorr_TOF 		);
        DataEff_corr -> afterNaF =  Correct_DataEff( (Basename + "2NaF"),	DataEff -> afterNaF,	LATcorr_NaF 		);
        DataEff_corr -> afterAgl =  Correct_DataEff( (Basename + "2Agl"),	DataEff -> afterAgl,	LATcorr_Agl 		);
	
	return;
}


void DatavsMC::Eval_DandMC_Eff(){
	DatavsMC::Eval_Corrected_DataEff();	
	
	DataEff_corr -> Eval_Efficiency();
	MCEff 	     -> Eval_Efficiency();	
	
	return;
}


void DatavsMC::Eval_Corrections(){

	DivideHisto( DataEff_corr -> effR   , MCEff -> effR  , Correction_R   );	
	DivideHisto( DataEff_corr -> effTOF , MCEff -> effTOF, Correction_TOF );
	DivideHisto( DataEff_corr -> effNaF , MCEff -> effNaF, Correction_NaF );
	DivideHisto( DataEff_corr -> effAgl , MCEff -> effAgl, Correction_Agl );

	return;
}



void DatavsMC::DivideHisto(TH1 *Histo1, TH1 *Histo2,TH1 * Correction){
    if(selections==1){
        if(mc_types ==1) {
            for(int iR=0;iR<Histo2->GetNbinsX();iR++){ 
		if(Histo2 -> GetBinContent(iR+1)<1&&Histo2 -> GetBinContent(iR+1)>0)
                Correction -> SetBinContent (iR+1,Histo1 -> GetBinContent(iR+1)/(float)Histo2 -> GetBinContent(iR+1));
                Correction -> SetBinError(iR+1,Histo1 -> GetBinError(iR+1));
            }
        }
        else{
            for(int mc_type=0;mc_type<mc_types;mc_type++){
		    for(int iR=0;iR<Histo2->GetNbinsX();iR++){
			    if(Histo2 -> GetBinContent(iR+1,mc_type+1)<1&&Histo2 -> GetBinContent(iR+1,mc_type+1)>0){
				    Correction -> SetBinContent (iR+1,mc_type+1,Histo1 -> GetBinContent(iR+1)/(float)Histo2 -> GetBinContent(iR+1,mc_type+1));
				    Correction -> SetBinError(iR+1,mc_type+1,Histo1 -> GetBinError(iR+1,mc_type+1)); }
			    else{
				    Correction -> SetBinContent (iR+1,mc_type+1,0);
				    Correction -> SetBinError(iR+1,mc_type+1,0); 
			    }	
		    }
		 
	    }
        }

    }
    else{
        if(mc_types == 1) {
            for(int iS=0;iS<selections;iS++)
                for(int iR=0;iR<Histo2->GetNbinsX();iR++){
                    Correction -> SetBinContent (iR+1,iS+1,Histo1 -> GetBinContent(iR+1,iS+1)/(float)Histo2 -> GetBinContent(iR+1,iS+1));
                    Correction -> SetBinError(iR+1,iS+1,Histo1 -> GetBinError(iR+1,iS+1));
                }
        }
        else{
            for(int iS=0;iS<selections;iS++)
                for(int mc_type=0;mc_type<mc_types;mc_type++){
                    for(int iR=0;iR<Histo1->GetNbinsX();iR++){
                        Correction -> SetBinContent (iR+1,mc_type+1,iS+1,Histo1 -> GetBinContent(iR+1,iS+1)/(float)Histo2 -> GetBinContent(iR+1,iS+1,mc_type+1));
                        Correction -> SetBinError(iR+1,mc_type+1,Histo1 -> GetBinError(iR+1,iS+1));
                    }
                }
        }

    }
    return; 
}


TH1 * Correct_DataEff(std::string histoname,TH1 * Histo, TH1 * LATcorr){
	TH1 * Histo_corr;
  
	int selections = Histo ->GetNbinsZ();
	int latzones   = Histo ->GetNbinsY();
	
	if(selections == 1){
		TH2F* temp = (TH2F *)Histo -> Clone();
		//Histo_corr = new TH1F("","",Histo-> GetNbinsX(),0,Histo-> GetNbinsX());
		for(int iR = 0;iR < Histo ->GetNbinsX();iR++) 
		     for(int lat =0; lat < latzones; lat++){	
				temp -> SetBinContent(iR+1,lat+1,Histo->GetBinContent(iR+1,lat+1)*LATcorr->GetBinContent(lat+1));
		}
		Histo_corr = (TH1F*) ProjectionXtoTH1F( temp, histoname,0,latzones );
	}
	else {
		TH3* temp = (TH3 *)Histo -> Clone();
		for(int iS = 0; iS < selections; iS++)	
			for(int iR = 0;iR < Histo ->GetNbinsX();iR++) 
				for(int lat =0; lat < latzones; lat++){	
					temp -> SetBinContent(iR+1,lat+1,iS+1,Histo->GetBinContent(iR+1,lat+1,iS+1)*LATcorr->GetBinContent(lat+1,iS+1));
				}
		Histo_corr = (TH2F*)((TH3*)temp) -> Project3D("zx") -> Clone();	
	}
	
	return Histo_corr;
}
