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
        
	TH1 * Correction_R  = NULL;
        TH1 * Correction_TOF= NULL;
        TH1 * Correction_NaF= NULL;
	TH1 * Correction_Agl= NULL;
	
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
	DatavsMC(TFile *file, std::string basename){
		MCEff   = new Efficiency(file,(basename + "_MC"  ).c_str());
		DataEff = new Efficiency(file,(basename + "_Data").c_str());

		DataEff_corr  = new Efficiency((basename + "_Data").c_str());

		latzones   = DataEff -> beforeR -> GetNbinsY();
		selections = DataEff -> beforeR -> GetNbinsZ();
	
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
	void Eval_Corrections();	

	TH1 * GetCorrection_R(int mc_type = -1){ 
		if(mc_type == -1) return Correction_R; else return (TH1F*)(((TH2*)Correction_R  )->ProjectionX("R DvsMC corr."  ,mc_type,mc_type))->Clone();
		};
	
	TH1 * GetCorrection_TOF(int mc_type = -1){ 
		if(mc_type == -1) return Correction_TOF; else return (TH1F*)(((TH2*)Correction_TOF)->ProjectionX("TOF DvsMC corr.",mc_type,mc_type))->Clone();
		};
	
	TH1 * GetCorrection_NaF(int mc_type = -1){ 
		if(mc_type == -1) return Correction_NaF; else return (TH1F*)(((TH2*)Correction_NaF)->ProjectionX("NaF DvsMC corr.",mc_type,mc_type))->Clone();
		};
	
	TH1 * GetCorrection_Agl(int mc_type = -1){ 
		if(mc_type == -1) return Correction_Agl; else return (TH1F*)(((TH2*)Correction_Agl)->ProjectionX("Agl DvsMC corr.",mc_type,mc_type))->Clone();		      };
};

void DatavsMC::Write(){
	MCEff 	-> Write();
	DataEff -> Write();
	return;
}

void DatavsMC::Eval_Corrected_DataEff(){

	if(selections == 1) {
		DataEff_corr -> beforeR   = (TH1F*)((TH1F *)((TH2F*)DataEff -> beforeR  ) -> ProjectionX((Basename + "1_R" ).c_str(),0,latzones)) -> Clone();
		DataEff_corr -> beforeTOF = (TH1F*)((TH1F *)((TH2F*)DataEff -> beforeTOF) -> ProjectionX((Basename + "1"   ).c_str(),0,latzones)) -> Clone();
		DataEff_corr -> beforeNaF = (TH1F*)((TH1F *)((TH2F*)DataEff -> beforeNaF) -> ProjectionX((Basename + "1NaF").c_str(),0,latzones)) -> Clone();
		DataEff_corr -> beforeAgl = (TH1F*)((TH1F *)((TH2F*)DataEff -> beforeAgl) -> ProjectionX((Basename + "1Agl").c_str(),0,latzones)) -> Clone();
	}

	else {
		DataEff_corr -> beforeR   = (TH2F*)((TH2F *)((TH3F*)DataEff -> beforeR  ) -> Project3D("xz")) -> Clone();
		DataEff_corr -> beforeTOF = (TH2F*)((TH2F *)((TH3F*)DataEff -> beforeTOF) -> Project3D("xz")) -> Clone();
		DataEff_corr -> beforeNaF = (TH2F*)((TH2F *)((TH3F*)DataEff -> beforeNaF) -> Project3D("xz")) -> Clone();
		DataEff_corr -> beforeAgl = (TH2F*)((TH2F *)((TH3F*)DataEff -> beforeAgl) -> Project3D("xz")) -> Clone();
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
		
		Correction_R  = DivideHisto( DataEff_corr -> effR   , MCEff -> effR    );	
                Correction_TOF=	DivideHisto( DataEff_corr -> effTOF , MCEff -> effTOF  );
	        Correction_NaF= DivideHisto( DataEff_corr -> effNaF , MCEff -> effNaF  );
                Correction_Agl= DivideHisto( DataEff_corr -> effAgl , MCEff -> effAgl  );

	return;
}



TH1 * DivideHisto(TH1 *Histo1, TH1 *Histo2){
	TH1 * result;
	if(Histo2->GetNbinsY()==1){
		result = (TH1F *)Histo1 -> Clone();		
		result -> Divide(Histo2);
	}
	else{
		result = (TH2F *)Histo2 -> Clone();
		for(int S = 0; S<Histo2->GetNbinsY();S++)
			for(int R=0;R<Histo1->GetNbinsX();R++){
			result -> SetBinContent (R+1,S+1,Histo1->GetBinContent(R+1)/Histo2->GetBinContent(R+1,S+1));
			}
	}
	return result;
}


TH1 * Correct_DataEff(std::string histoname,TH1 * Histo, TH1 * LATcorr){
	TH1 * Histo_corr;
	TH1 * temp = (TH1 *)Histo -> Clone();
  
	int selections = Histo ->GetNbinsZ();
	int latzones   = Histo ->GetNbinsY();
	
	if(selections == 1){
		//Histo_corr = new TH1F("","",Histo-> GetNbinsX(),0,Histo-> GetNbinsX());
		for(int R = 0;R < Histo ->GetNbinsX();R++) 
		     for(int lat =0; lat < latzones; lat++){	
				temp -> SetBinContent(R+1,lat+1,Histo->GetBinContent(R+1,lat+1)*LATcorr->GetBinContent(lat+1));
		}
		Histo_corr = (TH1F*)((TH1 *)((TH2*)temp ) -> ProjectionX(histoname.c_str(),0,latzones)) -> Clone();		
	}
	else {
		for(int S = 0; S < selections; S++)	
			for(int R = 0;R < Histo ->GetNbinsX();R++) 
				for(int lat =0; lat < latzones; lat++){	
					temp -> SetBinContent(R+1,lat+1,S+1,Histo->GetBinContent(R+1,lat+1,S+1)*LATcorr->GetBinContent(lat+1,S+1));
				}
		Histo_corr = (TH2F*)((TH3*)temp) -> Project3D("xz") -> Clone();	
	}
	
	return Histo_corr;
}


