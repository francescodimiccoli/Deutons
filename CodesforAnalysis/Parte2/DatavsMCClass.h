using namespace std;

TH1 * Correct_DataEff(std::string histoname,TH1 * Histo, TH1 * LATcorr);

class DatavsMC
{

public:
	Efficiency * MCEff;
	Efficiency * DataEff;

	//LAT corr.
	TH1 * LATcorr_R  = NULL; 	
	TH1 * LATcorr_TOF= NULL;
	TH1 * LATcorr_NaF= NULL;
	TH1 * LATcorr_Agl= NULL; 	
	
        Efficiency * DataEff_corr;
        
	std::string Basename;

	int latzones  = 0;
	int selections= 0;

	//creation constructors
	DatavsMC(std::string basename, int n){
		
		MCEff    = new Efficiency((basename + "_MC"  ).c_str());
		DataEff  = new Efficiency((basename + "_Data").c_str(),n);
		
		latzones = n;

		Basename = basename;
	} 
		
	DatavsMC(std::string basename, int n, int S){
		
		MCEff    = new Efficiency((basename + "_MC"  ).c_str());
		DataEff  = new Efficiency((basename + "_Data").c_str(),n);
		
		latzones   = n;
		selections = S;
	
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
	

};

void DatavsMC::Write(){
	MCEff 	-> Write();
	DataEff -> Write();
	return;
}

void DatavsMC::Eval_Corrected_DataEff(){

	DataEff_corr -> beforeR   = (TH1*)((TH1 *)((TH2*)DataEff -> beforeR  ) -> ProjectionX((Basename + "1_R" ).c_str(),0,latzones)) -> Clone();
	DataEff_corr -> beforeTOF = (TH1*)((TH1 *)((TH2*)DataEff -> beforeTOF) -> ProjectionX((Basename + "1"   ).c_str(),0,latzones)) -> Clone();
	DataEff_corr -> beforeNaF = (TH1*)((TH1 *)((TH2*)DataEff -> beforeNaF) -> ProjectionX((Basename + "1NaF").c_str(),0,latzones)) -> Clone();
	DataEff_corr -> beforeAgl = (TH1*)((TH1 *)((TH2*)DataEff -> beforeAgl) -> ProjectionX((Basename + "1Agl").c_str(),0,latzones)) -> Clone();

	if(!LATcorr_R   ) cout<<"ERROR: Lat. corr for R   histos not assigned"<<endl;
	if(!LATcorr_TOF ) cout<<"ERROR: Lat. corr for TOF histos not assigned"<<endl;
	if(!LATcorr_NaF ) cout<<"ERROR: Lat. corr for NaF histos not assigned"<<endl;
	if(!LATcorr_Agl ) cout<<"ERROR: Lat. corr for Agl histos not assigned"<<endl;
	
	DataEff_corr -> afterR   =  Correct_DataEff( (Basename + "2_R" ),	DataEff -> beforeR  ,	LATcorr_R   		);
        /*DataEff_corr -> afterTOF =  Correct_DataEff( (Basename + "2"   ),	DataEff -> beforeTOF,	LATcorr_TOF 		);
        DataEff_corr -> afterNaF =  Correct_DataEff( (Basename + "2NaF"),	DataEff -> beforeNaF,	LATcorr_NaF 		);
        DataEff_corr -> afterAgl =  Correct_DataEff( (Basename + "2Agl"),	DataEff -> beforeAgl,	LATcorr_Agl 		);
	*/
	return;
}


void DatavsMC::Eval_DandMC_Eff(){
	DatavsMC::Eval_Corrected_DataEff();	
	
	DataEff_corr -> Eval_Efficiency();
	MCEff 	     -> Eval_Efficiency();	
	
	return;
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
		Histo_corr = (TH1*)((TH1 *)((TH2*)temp ) -> ProjectionX(histoname.c_str(),0,latzones)) -> Clone();		
	}
	else {
		for(int S = 0; S < selections; S++)	
			for(int R = 0;R < Histo ->GetNbinsX();R++) 
				for(int lat =0; lat < latzones; lat++){	
					temp -> SetBinContent(R+1,lat+1,S+1,Histo->GetBinContent(R+1,lat+1,S+1)*LATcorr->GetBinContent(lat+1,S+1));
				}
		Histo_corr = (TH1*)((TH3*)temp) -> Project3D("xz") -> Clone();	
	}
	return Histo_corr;
}


