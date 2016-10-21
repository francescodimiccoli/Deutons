using namespace std;

TH1 * Correct_DataEff(std::string histoname,TH1 * Histo, TH1 * LATcorr);
TH1 * DivideHisto(TH1 *Histo1, TH1 *Histo2);

class DatavsMC
{

	private:
		//LAT corr.
		std::vector<TH1 *> LATcorr_R  ; 	
		std::vector<TH1 *> LATcorr_TOF;
		std::vector<TH1 *>  LATcorr_NaF;
		std::vector<TH1 *>  LATcorr_Agl; 	

		std::vector<TH1 *> Correction_R  ;
		std::vector<TH1 *> Correction_TOF;
		std::vector<TH1 *> Correction_NaF;
		std::vector<TH1 *> Correction_Agl;

		TH1 * Syst_R  ;	
		TH1 * Syst_TOF;
		TH1 * Syst_NaF;
		TH1 * Syst_Agl;

		TH1 * Abs_R  ;	
		TH1 * Abs_TOF;
		TH1 * Abs_NaF;
		TH1 * Abs_Agl;


		std::string Basename;

		int latzones  = 0;
		int selections= 0;
		int mc_types  = 0;
	public:
		std::vector<Efficiency *> MCEff;
		std::vector<Efficiency *> DataEff;

		std::vector<Efficiency *> DataEff_corr;

		//creation constructor

		DatavsMC(std::string basename, int n=1, int S=1, int mcs=1,int systerr=1){
			for(int l=0;l<systerr;l++){
				if(S==1){	
					if(mcs == 1){ 
						MCEff  .push_back(new Efficiency((basename + to_string(l) + "_MC"  ).c_str()));
						DataEff.push_back(new Efficiency((basename + to_string(l) + "_Data").c_str(),n));
					}
					else{
						MCEff  .push_back(new Efficiency((basename + to_string(l) +"_MC"  ).c_str(),mcs));  	
						DataEff.push_back(new Efficiency((basename + to_string(l) +"_Data").c_str(),n));
					}

				}
				else{	
					if(mcs == 1){ 
						MCEff  .push_back(new Efficiency((basename + to_string(l) +"_MC"  ).c_str(),S));
						DataEff.push_back(new Efficiency((basename + to_string(l) +"_Data").c_str(),n,S));
					}
					else{
						MCEff  .push_back(new Efficiency((basename + to_string(l) +"_MC"  ).c_str(),mcs,S));  	
						DataEff.push_back(new Efficiency((basename + to_string(l) +"_Data").c_str(),n,S));
					}

				}
			}

			latzones   = n;
			selections = S;
			mc_types   = mcs;	

			Basename = basename;

		}

		//reading constructors
		DatavsMC(TFile *file, std::string basename, int mcs=1,int systerr=1){	
			for(int l=0;l<systerr;l++){

				MCEff   .push_back(new Efficiency(file,(basename + to_string(l) + "_MC"  ).c_str()) );
				DataEff .push_back(new Efficiency(file,(basename + to_string(l) + "_Data").c_str()) );

				DataEff_corr.push_back(new Efficiency((basename  + to_string(l) + "_Data").c_str()) );

				Correction_R    .push_back(  (TH1 *) MCEff[l] -> beforeR   -> Clone());	
				Correction_TOF  .push_back(  (TH1 *) MCEff[l] -> beforeTOF -> Clone());
				Correction_NaF  .push_back(  (TH1 *) MCEff[l] -> beforeNaF -> Clone());
				Correction_Agl  .push_back(  (TH1 *) MCEff[l] -> beforeAgl -> Clone());

			}
			latzones   = DataEff[0] -> beforeR -> GetNbinsY();
			selections = DataEff[0] -> beforeR -> GetNbinsZ();
			mc_types   = mcs;

			Basename = basename;

		}	

		void Write();
		void Assign_LatCorr(TH1 * LatcorrR, TH1 * LatcorrTOF, TH1 * LatcorrNaF, TH1 * LatcorrAgl)
		{
			for(int l=0; l< MCEff.size();l++) {
	
				LATcorr_R   .push_back( (TH1*)LatcorrR   ->Clone());
				LATcorr_TOF .push_back( (TH1*)LatcorrTOF ->Clone());
				LATcorr_NaF .push_back( (TH1*)LatcorrNaF ->Clone());
				LATcorr_Agl .push_back( (TH1*)LatcorrAgl ->Clone());
			}

		};

		void Eval_Corrected_DataEff();
		void Eval_DandMC_Eff();
		void DivideHisto(TH1 *Histo1, TH1 *Histo2, TH1 * Correction);
		void Eval_Corrections();	
		void Eval_FittedCorrections();

		TH1 * GetCorrection_R(int l=0)  { return Correction_R[l]  ;};
		TH1 * GetCorrection_TOF(int l=0){ return Correction_TOF[l];};
		TH1 * GetCorrection_NaF(int l=0){ return Correction_NaF[l];};
		TH1 * GetCorrection_Agl(int l=0){ return Correction_Agl[l];}; 
};


void DatavsMC::Write(){
	if(MCEff.size()>1) cout<<"Write "<<MCEff.size()<<endl;

	for(int l=0; l< MCEff.size();l++) {
		cout<<l<<endl;
		MCEff[l] 	-> Write();
		DataEff[l]      -> Write();
		}	
	return;
	
}

void DatavsMC::Eval_Corrected_DataEff(){

	for(int l=0;l<DataEff.size();l++){
		cout<<DataEff.size()<<endl;
		if(selections == 1) {
			DataEff_corr[l] -> beforeR   = ProjectionXtoTH1F((TH2F*)DataEff[l] -> beforeR  , (Basename + to_string(l) + "1_R" ),0,latzones);
			DataEff_corr[l] -> beforeTOF = ProjectionXtoTH1F((TH2F*)DataEff[l] -> beforeTOF, (Basename + to_string(l) + "1"   ),0,latzones);
			DataEff_corr[l] -> beforeNaF = ProjectionXtoTH1F((TH2F*)DataEff[l] -> beforeNaF, (Basename + to_string(l) +"1NaF"),0,latzones);
			DataEff_corr[l] -> beforeAgl = ProjectionXtoTH1F((TH2F*)DataEff[l] -> beforeAgl, (Basename + to_string(l) +"1Agl"),0,latzones);
		}
		else {
			DataEff_corr[l] -> beforeR   = (TH2F*)((TH2F *)((TH3F*)DataEff[l] -> beforeR  ) -> Project3D("zx")) -> Clone();
			DataEff_corr[l] -> beforeTOF = (TH2F*)((TH2F *)((TH3F*)DataEff[l] -> beforeTOF) -> Project3D("zx")) -> Clone();
			DataEff_corr[l] -> beforeNaF = (TH2F*)((TH2F *)((TH3F*)DataEff[l] -> beforeNaF) -> Project3D("zx")) -> Clone();
			DataEff_corr[l] -> beforeAgl = (TH2F*)((TH2F *)((TH3F*)DataEff[l] -> beforeAgl) -> Project3D("zx")) -> Clone();
		}
		cout<<LATcorr_R.size() <<endl;
		if(!LATcorr_R[l]   ) cout<<"ERROR: Lat. corr for R   histos not assigned"<<endl;
		if(!LATcorr_TOF[l] ) cout<<"ERROR: Lat. corr for TOF histos not assigned"<<endl;
		if(!LATcorr_NaF[l] ) cout<<"ERROR: Lat. corr for NaF histos not assigned"<<endl;
		if(!LATcorr_Agl[l] ) cout<<"ERROR: Lat. corr for Agl histos not assigned"<<endl;
		cout<<DataEff.size()<<endl;
		DataEff_corr[l] -> afterR   =  Correct_DataEff( (Basename + to_string(l) +"2_R" ),	DataEff[l] -> afterR  ,	LATcorr_R[l]   		);
		DataEff_corr[l] -> afterTOF =  Correct_DataEff( (Basename + to_string(l) +"2"   ),	DataEff[l] -> afterTOF,	LATcorr_TOF[l] 		);
		DataEff_corr[l] -> afterNaF =  Correct_DataEff( (Basename + to_string(l) +"2NaF"),	DataEff[l] -> afterNaF,	LATcorr_NaF[l] 		);
		DataEff_corr[l] -> afterAgl =  Correct_DataEff( (Basename + to_string(l) +"2Agl"),	DataEff[l] -> afterAgl,	LATcorr_Agl[l] 		);
		cout<<DataEff.size()<<endl;
	}
	return;
}


void DatavsMC::Eval_DandMC_Eff(){

	DatavsMC::Eval_Corrected_DataEff();	

	for(int l=0; l< MCEff.size();l++) {
		DataEff_corr[l] -> Eval_Efficiency();
		MCEff[l] 	-> Eval_Efficiency();	
	}	
	return;
}


void DatavsMC::Eval_Corrections(){

	for(int l=0;l<DataEff.size();l++){

		DivideHisto( DataEff_corr[l] -> effR   , MCEff[l] -> effR  , Correction_R[l]   );	
		DivideHisto( DataEff_corr[l] -> effTOF , MCEff[l] -> effTOF, Correction_TOF[l] );
		DivideHisto( DataEff_corr[l] -> effNaF , MCEff[l] -> effNaF, Correction_NaF[l] );
		DivideHisto( DataEff_corr[l] -> effAgl , MCEff[l] -> effAgl, Correction_Agl[l] );
	}
	return;
}


void DatavsMC::Eval_FittedCorrections(){

	for(int l=0;l<DataEff.size();l++){


		DivideHisto( DataEff_corr[l] -> effR   , MCEff[l] -> effR  , Correction_R[l]   );
		DivideHisto( DataEff_corr[l] -> effTOF , MCEff[l] -> effTOF, Correction_TOF[l] );
		DivideHisto( DataEff_corr[l] -> effNaF , MCEff[l] -> effNaF, Correction_NaF[l] );
		DivideHisto( DataEff_corr[l] -> effAgl , MCEff[l] -> effAgl, Correction_Agl[l] );

		FitFunction * FitR   = new FitFunction( Correction_R[l]   ,4);
		FitFunction * FitTOF = new FitFunction( Correction_TOF[l] ,4);
		FitFunction * FitNaF = new FitFunction( Correction_NaF[l] ,4);
		FitFunction * FitAgl = new FitFunction( Correction_Agl[l] ,4);

		FitR  ->FitValues();
		FitTOF->FitValues();
		FitNaF->FitValues();
		FitAgl->FitValues();

		Correction_R[l]   =(TH1*)FitR  ->ReturnFittedValues();
		Correction_TOF[l] =(TH1*)FitTOF->ReturnFittedValues();
		Correction_NaF[l] =(TH1*)FitNaF->ReturnFittedValues();
		Correction_Agl[l] =(TH1*)FitAgl->ReturnFittedValues();
	}	

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




float SystErrorBin(int iR,std::vector<TH1 *> Correction,int Y=1,int Z=1){

	float syserr=0;
	for(int l=0;l<Correction.size();l++){
		syserr+=pow(Correction[10]->GetBinContent(iR+1,Y+1,Z+1)-Correction[l]->GetBinContent(iR+1,Y+1,Z+1) ,2)	
	}
	syserr/=Correction.size();
	syserr=pow(syserr,0.5);

	return syserr;
} 



float AbsErrorBin(int iR,std::vector<TH1 *> Correction,int Y=1,int Z=1){

	float abserr=0;
	abserr = abs(Correction[10]->GetBinContent(iR+1,Y+1,Z+1) - 1)/2;
	return abserr;

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
