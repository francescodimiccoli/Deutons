using namespace std;

TH1 * Correct_DataEff(std::string histoname,TH1 * Histo, TH1 * LATcorr,int s=0);
TH1 * DivideHisto(TH1 *Histo1, TH1 *Histo2);
void ScanCombinations(std::vector<std::vector<TH1 *>> SystError , int i, std::vector<int> comb,TH1 * SystCorr );
TH1 * ExtractError( TH1 * SystCorr, TH1 * Correction) ;
	
class DatavsMC
{

	public:
		//LAT corr.
		std::vector<TH1 *>  LATcorr_R  ; 	
		std::vector<TH1 *>  LATcorr_TOF;
		std::vector<TH1 *>  LATcorr_NaF;
		std::vector<TH1 *>  LATcorr_Agl; 	

		std::vector<TH1 *> Correction_R  ;
		std::vector<TH1 *> Correction_TOF;
		std::vector<TH1 *> Correction_NaF;
		std::vector<TH1 *> Correction_Agl;

		std::vector<std::vector<TH1 *>> SystError_R	; 
		std::vector<std::vector<TH1 *>> SystError_TOF	;	 
		std::vector<std::vector<TH1 *>> SystError_NaF	; 
		std::vector<std::vector<TH1 *>> SystError_Agl	; 
		

		TH1 * SystCorr_R  ;
		TH1 * SystCorr_TOF; 
		TH1 * SystCorr_NaF; 
		TH1 * SystCorr_Agl; 	

		TH1 * SystErr_R  =NULL;
                TH1 * SystErr_TOF=NULL;
                TH1 * SystErr_NaF=NULL;
                TH1 * SystErr_Agl=NULL;

		TH1 * StatErr_R  =NULL;
                TH1 * StatErr_TOF=NULL;
                TH1 * StatErr_NaF=NULL;
                TH1 * StatErr_Agl=NULL;


		std::string Basename;

		int latzones  = 0;
		int mc_types  = 0;
		
		std::vector<Efficiency *> MCEff;
		std::vector<Efficiency *> DataEff;

		std::vector<Efficiency *> DataEff_corr;


		//creation constructor

		DatavsMC(std::string basename, int n=1, int mcs=1,int systerr=1){
			for(int l=0;l<systerr;l++){
					if(mcs == 1){ 
						MCEff  .push_back(new Efficiency((basename + to_string(l) + "_MC"  ).c_str()));
						DataEff.push_back(new Efficiency((basename + to_string(l) + "_Data").c_str(),n));
					}
					else{
						MCEff  .push_back(new Efficiency((basename + to_string(l) +"_MC"  ).c_str(),mcs));  	
						DataEff.push_back(new Efficiency((basename + to_string(l) +"_Data").c_str(),n));
					}

			}

			latzones   = n;
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
			mc_types   = mcs;

			Basename = basename;

		}	

		
		//cloning constructor
		DatavsMC(DatavsMC * Cloned,std::string basename){
			MCEff .insert(MCEff.end(),Cloned->MCEff.begin(),Cloned->MCEff.end());
			DataEff .insert(DataEff.end(),Cloned->DataEff.begin(),Cloned->DataEff.end());
			
			DataEff_corr .insert(DataEff_corr.end(),Cloned->DataEff_corr.begin(),Cloned->DataEff_corr.end());
			
			Correction_R   .insert(Correction_R  .end(),Cloned-> Correction_R  .begin(),Cloned->Correction_R  .end());  	
			Correction_TOF .insert(Correction_TOF.end(),Cloned-> Correction_TOF.begin(),Cloned->Correction_TOF.end());  	
			Correction_NaF .insert(Correction_NaF.end(),Cloned-> Correction_NaF.begin(),Cloned->Correction_NaF.end());  	
			Correction_Agl .insert(Correction_Agl.end(),Cloned-> Correction_Agl.begin(),Cloned->Correction_Agl.end());  	
			
			
			SystError_R.  insert( SystError_R  .end(), Cloned->SystError_R  .begin(), Cloned->SystError_R  .end() );
        		SystError_TOF.insert( SystError_TOF.end(), Cloned->SystError_TOF.begin(), Cloned->SystError_TOF.end() );
        		SystError_NaF.insert( SystError_NaF.end(), Cloned->SystError_NaF.begin(), Cloned->SystError_NaF.end() );
        		SystError_Agl.insert( SystError_Agl.end(), Cloned->SystError_Agl.begin(), Cloned->SystError_Agl.end() );	

			Basename = basename;
			
			latzones = Cloned -> latzones;
			mc_types = Cloned -> mc_types;
				
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

		void Extract_SystError() ;
		void Eval_Corrected_DataEff(int sel=0);
		void Eval_DandMC_Eff(int s=0);
		void DivideHisto(TH1 *Histo1, TH1 *Histo2, TH1 * Correction);
		void Eval_Corrections();	
		void Eval_FittedCorrections();
		void Eval_SystError();
		void Initialize_SystError();

		//method for multiplication of different corrections
		void ComposeCorrection(DatavsMC * Factor,int i=10,int j=10);

		TH1 * GetCorrection_R(int l=0)  {EvalError(Correction_R[l]  ,SystErr_R  ,StatErr_R  );return Correction_R[l]  ;};
		TH1 * GetCorrection_TOF(int l=0){EvalError(Correction_TOF[l],SystErr_TOF,StatErr_TOF);return Correction_TOF[l];};
		TH1 * GetCorrection_NaF(int l=0){EvalError(Correction_NaF[l],SystErr_NaF,StatErr_NaF);return Correction_NaF[l];};
		TH1 * GetCorrection_Agl(int l=0){EvalError(Correction_Agl[l],SystErr_Agl,StatErr_Agl);return Correction_Agl[l];}; 

		TH1 * GetSystPlot_R()   {return SystCorr_R    ;};
		TH1 * GetSystPlot_TOF() {return SystCorr_TOF  ;};
		TH1 * GetSystPlot_NaF() {return SystCorr_NaF  ;};
		TH1 * GetSystPlot_Agl() {return SystCorr_Agl  ;};
		
		TH1 * GetSystErr_R()   {return SystErr_R    ;};
		TH1 * GetSystErr_TOF() {return SystErr_TOF  ;};
		TH1 * GetSystErr_NaF() {return SystErr_NaF  ;};
		TH1 * GetSystErr_Agl() {return SystErr_Agl  ;};
		
		TH1 * GetStatErr_R()   {return StatErr_R    ;};
		TH1 * GetStatErr_TOF() {return StatErr_TOF  ;};
		TH1 * GetStatErr_NaF() {return StatErr_NaF  ;};
		TH1 * GetStatErr_Agl() {return StatErr_Agl  ;};
		
		TH1 * GetMCEff_R(int l=0)   {return (TH1*)MCEff[l]->effR  ->Clone()     ;};
		TH1 * GetMCEff_TOF(int l=0) {return (TH1*)MCEff[l]->effTOF->Clone()     ;};
		TH1 * GetMCEff_NaF(int l=0) {return (TH1*)MCEff[l]->effNaF->Clone()     ;};
		TH1 * GetMCEff_Agl(int l=0) {return (TH1*)MCEff[l]->effAgl->Clone()     ;};
		
		TH1 * GetDataEff_R(int l=0)   {return (TH1*)DataEff_corr[l]->effR  ->Clone()   ;};
		TH1 * GetDataEff_TOF(int l=0) {return (TH1*)DataEff_corr[l]->effTOF->Clone()   ;};
		TH1 * GetDataEff_NaF(int l=0) {return (TH1*)DataEff_corr[l]->effNaF->Clone()   ;};
		TH1 * GetDataEff_Agl(int l=0) {return (TH1*)DataEff_corr[l]->effAgl->Clone()   ;};
		



};


void DatavsMC::Write(){
	if(MCEff.size()>1) cout<<"Write "<<MCEff.size()<<endl;

	for(int l=0; l< MCEff.size();l++) {
		MCEff[l] 	-> Write();
		DataEff[l]      -> Write();
		}	
	return;
	
}

void DatavsMC::ComposeCorrection(DatavsMC * Factor,int i,int j){

	
	SystError_R.  insert( SystError_R  .end(), Factor->SystError_R  .begin(), Factor->SystError_R  .end() );
        SystError_TOF.insert( SystError_TOF.end(), Factor->SystError_TOF.begin(), Factor->SystError_TOF.end() );
        SystError_NaF.insert( SystError_NaF.end(), Factor->SystError_NaF.begin(), Factor->SystError_NaF.end() );
        SystError_Agl.insert( SystError_Agl.end(), Factor->SystError_Agl.begin(), Factor->SystError_Agl.end() );


	TH1 * C_R   = (TH1 *) Correction_R  [i]   -> Clone();
	TH1 * C_TOF = (TH1 *) Correction_TOF[i]   -> Clone();
	TH1 * C_NaF = (TH1 *) Correction_NaF[i]   -> Clone();
	TH1 * C_Agl = (TH1 *) Correction_Agl[i]   -> Clone();

	C_R   -> Multiply(Factor->Correction_R[j]   );
	C_TOF -> Multiply(Factor->Correction_TOF[j] );
        C_NaF -> Multiply(Factor->Correction_NaF[j] );
        C_Agl -> Multiply(Factor->Correction_Agl[j] );

	Correction_R  [i] = (TH1 *)C_R  -> Clone();
        Correction_TOF[i] = (TH1 *)C_TOF-> Clone();
        Correction_NaF[i] = (TH1 *)C_NaF-> Clone();
        Correction_Agl[i] = (TH1 *)C_Agl-> Clone();

	return;
}



void DatavsMC::Eval_Corrected_DataEff(int sel){

	for(int l=0;l<DataEff.size();l++){
		DataEff_corr[l] -> beforeR   = ProjectionXtoTH1F((TH2F*)DataEff[l] -> beforeR  , (Basename + to_string(l) + "1_R" ),0,latzones);
		DataEff_corr[l] -> beforeTOF = ProjectionXtoTH1F((TH2F*)DataEff[l] -> beforeTOF, (Basename + to_string(l) + "1"   ),0,latzones);
		DataEff_corr[l] -> beforeNaF = ProjectionXtoTH1F((TH2F*)DataEff[l] -> beforeNaF, (Basename + to_string(l) +"1NaF"),0,latzones);
		DataEff_corr[l] -> beforeAgl = ProjectionXtoTH1F((TH2F*)DataEff[l] -> beforeAgl, (Basename + to_string(l) +"1Agl"),0,latzones);

		if(!LATcorr_R[l]   ) cout<<"ERROR: Lat. corr for R   histos not assigned"<<endl;
		if(!LATcorr_TOF[l] ) cout<<"ERROR: Lat. corr for TOF histos not assigned"<<endl;
		if(!LATcorr_NaF[l] ) cout<<"ERROR: Lat. corr for NaF histos not assigned"<<endl;
		if(!LATcorr_Agl[l] ) cout<<"ERROR: Lat. corr for Agl histos not assigned"<<endl;

		DataEff_corr[l] -> afterR   =  Correct_DataEff( (Basename + to_string(l) +"2_R" ),	DataEff[l] -> afterR  ,	LATcorr_R[l]  ,sel 		);
		DataEff_corr[l] -> afterTOF =  Correct_DataEff( (Basename + to_string(l) +"2"   ),	DataEff[l] -> afterTOF,	LATcorr_TOF[l],sel 		);
		DataEff_corr[l] -> afterNaF =  Correct_DataEff( (Basename + to_string(l) +"2NaF"),	DataEff[l] -> afterNaF,	LATcorr_NaF[l],sel 		);
		DataEff_corr[l] -> afterAgl =  Correct_DataEff( (Basename + to_string(l) +"2Agl"),	DataEff[l] -> afterAgl,	LATcorr_Agl[l],sel 		);
	}
	return;
}


void DatavsMC::Eval_DandMC_Eff(int sel){

	DatavsMC::Eval_Corrected_DataEff(sel);	

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

		StatErr_R  =SetErrors(Correction_R[0]);
                StatErr_TOF=SetErrors(Correction_TOF[0]);
                StatErr_NaF=SetErrors(Correction_NaF[0]);
                StatErr_Agl=SetErrors(Correction_Agl[0]);

                StatErr_R  ->SetName((Basename + "_statR"             ).c_str());
                StatErr_TOF->SetName((Basename + "_statTOF"           ).c_str());
                StatErr_NaF->SetName((Basename + "_statNaF"           ).c_str());
                StatErr_Agl->SetName((Basename + "_statAgl"           ).c_str());	




	return;
}


void DatavsMC::Eval_FittedCorrections(){

	for(int l=0;l<DataEff.size();l++){


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


		StatErr_R  =SetErrors(Correction_R[0]); 
                StatErr_TOF=SetErrors(Correction_TOF[0]);
		StatErr_NaF=SetErrors(Correction_NaF[0]);
                StatErr_Agl=SetErrors(Correction_Agl[0]);

		StatErr_R  ->SetName((Basename + "_statR"             ).c_str());
	        StatErr_TOF->SetName((Basename + "_statTOF"           ).c_str());
        	StatErr_NaF->SetName((Basename + "_statNaF"           ).c_str());
        	StatErr_Agl->SetName((Basename + "_statAgl"           ).c_str());


	return;
}




void DatavsMC::DivideHisto(TH1 *Histo1, TH1 *Histo2,TH1 * Correction){
		if(mc_types ==1) {
			for(int iR=0;iR<Histo2->GetNbinsX();iR++){ 
				if(Histo2 -> GetBinContent(iR+1)<1&&Histo2 -> GetBinContent(iR+1)>0)
					Correction -> SetBinContent (iR+1,Histo1 -> GetBinContent(iR+1)/(float)Histo2 -> GetBinContent(iR+1));
				Correction -> SetBinError(iR+1,pow(pow(Histo1 -> GetBinError(iR+1),2)+pow(Histo2 -> GetBinError(iR+1),2),0.5)*Correction->GetBinContent(iR+1));
			}
		}
		else{
			for(int mc_type=0;mc_type<mc_types;mc_type++){
				for(int iR=0;iR<Histo2->GetNbinsX();iR++){
					if(Histo2 -> GetBinContent(iR+1,mc_type+1)<1&&Histo2 -> GetBinContent(iR+1,mc_type+1)>0){
						Correction -> SetBinContent (iR+1,mc_type+1,Histo1 -> GetBinContent(iR+1)/(float)Histo2 -> GetBinContent(iR+1,mc_type+1));
						Correction -> SetBinError(iR+1,mc_type+1,pow(pow(Histo1 -> GetBinError(iR+1),2)+pow(Histo2 -> GetBinError(iR+1,mc_type+1),2),0.5)*Correction->GetBinContent(iR+1,mc_type+1)); }
					else{
						Correction -> SetBinContent (iR+1,mc_type+1,0);
						Correction -> SetBinError(iR+1,mc_type+1,0); 
					}	
				}

			}
		}

	return; 
}


void DatavsMC::Initialize_SystError(){

	if(SystError_R  .size()==0) SystError_R  .push_back(Correction_R  );	
	if(SystError_TOF.size()==0) SystError_TOF.push_back(Correction_TOF);	
	if(SystError_NaF.size()==0) SystError_NaF.push_back(Correction_NaF);	
	if(SystError_Agl.size()==0) SystError_Agl.push_back(Correction_Agl);	


	return;
}



void DatavsMC::Eval_SystError(){


	if(mc_types==1){
		SystCorr_R   = new TH2F((Basename + "SystCorr_R").c_str()  ,(Basename + "SystCorr_R").c_str()  ,Correction_R[0]    ->GetNbinsX(),0,Correction_R[0]    ->GetNbinsX(),400,0.8,1.6);
		SystCorr_TOF = new TH2F((Basename + "SystCorr_TOF").c_str(),(Basename + "SystCorr_TOF").c_str(),Correction_TOF[0]  ->GetNbinsX(),0,Correction_TOF[0]  ->GetNbinsX(),400,0.8,1.6);
		SystCorr_NaF = new TH2F((Basename + "SystCorr_NaF").c_str(),(Basename + "SystCorr_NaF").c_str(),Correction_NaF[0]  ->GetNbinsX(),0,Correction_NaF[0]  ->GetNbinsX(),400,0.8,1.6);
		SystCorr_Agl = new TH2F((Basename + "SystCorr_Agl").c_str(),(Basename + "SystCorr_Agl").c_str(),Correction_Agl[0]  ->GetNbinsX(),0,Correction_Agl[0]  ->GetNbinsX(),400,0.8,1.6);
	}


	else{
		SystCorr_R   = new TH3F((Basename + "SystCorr_R").c_str()  ,(Basename + "SystCorr_R").c_str()  ,Correction_R[0] ->GetNbinsX(),0,Correction_R[0]   ->GetNbinsX(),400,0.8,1.6,mc_types,0,mc_types);
		SystCorr_TOF = new TH3F((Basename + "SystCorr_TOF").c_str(),(Basename + "SystCorr_TOF").c_str(),Correction_TOF[0]->GetNbinsX(),0,Correction_TOF[0]->GetNbinsX(),400,0.8,1.6,mc_types,0,mc_types);
		SystCorr_NaF = new TH3F((Basename + "SystCorr_NaF").c_str(),(Basename + "SystCorr_NaF").c_str(),Correction_NaF[0]->GetNbinsX(),0,Correction_NaF[0]->GetNbinsX(),400,0.8,1.6,mc_types,0,mc_types);
		SystCorr_Agl = new TH3F((Basename + "SystCorr_Agl").c_str(),(Basename + "SystCorr_Agl").c_str(),Correction_Agl[0]->GetNbinsX(),0,Correction_Agl[0]->GetNbinsX(),400,0.8,1.6,mc_types,0,mc_types);

	}


	std::vector<int> comb_R;
	std::vector<int> comb_TOF;
	std::vector<int> comb_NaF;
	std::vector<int> comb_Agl;

	for(int i=0;i<SystError_R.size();i++) {
		 comb_R.push_back(0);
		 comb_TOF.push_back(0);
	 	 comb_NaF.push_back(0);
		 comb_Agl.push_back(0);
	}
	
	ScanCombinations(SystError_R  ,0,comb_R  ,SystCorr_R);
	ScanCombinations(SystError_TOF,0,comb_TOF,SystCorr_TOF);
	ScanCombinations(SystError_NaF,0,comb_NaF,SystCorr_NaF);
	ScanCombinations(SystError_Agl,0,comb_Agl,SystCorr_Agl);

	
	Extract_SystError();
	

	return;


}


void ScanCombinations(std::vector<std::vector<TH1 *>> SystError , int i, std::vector<int> comb,TH1 * SystCorr ){
		
		TH1 * TotalCorr;
		for(int c=0;c<SystError[i].size();c++){
			comb[i]=c;
			
			for (int z=0;z<comb.size();z++)
		
			TotalCorr = (TH1*)SystError[0][comb[0]]->Clone();
			for(int z=1;z<comb.size();z++) TotalCorr->Multiply(SystError[z][comb[z]]);
	
			if(SystCorr->GetNbinsZ()>1){		
				for(int l=0;l<SystCorr->GetNbinsX();l++)
					for(int mc_type=0;mc_type<SystCorr->GetNbinsZ();mc_type++)
						((TH3*)SystCorr)->Fill(l,TotalCorr->GetBinContent(l+1,mc_type),mc_type);	

			}
			else{
				for(int l=0;l<SystCorr->GetNbinsX();l++){
					((TH2*)SystCorr)->Fill(l,TotalCorr->GetBinContent(l+1));
				}

			} 

			if(i+1<SystError.size()) ScanCombinations(SystError,i+1,comb,SystCorr);			
		}
		return;
}





void DatavsMC::Extract_SystError(){

	SystErr_R   = ExtractError(SystCorr_R   , Correction_R[0]  );
	SystErr_TOF = ExtractError(SystCorr_TOF , Correction_TOF[0]);
	SystErr_NaF = ExtractError(SystCorr_NaF , Correction_NaF[0]);
	SystErr_Agl = ExtractError(SystCorr_Agl , Correction_Agl[0]);

	SystErr_R  ->SetName((Basename + "_systR"             ).c_str());
	SystErr_TOF->SetName((Basename + "_systTOF"           ).c_str());
        SystErr_NaF->SetName((Basename + "_systNaF"           ).c_str());
        SystErr_Agl->SetName((Basename + "_systAgl"           ).c_str());



			
	return;

}


TH1 * ExtractError( TH1 * SystCorr, TH1* Correction) {

        TH1* SystError =(TH1 *) Correction -> Clone();

	int mc_types=Correction->GetNbinsY();

        for (int iR=0;iR<SystCorr->GetNbinsX();iR++)    {
                for(int m=0;m<mc_types;m++){
			bool reverse=true;
                        TH1F* Slice = Extract_Bin(SystCorr, iR , m,reverse);
                        SystError->SetBinContent(iR+1,m+1, Slice->GetStdDev());
			cout<<Slice->GetEntries()<<" "<<Slice->GetStdDev()<<endl;
		}
        }

        return SystError;

}




TH1 * Correct_DataEff(std::string histoname,TH1 * Histo, TH1 * LATcorr,int sel){
	TH1 * Histo_corr;

	int selections = Histo ->GetNbinsZ();
	int latzones   = Histo ->GetNbinsY();

	if(selections == 1){
		TH2F* temp = (TH2F *)Histo -> Clone();
		//Histo_corr = new TH1F("","",Histo-> GetNbinsX(),0,Histo-> GetNbinsX());
		for(int iR = 0;iR < Histo ->GetNbinsX();iR++) 
			for(int lat =0; lat < latzones; lat++){	
				temp -> SetBinContent(iR+1,lat+1,Histo->GetBinContent(iR+1,lat+1)*LATcorr->GetBinContent(lat+1,sel+1));
			}
		Histo_corr = (TH1F*) ProjectionXtoTH1F( temp, histoname,0,latzones );
	}

	return Histo_corr;
}




