using namespace std;


#ifndef ACCEPTANCE_H
#define ACCEPTANCE_H

class ACCEPTANCE
{
private:

	//binning
	static const int nbinsr=43;
	static const int nbinsToF=18;
	static const int nbinsNaF=18;
	static const int nbinsAgl=18;

	float binsR[nbinsr+1];
	float binsBetaTOF[19];
	float binsBetaNaF[19];
	float binsBetaAgl[19];
	
public:
	//MC parameters
	float trigrate;
	float Rmin;
	float Rmax;	

	//sel. eff.
	TH1 * Efficiency_R  ;
	TH1 * Efficiency_TOF;
	TH1 * Efficiency_NaF;
	TH1 * Efficiency_Agl;

	//gen. events
	TH1 * before_R   ,  * after_R  ; 
	TH1 * before_TOF ,  * after_TOF; 
	TH1 * before_NaF ,  * after_NaF; 
	TH1 * before_Agl ,  * after_Agl; 

	//LAT corr.
	TH1 * LATcorr_R  ;
	TH1 * LATcorr_TOF;
	TH1 * LATcorr_NaF;
	TH1 * LATcorr_Agl;

	//Weighted LAT corr.
	TH1 * LATcorrW_R;
	TH1 * LATcorrW_TOF;
	TH1 * LATcorrW_NaF;
	TH1 * LATcorrW_Agl;

	//Acceptance

	TH1 * MCAcceptance_R   , * CorrectedAcceptance_R   ,  * Geomag_Acceptance_R   ,  * Gen_Acceptance_R    ;
	TH1 * MCAcceptance_TOF , * CorrectedAcceptance_TOF ,  * Geomag_Acceptance_TOF ,  * Gen_Acceptance_TOF  ;
	TH1 * MCAcceptance_NaF , * CorrectedAcceptance_NaF ,  * Geomag_Acceptance_NaF ,  * Gen_Acceptance_NaF  ;
	TH1 * MCAcceptance_Agl , * CorrectedAcceptance_Agl ,  * Geomag_Acceptance_Agl ,  * Gen_Acceptance_Agl  ;


	// reading constructor
	ACCEPTANCE(TFile * file , std::string dirname , std::string basename, std::string effname,std::string latcorrname , std::string wlatcorrname, int n)
	{
		after_TOF   = (TH1 *)file->Get((basename + "1"   ).c_str());
		after_NaF   = (TH1 *)file->Get((basename + "1NaF").c_str());
		after_Agl   = (TH1 *)file->Get((basename + "1Agl").c_str());
		after_R     = (TH1 *)file->Get((basename + "1_R" ).c_str());

		before_TOF  = new TH2F((basename + "Trig"   ).c_str(),(basename + "Trig"   ).c_str(),nbinsToF, 0, nbinsToF,  n, 0, n);
		before_NaF  = new TH2F((basename + "TrigNaF").c_str(),(basename + "TrigNaF").c_str(),nbinsNaF, 0, nbinsNaF,  n, 0, n);
		before_Agl  = new TH2F((basename + "TrigAgl").c_str(),(basename + "TrigAgl").c_str(),nbinsAgl, 0, nbinsAgl,  n, 0, n);
		before_R    = new TH2F((basename + "Trig_R" ).c_str(),(basename + "Trig_R" ).c_str(),nbinsr,   0, nbinsr,    n, 0, n);

		Efficiency_R  = (TH1 *)file->Get(("/" + dirname +"/" +effname + "_EffR"           ).c_str());
		Efficiency_TOF= (TH1 *)file->Get(("/" + dirname +"/" +effname + "_EffTOF"         ).c_str());
		Efficiency_NaF= (TH1 *)file->Get(("/" + dirname +"/" +effname + "_EffNaF"         ).c_str());
		Efficiency_Agl= (TH1 *)file->Get(("/" + dirname +"/" +effname + "_EffAgl"         ).c_str());

		LATcorr_R  = (TH1 *)file->Get(("/" + dirname +"/" +latcorrname + "_LATcorrR_fit"   ).c_str());
		LATcorr_TOF= (TH1 *)file->Get(("/" + dirname +"/" +latcorrname + "_LATcorrTOF_fit" ).c_str());
		LATcorr_NaF= (TH1 *)file->Get(("/" + dirname +"/" +latcorrname + "_LATcorrNaF_fit" ).c_str());
		LATcorr_Agl= (TH1 *)file->Get(("/" + dirname +"/" +latcorrname + "_LATcorrAgl_fit" ).c_str());

		LATcorrW_R  = (TH1 *)file->Get(("/" + dirname +"/" +wlatcorrname + "_R"         ).c_str());
		LATcorrW_TOF= (TH1 *)file->Get(("/" + dirname +"/" +wlatcorrname + "_TOF"       ).c_str());
		LATcorrW_NaF= (TH1 *)file->Get(("/" + dirname +"/" +wlatcorrname + "_NaF"       ).c_str());
		LATcorrW_Agl= (TH1 *)file->Get(("/" + dirname +"/" +wlatcorrname + "_Agl"       ).c_str());

	}

	void Set_MC_Par( float Trigrate, float rmin, float rmax); 
	void Set_Binning(bool deutons=1);
	TH1 * Triggerbin(int n , TH1 * after, float trigrate, float bins[]);
	void Eval_Gen_Acceptance( int n);
	void Eval_MC_Acceptance();
	TH1 * Geomag_Acceptance(int n, TH1* MCAcceptance, TH1* LATcorr);
	void Eval_Geomag_Acceptance(int n);
	TH1 * Corrected_Acceptance(int n, TH1* MCAcceptance, TH1 *LATcorrW);
	void Eval_Corrected_Acceptance(int n);	

	void ApplyDvsMCcorrection(int n,int S, TH1* Correction, TH1* Corrected_Acc, TH1* Geom_Acc);
	
	void Apply_DvsMCcorrection_R  (TH1 * R_Correction  ,int n=1,int S=1){ if (LATcorrW_R   )  ApplyDvsMCcorrection(n,S,R_Correction  , CorrectedAcceptance_R   , Geomag_Acceptance_R   ); return;}
        void Apply_DvsMCcorrection_TOF(TH1 * TOF_Correction,int n=1,int S=1){ if (LATcorrW_TOF )  ApplyDvsMCcorrection(n,S,TOF_Correction, CorrectedAcceptance_TOF , Geomag_Acceptance_TOF ); return;}
        void Apply_DvsMCcorrection_NaF(TH1 * NaF_Correction,int n=1,int S=1){ if (LATcorrW_NaF )  ApplyDvsMCcorrection(n,S,NaF_Correction, CorrectedAcceptance_NaF , Geomag_Acceptance_NaF ); return;}
        void Apply_DvsMCcorrection_Agl(TH1 * Agl_Correction,int n=1,int S=1){ if (LATcorrW_Agl )  ApplyDvsMCcorrection(n,S,Agl_Correction, CorrectedAcceptance_Agl , Geomag_Acceptance_Agl ); return;}
 		

	void ApplyGlobalFactor(float factor, float error);
};

void ACCEPTANCE::Set_MC_Par( float Trigrate, float rmin, float rmax){
	trigrate = Trigrate;
	Rmin = rmin;
	Rmax = rmax;
	return;
}



void ACCEPTANCE::Set_Binning(bool deutons){

	int nbins_beta =   after_TOF -> GetNbinsX() + 1; //TOF,NaF,Agl: same number of bins

	if (deutons) {

		for(int i=0;i<nbins_beta; i++) {
			binsBetaTOF[i] = ToFDB.MomBin(i); 
			binsBetaNaF[i] = NaFDB.MomBin(i); 
			binsBetaAgl[i] = AglDB.MomBin(i); 	
		}

	} else { // protons

		for(int i=0;i<nbins_beta; i++) {
			binsBetaTOF[i] = ToFPB.MomBin(i); 
			binsBetaNaF[i] = NaFPB.MomBin(i); 
			binsBetaAgl[i] = AglPB.MomBin(i); 	
		}

	}


}


TH1 * ACCEPTANCE::Triggerbin(int n , TH1 * after , float trigrate, float bins[]){
	//estimation of gen. particles in every bin

	float eventiprova[n]; //counts total number of events
	TH1 * triggerbin = (TH1 *) after -> Clone();	
	float gen_bin=0;
	if(n>1)	
		for(int m=0;m<n;m++){
			eventiprova[m]   = ((TH2 *) after_R ) -> Integral(0,(int)after_R  ->GetNbinsX(),m+1,m+1);
			for(int i=0; i<(int)after ->GetNbinsX();i++ ){
				gen_bin = eventiprova[m]*(pow(trigrate,-1))*(log(bins[i+1])-log(bins[i]))/(log(Rmax)-log(Rmin));
				triggerbin -> SetBinContent(i+1,m+1,gen_bin);
				triggerbin -> SetBinError(i+1,m+1,pow(gen_bin,0.5));
			}
		}

	else
		for(int m=0;m<n;m++){
			eventiprova[m]   = ((TH1 *) after_R  ) -> Integral(0,(int)after_R   ->GetNbinsX());
			for(int i=0; i<(int)after ->GetNbinsX();i++ ){
				gen_bin = eventiprova[m]*(pow(trigrate,-1))*(log(bins[i+1])-log(bins[i]))/(log(Rmax)-log(Rmin));
				triggerbin -> SetBinContent(i+1,gen_bin);
				triggerbin -> SetBinError(i+1,pow(gen_bin,0.5));
			}
		}
	return triggerbin;
}



void ACCEPTANCE::Eval_Gen_Acceptance(int n){
	before_R    = ACCEPTANCE::Triggerbin( n, after_R   ,trigrate , PRB.RigBins().data());	
	before_TOF 	= ACCEPTANCE::Triggerbin( n, after_TOF ,trigrate , binsBetaTOF);
	before_NaF 	= ACCEPTANCE::Triggerbin( n, after_NaF ,trigrate , binsBetaNaF);	
	before_Agl 	= ACCEPTANCE::Triggerbin( n, after_Agl ,trigrate , binsBetaAgl);

	if(after_R)	 Gen_Acceptance_R     = (TH1 *)after_R		->Clone();
	if(after_TOF)	 Gen_Acceptance_TOF   = (TH1 *)after_TOF	->Clone();
	if(after_NaF)  	 Gen_Acceptance_NaF   = (TH1 *)after_NaF	->Clone();
	if(after_Agl) 	 Gen_Acceptance_Agl   = (TH1 *)after_Agl	->Clone();

	Gen_Acceptance_R    -> Divide ( before_R	) ;
	Gen_Acceptance_TOF  -> Divide ( before_TOF	) ;
	Gen_Acceptance_NaF  -> Divide ( before_NaF	) ;
	Gen_Acceptance_Agl  -> Divide ( before_Agl	) ;

	//geometric acc. factor
	Gen_Acceptance_R    -> Scale ( 47.78 ) ;
	Gen_Acceptance_TOF  -> Scale ( 47.78 ) ;
	Gen_Acceptance_NaF  -> Scale ( 47.78 ) ;
	Gen_Acceptance_Agl  -> Scale ( 47.78 ) ;

} 


void ACCEPTANCE::Eval_MC_Acceptance(){

	if(Gen_Acceptance_R  )  MCAcceptance_R   =(TH1 *) Gen_Acceptance_R  ->Clone();
	if(Gen_Acceptance_TOF)  MCAcceptance_TOF =(TH1 *) Gen_Acceptance_TOF->Clone();
	if(Gen_Acceptance_NaF)  MCAcceptance_NaF =(TH1 *) Gen_Acceptance_NaF->Clone();
	if(Gen_Acceptance_Agl)  MCAcceptance_Agl =(TH1 *) Gen_Acceptance_Agl->Clone();

	if(	Efficiency_R  	) MCAcceptance_R  -> Multiply (	Efficiency_R  	);
	if(	Efficiency_TOF	) MCAcceptance_TOF-> Multiply (	Efficiency_TOF	);
	if(	Efficiency_NaF	) MCAcceptance_NaF-> Multiply (	Efficiency_NaF	);
	if(	Efficiency_Agl	) MCAcceptance_Agl-> Multiply (	Efficiency_Agl	);

}

TH1 * ACCEPTANCE::Geomag_Acceptance(int n, TH1* MCAcceptance, TH1* LATcorr){
	TH1 * Geomag_Acceptance;
	float error=0;
	if(n>1){
		Geomag_Acceptance = new TH3F ("","",MCAcceptance->GetNbinsX(),0,MCAcceptance->GetNbinsX(),LATcorr->GetNbinsX(),0,LATcorr->GetNbinsX(),n,0,n);		      
		for(int m=0;m<n;m++)
			for(int iR=0;iR<MCAcceptance->GetNbinsX();iR++)
				for(int lat=1;lat<LATcorr->GetNbinsX();lat++){
					Geomag_Acceptance->SetBinContent(iR+1,lat+1,m+1,MCAcceptance->GetBinContent(iR+1,m+1)/LATcorr->GetBinContent(lat+1));
					error=pow(MCAcceptance->GetBinContent(iR+1,m+1)*LATcorr->GetBinError(lat+1),2);
				        error=pow(MCAcceptance->GetBinError(iR+1,m+1)*LATcorr->GetBinContent(lat+1),2);
					error=pow(error,0.5);
					Geomag_Acceptance->SetBinError(iR+1,lat+1,m+1,error);				
					}			
	}
	else{
		Geomag_Acceptance = new TH2F ("","",MCAcceptance->GetNbinsX(),0,MCAcceptance->GetNbinsX(),LATcorr->GetNbinsX(),0,LATcorr->GetNbinsX());
		for(int iR=0;iR<MCAcceptance->GetNbinsX();iR++)
			for(int lat=1;lat<LATcorr->GetNbinsX();lat++){
				Geomag_Acceptance->SetBinContent(iR+1,lat+1,MCAcceptance->GetBinContent(iR+1)/LATcorr->GetBinContent(lat+1));
                                error=pow(MCAcceptance->GetBinContent(iR+1)*LATcorr->GetBinError(lat+1),2);
				error=pow(MCAcceptance->GetBinError(iR+1)*LATcorr->GetBinContent(lat+1),2);
                                error=pow(error,0.5);
				Geomag_Acceptance->SetBinError(iR+1,lat+1,error);
				}
	}

	return Geomag_Acceptance; 
}

void ACCEPTANCE::Eval_Geomag_Acceptance(int n){

	if(LATcorr_R	)   Geomag_Acceptance_R   = ACCEPTANCE::Geomag_Acceptance(n,MCAcceptance_R	, LATcorr_R	);
	if(LATcorr_TOF	)   Geomag_Acceptance_TOF = ACCEPTANCE::Geomag_Acceptance(n,MCAcceptance_TOF	, LATcorr_TOF	);
	if(LATcorr_NaF	)   Geomag_Acceptance_NaF = ACCEPTANCE::Geomag_Acceptance(n,MCAcceptance_NaF	, LATcorr_NaF	);
	if(LATcorr_Agl	)   Geomag_Acceptance_Agl = ACCEPTANCE::Geomag_Acceptance(n,MCAcceptance_Agl	, LATcorr_Agl	);

	return;
}



TH1 * ACCEPTANCE::Corrected_Acceptance(int n, TH1* MCAcceptance, TH1 *LATcorrW){
	TH1 * Corrected_Acceptance = (TH1 *) MCAcceptance -> Clone();
	float error=0;
	if(n>1){
		for(int m=0;m<n;m++)
			for(int iR=0;iR<MCAcceptance->GetNbinsX();iR++)
				if(MCAcceptance->GetBinContent(iR+1,m+1)>0&&LATcorrW->GetBinContent(iR)>0){
					Corrected_Acceptance -> SetBinContent(iR+1,m+1,MCAcceptance->GetBinContent(iR+1,m+1)/*LATcorrW->GetBinContent(iR+1)*/);
					error=pow(MCAcceptance->GetBinContent(iR+1,m+1)*LATcorrW->GetBinError(iR+1),2);
                                        error=pow(MCAcceptance->GetBinError(iR+1,m+1)*LATcorrW->GetBinContent(iR+1),2);
                                        error=pow(error,0.5);
                                        Corrected_Acceptance->SetBinError(iR+1,m+1,error);
					}
	}	

	else{
		for(int iR=0;iR<MCAcceptance->GetNbinsX();iR++){
			if(MCAcceptance->GetBinContent(iR+1)>0&&LATcorrW->GetBinContent(iR)>0){
				Corrected_Acceptance -> SetBinContent(iR+1,MCAcceptance->GetBinContent(iR+1)/*LATcorrW->GetBinContent(iR+1)*/);}
				error=pow(MCAcceptance->GetBinContent(iR+1)*LATcorrW->GetBinError(iR+1),2);
                                error=pow(MCAcceptance->GetBinError(iR+1)*LATcorrW->GetBinContent(iR+1),2);
                                error=pow(error,0.5);
                                Corrected_Acceptance->SetBinError(iR+1,error);
				
			}
	}
	return Corrected_Acceptance;
}


void ACCEPTANCE::Eval_Corrected_Acceptance(int n){
	if(LATcorrW_R   )  CorrectedAcceptance_R    = ACCEPTANCE::Corrected_Acceptance(n, MCAcceptance_R  , LATcorrW_R  );
	if(LATcorrW_TOF )  CorrectedAcceptance_TOF  = ACCEPTANCE::Corrected_Acceptance(n, MCAcceptance_TOF, LATcorrW_TOF);
	if(LATcorrW_NaF )  CorrectedAcceptance_NaF  = ACCEPTANCE::Corrected_Acceptance(n, MCAcceptance_NaF, LATcorrW_NaF);
	if(LATcorrW_Agl )  CorrectedAcceptance_Agl  = ACCEPTANCE::Corrected_Acceptance(n, MCAcceptance_Agl, LATcorrW_Agl);

	return;

}


void GlobalFactor(TH1* Acceptance,float factor, float error){
		Acceptance -> Scale(factor);
		float Error = 0;
		for(int m=0;m<Acceptance->GetNbinsY();m++)
                        for(int R=0;R<Acceptance->GetNbinsX();R++)
				for(int R=0;R<Acceptance->GetNbinsX();R++){
				Error = pow(error/factor,2);
				Error +=pow( Acceptance -> GetBinError(R+1,m+1) /  Acceptance -> GetBinContent(R+1,m+1),2);
				Error = pow(Error,0.5)*Acceptance -> GetBinContent(R+1,m+1);
				Acceptance -> SetBinError(R+1,m+1,Error);
			}
		return;
}


void ACCEPTANCE::ApplyGlobalFactor(float factor, float error){
	
	if(LATcorrW_R   )  GlobalFactor( CorrectedAcceptance_R  ,factor,error);
        if(LATcorrW_TOF )  GlobalFactor( CorrectedAcceptance_TOF,factor,error);
        if(LATcorrW_NaF )  GlobalFactor( CorrectedAcceptance_NaF,factor,error);
        if(LATcorrW_Agl )  GlobalFactor( CorrectedAcceptance_Agl,factor,error);

	return;
}




void ACCEPTANCE::ApplyDvsMCcorrection(int n, int S, TH1* Correction, TH1* Corrected_Acc, TH1* Geom_Acc){
	float error=0;
	//correct Geom. acceptance
	if(n>1){
		for(int m=0;m<n;m++)
			for(int R=0;R<Geom_Acc->GetNbinsX();R++)
				for(int lat=1;lat<Geom_Acc->GetNbinsY();lat++){
					for(int sel =0; sel < S;sel++){
						Geom_Acc->SetBinContent(R+1,lat+1,m+1,Geom_Acc->GetBinContent(R+1,lat+1,m+1)*Correction->GetBinContent(R+1,m+1,sel+1));
						error=pow(Geom_Acc->GetBinContent(R+1,lat+1,m+1)*Correction->GetBinError(R+1,m+1,sel+1),2);
						error+=pow(Geom_Acc->GetBinError(R+1,lat+1,m+1)*Correction->GetBinContent(R+1,m+1,sel+1),2);
						error=pow(error,0.5);
						Geom_Acc->SetBinError(R+1,lat+1,m+1,error);
					}
				}
	}
	else {
		for(int R=0;R<Geom_Acc->GetNbinsX();R++)
			for(int lat=1;lat<Geom_Acc->GetNbinsY();lat++){
				for(int sel =0; sel < S;sel++){
					Geom_Acc->SetBinContent(R+1,lat+1,Geom_Acc->GetBinContent(R+1,lat+1)*Correction->GetBinContent(R+1,sel+1));
					error=pow(Geom_Acc->GetBinContent(R+1,lat+1)*Correction->GetBinError(R+1,sel+1),2);
					error+=pow(Geom_Acc->GetBinError(R+1,lat+1)*Correction->GetBinContent(R+1,sel+1),2);
					error=pow(error,0.5);
					Geom_Acc->SetBinError(R+1,lat+1,error);
				}
			}

	}

	//correct total acceptance
	if(S == 1) {
		Corrected_Acc    -> Multiply( Correction   );	
	}

	else {
		if(n>1){
			for(int m=0;m<n;m++)
				for(int R=0;R<Corrected_Acc->GetNbinsX();R++)
					for(int sel =0; sel < S; sel++){
						Corrected_Acc -> SetBinContent(R+1,m+1, Corrected_Acc -> GetBinContent(R+1,m+1)*Correction -> GetBinContent(R+1,m+1,sel+1) );
						error=pow(Corrected_Acc -> GetBinContent(R+1,m+1)*Correction->GetBinError(R+1,m+1,sel+1),2);
						error+=pow(Corrected_Acc -> GetBinError(R+1,m+1)*Correction->GetBinContent(R+1,m+1,sel+1),2);
						error=pow(error,0.5);
						Corrected_Acc -> SetBinError(R+1,m+1,sel+1,error);

					}			
		}
		else{
			cout<<Corrected_Acc -> ClassName()<<endl;
			for(int R=0;R<Corrected_Acc->GetNbinsX();R++)
				for(int sel =0; sel < S; sel++){
					Corrected_Acc -> SetBinContent(R+1,Corrected_Acc -> GetBinContent(R+1)*Correction -> GetBinContent(R+1,sel+1) );
					error=pow(Corrected_Acc -> GetBinContent(R+1)*Correction->GetBinError(R+1,sel+1),2);
					error+=pow(Corrected_Acc -> GetBinError(R+1)*Correction->GetBinContent(R+1,sel+1),2);
					error=pow(error,0.5);
					Corrected_Acc -> SetBinError(R+1,error);
				}

		}
	}

}
#endif
