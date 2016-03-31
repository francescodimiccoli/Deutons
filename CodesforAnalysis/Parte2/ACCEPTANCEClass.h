using namespace std;


class ACCEPTANCE
{
private:

	//binning
	float binsR[44];
	float binsBetaTOF[19];
	float binsBetaNaF[19];
	float binsBetaAgl[19];
public:
	//triggered/generated eff.
	float trigrate;
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
	
	TH1 * MCAcceptance_R   , * CorrectedAcceptance_R   ,  * Geom_Acceptance_R   ,  * Gen_Acceptance_R    ;
	TH1 * MCAcceptance_TOF , * CorrectedAcceptance_TOF ,  * Geom_Acceptance_TOF ,  * Gen_Acceptance_TOF  ;
	TH1 * MCAcceptance_NaF , * CorrectedAcceptance_NaF ,  * Geom_Acceptance_NaF ,  * Gen_Acceptance_NaF  ;
	TH1 * MCAcceptance_Agl , * CorrectedAcceptance_Agl ,  * Geom_Acceptance_Agl ,  * Gen_Acceptance_Agl  ;

	
	// reading constructor

	ACCEPTANCE(TFile * file , float bin[] , float binBetaTOF[], float binBetaNaF[], float binBetaAgl[],  float Trigrate, std::string dirname , std::string basename, std::string effname,std::string latcorrname , std::string wlatcorrname, int n)
		{
		after_TOF   = (TH1 *)file->Get((basename + "1"   ).c_str());
		after_NaF   = (TH1 *)file->Get((basename + "1NaF").c_str());
		after_Agl   = (TH1 *)file->Get((basename + "1Agl").c_str());
		after_R     = (TH1 *)file->Get((basename + "1_R" ).c_str());

		before_TOF  = new TH2F((basename + "1"   ).c_str(),(basename + "1"   ).c_str(),18,0,18, n, 0 ,n);
		before_NaF  = new TH2F((basename + "1NaF").c_str(),(basename + "1NaF").c_str(),18,0,18, n, 0 ,n);
		before_Agl  = new TH2F((basename + "1Agl").c_str(),(basename + "1Agl").c_str(),18,0,18, n, 0 ,n);
		before_R    = new TH2F((basename + "1_R" ).c_str(),(basename + "1_R" ).c_str(),43,0,43, n, 0 ,n);

		Efficiency_R  = (TH1 *)file->Get(("/" + dirname +"/" +effname + "TOF"         ).c_str());
                Efficiency_TOF= (TH1 *)file->Get(("/" + dirname +"/" +effname + "TOF"         ).c_str());
		Efficiency_NaF= (TH1 *)file->Get(("/" + dirname +"/" +effname + "NaF"         ).c_str());
                Efficiency_Agl= (TH1 *)file->Get(("/" + dirname +"/" +effname + "Agl"         ).c_str());

		LATcorr_R  = (TH1 *)file->Get(("/" + dirname +"/" +latcorrname + "TOF"         ).c_str());
		LATcorr_TOF= (TH1 *)file->Get(("/" + dirname +"/" +latcorrname + "TOF"         ).c_str());
		LATcorr_NaF= (TH1 *)file->Get(("/" + dirname +"/" +latcorrname + "NaF"         ).c_str());
		LATcorr_Agl= (TH1 *)file->Get(("/" + dirname +"/" +latcorrname + "Agl"         ).c_str());

		LATcorrW_R  = (TH1 *)file->Get(("/" + dirname +"/" +wlatcorrname + "_R"         ).c_str());
		LATcorrW_TOF= (TH1 *)file->Get(("/" + dirname +"/" +wlatcorrname + "_TOF"       ).c_str());
		LATcorrW_NaF= (TH1 *)file->Get(("/" + dirname +"/" +wlatcorrname + "_NaF"       ).c_str());
		LATcorrW_Agl= (TH1 *)file->Get(("/" + dirname +"/" +wlatcorrname + "_Agl"       ).c_str());
		
		for(int i=0;i<44; i++) {binsR[i] = bin[i];cout<<bin[i]<<endl;}
		for(int i=0;i<19; i++) {
		
			binsBetaTOF[i] = binBetaTOF[i]; 
                        binsBetaNaF[i] = binBetaNaF[i]; 
			binsBetaAgl[i] = binBetaAgl[i]; 	
			cout<<binsBetaTOF[i]<<" "<<binsBetaNaF[i]<<" "<<binsBetaAgl[i]<<endl;
		}
		trigrate = Trigrate;
	}


TH1 * Triggerbin(int n , TH1 * after, float trigrate, float bins[]);
void Eval_Gen_Acceptance( int n);
void Eval_MC_Acceptance();

};


TH1 * ACCEPTANCE::Triggerbin(int n , TH1 * after , float trigrate, float bins[]){
	//estimation of gen. particles in every bin
	float eventiprova[n];
	TH1 * triggerbin = (TH1 *) after -> Clone();
	float gen_bin=0;
	if(n>1)	
		for(int m=0;m<n;m++){
			eventiprova[m]   = ((TH2 *) after  ) -> Integral(0,(int)after   ->GetNbinsX(),m+1,m+1);
			for(int i=0; i<(int)after ->GetNbinsX();i++ ){
				gen_bin = eventiprova[m]*(pow(trigrate,-1))*(log(bins[i+1])-log(bins[i]))/(log(bins[(int)after ->GetNbinsX()])-log(bins[0]));
				triggerbin -> SetBinContent(i+1,m+1,gen_bin);
			}
		}

	else
		for(int m=0;m<n;m++){
			eventiprova[m]   = ((TH1 *) after  ) -> Integral(0,(int)after   ->GetNbinsX());
			for(int i=0; i<(int)after ->GetNbinsX();i++ ){
				gen_bin = eventiprova[m]*(pow(trigrate,-1))*(log(bins[i+1])-log(bins[i]))/(log(bins[(int)after ->GetNbinsX()])-log(bins[0]));
				triggerbin -> SetBinContent(i+1,gen_bin);
			}
		}
	return triggerbin;
}



void ACCEPTANCE::Eval_Gen_Acceptance(int n){
	before_R	= ACCEPTANCE::Triggerbin( n, after_R   ,trigrate , binsR);	
	before_TOF	= ACCEPTANCE::Triggerbin( n, after_TOF ,trigrate , binsBetaTOF);
	before_NaF	= ACCEPTANCE::Triggerbin( n, after_NaF ,trigrate , binsBetaNaF);	
	before_Agl	= ACCEPTANCE::Triggerbin( n, after_Agl ,trigrate , binsBetaAgl);

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

Gen_Acceptance_R  -> Multiply (	Efficiency_R  	);
Gen_Acceptance_TOF-> Multiply (	Efficiency_TOF	);
Gen_Acceptance_NaF-> Multiply (	Efficiency_NaF	);
Gen_Acceptance_Agl-> Multiply (	Efficiency_Agl	);

}


