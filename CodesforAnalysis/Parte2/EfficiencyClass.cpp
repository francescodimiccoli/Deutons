using namespace std;


class Efficiency
{
public:
    //counts
    TH1 * beforeR,   * afterR;
    TH1 * beforeTOF, * afterTOF;
    TH1 * beforeNaF, * afterNaF;
    TH1 * beforeAgl, * afterAgl;

    //eff	 	
    TH1 * effR;   
    TH1 * effTOF;
    TH1 * effNaF;
    TH1 * effAgl;
   
    //eff_fit	
    TH1 * effR_fit;
    TH1 * effTOF_fit; 	
    TH1 * effNaF_fit;
    TH1 * effAgl_fit;


    //stat. errors	 	
    TH1 * err_statR  ;   
    TH1 * err_statTOF;
    TH1 * err_statNaF;
    TH1 * err_statAgl;
  
    //syst errors	 	
    TH1 * err_systR;   
    TH1 * err_systTOF;
    TH1 * err_systNaF;
    TH1 * err_systAgl;
   
 
	
 
    //name
    std::string name;				 
    //  Creation constructors:
   
     Efficiency(std::string basename){
        beforeTOF = new TH1F((basename + "1"   ).c_str(),(basename + "1"   ).c_str(),nbinsToF,0,nbinsToF);
        afterTOF  = new TH1F((basename + "2"   ).c_str(),(basename + "2"   ).c_str(),nbinsToF,0,nbinsToF);
        beforeNaF = new TH1F((basename + "1NaF").c_str(),(basename + "1NaF").c_str(),nbinsNaF,0,nbinsNaF);
        afterNaF  = new TH1F((basename + "2NaF").c_str(),(basename + "2NaF").c_str(),nbinsNaF,0,nbinsNaF);
        beforeAgl = new TH1F((basename + "1Agl").c_str(),(basename + "1Agl").c_str(),nbinsAgl,0,nbinsAgl);
        afterAgl  = new TH1F((basename + "2Agl").c_str(),(basename + "2Agl").c_str(),nbinsAgl,0,nbinsAgl);
        beforeR   = new TH1F((basename + "1_R" ).c_str(),(basename + "1_R" ).c_str(),nbinsr,  0,nbinsr);
        afterR    = new TH1F((basename + "2_R" ).c_str(),(basename + "2_R" ).c_str(),nbinsr,  0,nbinsr);
        name = basename;
    }

    Efficiency(std::string basename, int n){
        beforeTOF = new TH2F((basename + "1"   ).c_str(),(basename + "1"   ).c_str(),nbinsToF,0,nbinsToF, n, 0 ,n);
        afterTOF  = new TH2F((basename + "2"   ).c_str(),(basename + "2"   ).c_str(),nbinsToF,0,nbinsToF, n, 0 ,n);
        beforeNaF = new TH2F((basename + "1NaF").c_str(),(basename + "1NaF").c_str(),nbinsNaF,0,nbinsNaF, n, 0 ,n);
        afterNaF  = new TH2F((basename + "2NaF").c_str(),(basename + "2NaF").c_str(),nbinsNaF,0,nbinsNaF, n, 0 ,n);
        beforeAgl = new TH2F((basename + "1Agl").c_str(),(basename + "1Agl").c_str(),nbinsAgl,0,nbinsAgl, n, 0 ,n);
        afterAgl  = new TH2F((basename + "2Agl").c_str(),(basename + "2Agl").c_str(),nbinsAgl,0,nbinsAgl, n, 0 ,n);
        beforeR   = new TH2F((basename + "1_R" ).c_str(),(basename + "1_R" ).c_str(),nbinsr,  0,nbinsr, n, 0 ,n);
        afterR    = new TH2F((basename + "2_R" ).c_str(),(basename + "2_R" ).c_str(),nbinsr,  0,nbinsr, n, 0 ,n);
   	name = basename; 
   }

   Efficiency(std::string basename, int n, int m){
        beforeTOF = new TH3F((basename + "1"   ).c_str(),(basename + "1"   ).c_str(),nbinsToF,0,nbinsToF, n, 0 ,n, m, 0,m);
        afterTOF  = new TH3F((basename + "2"   ).c_str(),(basename + "2"   ).c_str(),nbinsToF,0,nbinsToF, n, 0 ,n, m, 0,m);
        beforeNaF = new TH3F((basename + "1NaF").c_str(),(basename + "1NaF").c_str(),nbinsNaF,0,nbinsNaF, n, 0 ,n, m, 0,m);
        afterNaF  = new TH3F((basename + "2NaF").c_str(),(basename + "2NaF").c_str(),nbinsNaF,0,nbinsNaF, n, 0 ,n, m, 0,m);
        beforeAgl = new TH3F((basename + "1Agl").c_str(),(basename + "1Agl").c_str(),nbinsAgl,0,nbinsAgl, n, 0 ,n, m, 0,m);
        afterAgl  = new TH3F((basename + "2Agl").c_str(),(basename + "2Agl").c_str(),nbinsAgl,0,nbinsAgl, n, 0 ,n, m, 0,m);
        beforeR   = new TH3F((basename + "1_R" ).c_str(),(basename + "1_R" ).c_str(),nbinsr,  0,nbinsr, n, 0 ,n, m, 0,m);
        afterR    = new TH3F((basename + "2_R" ).c_str(),(basename + "2_R" ).c_str(),nbinsr,  0,nbinsr, n, 0 ,n, m, 0,m);
   	name = basename; 
   }

  //   Reading constructors

    Efficiency(TFile * file, std::string basename){
        if (!file->IsOpen()) file->Open("READ");
        beforeTOF = (TH1 *)file->Get((basename + "1"   ).c_str());
        afterTOF  = (TH1 *)file->Get((basename + "2"   ).c_str());
        beforeNaF = (TH1 *)file->Get((basename + "1NaF").c_str());
        afterNaF  = (TH1 *)file->Get((basename + "2NaF").c_str());
        beforeAgl = (TH1 *)file->Get((basename + "1Agl").c_str());
        afterAgl  = (TH1 *)file->Get((basename + "2Agl").c_str());
        beforeR   = (TH1 *)file->Get((basename + "1_R" ).c_str());
        afterR    = (TH1 *)file->Get((basename + "2_R" ).c_str());
    	name = basename;   
    }


   Efficiency(TFile * file, std::string basename,std::string dirname){
        if (!file->IsOpen()) file->Open("READ");
        beforeTOF = (TH1 *)file->Get((basename + "1"   ).c_str());
        afterTOF  = (TH1 *)file->Get((basename + "2"   ).c_str());
        beforeNaF = (TH1 *)file->Get((basename + "1NaF").c_str());
        afterNaF  = (TH1 *)file->Get((basename + "2NaF").c_str());
        beforeAgl = (TH1 *)file->Get((basename + "1Agl").c_str());
        afterAgl  = (TH1 *)file->Get((basename + "2Agl").c_str());
        beforeR   = (TH1 *)file->Get((basename + "1_R" ).c_str());
        afterR    = (TH1 *)file->Get((basename + "2_R" ).c_str());

	effR	  = (TH1 *)file->Get(("/" + dirname +"/" +basename + "_EffR"             ).c_str());
        effTOF    = (TH1 *)file->Get(("/" + dirname +"/" +basename + "_EffTOF"           ).c_str());
        effNaF    = (TH1 *)file->Get(("/" + dirname +"/" +basename + "_EffNaF"           ).c_str());
        effAgl    = (TH1 *)file->Get(("/" + dirname +"/" +basename + "_EffAgl"           ).c_str());

	effR_fit      = (TH1 *)file->Get(("/" + dirname +"/" +basename + "_EffFitR"             ).c_str());
        effTOF_fit    = (TH1 *)file->Get(("/" + dirname +"/" +basename + "_EffFitTOF"           ).c_str());
        effNaF_fit    = (TH1 *)file->Get(("/" + dirname +"/" +basename + "_EffFitNaF"           ).c_str());
        effAgl_fit    = (TH1 *)file->Get(("/" + dirname +"/" +basename + "_EffFitAgl"           ).c_str());

	err_statR     = (TH1 *)file->Get(("/" + dirname +"/" +basename + "_statR"             ).c_str());
        err_statTOF   = (TH1 *)file->Get(("/" + dirname +"/" +basename + "_statTOF"           ).c_str());
        err_statNaF   = (TH1 *)file->Get(("/" + dirname +"/" +basename + "_statNaF"           ).c_str());
        err_statAgl   = (TH1 *)file->Get(("/" + dirname +"/" +basename + "_statAgl"           ).c_str());

	err_systR     = (TH1 *)file->Get(("/" + dirname +"/" +basename + "_systR"             ).c_str());
        err_systTOF   = (TH1 *)file->Get(("/" + dirname +"/" +basename + "_systTOF"           ).c_str());
        err_systNaF   = (TH1 *)file->Get(("/" + dirname +"/" +basename + "_systNaF"           ).c_str());
        err_systAgl   = (TH1 *)file->Get(("/" + dirname +"/" +basename + "_systAgl"           ).c_str());

	
	
        name = basename;
    }
	

    //  Cloning constructor
    
    Efficiency( Efficiency * Clone, std::string basename){
		   
	beforeTOF = (TH1 *) Clone -> beforeTOF ->Clone(); 
        afterTOF  = (TH1 *) Clone -> afterTOF  ->Clone(); 
        beforeNaF = (TH1 *) Clone -> beforeNaF ->Clone(); 
        afterNaF  = (TH1 *) Clone -> afterNaF  ->Clone(); 
        beforeAgl = (TH1 *) Clone -> beforeAgl ->Clone(); 
        afterAgl  = (TH1 *) Clone -> afterAgl  ->Clone(); 
        beforeR   = (TH1 *) Clone -> beforeR   ->Clone(); 
        afterR    = (TH1 *) Clone -> afterR    ->Clone(); 
                                               
        effR	  = (TH1 *) Clone -> effR      ->Clone();		  
        effTOF    = (TH1 *) Clone -> effTOF    ->Clone(); 
        effNaF    = (TH1 *) Clone -> effNaF    ->Clone(); 
        effAgl    = (TH1 *) Clone -> effAgl    ->Clone(); 

	effR	  ->SetName((basename + "_EffR"             ).c_str());		  
        effTOF    ->SetName((basename + "_EffTOF"           ).c_str());	
        effNaF    ->SetName((basename + "_EffNaF"           ).c_str());
        effAgl    ->SetName((basename + "_EffAgl"           ).c_str());


	name = basename;
 
    }		    
   

    void Write();
    void UpdateErrorbars();	
    void Eval_Efficiency();	
    void Eval_FittedEfficiency();		

    void Compose_Efficiency(Efficiency * Factor);
   			
};


void Efficiency::Write()
{
	if(afterR)	afterR->Write(); 
	if(beforeR)	beforeR->Write();	   
	if(beforeTOF)	beforeTOF->Write(); 
	if(beforeNaF)	beforeNaF->Write();
	if(beforeAgl)	beforeAgl->Write();
	if(afterTOF)	afterTOF->Write();
	if(afterNaF)	afterNaF->Write();
	if(afterAgl)	afterAgl->Write();
}

void Efficiency::UpdateErrorbars()
{

   if(afterR)	 afterR   ->Sumw2();
   if(beforeR)	 beforeR  ->Sumw2(); 
   if(beforeTOF) beforeTOF->Sumw2(); 
   if(beforeNaF) beforeNaF->Sumw2();
   if(beforeAgl) beforeAgl->Sumw2();
   if(afterTOF)	 afterTOF->Sumw2();
   if(afterNaF)	 afterNaF->Sumw2();
   if(afterAgl)	 afterAgl->Sumw2();

}

 
void Efficiency::Eval_Efficiency(){
	
	if(afterR)	 effR     = (TH1 *)afterR	->Clone();
	if(afterTOF)	 effTOF   = (TH1 *)afterTOF	->Clone();	
	if(afterNaF)  	 effNaF   = (TH1 *)afterNaF	->Clone();
	if(afterAgl) 	 effAgl   = (TH1 *)afterAgl	->Clone();
	
	Efficiency::UpdateErrorbars();
	if(effR)   	effR   ->Divide( beforeR   );	
        if(effTOF) 	effTOF ->Divide( beforeTOF );
        if(effNaF) 	effNaF ->Divide( beforeNaF );
        if(effAgl) 	effAgl ->Divide( beforeAgl );

	if(effR)  	effR   ->SetName((name	+ "_EffR"  ).c_str()); 
	if(effTOF)	effTOF ->SetName((name	+ "_EffTOF").c_str());
	if(effNaF)	effNaF ->SetName((name	+ "_EffNaF").c_str());
	if(effAgl)	effAgl ->SetName((name	+ "_EffAgl").c_str());

	if(effR)    err_statR  =SetErrors(effR   ); 
	if(effTOF)  err_statTOF=SetErrors(effTOF );
	if(effNaF)  err_statNaF=SetErrors(effNaF );
	if(effAgl)  err_statAgl=SetErrors(effAgl );

	if(effR)    err_statR  ->SetName((name + "_statR"             ).c_str()); 
	if(effTOF)  err_statTOF->SetName((name + "_statTOF"           ).c_str());
	if(effNaF)  err_statNaF->SetName((name + "_statNaF"           ).c_str());
	if(effAgl)  err_statAgl->SetName((name + "_statAgl"           ).c_str());

	return;
}


void Efficiency::Compose_Efficiency(Efficiency * Factor){
	
	 effR  ->Multiply((TH1*)Factor->effR  ->Clone() );
         effTOF->Multiply((TH1*)Factor->effTOF->Clone() );
         effNaF->Multiply((TH1*)Factor->effNaF->Clone() );
         effAgl->Multiply((TH1*)Factor->effAgl->Clone() );

	if(effR)    err_statR  =SetErrors(effR   ); 
	if(effTOF)  err_statTOF=SetErrors(effTOF );
	if(effNaF)  err_statNaF=SetErrors(effNaF );
	if(effAgl)  err_statAgl=SetErrors(effAgl );

	if(effR)    err_statR  ->SetName((name + "_statR"             ).c_str()); 
	if(effTOF)  err_statTOF->SetName((name + "_statTOF"           ).c_str());
	if(effNaF)  err_statNaF->SetName((name + "_statNaF"           ).c_str());
	if(effAgl)  err_statAgl->SetName((name + "_statAgl"           ).c_str());

	return;
}

void Efficiency::Eval_FittedEfficiency(){

	cout<<effR<<" "<<effTOF<<" "<<effNaF<<" "<<effAgl<<endl;	
	
	if(effR) {
		FitFunction * FitR = new FitFunction( effR,1);	
		FitR->FitValues(); 
		effR_fit =FitR-> ReturnFittedValues(); 
	}

	if(effTOF) {
		FitFunction * FitTOF = new FitFunction( effTOF,1);	
		FitTOF->FitValues(); 
		effTOF_fit =FitTOF-> ReturnFittedValues(); 
	}
	if(effNaF) {
		FitFunction * FitNaF = new FitFunction( effNaF,1);	
		FitNaF->FitValues(); 
		effNaF_fit =FitNaF-> ReturnFittedValues(); 
	}
	
	if(effAgl) {
		FitFunction * FitAgl = new FitFunction( effAgl,1);	
		FitAgl->FitValues(); 
		effAgl_fit =FitAgl-> ReturnFittedValues(); 
	}

	if(effR_fit)	effR_fit   ->SetName((name	+ "_EffFitR"  ).c_str());
	if(effTOF_fit)	effTOF_fit ->SetName((name	+ "_EffFitTOF").c_str());
        if(effNaF_fit)	effNaF_fit ->SetName((name	+ "_EffFitNaF").c_str());
        if(effAgl_fit)	effAgl_fit ->SetName((name	+ "_EffFitAgl").c_str());

	if(effR)   err_systR  = SetErrors(effR_fit   );
        if(effTOF) err_systTOF= SetErrors(effTOF_fit );
        if(effNaF) err_systNaF= SetErrors(effNaF_fit );
        if(effAgl) err_systAgl= SetErrors(effAgl_fit );

	if(effR)    err_systR  ->SetName((name + "_systR"             ).c_str()); 
	if(effTOF)  err_systTOF->SetName((name + "_systTOF"           ).c_str());
	if(effNaF)  err_systNaF->SetName((name + "_systNaF"           ).c_str());
	if(effAgl)  err_systAgl->SetName((name + "_systAgl"           ).c_str());

	if(effR)   EvalError(effR_fit   ,err_statR   ,err_systR   );
        if(effTOF) EvalError(effTOF_fit ,err_statTOF ,err_systTOF );
        if(effNaF) EvalError(effNaF_fit ,err_statNaF ,err_systNaF );
        if(effAgl) EvalError(effAgl_fit ,err_statAgl ,err_systAgl );



}



