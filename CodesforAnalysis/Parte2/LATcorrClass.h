using namespace std;

class LATcorr
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
    
    //LAT corr
    TH1 * LATcorrR;
    TH1 * LATcorrTOF;
    TH1 * LATcorrNaF;
    TH1 * LATcorrAgl;

    //Fit
    TH1 * LATcorrR_fit;
    TH1 * LATcorrTOF_fit;
    TH1 * LATcorrNaF_fit;
    TH1 * LATcorrAgl_fit;

    //name
    std::string name;				 
    
    //  Creation constructors:

    LATcorr(std::string basename, int n){
        beforeTOF = new TH2F((basename + "1"   ).c_str(),(basename + "1"   ).c_str(),43,0,43, n, 0 ,n);
        afterTOF  = new TH2F((basename + "2"   ).c_str(),(basename + "2"   ).c_str(),43,0,43, n, 0 ,n);
        beforeNaF = new TH2F((basename + "1NaF").c_str(),(basename + "1NaF").c_str(),43,0,43, n, 0 ,n);
        afterNaF  = new TH2F((basename + "2NaF").c_str(),(basename + "2NaF").c_str(),43,0,43, n, 0 ,n);
        beforeAgl = new TH2F((basename + "1Agl").c_str(),(basename + "1Agl").c_str(),43,0,43, n, 0 ,n);
        afterAgl  = new TH2F((basename + "2Agl").c_str(),(basename + "2Agl").c_str(),43,0,43, n, 0 ,n);
        beforeR   = new TH2F((basename + "1_R" ).c_str(),(basename + "1_R" ).c_str(),43,0,43, n, 0 ,n);
        afterR    = new TH2F((basename + "2_R" ).c_str(),(basename + "2_R" ).c_str(),43,0,43, n, 0 ,n);
	name = basename; 
   }

  //   Reading constructors
    //standard	
    LATcorr(TFile * file, std::string basename){
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
    //read all results
    LATcorr(TFile * file, std::string basename,std::string dirname){
        beforeTOF = (TH1 *)file->Get((basename + "1"   ).c_str());
        afterTOF  = (TH1 *)file->Get((basename + "2"   ).c_str());
        beforeNaF = (TH1 *)file->Get((basename + "1NaF").c_str());
        afterNaF  = (TH1 *)file->Get((basename + "2NaF").c_str());
        beforeAgl = (TH1 *)file->Get((basename + "1Agl").c_str());
        afterAgl  = (TH1 *)file->Get((basename + "2Agl").c_str());
        beforeR   = (TH1 *)file->Get((basename + "1_R" ).c_str());
        afterR    = (TH1 *)file->Get((basename + "2_R" ).c_str());
	LATcorrR  	= (TH1 *)file->Get(("/" + dirname +"/" + basename + "_LATcorrR 	 	 "  ).c_str());
        LATcorrTOF 	= (TH1 *)file->Get(("/" + dirname +"/" +basename  + "_LATcorrTOF	 "  ).c_str());
        LATcorrNaF 	= (TH1 *)file->Get(("/" + dirname +"/" +basename  + "_LATcorrNaF	 "  ).c_str());
        LATcorrAgl 	= (TH1 *)file->Get(("/" + dirname +"/" +basename  + "_LATcorrAgl	 "  ).c_str());
	LATcorrR_fit    = (TH1 *)file->Get(("/" + dirname +"/" +basename  + "_LATcorrR_fit   "  ).c_str());
        LATcorrTOF_fit  = (TH1 *)file->Get(("/" + dirname +"/" +basename  + "_LATcorrTOF_fit "  ).c_str());
        LATcorrNaF_fit  = (TH1 *)file->Get(("/" + dirname +"/" +basename  + "_LATcorrNaF_fit "  ).c_str());
        LATcorrAgl_fit  = (TH1 *)file->Get(("/" + dirname +"/" +basename  + "_LATcorrAgl_fit "  ).c_str());
	name = basename;   
    }	

    void Write();
    void UpdateErrorbars();	
    void Eval_Efficiency();		
    void Eval_LATcorr();	
    void Fit_LATcorr();	
};


void LATcorr::Write()
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

void LATcorr::UpdateErrorbars()
{

   if(afterR)	 afterR->Sumw2();
   if(beforeR)	 beforeR->Sumw2(); 
   if(beforeTOF) beforeTOF->Sumw2(); 
   if(beforeNaF) beforeNaF->Sumw2();
   if(beforeAgl) beforeAgl->Sumw2();
   if(afterTOF)	 afterTOF->Sumw2();
   if(afterNaF)	 afterNaF->Sumw2();
   if(afterAgl)	 afterAgl->Sumw2();

}
 
void LATcorr::Eval_Efficiency(){
	
	if(afterR)	 effR     = (TH1 *)afterR	->Clone();
	if(afterTOF)	 effTOF   = (TH1 *)afterTOF	->Clone();	
	if(afterNaF)  	 effNaF   = (TH1 *)afterNaF	->Clone();
	if(afterAgl) 	 effAgl   = (TH1 *)afterAgl	->Clone();
	
	LATcorr::UpdateErrorbars();
	if(effR)   	effR   ->Divide( beforeR   );	
        if(effTOF) 	effTOF ->Divide( beforeTOF );
        if(effNaF) 	effNaF ->Divide( beforeNaF );
        if(effAgl) 	effAgl ->Divide( beforeAgl );

	if(effR)  	effR   ->SetName((name	+ "_EffR"  ).c_str()); 
        if(effTOF)	effTOF ->SetName((name	+ "_EffTOF").c_str());
        if(effNaF)	effNaF ->SetName((name	+ "_EffNaF").c_str());
        if(effAgl)	effAgl ->SetName((name	+ "_EffAgl").c_str());
}


void  CalcLATcorr(TH2 * before,TH2 * after, TH1 * LATcorr){
		float HEeff_before[11];
                float HEeff_after[11];

                for(int i=1;i<11;i++) {
                                HEeff_before[i] = before-> Integral(30,43,i+1,i+1);
                                HEeff_after[i]  = after -> Integral(30,43,i+1,i+1);
                        } 
		for(int i=1;i<11;i++){
			LATcorr -> SetBinContent(i+1,(HEeff_after[1]/HEeff_before[1])/(HEeff_after[i]/HEeff_before[i]));
			LATcorr -> SetBinError(i+1,pow(HEeff_after[i],-0.5)*LATcorr -> GetBinContent(i+1));
                }
}


void LATcorr::Eval_LATcorr(){
	if(effR){
                LATcorrR = new TH1F((name  + "_LATcorrR"  ).c_str(),(name  + "_LATcorrR"  ).c_str(),11,0,11);
		CalcLATcorr( ((TH2 *)beforeR),((TH2 *)afterR),((TH1 *)LATcorrR));
        }

	if(effTOF){
		LATcorrTOF = new TH1F((name  + "_LATcorrTOF"  ).c_str(),(name  + "_LATcorrTOF"  ).c_str(),11,0,11);
		CalcLATcorr( ((TH2 *)beforeTOF),((TH2 *)afterTOF),LATcorrTOF);
	}
	if(effNaF){
                LATcorrNaF = new TH1F((name  + "_LATcorrNaF"  ).c_str(),(name  + "_LATcorrNaF"  ).c_str(),11,0,11); 
		CalcLATcorr( ((TH2 *)beforeNaF),((TH2 *)afterNaF),LATcorrNaF);
	}
	if(effAgl){
                LATcorrAgl = new TH1F((name  + "_LATcorrAgl"  ).c_str(),(name  + "_LATcorrAgl"  ).c_str(),11,0,11); 
		CalcLATcorr( ((TH2 *)beforeAgl),((TH2 *)afterAgl),LATcorrAgl);
	}
	
}


void LATcorr::Fit_LATcorr(){
	if(LATcorrR){
		TF1 * Fitcorr = new TF1("Fitcorr","pol3");	
		LATcorrR -> Fit("Fitcorr");
		LATcorrR_fit = new TH1F((name  + "_LATcorrR_fit"  ).c_str(),(name  + "_LATcorrR_fit"  ).c_str(),11,0,11);
		for(int i=1;i<11;i++)  LATcorrR_fit -> SetBinContent(i+1,Fitcorr->Eval(i));
		for(int i=1;i<11;i++) LATcorrR_fit -> SetBinError(i+1,FitError(LATcorrR_fit,LATcorrR,11,3));		
	}

	if(LATcorrTOF){
                TF1 * Fitcorr = new TF1("Fitcorr","pol3");
                LATcorrTOF -> Fit("Fitcorr");
                LATcorrTOF_fit = new TH1F((name  + "_LATcorrTOF_fit"  ).c_str(),(name  + "_LATcorrTOF_fit"  ).c_str(),11,0,11);
                for(int i=1;i<11;i++)  LATcorrTOF_fit -> SetBinContent(i+1,Fitcorr->Eval(i));
		for(int i=1;i<11;i++) LATcorrTOF_fit -> SetBinError(i+1, FitError(LATcorrTOF_fit,LATcorrTOF,11,3));
        }
	if(LATcorrNaF){
                TF1 * Fitcorr = new TF1("Fitcorr","pol3");
                LATcorrNaF -> Fit("Fitcorr");
                LATcorrNaF_fit = new TH1F((name  + "_LATcorrNaF_fit"  ).c_str(),(name  + "_LATcorrNaF_fit"  ).c_str(),11,0,11);
                for(int i=1;i<11;i++)  LATcorrNaF_fit -> SetBinContent(i+1,Fitcorr->Eval(i));
		for(int i=1;i<11;i++) LATcorrNaF_fit -> SetBinError(i+1,FitError(LATcorrNaF_fit,LATcorrNaF,11,3));
        }
	if(LATcorrAgl){
                TF1 * Fitcorr = new TF1("Fitcorr","pol3");
                LATcorrAgl -> Fit("Fitcorr");
                LATcorrAgl_fit = new TH1F((name  + "_LATcorrAgl_fit"  ).c_str(),(name  + "_LATcorrAgl_fit"  ).c_str(),11,0,11);
                for(int i=1;i<11;i++)  LATcorrAgl_fit -> SetBinContent(i+1,Fitcorr->Eval(i));
		for(int i=1;i<11;i++) LATcorrAgl_fit -> SetBinError(i+1,FitError(LATcorrAgl_fit,LATcorrAgl,11,3));
        }

		
}
