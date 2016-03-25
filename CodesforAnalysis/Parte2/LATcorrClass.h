using namespace std;

void  CalcLATcorr(TH1 * before,TH1 * after, TH1 * LATcorr,int n);
void  FitLATcorr (TH1 * LATcorr,TH1 * LATcorr_fit,int n);

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
    LATcorr(std::string basename){
        beforeTOF  	= new TH2F((basename  + "1"   		  ).c_str(),(basename  + "1"   		  ).c_str(),43,0,43, 11,0,11);
        afterTOF   	= new TH2F((basename  + "2"   		  ).c_str(),(basename  + "2"   		  ).c_str(),43,0,43, 11,0,11);
        beforeNaF  	= new TH2F((basename  + "1NaF"		  ).c_str(),(basename  + "1NaF"		  ).c_str(),43,0,43, 11,0,11);
        afterNaF   	= new TH2F((basename  + "2NaF"		  ).c_str(),(basename  + "2NaF"		  ).c_str(),43,0,43, 11,0,11);
        beforeAgl  	= new TH2F((basename  + "1Agl"		  ).c_str(),(basename  + "1Agl"		  ).c_str(),43,0,43, 11,0,11);
        afterAgl   	= new TH2F((basename  + "2Agl"		  ).c_str(),(basename  + "2Agl"		  ).c_str(),43,0,43, 11,0,11);
        beforeR    	= new TH2F((basename  + "1_R" 		  ).c_str(),(basename  + "1_R" 		  ).c_str(),43,0,43, 11,0,11);
        afterR     	= new TH2F((basename  + "2_R" 		  ).c_str(),(basename  + "2_R" 		  ).c_str(),43,0,43, 11,0,11);
	LATcorrR   	= new TH1F((basename  + "_LATcorrR"       ).c_str(),(basename  + "_LATcorrR"      ).c_str(),11,0,11);
        LATcorrTOF 	= new TH1F((basename  + "_LATcorrTOF"     ).c_str(),(basename  + "_LATcorrTOF"    ).c_str(),11,0,11);
        LATcorrNaF 	= new TH1F((basename  + "_LATcorrNaF"     ).c_str(),(basename  + "_LATcorrNaF"    ).c_str(),11,0,11);
        LATcorrAgl 	= new TH1F((basename  + "_LATcorrAgl"     ).c_str(),(basename  + "_LATcorrAgl"    ).c_str(),11,0,11);
        LATcorrR_fit    = new TH1F((basename  + "_LATcorrR_fit"   ).c_str(),(basename  + "_LATcorrR_fit"  ).c_str(),11,0,11);
        LATcorrTOF_fit  = new TH1F((basename  + "_LATcorrTOF_fit" ).c_str(),(basename  + "_LATcorrTOF_fit").c_str(),11,0,11);
        LATcorrNaF_fit  = new TH1F((basename  + "_LATcorrNaF_fit" ).c_str(),(basename  + "_LATcorrNaF_fit").c_str(),11,0,11);
        LATcorrAgl_fit  = new TH1F((basename  + "_LATcorrAgl_fit" ).c_str(),(basename  + "_LATcorrAgl_fit").c_str(),11,0,11);
	name = basename; 
    }

    LATcorr(std::string basename, int n){
        beforeTOF 	= new TH3F((basename  + "1"   		  ).c_str(),(basename  + "1"   		  ).c_str(),43,0,43, 11,0,11, n, 0 ,n);
        afterTOF  	= new TH3F((basename  + "2"   		  ).c_str(),(basename  + "2"   		  ).c_str(),43,0,43, 11,0,11, n, 0 ,n);
        beforeNaF 	= new TH3F((basename  + "1NaF"		  ).c_str(),(basename  + "1NaF"		  ).c_str(),43,0,43, 11,0,11, n, 0 ,n);
        afterNaF  	= new TH3F((basename  + "2NaF"		  ).c_str(),(basename  + "2NaF"		  ).c_str(),43,0,43, 11,0,11, n, 0 ,n);
        beforeAgl 	= new TH3F((basename  + "1Agl"		  ).c_str(),(basename  + "1Agl"		  ).c_str(),43,0,43, 11,0,11, n, 0 ,n);
        afterAgl  	= new TH3F((basename  + "2Agl"		  ).c_str(),(basename  + "2Agl"		  ).c_str(),43,0,43, 11,0,11, n, 0 ,n);
        beforeR   	= new TH3F((basename  + "1_R" 		  ).c_str(),(basename  + "1_R" 		  ).c_str(),43,0,43, 11,0,11, n, 0 ,n);
        afterR    	= new TH3F((basename  + "2_R" 		  ).c_str(),(basename  + "2_R" 		  ).c_str(),43,0,43, 11,0,11, n, 0 ,n);
	LATcorrR   	= new TH2F((basename  + "_LATcorrR"    	  ).c_str(),(basename  + "_LATcorrR"      ).c_str(),11,0,11,n,0,n);
	LATcorrTOF 	= new TH2F((basename  + "_LATcorrTOF"  	  ).c_str(),(basename  + "_LATcorrTOF"    ).c_str(),11,0,11,n,0,n);
	LATcorrNaF 	= new TH2F((basename  + "_LATcorrNaF"  	  ).c_str(),(basename  + "_LATcorrNaF"    ).c_str(),11,0,11,n,0,n);
	LATcorrAgl 	= new TH2F((basename  + "_LATcorrAgl"  	  ).c_str(),(basename  + "_LATcorrAgl"    ).c_str(),11,0,11,n,0,n);
	LATcorrR_fit   	= new TH2F((basename  + "_LATcorrR_fit"   ).c_str(),(basename  + "_LATcorrR_fit"  ).c_str(),11,0,11,n,0,n);
	LATcorrTOF_fit	= new TH2F((basename  + "_LATcorrTOF_fit" ).c_str(),(basename  + "_LATcorrTOF_fit").c_str(),11,0,11,n,0,n);
	LATcorrNaF_fit  = new TH2F((basename  + "_LATcorrNaF_fit" ).c_str(),(basename  + "_LATcorrNaF_fit").c_str(),11,0,11,n,0,n);
	LATcorrAgl_fit 	= new TH2F((basename  + "_LATcorrAgl_fit" ).c_str(),(basename  + "_LATcorrAgl_fit").c_str(),11,0,11,n,0,n);
	
	name = basename; 
   }

  //   Reading constructors
    //standard	
    
   LATcorr(TFile * file, std::string basename){
        beforeTOF 	= (TH1 *)file->Get((basename + "1"   ).c_str());
        afterTOF  	= (TH1 *)file->Get((basename + "2"   ).c_str());
        beforeNaF 	= (TH1 *)file->Get((basename + "1NaF").c_str());
        afterNaF  	= (TH1 *)file->Get((basename + "2NaF").c_str());
        beforeAgl 	= (TH1 *)file->Get((basename + "1Agl").c_str());
        afterAgl  	= (TH1 *)file->Get((basename + "2Agl").c_str());
        beforeR   	= (TH1 *)file->Get((basename + "1_R" ).c_str());
        afterR    	= (TH1 *)file->Get((basename + "2_R" ).c_str());
	LATcorrR   	= new TH1F((basename  + "_LATcorrR"	    ).c_str(),(basename  + "_LATcorrR"    	).c_str(),11,0,11);
	LATcorrTOF 	= new TH1F((basename  + "_LATcorrTOF"	    ).c_str(),(basename  + "_LATcorrTOF"  	).c_str(),11,0,11);
	LATcorrNaF 	= new TH1F((basename  + "_LATcorrNaF"       ).c_str(),(basename  + "_LATcorrNaF"  	).c_str(),11,0,11);
	LATcorrAgl 	= new TH1F((basename  + "_LATcorrAgl"  	    ).c_str(),(basename  + "_LATcorrAgl"  	).c_str(),11,0,11);
	LATcorrR_fit    = new TH1F((basename  + "_LATcorrR_fit"     ).c_str(),(basename  + "_LATcorrR_fit"    	).c_str(),11,0,11);
	LATcorrTOF_fit	= new TH1F((basename  + "_LATcorrTOF_fit"   ).c_str(),(basename  + "_LATcorrTOF_fit"  	).c_str(),11,0,11);
	LATcorrNaF_fit  = new TH1F((basename  + "_LATcorrNaF_fit"   ).c_str(),(basename  + "_LATcorrNaF_fit"  	).c_str(),11,0,11);
	LATcorrAgl_fit	= new TH1F((basename  + "_LATcorrAgl_fit"   ).c_str(),(basename  + "_LATcorrAgl_fit"  	).c_str(),11,0,11);
	name = basename; 
  
    }

   LATcorr(TFile * file, std::string basename,int n){
         beforeTOF 	= (TH1 *)file->Get((basename + "1"   ).c_str());
        afterTOF  	= (TH1 *)file->Get((basename + "2"   ).c_str());
        beforeNaF 	= (TH1 *)file->Get((basename + "1NaF").c_str());
        afterNaF  	= (TH1 *)file->Get((basename + "2NaF").c_str());
        beforeAgl 	= (TH1 *)file->Get((basename + "1Agl").c_str());
        afterAgl  	= (TH1 *)file->Get((basename + "2Agl").c_str());
        beforeR   	= (TH1 *)file->Get((basename + "1_R" ).c_str());
        afterR    	= (TH1 *)file->Get((basename + "2_R" ).c_str());
	LATcorrR   	= new TH2F((basename  + "_LATcorrR"	    ).c_str(),(basename  + "_LATcorrR"    	).c_str(),11,0,11,n,0,n);
	LATcorrTOF 	= new TH2F((basename  + "_LATcorrTOF"	    ).c_str(),(basename  + "_LATcorrTOF"  	).c_str(),11,0,11,n,0,n);
	LATcorrNaF 	= new TH2F((basename  + "_LATcorrNaF"       ).c_str(),(basename  + "_LATcorrNaF"  	).c_str(),11,0,11,n,0,n);
	LATcorrAgl 	= new TH2F((basename  + "_LATcorrAgl"  	    ).c_str(),(basename  + "_LATcorrAgl"  	).c_str(),11,0,11,n,0,n);
	LATcorrR_fit    = new TH2F((basename  + "_LATcorrR_fit"     ).c_str(),(basename  + "_LATcorrR_fit"    	).c_str(),11,0,11,n,0,n);
	LATcorrTOF_fit	= new TH2F((basename  + "_LATcorrTOF_fit"   ).c_str(),(basename  + "_LATcorrTOF_fit"  	).c_str(),11,0,11,n,0,n);
	LATcorrNaF_fit  = new TH2F((basename  + "_LATcorrNaF_fit"   ).c_str(),(basename  + "_LATcorrNaF_fit"  	).c_str(),11,0,11,n,0,n);
	LATcorrAgl_fit	= new TH2F((basename  + "_LATcorrAgl_fit"   ).c_str(),(basename  + "_LATcorrAgl_fit"  	).c_str(),11,0,11,n,0,n);
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
	LATcorrR_fit    = (TH1 *)file->Get(("/" + dirname +"/" +basename  + "_LATcorrR_fit   	 "  ).c_str());
        LATcorrTOF_fit  = (TH1 *)file->Get(("/" + dirname +"/" +basename  + "_LATcorrTOF_fit 	 "  ).c_str());
        LATcorrNaF_fit  = (TH1 *)file->Get(("/" + dirname +"/" +basename  + "_LATcorrNaF_fit 	 "  ).c_str());
        LATcorrAgl_fit  = (TH1 *)file->Get(("/" + dirname +"/" +basename  + "_LATcorrAgl_fit  	 "  ).c_str());
	name = basename;   
    }	

    void Write();
    void UpdateErrorbars();	
    void Eval_Efficiency();		
    void Eval_LATcorr(int n);	
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

void LATcorr::Eval_LATcorr(int n){
	if(effR->GetEntries()>0){
		CalcLATcorr( ((TH2 *)beforeR),((TH2 *)afterR),LATcorrR,n);
		FitLATcorr(LATcorrR,LATcorrR_fit,n);
        }

	if(effTOF->GetEntries()>0){
		CalcLATcorr( ((TH2 *)beforeTOF),((TH2 *)afterTOF),((TH2 *)LATcorrTOF),n);
		FitLATcorr(LATcorrTOF,LATcorrTOF_fit,n);
	}
	if(effNaF->GetEntries()>0){
		CalcLATcorr( ((TH2 *)beforeNaF),((TH2 *)afterNaF),LATcorrNaF,n);
		FitLATcorr(LATcorrNaF,LATcorrNaF_fit,n);
	}
	if(effAgl->GetEntries()>0){
		CalcLATcorr( ((TH2 *)beforeAgl),((TH2 *)afterAgl),LATcorrAgl,n);
		FitLATcorr(LATcorrAgl,LATcorrAgl_fit,n);
	}
}




// Functions



void  CalcLATcorr(TH1 * before,TH1 * after, TH1 * LATcorr,int n){
	if(n>1){
		for(int m=0;m<n;m++){
			float HEeff_before[11];
			float HEeff_after[11];

			for(int i=1;i<11;i++) {
				HEeff_before[i] = ((TH3 *)before)-> Integral(30,43,i+1,i+1,m+1,m+1);
				HEeff_after[i]  = ((TH3 *)after )-> Integral(30,43,i+1,i+1,m+1,m+1);
			} 
			for(int i=1;i<11;i++){
				LATcorr -> SetBinContent(i+1,m+1,(HEeff_after[1]/HEeff_before[1])/(HEeff_after[i]/HEeff_before[i]));
				LATcorr -> SetBinError(i+1,m+1,pow(HEeff_after[i],-0.5)*LATcorr -> GetBinContent(i+1,m+1));
			}
		}
	}
	else{
		float HEeff_before[11];
		float HEeff_after[11];

		for(int i=1;i<11;i++) {
			HEeff_before[i] = ((TH2 *)before)-> Integral(30,43,i+1,i+1);
			HEeff_after[i]  = ((TH2 *)after )-> Integral(30,43,i+1,i+1);
		}
		for(int i=1;i<11;i++){
			  LATcorr -> SetBinContent(i+1,1,(HEeff_after[1]/HEeff_before[1])/(HEeff_after[i]/HEeff_before[i]));
			  LATcorr -> SetBinError(i+1,1,pow(HEeff_after[i],-0.5)*LATcorr -> GetBinContent(i+1));
		}
	}	
}


void FitLATcorr( TH1 * LATcorr,TH1 * LATcorr_fit,int n){
	if(n>1){
		for(int m=0;m<n;m++){
			TH1F * latcorr     = new TH1F("","",11,0,11);
			TH1F * latcorr_fit = new TH1F("","",11,0,11);

			for(int i=1;i<11;i++) {
				latcorr -> SetBinContent(i+1,((TH2 *)LATcorr)->GetBinContent(i+1,m+1));
				latcorr -> SetBinError  (i+1,((TH2 *)LATcorr)->GetBinError(i+1,m+1));
			}
			TF1 * Fitcorr = new TF1("Fitcorr","pol3");
			latcorr -> Fit ("Fitcorr");
			for(int i=1;i<11;i++)
				 latcorr_fit ->  SetBinContent(i+1,Fitcorr->Eval(i+0.5));
			
			for(int i=1;i<11;i++)  { 
				LATcorr_fit -> SetBinContent(i+1,m+1,Fitcorr->Eval(i+0.5));
				LATcorr_fit -> SetBinError(i+1,m+1,FitError(latcorr_fit,latcorr,11,3));			
			}
		}
	}
	else{
		TH1F * latcorr = ((TH1F *)LATcorr);
		TF1 * Fitcorr = new TF1("Fitcorr","pol3");
		latcorr -> Fit ("Fitcorr");
		for(int i=1;i<11;i++) {
			LATcorr_fit -> SetBinContent(i+1,Fitcorr->Eval(i+0.5));
                	LATcorr_fit -> SetBinError(i+1,FitError(((TH1F *)LATcorr_fit),latcorr,11,3));		
		}
	}
}
