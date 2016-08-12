using namespace std;

void  CalcLATcorr(TH1 * before,TH1 * after, TH1 * LATcorr,int n);
void  FitLATcorr (TH1 * LATcorr,TH1 * LATcorr_fit,int n);

class LATcorr
{
public:
    //counts
    TH1 * beforeR,   * afterR;

    //eff	 	
    TH1 * effR;   
    
    //LAT corr
    TH1 * LATcorrR;

    //Fit
    TH1 * LATcorrR_fit;

    //name
    std::string name;				 
    
    //  Creation constructors:
    LATcorr(std::string basename){
        beforeR    	= new TH2F((basename  + "1" 		  ).c_str(),(basename  + "1" 		  ).c_str(),nbinsr,0,nbinsr, 11,0,11);
        afterR     	= new TH2F((basename  + "2" 		  ).c_str(),(basename  + "2" 		  ).c_str(),nbinsr,0,nbinsr, 11,0,11);
	LATcorrR   	= new TH1F((basename  + "_LATcorr"       ).c_str(),(basename  + "_LATcorr"      ).c_str(),11,0,11);
        LATcorrR_fit    = new TH1F((basename  + "_LATcorr_fit"   ).c_str(),(basename  + "_LATcorr_fit"  ).c_str(),11,0,11);
	name = basename; 
    }

    LATcorr(std::string basename, int n){
        beforeR   	= new TH3F((basename  + "1" 		  ).c_str(),(basename  + "1" 		  ).c_str(),nbinsr,0,nbinsr, 11,0,11, n, 0 ,n);
        afterR    	= new TH3F((basename  + "2" 		  ).c_str(),(basename  + "2" 		  ).c_str(),nbinsr,0,nbinsr, 11,0,11, n, 0 ,n);
	LATcorrR   	= new TH2F((basename  + "_LATcorr"    	  ).c_str(),(basename  + "_LATcorr"      ).c_str(),11,0,11,n,0,n);
	LATcorrR_fit   	= new TH2F((basename  + "_LATcorr_fit"   ).c_str(),(basename  + "_LATcorr_fit"  ).c_str(),11,0,11,n,0,n);
	
	name = basename; 
   }

  //   Reading constructors
    //standard	
    
   LATcorr(TFile * file, std::string basename){
        beforeR   	= (TH1 *)file->Get((basename + "1" ).c_str());
        afterR    	= (TH1 *)file->Get((basename + "2" ).c_str());
	LATcorrR   	= new TH1F((basename  + "_LATcorr"	    ).c_str(),(basename  + "_LATcorr"    	).c_str(),11,0,11);
	LATcorrR_fit    = new TH1F((basename  + "_LATcorr_fit"     ).c_str(),(basename  + "_LATcorr_fit"    	).c_str(),11,0,11);
	name = basename; 
  
    }

   LATcorr(TFile * file, std::string basename,int n){
        beforeR   	= (TH1 *)file->Get((basename + "1" ).c_str());
        afterR    	= (TH1 *)file->Get((basename + "2" ).c_str());
	LATcorrR   	= new TH2F((basename  + "_LATcorr"	    ).c_str(),(basename  + "_LATcorr"    	).c_str(),11,0,11,n,0,n);
	LATcorrR_fit    = new TH2F((basename  + "_LATcorr_fit"     ).c_str(),(basename  + "_LATcorr_fit"    	).c_str(),11,0,11,n,0,n);
	name = basename; 
  
    }
    //read all results
    LATcorr(TFile * file, std::string basename,std::string dirname){
        beforeR   = (TH1 *)file->Get((basename + "1" ).c_str());
        afterR    = (TH1 *)file->Get((basename + "2" ).c_str());
	LATcorrR  	= (TH1 *)file->Get(("/" + dirname +"/" + basename + "_LATcorr"  			).c_str());
	LATcorrR_fit    = (TH1 *)file->Get(("/" + dirname +"/" +basename  + "_LATcorr_fit"  			).c_str());
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
}

void LATcorr::UpdateErrorbars()
{

   if(afterR)	 afterR->Sumw2();
   if(beforeR)	 beforeR->Sumw2(); 

}
 
void LATcorr::Eval_Efficiency(){
	
	if(afterR)	 effR     = (TH1 *)afterR	->Clone();
	LATcorr::UpdateErrorbars();
	if(effR)   	effR   ->Divide( beforeR   );	
	if(effR)  	effR   ->SetName((name	+ "_Eff"  ).c_str()); 
}

void LATcorr::Eval_LATcorr(int n){
	if(effR->GetEntries()>0){
		CalcLATcorr( ((TH2 *)beforeR),((TH2 *)afterR),LATcorrR,n);
		FitLATcorr(LATcorrR,LATcorrR_fit,n);
        }
}




// Functions



void  CalcLATcorr(TH1 * before,TH1 * after, TH1 * LATcorr,int n){
	if(n>1){
		for(int m=0;m<n;m++){
			float HEeff_before[11];
			float HEeff_after[11];

			for(int i=1;i<11;i++) {
				HEeff_before[i] = ((TH3 *)before)-> Integral(30,nbinsr,i+1,i+1,m+1,m+1);
				HEeff_after[i]  = ((TH3 *)after )-> Integral(30,nbinsr,i+1,i+1,m+1,m+1);
			} 
			for(int i=1;i<11;i++){
			    if(HEeff_before[i]>0){	
				LATcorr -> SetBinContent(i+1,m+1,(HEeff_after[1]/HEeff_before[1])/(HEeff_after[i]/HEeff_before[i]));
				LATcorr -> SetBinError(i+1,m+1,pow(HEeff_after[i],-0.5)*LATcorr -> GetBinContent(i+1,m+1));
				}
			    else{
				LATcorr -> SetBinContent(i+1,m+1,1);
                                LATcorr -> SetBinError(i+1,m+1,0.5);
    
				}	
			   }
		}
	}
	else{
		float HEeff_before[11];
		float HEeff_after[11];

		for(int i=1;i<11;i++) {
			HEeff_before[i] = ((TH2 *)before)-> Integral(30,nbinsr,i+1,i+1);
			HEeff_after[i]  = ((TH2 *)after )-> Integral(30,nbinsr,i+1,i+1);
		}
		for(int i=1;i<11;i++){
		       if(HEeff_before[i]>0){	  
			  LATcorr -> SetBinContent(i+1,1,(HEeff_after[1]/HEeff_before[1])/(HEeff_after[i]/HEeff_before[i]));
			  LATcorr -> SetBinError(i+1,1,pow(HEeff_after[i],-0.5)*LATcorr -> GetBinContent(i+1));
			}
			else{
			  LATcorr -> SetBinContent(i+1,1,1);
                          LATcorr -> SetBinError(i+1,1,0.5);	
			}
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
			string fitname = LATcorr->GetName() ;
			TF1 * Fitcorr = new TF1(fitname.c_str(),"pol4");
			latcorr -> Fit (fitname.c_str());
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
		string fitname = LATcorr->GetName() ;
		TF1 * Fitcorr = new TF1(fitname.c_str(),"pol4");
		latcorr -> Fit (fitname.c_str());
		for(int i=1;i<11;i++) {
			cout<<Fitcorr->Eval(i+0.5)<<endl;
			LATcorr_fit -> SetBinContent(i+1,Fitcorr->Eval(i+0.5));}
                for(int i=1;i<11;i++)	
			LATcorr_fit -> SetBinError(i+1,FitError(((TH1F *)LATcorr_fit),latcorr,11,3));
		
	}
}
