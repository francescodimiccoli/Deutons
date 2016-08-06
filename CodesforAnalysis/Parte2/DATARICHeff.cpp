#include "PlottingFunctions/DATARICHeff_Plot.h"

using namespace std;

LATcorr * LATrichDATA_NaF   = new LATcorr("LATrichDATA_NaF");
LATcorr * LATrichDATA_Agl   = new LATcorr("LATrichDATA_Agl");


void DATARICHeff_Fill(int zona) {

	//cuts
	if(!((Tup.R>Rcut[zona]&&zona<10)||(zona==10)))  return;
	if(!(Distcut&&Likcut)) return;

	int Kbin=PRB.GetRBin(Tup.R);

	LATrichDATA_NaF -> beforeR -> Fill(Kbin,zona);
	LATrichDATA_Agl	-> beforeR -> Fill(Kbin,zona);
	
	if (cmask.isFromNaF()) LATrichDATA_NaF -> afterR -> Fill(Kbin,zona); 
	if (cmask.isFromAgl()) LATrichDATA_Agl -> afterR -> Fill(Kbin,zona); 

   	return;
}

void DATARICHeff_Write() {
	LATrichDATA_NaF -> Write();
	LATrichDATA_Agl -> Write();  
	return;
}




void DATARICHeff(string filename) {

   cout<<"****************************** DATA RICH SEL. EFFICIENCIES **************************************"<<endl;

	cout<<"*** Reading  P1 file ****"<<endl;
        TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");

   LATcorr * LATrichDATA_NaF = new LATcorr(inputHistoFile,"LATrichDATA_NaF");
   LATcorr * LATrichDATA_Agl = new LATcorr(inputHistoFile,"LATrichDATA_Agl");


   cout<<"****************************** DATA RICH SEL. EFFICIENCIES **************************************"<<endl;


   LATrichDATA_NaF -> Eval_Efficiency();
   LATrichDATA_Agl -> Eval_Efficiency();


   TH2F *LATrichDATANaF = (TH2F *)  LATrichDATA_NaF  -> effR -> Clone();
   TH2F *LATrichDATAAgl = (TH2F *)  LATrichDATA_Agl  -> effR -> Clone();


   cout<<"****************************** LAT. Eff. CORRECTION *************************************************"<<endl;

   LATrichDATA_NaF ->  Eval_LATcorr(1);
   LATrichDATA_Agl ->  Eval_LATcorr(1);

   TH2F *LATrichcorr_NaF  =	(TH2F *) LATrichDATA_NaF  -> LATcorrR -> Clone();
   TH2F *LATrichcorr_Agl  =	(TH2F *) LATrichDATA_Agl  -> LATcorrR -> Clone();

   TH1F *LATrichcorr_NaF_fit  	= (TH1F *) LATrichDATA_NaF   -> LATcorrR_fit-> Clone();
   TH1F *LATrichcorr_Agl_fit	= (TH1F *) LATrichDATA_Agl   -> LATcorrR_fit-> Clone();

	
	finalHistos.Add(LATrichDATANaF 		);	
        finalHistos.Add(LATrichDATAAgl 		);
                             
        finalHistos.Add(LATrichcorr_NaF		);
        finalHistos.Add(LATrichcorr_Agl		);
             
        finalHistos.Add(LATrichcorr_NaF_fit	);
        finalHistos.Add(LATrichcorr_Agl_fit	);

	finalHistos.writeObjsInFolder("Results");

        cout<<"*** Plotting ...  ****"<<endl;
	
	DATARICHeff_Plot(LATrichDATANaF 	 ,  
                         LATrichDATAAgl 	 ,  
                              
                         LATrichcorr_NaF	 ,  
                         LATrichcorr_Agl	 ,  
                                            
                         LATrichcorr_NaF_fit	,
                         LATrichcorr_Agl_fit
	);


	return;
}

