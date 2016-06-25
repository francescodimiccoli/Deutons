#include "PlottingFunctions/MCQualeff_Plot.h"

using namespace std;

// eff. likelihood sel.
Efficiency * EffLikMCP  = new Efficiency ("EffLikMCP");
Efficiency * EffLikMCD  = new Efficiency ("EffLikMCD",6);
// eff. distance sel.
Efficiency * EffDistMCP = new Efficiency ("EffDistMCP");
Efficiency * EffDistMCD = new Efficiency ("EffDistMCD",6);

void MCQualeff_Fill() {

	int Kbin;
	if(Massa_gen<1) {
		//R bins
		Kbin=PRB.GetRBin(Tup.R);

		EffLikMCP->beforeR->Fill(Kbin);
		EffDistMCP->beforeR->Fill(Kbin);
		if(Tup.Dist5D_P<6)          EffDistMCP->afterR->Fill(Kbin);
		if(Tup.Dist5D_P<6&&Likcut)	EffLikMCP->afterR->Fill(Kbin);


		//Beta bins
		Kbin=ToFPB.GetRBin(RUsed);

		EffLikMCP ->beforeTOF->Fill(Kbin);
		EffDistMCP->beforeTOF->Fill(Kbin);
		if(Distcut)             EffDistMCP ->afterTOF ->Fill(Kbin);
		if(Distcut&&Likcut)	EffLikMCP->afterTOF ->Fill(Kbin);


		if(cmask.isFromNaF()) {
			Kbin=NaFPB.GetRBin(RUsed);
			EffLikMCP  ->beforeNaF->Fill(Kbin);
			EffDistMCP ->beforeNaF->Fill(Kbin);
			if(Distcut)             EffDistMCP  ->afterNaF ->Fill(Kbin);
			if(Distcut&&Likcut)	EffLikMCP ->afterNaF ->Fill(Kbin);
		}

		if(cmask.isFromAgl()) {
			Kbin=AglPB.GetRBin(RUsed);
			EffLikMCP->beforeAgl->Fill(Kbin);
			EffDistMCP->beforeAgl->Fill(Kbin);
			if(Distcut) 		EffDistMCP->afterAgl->Fill(Kbin);
			if(Distcut&&Likcut)	EffLikMCP->afterAgl->Fill(Kbin);
		}


	}
	if(Massa_gen<2&&Massa_gen>1) {
		//R bins
		Kbin=PRB.GetRBin(Tup.R);

		FillBinMGen(EffLikMCD ->beforeR, Kbin);
		FillBinMGen(EffDistMCD->beforeR, Kbin);
		if(Distcut)             FillBinMGen(EffDistMCD ->afterR,  Kbin);
		if(Tup.Dist5D<6&&Likcut)	FillBinMGen(EffLikMCD->afterR,  Kbin);


		//Beta bins
		Kbin=ToFDB.GetRBin(RUsed);
		FillBinMGen(EffLikMCD ->beforeTOF, Kbin);
		FillBinMGen(EffDistMCD->beforeTOF, Kbin);
		if(Distcut)             FillBinMGen(EffDistMCD ->afterTOF , Kbin);
		if(Distcut&&Likcut)	FillBinMGen(EffLikMCD->afterTOF , Kbin);

		if(cmask.isFromNaF()) {
			Kbin=NaFDB.GetRBin(RUsed);
			FillBinMGen(EffLikMCD ->beforeNaF, Kbin);
			FillBinMGen(EffDistMCD->beforeNaF, Kbin);
			if(Distcut)             FillBinMGen(EffDistMCD ->afterNaF , Kbin);
			if(Distcut&&Likcut)	FillBinMGen(EffLikMCD->afterNaF , Kbin);
		}

		if(cmask.isFromAgl()) {
			Kbin=AglDB.GetRBin(RUsed);
			FillBinMGen(EffLikMCD ->beforeAgl, Kbin);
			FillBinMGen(EffDistMCD->beforeAgl, Kbin);
			if(Distcut)          FillBinMGen(EffDistMCD ->afterAgl , Kbin);
			if(Distcut&&Likcut)	FillBinMGen(EffLikMCD->afterAgl , Kbin);
		}
	}

return;
}




void MCQualeff_Write() {
   EffLikMCP  -> Write();
   EffLikMCD  -> Write();
   EffDistMCP -> Write();
   EffDistMCD -> Write();
   return;
}


void MCQualeff(string filename) {

   cout<<"******* MC QUALITY SEL. EFFICIENCIES ********"<<endl;
    
   cout<<"*** Reading  P1 file ****"<<endl;
        TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");

   // eff. likelihood sel.
   Efficiency * EffLikMCP  = new Efficiency (inputHistoFile,"EffLikMCP");
   Efficiency * EffLikMCD  = new Efficiency (inputHistoFile,"EffLikMCD");
   // eff. distance sel.
   Efficiency * EffDistMCP = new Efficiency (inputHistoFile,"EffDistMCP");
   Efficiency * EffDistMCD = new Efficiency (inputHistoFile,"EffDistMCD");


   cout<<"******* MC QUALITY SEL. EFFICIENCIES ********"<<endl;


   EffLikMCP ->Eval_Efficiency();
   EffLikMCD ->Eval_Efficiency();

   EffDistMCP->Eval_Efficiency();
   EffDistMCD->Eval_Efficiency();


   TH1F * EffMCLikP_TH1F 		=(TH1F *)EffLikMCP ->effR  ->Clone();
   TH2F * EffMCLikD_TH2F 		=(TH2F *)EffLikMCD ->effR  ->Clone();
   TH1F * EffMCLikP_Beta_TH1F 	=(TH1F *)EffLikMCP ->effTOF->Clone();
   TH2F * EffMCLikD_Beta_TH2F 	=(TH2F *)EffLikMCD ->effTOF->Clone();
   TH1F * EffMCLikP_BetaNaF_TH1F 	=(TH1F *)EffLikMCP ->effNaF->Clone();
   TH2F * EffMCLikD_BetaNaF_TH2F 	=(TH2F *)EffLikMCD ->effNaF->Clone();
   TH1F * EffMCLikP_BetaAgl_TH1F 	=(TH1F *)EffLikMCP ->effAgl->Clone();
   TH2F * EffMCLikD_BetaAgl_TH2F 	=(TH2F *)EffLikMCD ->effAgl->Clone();

   TH1F * EffMCDistP_TH1F 		=(TH1F *)EffDistMCP->effR  ->Clone();
   TH2F * EffMCDistD_TH2F 		=(TH2F *)EffDistMCD->effR  ->Clone();
   TH1F * EffMCDistP_Beta_TH1F 	=(TH1F *)EffDistMCP->effTOF->Clone();
   TH2F * EffMCDistD_Beta_TH2F 	=(TH2F *)EffDistMCD->effTOF->Clone();
   TH1F * EffMCDistP_BetaNaF_TH1F  =(TH1F *)EffDistMCP->effNaF->Clone();
   TH2F * EffMCDistD_BetaNaF_TH2F  =(TH2F *)EffDistMCD->effNaF->Clone();
   TH1F * EffMCDistP_BetaAgl_TH1F  =(TH1F *)EffDistMCP->effAgl->Clone();
   TH2F * EffMCDistD_BetaAgl_TH2F  =(TH2F *)EffDistMCD->effAgl->Clone();


   finalHistos.Add(EffMCLikP_TH1F 	    );
   finalHistos.Add(EffMCLikD_TH2F 	  );
   finalHistos.Add(EffMCLikP_Beta_TH1F    );
   finalHistos.Add(EffMCLikD_Beta_TH2F    );
   finalHistos.Add(EffMCLikP_BetaNaF_TH1F );
   finalHistos.Add(EffMCLikD_BetaNaF_TH2F );
   finalHistos.Add(EffMCLikP_BetaAgl_TH1F );
   finalHistos.Add(EffMCLikD_BetaAgl_TH2F );
                          
   finalHistos.Add(EffMCDistP_TH1F 	  );
   finalHistos.Add(EffMCDistD_TH2F 	  );
   finalHistos.Add(EffMCDistP_Beta_TH1F   );
   finalHistos.Add(EffMCDistD_Beta_TH2F   );
   finalHistos.Add(EffMCDistP_BetaNaF_TH1F);
   finalHistos.Add(EffMCDistD_BetaNaF_TH2F);
   finalHistos.Add(EffMCDistP_BetaAgl_TH1F);
   finalHistos.Add(EffMCDistD_BetaAgl_TH2F);
   finalHistos.writeObjsInFolder("Results");

   cout<<"*** Plotting ...  ****"<<endl; 

   MCQualeff_Plot(	

	EffMCLikP_TH1F 	  	,
        EffMCLikD_TH2F 	   	,
        EffMCLikP_Beta_TH1F	,    
        EffMCLikD_Beta_TH2F    ,  
        EffMCLikP_BetaNaF_TH1F , 
        EffMCLikD_BetaNaF_TH2F ,
        EffMCLikP_BetaAgl_TH1F ,
        EffMCLikD_BetaAgl_TH2F ,
               
        EffMCDistP_TH1F 	,  
        EffMCDistD_TH2F 	, 
        EffMCDistP_Beta_TH1F   ,
        EffMCDistD_Beta_TH2F   ,
        EffMCDistP_BetaNaF_TH1F,
        EffMCDistD_BetaNaF_TH2F,
        EffMCDistP_BetaAgl_TH1F,
        EffMCDistD_BetaAgl_TH2F
	);
	
   return;
}

