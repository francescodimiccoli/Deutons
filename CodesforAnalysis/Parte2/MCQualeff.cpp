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
	
	float mass=0;	
	if(Massa_gen<1) {
		//R bins
		Kbin=PRB.GetRBin(Tup.R);

		EffLikMCP->beforeR->Fill(Kbin,Tup.mcweight);
		EffDistMCP->beforeR->Fill(Kbin,Tup.mcweight);
		if(Likcut)          EffLikMCP->afterR->Fill(Kbin,Tup.mcweight);
		if(Distcut&&Likcut)	EffDistMCP->afterR->Fill(Kbin,Tup.mcweight);


		//Beta bins
		Kbin=ToFPB.GetBin(RUsed);

		EffLikMCP ->beforeTOF->Fill(Kbin,Tup.mcweight);
		EffDistMCP->beforeTOF->Fill(Kbin,Tup.mcweight);
		if(Likcut)             EffLikMCP ->afterTOF ->Fill(Kbin,Tup.mcweight);
		if(Distcut&&Likcut)	EffDistMCP->afterTOF ->Fill(Kbin,Tup.mcweight);


		if(cmask.isFromNaF()) {
			Kbin=NaFPB.GetBin(RUsed);
			EffLikMCP  ->beforeNaF->Fill(Kbin,Tup.mcweight);
			EffDistMCP ->beforeNaF->Fill(Kbin,Tup.mcweight);
			if(Likcut)             EffLikMCP  ->afterNaF ->Fill(Kbin,Tup.mcweight);
			if(Distcut&&Likcut)	EffDistMCP ->afterNaF ->Fill(Kbin,Tup.mcweight);
		}

		if(cmask.isFromAgl()) {
			Kbin=AglPB.GetBin(RUsed);
			EffLikMCP->beforeAgl->Fill(Kbin,Tup.mcweight);
			EffDistMCP->beforeAgl->Fill(Kbin,Tup.mcweight);
			if(Likcut) 		EffLikMCP->afterAgl->Fill(Kbin,Tup.mcweight);
			if(Distcut&&Likcut)	EffDistMCP->afterAgl->Fill(Kbin,Tup.mcweight);
		}


	}
	if(Massa_gen<2&&Massa_gen>1) {
		//R bins
		Kbin=PRB.GetRBin(Tup.R);

		FillBinMGen(EffLikMCD ->beforeR, Kbin);
		FillBinMGen(EffDistMCD->beforeR, Kbin);
		if(Likcut)             FillBinMGen(EffLikMCD ->afterR,  Kbin);
		if(Distcut&&Likcut)	FillBinMGen(EffDistMCD->afterR,  Kbin);


		//Beta bins
		Kbin=ToFDB.GetBin(RUsed);
		mass = ((Tup.R/Tup.Beta)*pow((1-pow(Tup.Beta,2)),0.5));
		FillBinMGen(EffLikMCD ->beforeTOF, Kbin);
		FillBinMGen(EffDistMCD->beforeTOF, Kbin);
		if(Likcut)             FillBinMGen(EffLikMCD ->afterTOF , Kbin);
		if(Distcut&&Likcut&&mass>=0 && mass<3)	FillBinMGen(EffDistMCD->afterTOF , Kbin);

		if(cmask.isFromNaF()) {
			Kbin=NaFDB.GetBin(RUsed);
			FillBinMGen(EffLikMCD ->beforeNaF, Kbin);
			FillBinMGen(EffDistMCD->beforeNaF, Kbin);
			mass = ((Tup.R/Tup.BetaRICH)*pow((1-pow(Tup.BetaRICH,2)),0.5));
			if(Likcut)             FillBinMGen(EffLikMCD ->afterNaF , Kbin);
			if(Distcut&&Likcut&&mass>=0&&mass<3)	FillBinMGen(EffDistMCD->afterNaF , Kbin);
		}

		if(cmask.isFromAgl()) {
			Kbin=AglDB.GetBin(RUsed);
			FillBinMGen(EffLikMCD ->beforeAgl, Kbin);
			FillBinMGen(EffDistMCD->beforeAgl, Kbin);
			mass = ((Tup.R/Tup.BetaRICH)*pow((1-pow(Tup.BetaRICH,2)),0.5));
			if(Likcut)          FillBinMGen(EffLikMCD ->afterAgl , Kbin);
			if(Distcut&&Likcut && mass>=0 &&mass <3 )FillBinMGen(EffDistMCD->afterAgl , Kbin);
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
   TH1F * EffMCLikP_Beta_TH1F 		=(TH1F *)EffLikMCP ->effTOF->Clone();
   TH2F * EffMCLikD_Beta_TH2F 		=(TH2F *)EffLikMCD ->effTOF->Clone();
   TH1F * EffMCLikP_BetaNaF_TH1F 	=(TH1F *)EffLikMCP ->effNaF->Clone();
   TH2F * EffMCLikD_BetaNaF_TH2F 	=(TH2F *)EffLikMCD ->effNaF->Clone();
   TH1F * EffMCLikP_BetaAgl_TH1F 	=(TH1F *)EffLikMCP ->effAgl->Clone();
   TH2F * EffMCLikD_BetaAgl_TH2F 	=(TH2F *)EffLikMCD ->effAgl->Clone();

   TH1F * EffMCDistP_TH1F 		=(TH1F *)EffDistMCP->effR  ->Clone();
   TH2F * EffMCDistD_TH2F 		=(TH2F *)EffDistMCD->effR  ->Clone();
   TH1F * EffMCDistP_Beta_TH1F 	        =(TH1F *)EffDistMCP->effTOF->Clone();
   TH2F * EffMCDistD_Beta_TH2F 	        =(TH2F *)EffDistMCD->effTOF->Clone();
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

