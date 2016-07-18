#include "PlottingFunctions/MCUnbiaseff_Plot.h"

using namespace std;


Efficiency * EffUnbiasMCP;
Efficiency * EffUnbiasMCD;


void MCUnbiaseff_Fill() {

   if(!cmask.isPreselected()||Tup.Beta_pre<=0||Tup.R_pre<=0) return;
   //if(!(Tup.EdepTrack<EdepTrackbeta->Eval(Tup.Beta_pre)+0.2&&Tup.EdepTrack>EdepTrackbeta->Eval(Tup.Beta_pre)-0.2)) return;
   if(!Herejcut) return;
   int Kbin;
    	   
    if(Massa_gen<1&&Massa_gen>0.5) {
      //R bins
      Kbin=PRB.GetRBin(Tup.R_pre);
      if(trgpatt.IsPhysical()) EffUnbiasMCP->beforeR->Fill(Kbin,Tup.mcweight);
      if(trgpatt.IsUnbias()  ) EffUnbiasMCP->beforeR->Fill(Kbin,Tup.mcweight);
      if(trgpatt.IsPhysical()) EffUnbiasMCP->afterR->Fill(Kbin,Tup.mcweight);

      //Beta bins
      Kbin=ToFPB.GetBin(Tup.R_pre);
      if(trgpatt.IsPhysical()) EffUnbiasMCP->beforeTOF->Fill(Kbin,Tup.mcweight);
      if(trgpatt.IsUnbias()  ) EffUnbiasMCP->beforeTOF->Fill(Kbin,Tup.mcweight);
      if(trgpatt.IsPhysical()) EffUnbiasMCP->afterTOF ->Fill(Kbin,Tup.mcweight);  
   }

   if(Massa_gen>1&&Massa_gen<2) {
      //R bins
      Kbin=PRB.GetRBin(fabs(Tup.Momento_gen));
      FillBinMGen(EffUnbiasMCD->beforeR, Kbin);
      if(trgpatt.IsPhysical()) FillBinMGen(EffUnbiasMCD->afterR , Kbin);
      
      //Beta bins
      Kbin=ToFDB.GetBin(Tup.Momento_gen);
      FillBinMGen(EffUnbiasMCD->beforeTOF, Kbin);
      if(trgpatt.IsPhysical()) FillBinMGen(EffUnbiasMCD->afterTOF , Kbin);
      
   }

   return;
}


void MCUnbiaseff_Write() {
   EffUnbiasMCP->Write();
   EffUnbiasMCD->Write();
   return;
}



void MCUnbiaseff(string filename) {

	cout<<"**** MC Unbias TRIGGER EFF. ****"<<endl;

	cout<<"*** Reading  P1 file ****"<<endl;
	TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");

	EffUnbiasMCP = new Efficiency(inputHistoFile, "EffUnbiasMCP");
	EffUnbiasMCD = new Efficiency(inputHistoFile, "EffUnbiasMCD");

	cout<<"**** MC Unbias TRIGGER EFF. ****"<<endl;

	EffUnbiasMCP -> Eval_Efficiency();
	EffUnbiasMCD -> Eval_Efficiency();

	TH1F *EffUnbMCP_R_TH1F = (TH1F*)  EffUnbiasMCP->effR   ->Clone();
	TH1F *EffUnbMCP_TH1F   = (TH1F*)  EffUnbiasMCP->effTOF->Clone();
	TH2F *EffUnbMCD_R_TH2F = (TH2F*)  EffUnbiasMCD->effR   ->Clone();
	TH2F *EffUnbMCD_TH2F   = (TH2F*)  EffUnbiasMCD->effTOF->Clone();


	finalHistos.Add(EffUnbMCP_R_TH1F ); 
	finalHistos.Add(EffUnbMCP_TH1F   );
	finalHistos.Add(EffUnbMCD_R_TH2F );
	finalHistos.Add(EffUnbMCD_TH2F   );

	finalHistos.writeObjsInFolder("Results");

	cout<<"*** Plotting ...  ****"<<endl;

	MCUnbiaseff_Plot( EffUnbMCP_R_TH1F, 
			EffUnbMCP_TH1F   ,
			EffUnbMCD_R_TH2F ,
			EffUnbMCD_TH2F  ); 

	return;	
}

