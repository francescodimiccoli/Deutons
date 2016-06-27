#include "PlottingFunctions/MCTrackeff_Plot.h"

Efficiency * EffTriggMCP = new Efficiency ("EffTriggMCP");
Efficiency * EffTriggMCD = new Efficiency ("EffTriggMCD",6);

Efficiency * EffTrackMCP = new Efficiency ("EffTrackMCP");
Efficiency * EffTrackMCD = new Efficiency ("EffTrackMCD",6);

Efficiency * EffTOFMCP 	 = new Efficiency ("EffTOFMCP");
Efficiency * EffTOFMCD   = new Efficiency ("EffTOFMCD",6);


void MCTrackeff_Fill ()
{
   int Kbin;

   if (Massa_gen<1&&Massa_gen>0.5) {
      //R bins
      Kbin=PRB.GetRBin (fabs (Tup.Momento_gen) );
      EffTriggMCP->beforeR->Fill (Kbin);
      if ( cmask.isMinimumBiasTrigger() ) {
         EffTriggMCP->afterR->Fill (Kbin);
         EffTOFMCP->beforeR->Fill (Kbin);
      }
      if ( cmask.isMinimumBiasTrigger() && cmask.isMinimumBiasToF3or4Layers() && Tup.Beta_pre>0) EffTOFMCP->afterR->Fill (Kbin);

      if (Tup.EdepTOFU < EdepTOFbeta->Eval (Tup.Beta_pre)+1
            && Tup.EdepTOFU > EdepTOFbeta->Eval (Tup.Beta_pre)-1
            && cmask.isMinimumBiasTrigger() && cmask.isMinimumBiasToF3or4Layers()
            && Tup.Beta_pre>0) {
         EffTrackMCP->beforeR->Fill (Kbin);
         if ( cmask.isMinimumBiasTrigger() && cmask.isMinimumBiasToF3or4Layers() && cmask.isMinimumBiasTracker())
            EffTrackMCP->afterR->Fill (Kbin );
      }

      //Beta bins
      Kbin=ToFPB.GetRBin (Tup.Momento_gen);
      EffTriggMCP->beforeTOF->Fill (Kbin);
      if ( cmask.isMinimumBiasTrigger() ) {
         EffTriggMCP->afterTOF->Fill (Kbin);
         EffTOFMCP->beforeTOF->Fill (Kbin);
      }
      if ( cmask.isMinimumBiasTrigger() && cmask.isMinimumBiasToF3or4Layers() &&Tup.Beta_pre>0)
         EffTOFMCP->afterTOF->Fill (Kbin);

      if (Tup.EdepTOFU < EdepTOFbeta->Eval (Tup.Beta_pre)+1
            && Tup.EdepTOFU > EdepTOFbeta->Eval (Tup.Beta_pre)-1
      &&  cmask.isMinimumBiasTrigger() && cmask.isMinimumBiasToF3or4Layers() && Tup.Beta_pre > 0) {
      Kbin=ToFPB.GetRBin (Tup.Momento_gen);
         EffTrackMCP->beforeTOF->Fill (Kbin);
         if ( cmask.isMinimumBiasTrigger() && cmask.isMinimumBiasToF3or4Layers() && cmask.isMinimumBiasTracker())
            EffTrackMCP->afterTOF->Fill (Kbin);
      }
   }

   if (Massa_gen>1&&Massa_gen<2) {
      //R bins
      Kbin=PRB.GetRBin (Tup.Momento_gen);
      FillBinMGen (EffTriggMCD->beforeR, Kbin);
      if ( cmask.isMinimumBiasTrigger() ) {
         FillBinMGen (EffTriggMCD->afterR , Kbin);
         FillBinMGen (EffTOFMCD  ->beforeR, Kbin);
      }
      if ( cmask.isMinimumBiasTrigger() && cmask.isMinimumBiasToF3or4Layers() &&Tup.Beta_pre>0) FillBinMGen (EffTOFMCD  ->afterR , Kbin);


      if (Tup.EdepTOFU < EdepTOFbeta->Eval (Tup.Beta_pre)+1
            && Tup.EdepTOFU > EdepTOFbeta->Eval (Tup.Beta_pre)-1
      &&  cmask.isMinimumBiasTrigger() && cmask.isMinimumBiasToF3or4Layers() && Tup.Beta_pre > 0) {
      Kbin=PRB.GetRBin (Tup.Momento_gen);
         FillBinMGen (EffTrackMCD->beforeR, Kbin);
         if ( cmask.isMinimumBiasTrigger() && cmask.isMinimumBiasToF3or4Layers() && cmask.isMinimumBiasTracker()&&Tup.R_pre>0) FillBinMGen (EffTrackMCD->afterR , Kbin);

      }
      //Beta bins
      Kbin=ToFDB.GetRBin (Tup.Momento_gen);
           FillBinMGen (EffTriggMCD->beforeTOF, Kbin);
      if ( cmask.isMinimumBiasTrigger() )   {
      FillBinMGen (EffTriggMCD->afterTOF , Kbin);
         FillBinMGen (EffTOFMCD  ->beforeTOF, Kbin);
      }
      if ( cmask.isMinimumBiasTrigger() && cmask.isMinimumBiasToF3or4Layers() &&Tup.Beta_pre>0) FillBinMGen (EffTOFMCD  ->afterTOF , Kbin);

      if (Tup.EdepTOFU < EdepTOFbeta->Eval (Tup.Beta_pre)+1 && Tup.EdepTOFU > EdepTOFbeta->Eval (Tup.Beta_pre)-1 && ( (int) Tup.Cutmask&3 ) == 3 && Tup.Beta_pre > 0) {
         FillBinMGen (EffTrackMCD->beforeTOF, Kbin);
            if ( cmask.isMinimumBiasTrigger() && cmask.isMinimumBiasToF3or4Layers() && cmask.isMinimumBiasTracker()&&Tup.R_pre!=0) FillBinMGen (EffTrackMCD->afterTOF , Kbin);
         }

   }

   return;
}



void MCTrackeff_Write()
{
   EffTriggMCP->Write();
   EffTriggMCD->Write();
   EffTrackMCP->Write();
   EffTrackMCD->Write();
   EffTOFMCP  ->Write();
   EffTOFMCD  ->Write();

   return;
}



void MCTrackeff (string filename)
{

  cout<<"**** MC BASIC SEL. EFFICIENCIES ****"<<endl;

  cout<<"*** Reading  P1 file ****"<<endl;
  TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");



   Efficiency * EffTriggMCP = new Efficiency (inputHistoFile,"EffTriggMCP");
   Efficiency * EffTriggMCD = new Efficiency (inputHistoFile,"EffTriggMCD");

   Efficiency * EffTrackMCP = new Efficiency (inputHistoFile,"EffTrackMCP");
   Efficiency * EffTrackMCD = new Efficiency (inputHistoFile,"EffTrackMCD");

   Efficiency * EffTOFMCP 	 = new Efficiency (inputHistoFile,"EffTOFMCP");
   Efficiency * EffTOFMCD   = new Efficiency (inputHistoFile,"EffTOFMCD");


   cout<<"**** MC BASIC SEL. EFFICIENCIES ****"<<endl;


   EffTriggMCP  -> Eval_Efficiency();
   EffTriggMCD  -> Eval_Efficiency();

   EffTrackMCP  -> Eval_Efficiency();
   EffTrackMCD  -> Eval_Efficiency();

   EffTOFMCP    -> Eval_Efficiency();
   EffTOFMCD    -> Eval_Efficiency();


   TH1F * EffTriggerMCP_R_TH1F 	= (TH1F *) EffTriggMCP  -> effR    ->	Clone();
   TH1F * EffTriggerMCP_TH1F       = (TH1F *) EffTriggMCP  -> effTOF  ->	Clone();
   TH2F * EffTriggerMCD_R_TH2F 	= (TH2F *) EffTriggMCD  -> effR    ->	Clone();
   TH2F * EffTriggerMCD_TH2F       = (TH2F *) EffTriggMCD  -> effTOF  ->	Clone();

   TH1F * EffTrackerMCP_R_TH1F 	= (TH1F *) EffTrackMCP  -> effR    ->	Clone();
   TH1F * EffTrackerMCP_TH1F       = (TH1F *) EffTrackMCP  -> effTOF  ->	Clone();
   TH2F * EffTrackerMCD_R_TH2F 	= (TH2F *) EffTrackMCD  -> effR    ->	Clone();
   TH2F * EffTrackerMCD_TH2F       = (TH2F *) EffTrackMCD  -> effTOF  ->	Clone();

   TH1F * EffTOF_MCP_R_TH1F        = (TH1F *) EffTOFMCP  	-> effR   ->	Clone();
   TH1F * EffTOF_MCP_TH1F	        = (TH1F *) EffTOFMCP  	-> effTOF ->	Clone();
   TH2F * EffTOF_MCD_R_TH2F        = (TH2F *) EffTOFMCD  	-> effR   ->	Clone();
   TH2F * EffTOF_MCD_TH2F	        = (TH2F *) EffTOFMCD  	-> effTOF ->	Clone();






    finalHistos.Add(EffTriggerMCP_R_TH1F 	);	
    finalHistos.Add(EffTriggerMCP_TH1F 		);
    finalHistos.Add(EffTriggerMCD_R_TH2F	);
    finalHistos.Add(EffTriggerMCD_TH2F		);
    finalHistos.Add(EffTrackerMCP_R_TH1F	);
    finalHistos.Add(EffTrackerMCP_TH1F		);	
    finalHistos.Add(EffTrackerMCD_R_TH2F	);
    finalHistos.Add(EffTrackerMCD_TH2F		);	
    finalHistos.Add(EffTOF_MCP_R_TH1F		);		
    finalHistos.Add(EffTOF_MCP_TH1F		);	
    finalHistos.Add(EffTOF_MCD_R_TH2F		);	
    finalHistos.Add(EffTOF_MCD_TH2F		);	
    finalHistos.writeObjsInFolder("Results");
   	


	cout<<"*** Plotting ...  ****"<<endl;

	MCTrackeff_Plot(EffTriggerMCP_R_TH1F,
                        EffTriggerMCP_TH1F 	,
                        EffTriggerMCD_R_TH2F,
                        EffTriggerMCD_TH2F,	
                        EffTrackerMCP_R_TH1F,
                        EffTrackerMCP_TH1F,	
                        EffTrackerMCD_R_TH2F,
                        EffTrackerMCD_TH2F,	
                        EffTOF_MCP_R_TH1F,	
                        EffTOF_MCP_TH1F	,
                        EffTOF_MCD_R_TH2F,	
                        EffTOF_MCD_TH2F	
	);


	return;
}


