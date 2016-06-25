#include "PlottingFunctions/MCpreeff_Plot.h"

using namespace std;



Efficiency * EffpreselMCP = new Efficiency("EffpreselMCP");
Efficiency * EffpreselMCD = new Efficiency("EffpreselMCD", 6);


void MCpreseff_Fill() {


   if(Massa_gen<1&&Massa_gen>0.5) {
      //R bins
      EffpreselMCP->beforeR->Fill(PRB.GetRBin(fabs(Tup.Momento_gen)));
      if(Tup.Unbias==0&&cmask.isPreselected()&&Tup.Beta_pre>0&&Tup.R_pre>0) EffpreselMCP->afterR->Fill(PRB.GetRBin(fabs(Tup.R_pre)));
         
      // Beta bins
 ToFPB.RFill(EffpreselMCP->beforeTOF, Tup.Momento_gen);
 NaFPB.RFill(EffpreselMCP->beforeNaF, Tup.Momento_gen);
 AglPB.RFill(EffpreselMCP->beforeAgl, Tup.Momento_gen);

      if(Tup.Unbias==0 && cmask.isPreselected() && Tup.Beta_pre>0 && Tup.R_pre>0)
      {
                                     ToFPB.RFill(EffpreselMCP->afterTOF, RUsed);
         if(cmask.isFromNaF()) NaFPB.RFill(EffpreselMCP->afterNaF, RUsed);
         if(cmask.isFromAgl()  ) AglPB.RFill(EffpreselMCP->afterAgl, RUsed);
      }

   }

   if(Massa_gen>1&&Massa_gen<2) {
      // R bins      
      FillBinMGen(EffpreselMCD->beforeR, PRB.GetRBin(fabs(Tup.Momento_gen)));
      if(cmask.isPreselected()&&Tup.Beta_pre>0&&Tup.Unbias==0&&Tup.R_pre>0)
         FillBinMGen(EffpreselMCD->afterR, PRB.GetRBin(fabs(Tup.R_pre)));

      // Beta bins

         FillBinMGen(EffpreselMCD->beforeTOF, ToFDB.GetRBin(Tup.Momento_gen) );
         FillBinMGen(EffpreselMCD->beforeNaF, NaFDB.GetRBin(Tup.Momento_gen) );
         FillBinMGen(EffpreselMCD->beforeAgl, AglDB.GetRBin(Tup.Momento_gen) );

         if(cmask.isPreselected() && Tup.Beta_pre>0 && Tup.Unbias==0 && Tup.R_pre>0)
         {
                                           FillBinMGen(EffpreselMCD->afterTOF, ToFDB.GetRBin(RUsed));
            if(cmask.isFromNaF()) FillBinMGen(EffpreselMCD->afterNaF, NaFDB.GetRBin(RUsed));
            if(cmask.isFromAgl()) FillBinMGen(EffpreselMCD->afterAgl, AglDB.GetRBin(RUsed));
         }
      
   }
   return;
}


void MCpreeff_Write() {
   EffpreselMCP->Write();
   EffpreselMCD->Write();
   return;
}



void MCpreeff(string filename) {
  
   cout<<"**** MC PRESELECTIONS EFFICIENCY (FULL SET) ****"<<endl;
   
   cout<<"*** Reading  P1 file ****"<<endl;
   TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");

   Efficiency * EffpreselMCP = new Efficiency(inputHistoFile, "EffpreselMCP");
   Efficiency * EffpreselMCD = new Efficiency(inputHistoFile, "EffpreselMCD");

   string tagli[10]= {"Trigger","3of4 TOF","TRD Segments","Rigidity exists","Chi^2 R","Matching TOF","Matching TRD","In TRD Accept.","1 Particle","1 Tr. Track"};
   string nome;

   cout<<"**** MC PRESELECTIONS EFFICIENCY (FULL SET) ****"<<endl;

   EffpreselMCP -> Eval_Efficiency();
   EffpreselMCD -> Eval_Efficiency();

   TH1F * EffPreMCP_R_TH1F  =  (TH1F *)EffpreselMCP->effR	->Clone();
   TH1F * EffPreMCP_TH1F    =  (TH1F *)EffpreselMCP->effTOF->Clone();
   TH1F * EffPreMCPNaF_TH1F =  (TH1F *)EffpreselMCP->effNaF->Clone();
   TH1F * EffPreMCPAgl_TH1F =  (TH1F *)EffpreselMCP->effAgl->Clone();
   TH2F * EffPreMCD_R_TH2F  =  (TH2F *)EffpreselMCD->effR  ->Clone();
   TH2F * EffPreMCD_TH2F    =  (TH2F *)EffpreselMCD->effTOF->Clone();
   TH2F * EffPreMCDNaF_TH2F =  (TH2F *)EffpreselMCD->effNaF->Clone();
   TH2F * EffPreMCDAgl_TH2F =  (TH2F *)EffpreselMCD->effAgl->Clone();

	   finalHistos.Add(EffPreMCP_R_TH1F  	);
	   finalHistos.Add(EffPreMCP_TH1F       );
	   finalHistos.Add(EffPreMCPNaF_TH1F    );
	   finalHistos.Add(EffPreMCPAgl_TH1F 	);
	   finalHistos.Add(EffPreMCD_R_TH2F   	);
	   finalHistos.Add(EffPreMCD_TH2F    	);
	   finalHistos.Add(EffPreMCDNaF_TH2F 	);
	   finalHistos.Add(EffPreMCDAgl_TH2F 	);
	   finalHistos.writeObjsInFolder("Results");

	
	cout<<"*** Plotting ...  ****"<<endl;

	MCpreeff_Plot(

			EffPreMCP_R_TH1F , 
			EffPreMCP_TH1F   , 
			EffPreMCPNaF_TH1F, 
			EffPreMCPAgl_TH1F, 
			EffPreMCD_R_TH2F , 
			EffPreMCD_TH2F   , 
			EffPreMCDNaF_TH2F, 
			EffPreMCDAgl_TH2F 

		     );	


	return;

	
}


