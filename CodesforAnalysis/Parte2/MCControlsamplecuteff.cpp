#include "PlottingFunctions/MCControlsamplecuteff_Plot.h"

using namespace std;

// eff. likelihood sel.
Efficiency * EffCSCMCP  = new Efficiency ("EffCSCMCP");
// eff. distance sel.

void MCControlsamplecuteff_Fill() {

	if(!trgpatt.IsPhysical()) return;
	if(Tup.Beta<=0||Tup.R<=0) return;

	int Kbin;
	
	if(Massa_gen<1) {
		//R bins
		Kbin=PRB.GetRBin(Tup.R);

		EffCSCMCP->beforeR->Fill(Kbin,Tup.mcweight);
		if(Herejcut && ProtonsMassWindow) EffCSCMCP->afterR->Fill(Kbin,Tup.mcweight);
	}

return;
}




void MCControlsamplecuteff_Write() {
   EffCSCMCP -> Write();
   return;
}


void MCControlsamplecuteff(string filename) {

   cout<<"******* MC Controlsamplecut SEL. EFFICIENCIES ********"<<endl;
    
   cout<<"*** Reading  P1 file ****"<<endl;
        TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");

   // eff. control sampel cuts
   Efficiency * EffCSCMCP  = new Efficiency (inputHistoFile,"EffCSCMCP");
   // eff. preselections
   Efficiency * EffpreselMCP = new Efficiency(inputHistoFile, "EffpreselMCP");


   cout<<"******* MC Controlsamplecut SEL. EFFICIENCIES ********"<<endl;


   EffCSCMCP ->Eval_Efficiency();
   EffpreselMCP ->Eval_Efficiency();

   TH1F * EffCSCMCP_TH1F 		=(TH1F *)EffCSCMCP ->effR  ->Clone();
   TH1F * EffCSCFullsetMCP_TH1F         =(TH1F *)EffpreselMCP ->effR  ->Clone();

   EffCSCFullsetMCP_TH1F -> Multiply(EffCSCMCP_TH1F);

   EffCSCFullsetMCP_TH1F -> SetName("EffCSCFullsetMCP_EffR");

   finalHistos.Add(EffCSCMCP_TH1F 	    );
   finalHistos.Add(EffCSCFullsetMCP_TH1F    ); 
   finalHistos.writeObjsInFolder("Results");

   cout<<"*** Plotting ...  ****"<<endl; 

   MCControlsamplecuteff_Plot(	

	EffCSCMCP_TH1F 	  	
	);
	
   return;
}

