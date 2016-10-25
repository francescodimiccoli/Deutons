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

   
   Efficiency * EffCSCFullsetMCP = new Efficiency(EffpreselMCP,"EffCSCFullsetMCP"); 

   EffCSCFullsetMCP->Compose_Efficiency(EffCSCMCP);	

   EffCSCFullsetMCP->Eval_FittedEfficiency();

   
   finalHistos.Add(EffCSCFullsetMCP->effR  ); 
   finalHistos.Add(EffCSCFullsetMCP->effR_fit);
   finalHistos.writeObjsInFolder("Results");

   cout<<"*** Plotting ...  ****"<<endl; 

   MCControlsamplecuteff_Plot(	

	EffCSCMCP ->effR);
	
   return;
}

