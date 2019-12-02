#include <bitset>
#include "TROOT.h"
#include "TNtuple.h"
#include <TSpline.h>
#include "../include/binning.h"
#include "TFile.h"
#include "TH2.h"
#include "TF2.h"
#include <TVector3.h>
#include "TMath.h"
#include <TFile.h>
#include "TFile.h"
#include "TH2.h"
#include "TF2.h"
#include <TVector3.h>
#include "TMath.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TKey.h"
#include "../include/Globals.h"
#include "TFractionFitter.h"

#include "../include/Variables.hpp"
#include "../include/Cuts.h"
#include "../include/filesaver.h"

#include "../include/Efficiency.h"

#include "../include/PlottingFunctions.h"

#include "../include/AllRangesEfficiency.h"
#include "../include/Acceptance.h"


int main(int argc, char * argv[]){

        cout<<"****************************** FILES OPENING ***************************************"<<endl;

        string INPUT(argv[1]);
        string OUTPUT(argv[2]);

        FileSaver finalHistos;
        FileSaver Plots;

        finalHistos.setName(INPUT.c_str());
        Plots.setName(OUTPUT.c_str());


        bool checkfile = finalHistos.CheckFile();


        cout<<"****************************** BINS ***************************************"<<endl;
	SetUpTOIBinning();

        cout<<"****************************** PLOTTING EFFICIENCIES ***************************************"<<endl;


	// Efficiency
        Efficiency * Cascade1 = new Efficiency(finalHistos,"Cascade1","Cascade1",PRB);
        Efficiency * Cascade2 = new Efficiency(finalHistos,"Cascade2","Cascade2",PRB);
        Efficiency * Cascade3 = new Efficiency(finalHistos,"Cascade3","Cascade3",PRB);
        Efficiency * Cascade4 = new Efficiency(finalHistos,"Cascade4","Cascade4",PRB);
        Efficiency * Cascade5 = new Efficiency(finalHistos,"Cascade5","Cascade5",PRB);
        Efficiency * Cascade6 = new Efficiency(finalHistos,"Cascade6","Cascade6",PRB);
        Efficiency * Cascade7 = new Efficiency(finalHistos,"Cascade7","Cascade7",PRB);
        Efficiency * Cascade8 = new Efficiency(finalHistos,"Cascade8","Cascade8",PRB);


	TCanvas * c1 = new TCanvas("Efficiency Cascade"); 
        c1->SetCanvasSize(2000,1500);
	gPad->SetLogy();
	gPad->SetLogx();

	TH1F * Cascadehistos[8];
	Cascadehistos[0] = (TH1F *) Cascade1 ->GetEfficiency();
	Cascadehistos[1] = (TH1F *) Cascade2 ->GetEfficiency();
	Cascadehistos[2] = (TH1F *) Cascade3 ->GetEfficiency();
	Cascadehistos[3] = (TH1F *) Cascade4 ->GetEfficiency();
	Cascadehistos[4] = (TH1F *) Cascade5 ->GetEfficiency();
	Cascadehistos[5] = (TH1F *) Cascade6 ->GetEfficiency();
	Cascadehistos[6] = (TH1F *) Cascade7 ->GetEfficiency();
	Cascadehistos[7] = (TH1F *) Cascade8 ->GetEfficiency();

	//for(int i=0;i<8;i++) Cascadehistos[i]->Divide(Cascadehistos[0]);

	PlotTH1FintoGraph(gPad,PRB, (TH1F*)Cascadehistos[0],"Kinetic Energy [GeV/n.]", "Efficiency",1,true,"same",0.1,50,5e-4,1.1,"MinimumBias",8);
	PlotTH1FintoGraph(gPad,PRB, (TH1F*)Cascadehistos[1],"Kinetic Energy [GeV/n.]", "Efficiency",2,true,"same",0.1,50,5e-4,1.1,"LooseQ1",8);
	PlotTH1FintoGraph(gPad,PRB, (TH1F*)Cascadehistos[2],"Kinetic Energy [GeV/n.]", "Efficiency",3,true,"same",0.1,50,5e-4,1.1,"Cleaning",8);
	PlotTH1FintoGraph(gPad,PRB, (TH1F*)Cascadehistos[3],"Kinetic Energy [GeV/n.]", "Efficiency",4,true,"same",0.1,50,5e-4,1.1,"Good Time",8);
	PlotTH1FintoGraph(gPad,PRB, (TH1F*)Cascadehistos[4],"Kinetic Energy [GeV/n.]", "Efficiency",5,true,"same",0.1,50,5e-4,1.1,"CIEMAT NaF",8);
	PlotTH1FintoGraph(gPad,PRB, (TH1F*)Cascadehistos[5],"Kinetic Energy [GeV/n.]", "Efficiency",6,true,"same",0.1,50,5e-4,1.1,"CIEMAT Agl",8);
	PlotTH1FintoGraph(gPad,PRB, (TH1F*)Cascadehistos[6],"Kinetic Energy [GeV/n.]", "Efficiency",7,true,"same",0.1,50,5e-4,1.1,"BDT NaF",8);
	PlotTH1FintoGraph(gPad,PRB, (TH1F*)Cascadehistos[7],"Kinetic Energy [GeV/n.]", "Efficiency",8,true,"same",0.1,50,5e-4,1.1,"BDT Agl",8);

	Plots.Add(c1);
	Plots.writeObjsInFolder("Efficiency Cascade");

	//Acceptance
	TCanvas * c2 = new TCanvas("Effective Acceptance Ekin"); 
        c2->SetCanvasSize(2000,1500);
	gPad->SetLogy();
	gPad->SetLogx();
	
	Acceptance * Acceptance_HE     = new Acceptance(finalHistos,"Acceptance_HE"	,"Acceptance",PRB);
	Acceptance * Acceptance_L1HE   = new Acceptance(finalHistos,"Acceptance_L1HE"	,"Acceptance",PRB);
	Acceptance * Acceptance_QualHE = new Acceptance(finalHistos,"Acceptance_QualHE" ,"Acceptance",PRB);

	Acceptance * Acceptance_PTOF = new Acceptance(finalHistos,"Acceptance_PTOF","Acceptance",Global.GetToFPBins());
        Acceptance * Acceptance_PNaF = new Acceptance(finalHistos,"Acceptance_PNaF","Acceptance",Global.GetNaFPBins());
        Acceptance * Acceptance_PAgl = new Acceptance(finalHistos,"Acceptance_PAgl","Acceptance",Global.GetAglPBins());

        Acceptance * Acceptance_DTOF = new Acceptance(finalHistos,"Acceptance_DTOF","Acceptance",Global.GetToFDBins());
        Acceptance * Acceptance_DNaF = new Acceptance(finalHistos,"Acceptance_DNaF","Acceptance",Global.GetNaFDBins());
        Acceptance * Acceptance_DAgl = new Acceptance(finalHistos,"Acceptance_DAgl","Acceptance",Global.GetAglDBins());

	PlotTH1FintoGraph(gPad,Global.GetToFPBins(),   Acceptance_PTOF->GetEffAccMC(), "Ekin [GeV/n]", "Acceptance [m^2 sr]",2,true,"Csame",0.1,50,1e-4,3,"Proton",8);	
	PlotTH1FintoGraph(gPad,Global.GetNaFPBins(),   Acceptance_PNaF->GetEffAccMC(), "Ekin [GeV/n]", "Acceptance [m^2 sr]",2,true,"Csame",0.1,50,1e-4,3,"",8,true);	
	PlotTH1FintoGraph(gPad,Global.GetAglPBins(),   Acceptance_PAgl->GetEffAccMC(), "Ekin [GeV/n]", "Acceptance [m^2 sr]",2,true,"Csame",0.1,50,1e-4,3,"",8,true);	
	
	PlotTH1FintoGraph(gPad,Global.GetToFDBins(),   Acceptance_DTOF->GetEffAccMC(), "Ekin [GeV/n]", "Acceptance [m^2 sr]",4,true,"Csame",0.1,50,1e-4,3,"Deuton",8);	
	PlotTH1FintoGraph(gPad,Global.GetNaFDBins(),   Acceptance_DNaF->GetEffAccMC(), "Ekin [GeV/n]", "Acceptance [m^2 sr]",4,true,"Csame",0.1,50,1e-4,3,"",8,true);	
	PlotTH1FintoGraph(gPad,Global.GetAglDBins(),   Acceptance_DAgl->GetEffAccMC(), "Ekin [GeV/n]", "Acceptance [m^2 sr]",4,true,"Csame",0.1,50,1e-4,3,"",8,true);	

	PlotTH1FintoGraph(gPad,Global.GetToFPBins(),   Acceptance_PTOF->GetEffAcc(), "Ekin [GeV/n]", "Acceptance [m^2 sr]",2,true,"Csame",0.1,50,1e-4,3,"Proton",8);	
	PlotTH1FintoGraph(gPad,Global.GetNaFPBins(),   Acceptance_PNaF->GetEffAcc(), "Ekin [GeV/n]", "Acceptance [m^2 sr]",2,true,"Csame",0.1,50,1e-4,3,"",8,true);	
	PlotTH1FintoGraph(gPad,Global.GetAglPBins(),   Acceptance_PAgl->GetEffAcc(), "Ekin [GeV/n]", "Acceptance [m^2 sr]",2,true,"Csame",0.1,50,1e-4,3,"",8,true);	
	
	PlotTH1FintoGraph(gPad,Global.GetToFDBins(),   Acceptance_DTOF->GetEffAcc(), "Ekin [GeV/n]", "Acceptance [m^2 sr]",4,true,"Csame",0.1,50,1e-4,3,"Deuton",8);	
	PlotTH1FintoGraph(gPad,Global.GetNaFDBins(),   Acceptance_DNaF->GetEffAcc(), "Ekin [GeV/n]", "Acceptance [m^2 sr]",4,true,"Csame",0.1,50,1e-4,3,"",8,true);	
	PlotTH1FintoGraph(gPad,Global.GetAglDBins(),   Acceptance_DAgl->GetEffAcc(), "Ekin [GeV/n]", "Acceptance [m^2 sr]",4,true,"Csame",0.1,50,1e-4,3,"",8,true);	

	PlotTH1FintoGraph(gPad,PRB,   Acceptance_L1HE  ->GetEffAccMC(), "Ekin [GeV/n]", "Acceptance [m^2 sr]",1,true,"Csame",0.1,50,1e-4,3,"",8,true);	
	PlotTH1FintoGraph(gPad,PRB,   Acceptance_QualHE->GetEffAccMC(), "Ekin [GeV/n]", "Acceptance [m^2 sr]",1,true,"Csame",0.1,50,1e-4,3,"",8,true);	


	Plots.Add(c2);
	Plots.writeObjsInFolder("Effective Acceptance");

	TCanvas * c3 = new TCanvas("Effective Acceptance R"); 
        c3->SetCanvasSize(2000,1500);
	gPad->SetLogy();
	gPad->SetLogx();

	SetUpRigTOIBinning();
	PlotTH1FintoGraph(gPad,Global.GetToFPBins(),   Acceptance_PTOF->GetEffAccMC(), "R [GV]", "Acceptance [m^2 sr]",2,false,"Csame",0.1,50,1e-4,3,"Proton",8);	
	PlotTH1FintoGraph(gPad,Global.GetNaFPBins(),   Acceptance_PNaF->GetEffAccMC(), "R [GV]", "Acceptance [m^2 sr]",2,false,"Csame",0.1,50,1e-4,3,"",8,true);	
	PlotTH1FintoGraph(gPad,Global.GetAglPBins(),   Acceptance_PAgl->GetEffAccMC(), "R [GV]", "Acceptance [m^2 sr]",2,false,"Csame",0.1,50,1e-4,3,"",8,true);	
	
	PlotTH1FintoGraph(gPad,Global.GetToFDBins(),   Acceptance_DTOF->GetEffAccMC(), "R [GV]", "Acceptance [m^2 sr]",4,false,"Csame",0.1,50,1e-4,3,"Deuton",8);	
	PlotTH1FintoGraph(gPad,Global.GetNaFDBins(),   Acceptance_DNaF->GetEffAccMC(), "R [GV]", "Acceptance [m^2 sr]",4,false,"Csame",0.1,50,1e-4,3,"",8,true);	
	PlotTH1FintoGraph(gPad,Global.GetAglDBins(),   Acceptance_DAgl->GetEffAccMC(), "R [GV]", "Acceptance [m^2 sr]",4,false,"Csame",0.1,50,1e-4,3,"",8,true);	

	PlotTH1FintoGraph(gPad,Global.GetToFPBins(),   Acceptance_PTOF->GetEffAcc(), "R [GV]", "Acceptance [m^2 sr]",2,false,"Csame",0.1,50,1e-4,3,"Proton",8);	
	PlotTH1FintoGraph(gPad,Global.GetNaFPBins(),   Acceptance_PNaF->GetEffAcc(), "R [GV]", "Acceptance [m^2 sr]",2,false,"Csame",0.1,50,1e-4,3,"",8,true);	
	PlotTH1FintoGraph(gPad,Global.GetAglPBins(),   Acceptance_PAgl->GetEffAcc(), "R [GV]", "Acceptance [m^2 sr]",2,false,"Csame",0.1,50,1e-4,3,"",8,true);	
	
	PlotTH1FintoGraph(gPad,Global.GetToFDBins(),   Acceptance_DTOF->GetEffAcc(), "R [GV]", "Acceptance [m^2 sr]",4,false,"Csame",0.1,50,1e-4,3,"Deuton",8);	
	PlotTH1FintoGraph(gPad,Global.GetNaFDBins(),   Acceptance_DNaF->GetEffAcc(), "R [GV]", "Acceptance [m^2 sr]",4,false,"Csame",0.1,50,1e-4,3,"",8,true);	
	PlotTH1FintoGraph(gPad,Global.GetAglDBins(),   Acceptance_DAgl->GetEffAcc(), "R [GV]", "Acceptance [m^2 sr]",4,false,"Csame",0.1,50,1e-4,3,"",8,true);	

	PlotTH1FintoGraph(gPad,PRB,   Acceptance_L1HE->GetEffAccMC(), "R [GV]", "Acceptance [m^2 sr]",1,false,"Csame",0.1,50,1e-4,3,"",8,true);	
	PlotTH1FintoGraph(gPad,PRB,   Acceptance_QualHE->GetEffAccMC(), "R [GV]", "Acceptance [m^2 sr]",1,false,"Csame",0.1,50,1e-4,3,"",8,true);	


	Plots.Add(c3);
	Plots.writeObjsInFolder("Effective Acceptance");


	return 0;
}

