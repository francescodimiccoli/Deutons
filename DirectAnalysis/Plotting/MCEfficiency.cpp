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
#include "../include/GlobalBinning.h"
#include "TKey.h"
#include "TFractionFitter.h"

#include "../include/Variables.hpp"
#include "../include/Cuts.h"
#include "../include/filesaver.h"

#include "../include/FitError.h"
#include "../include/Resolution.h"
#include "../include/Efficiency.h"

#include "../include/PlottingFunctions.h"

#include "../include/AllRangesEfficiency.h"


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

        SetBins();

	PRB.Print();

        cout<<"**TOF**"<<endl;
        ToFDB.Print();

        cout<<"**NaF**"<<endl;
        NaFDB.Print();

        cout<<"**Agl**"<<endl;
        AglDB.Print();

        ToFDB.UseBetaEdges();
        NaFDB.UseBetaEdges();
        AglDB.UseBetaEdges();

        PRB.UseREdges();


        cout<<endl;

        cout<<"****************************** PLOTTING EFFICIENCIES ***************************************"<<endl;

	AllRangesEfficiency * Preselections_P = new AllRangesEfficiency(finalHistos,"PresEff_P","PreselectionEfficiency");
	AllRangesEfficiency * Quality_P       = new AllRangesEfficiency(finalHistos,"QualEff_P","QualityEfficiency"     );
        AllRangesEfficiency * FullSet_P      = new AllRangesEfficiency(finalHistos,"FullsetEff_P","FullsetEfficiency"   );


	AllRangesEfficiency * Preselections_D = new AllRangesEfficiency(finalHistos,"PresEff_D","PreselectionEfficiency");
	AllRangesEfficiency * Quality_D       = new AllRangesEfficiency(finalHistos,"QualEff_D","QualityEfficiency"     );
        AllRangesEfficiency * FullSet_D       = new AllRangesEfficiency(finalHistos,"FullsetEff_D","FullsetEfficiency"  );

	AllRangesEfficiency * RICH_P = new AllRangesEfficiency(finalHistos,"RICHEff_P","RICHEfficiency");
	AllRangesEfficiency * RICH_D = new AllRangesEfficiency(finalHistos,"RICHEff_D","RICHEfficiency");




	TCanvas * c1 = new TCanvas("Efficiency Cascade"); 
        c1->SetCanvasSize(2000,1500);

	// Fix Legend	
	PlotTH1FintoGraph(gPad,ToFDB, (TH1F*)Preselections_P->EffTOF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",1,true,"Psame",0.1,10,0,1,"TOF range",8);
	PlotTH1FintoGraph(gPad,NaFDB, (TH1F*)Preselections_P->EffNaF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",1,true,"Psame",0.1,10,0,1,"NaF range",22);
	PlotTH1FintoGraph(gPad,AglDB, (TH1F*)Preselections_P->EffAgl->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",1,true,"Psame",0.1,10,0,1,"Agl range",29);

	PlotTH1FintoGraph(gPad,ToFDB, (TH1F*)Preselections_P->EffTOF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",1,true,"Psame",0.1,10,0,1,"Basic + Clean-Event",8);
	PlotTH1FintoGraph(gPad,ToFDB, (TH1F*)FullSet_P->EffTOF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",1,true,"Psame",0.1,10,0,1,"Full Set (Basic + Clean-Event + Quality)",24);	
	
	PlotTH1FintoGraph(gPad,ToFDB, (TH1F*)Preselections_D->EffTOF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",4,true,"Psame",0.1,10,0,1,"Deutons",8);
	PlotTH1FintoGraph(gPad,ToFDB, (TH1F*)Preselections_P->EffTOF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",2,true,"Psame",0.1,10,0,1,"Protons",8);
	

	//
		

	PlotTH1FintoGraph(gPad,ToFDB, (TH1F*)Preselections_P->EffTOF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",2,true,"Psame",0.1,10,0,1,"TOF range",8,true);
	PlotTH1FintoGraph(gPad,NaFDB, (TH1F*)Preselections_P->EffNaF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",2,true,"Psame",0.1,10,0,1,"NaF range",22,true);
	PlotTH1FintoGraph(gPad,AglDB, (TH1F*)Preselections_P->EffAgl->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",2,true,"Psame",0.1,10,0,1,"Agl range",29,true);

	PlotTH1FintoGraph(gPad,ToFDB, (TH1F*)Preselections_D->EffTOF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",4,true,"Psame",0.1,10,0,1,"Deutons TOF",8,true);
	PlotTH1FintoGraph(gPad,NaFDB, (TH1F*)Preselections_D->EffNaF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",4,true,"Psame",0.1,10,0,1,"Deutons NaF",22,true);
	PlotTH1FintoGraph(gPad,AglDB, (TH1F*)Preselections_D->EffAgl->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",4,true,"Psame",0.1,10,0,1,"Deutons Agl",29,true);
	

	PlotTH1FintoGraph(gPad,ToFDB, (TH1F*)FullSet_P->EffTOF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",2,true,"Psame",0.1,10,0,1,"Protons TOF",24,true);
	PlotTH1FintoGraph(gPad,NaFDB, (TH1F*)FullSet_P->EffNaF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",2,true,"Psame",0.1,10,0,1,"Protons NaF",26,true);
	PlotTH1FintoGraph(gPad,AglDB, (TH1F*)FullSet_P->EffAgl->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",2,true,"Psame",0.1,10,0,1,"Protons Agl",30,true);

	PlotTH1FintoGraph(gPad,ToFDB, (TH1F*)FullSet_D->EffTOF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",4,true,"Psame",0.1,10,0,1,"Deutons TOF",24,true);
	PlotTH1FintoGraph(gPad,NaFDB, (TH1F*)FullSet_D->EffNaF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",4,true,"Psame",0.1,10,0,1,"Deutons NaF",26,true);
	PlotTH1FintoGraph(gPad,AglDB, (TH1F*)FullSet_D->EffAgl->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",4,true,"Psame",0.1,10,0,1,"Deutons Agl",30,true);

		
	Plots.Add(c1);
	Plots.writeObjsInFolder("Efficiencies");

	TCanvas * c2 = new TCanvas(" Quality Efficiency"); 
        c2->SetCanvasSize(2000,1500);

	//fix legend
	PlotTH1FintoGraph(gPad,ToFDB, (TH1F*)Quality_P->EffTOF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",1,true,"Psame",0.1,10,0,1,"TOF range",8);
	PlotTH1FintoGraph(gPad,NaFDB, (TH1F*)Quality_P->EffNaF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",1,true,"Psame",0.1,10,0,1,"NaF range",22);
	PlotTH1FintoGraph(gPad,AglDB, (TH1F*)Quality_P->EffAgl->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",1,true,"Psame",0.1,10,0,1,"Agl range",29);

	PlotTH1FintoGraph(gPad,ToFDB, (TH1F*)Quality_P->EffTOF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",2,true,"Psame",0.1,10,0,1,"Protons",8);
	PlotTH1FintoGraph(gPad,ToFDB, (TH1F*)Quality_D->EffTOF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",4,true,"Psame",0.1,10,0,1,"Deutons",8);
	//	

	PlotTH1FintoGraph(gPad,ToFDB, (TH1F*)Quality_P->EffTOF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",2,true,"Psame",0.1,10,0,1,"TOF range",8,true);
	PlotTH1FintoGraph(gPad,NaFDB, (TH1F*)Quality_P->EffNaF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",2,true,"Psame",0.1,10,0,1,"NaF range",22,true);
	PlotTH1FintoGraph(gPad,AglDB, (TH1F*)Quality_P->EffAgl->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",2,true,"Psame",0.1,10,0,1,"Agl range",29,true);

	PlotTH1FintoGraph(gPad,ToFDB, (TH1F*)Quality_D->EffTOF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",4,true,"Psame",0.1,10,0,1,"Deutons TOF",8,true);
	PlotTH1FintoGraph(gPad,NaFDB, (TH1F*)Quality_D->EffNaF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",4,true,"Psame",0.1,10,0,1,"Deutons NaF",22,true);
	PlotTH1FintoGraph(gPad,AglDB, (TH1F*)Quality_D->EffAgl->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",4,true,"Psame",0.1,10,0,1,"Deutons Agl",29,true);
	
	Plots.Add(c2);
	Plots.writeObjsInFolder("Efficiencies");

	TCanvas * c3 = new TCanvas(" RICH rec. Efficiency"); 
        c3->SetCanvasSize(2000,1500);

	//fix legend
	PlotTH1FintoGraph(gPad,NaFDB, (TH1F*)RICH_P->EffNaF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",1,true,"Psame",0.5,10,0,0.9,"NaF range",22);
	PlotTH1FintoGraph(gPad,AglDB, (TH1F*)RICH_P->EffAgl->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",1,true,"Psame",0.5,10,0,0.9,"Agl range",29);

	PlotTH1FintoGraph(gPad,ToFDB, (TH1F*)RICH_P->EffTOF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",2,true,"Psame",0.5,10,0,0.9,"Protons",8);
	PlotTH1FintoGraph(gPad,ToFDB, (TH1F*)RICH_D->EffTOF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",4,true,"Psame",0.5,10,0,0.9,"Deutons",8);
	//

	PlotTH1FintoGraph(gPad,NaFDB, (TH1F*)RICH_P->EffNaF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",2,true,"Psame",0.5,10,0,0.9,"NaF range",22,true);
	PlotTH1FintoGraph(gPad,AglDB, (TH1F*)RICH_P->EffAgl->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",2,true,"Psame",0.5,10,0,0.9,"Agl range",29,true);

	PlotTH1FintoGraph(gPad,NaFDB, (TH1F*)RICH_D->EffNaF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",4,true,"Psame",0.5,10,0,0.9,"Deutons NaF",22,true);
	PlotTH1FintoGraph(gPad,AglDB, (TH1F*)RICH_D->EffAgl->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",4,true,"Psame",0.5,10,0,0.9,"Deutons Agl",29,true);
	
	Plots.Add(c3);
	Plots.writeObjsInFolder("Efficiencies");



	return 0;
}
