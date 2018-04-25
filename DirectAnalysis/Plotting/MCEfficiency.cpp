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

#include "../include/FitError.h"
#include "../include/Efficiency.h"

#include "../include/PlottingFunctions.h"

#include "../include/AllRangesEfficiency.h"


void DrawEfficiencyComparison( AllRangesEfficiency * EffP,AllRangesEfficiency * EffD );


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

	AllRangesEfficiency * Trigger_P       = new AllRangesEfficiency(finalHistos,"Trigger_P","Trigger");
	AllRangesEfficiency * Fragmentation_P = new AllRangesEfficiency(finalHistos,"Fragmentation_P","Fragmentation");
	AllRangesEfficiency * Preselections_P = new AllRangesEfficiency(finalHistos,"PresEff_P","PreselectionEfficiency");
	AllRangesEfficiency * Quality_P       = new AllRangesEfficiency(finalHistos,"QualEff_P","QualityEfficiency"     );
        AllRangesEfficiency * FullSet_P      = new AllRangesEfficiency(finalHistos,"FullsetEff_P","FullsetEfficiency"   );
	AllRangesEfficiency * FullSetTOT_P      = new AllRangesEfficiency(finalHistos,"FullsetTOTEff_P","FullsetTOTEfficiency"   );


	AllRangesEfficiency * Trigger_D       = new AllRangesEfficiency(finalHistos,"Trigger_D","Trigger");
	AllRangesEfficiency * Fragmentation_D = new AllRangesEfficiency(finalHistos,"Fragmentation_D","Fragmentation");
	AllRangesEfficiency * Preselections_D = new AllRangesEfficiency(finalHistos,"PresEff_D","PreselectionEfficiency");
	AllRangesEfficiency * Quality_D       = new AllRangesEfficiency(finalHistos,"QualEff_D","QualityEfficiency"     );
        AllRangesEfficiency * FullSet_D       = new AllRangesEfficiency(finalHistos,"FullsetEff_D","FullsetEfficiency"  );
	AllRangesEfficiency * FullSetTOT_D    = new AllRangesEfficiency(finalHistos,"FullsetTOTEff_D","FullsetTOTEfficiency"   );


	AllRangesEfficiency * RICH_P = new AllRangesEfficiency(finalHistos,"RICHEff_P","RICHEfficiency");
	AllRangesEfficiency * RICH_D = new AllRangesEfficiency(finalHistos,"RICHEff_D","RICHEfficiency");

	AllRangesEfficiency * RICH_PQual = new AllRangesEfficiency(finalHistos,"RICHEff_PQual","RICHQualEfficiency");
	AllRangesEfficiency * RICH_DQual = new AllRangesEfficiency(finalHistos,"RICHEff_DQual","RICHQualEfficiency");




	TCanvas * c1 = new TCanvas("Efficiency Cascade"); 
        c1->SetCanvasSize(2000,1500);
	gPad->SetLogy();
	gPad->SetLogx();

	// Fix Legend	
	PlotTH1FintoGraph(gPad,ToFDB, (TH1F*)Preselections_P->EffTOF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",1,true,"Psame",0.1,10,1e-6,1,"TOF range",8);
	PlotTH1FintoGraph(gPad,NaFDB, (TH1F*)Preselections_P->EffNaF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",1,true,"Psame",0.1,10,1e-6,1,"NaF range",22);
	PlotTH1FintoGraph(gPad,AglDB, (TH1F*)Preselections_P->EffAgl->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",1,true,"Psame",0.1,10,1e-6,1,"Agl range",29);

	PlotTH1FintoGraph(gPad,ToFDB, (TH1F*)Preselections_P->EffTOF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",1,true,"Psame",0.1,10,1e-6,1,"Basic + Clean-Event",8);
	PlotTH1FintoGraph(gPad,ToFDB, (TH1F*)FullSet_P->EffTOF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",1,true,"Psame",0.1,10,1e-6,1,"Full Set (Basic + Clean-Event + Quality)",24);	
	
	PlotTH1FintoGraph(gPad,ToFDB, (TH1F*)Preselections_D->EffTOF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",4,true,"Psame",0.1,10,1e-6,1,"Deutons",8);
	PlotTH1FintoGraph(gPad,ToFDB, (TH1F*)Preselections_P->EffTOF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",2,true,"Psame",0.1,10,1e-6,1,"Protons",8);
	

	//
		

	PlotTH1FintoGraph(gPad,ToFDB, (TH1F*)Preselections_P->EffTOF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",2,true,"Psame",0.1,10,1e-6,1,"TOF range",8,true);
	PlotTH1FintoGraph(gPad,NaFDB, (TH1F*)Preselections_P->EffNaF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",2,true,"Psame",0.1,10,1e-6,1,"NaF range",22,true);
	PlotTH1FintoGraph(gPad,AglDB, (TH1F*)Preselections_P->EffAgl->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",2,true,"Psame",0.1,10,1e-6,1,"Agl range",29,true);

	PlotTH1FintoGraph(gPad,ToFDB, (TH1F*)Preselections_D->EffTOF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",4,true,"Psame",0.1,10,1e-6,1,"Deutons TOF",8,true);
	PlotTH1FintoGraph(gPad,NaFDB, (TH1F*)Preselections_D->EffNaF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",4,true,"Psame",0.1,10,1e-6,1,"Deutons NaF",22,true);
	PlotTH1FintoGraph(gPad,AglDB, (TH1F*)Preselections_D->EffAgl->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",4,true,"Psame",0.1,10,1e-6,1,"Deutons Agl",29,true);
	

	PlotTH1FintoGraph(gPad,ToFDB, (TH1F*)FullSet_P->EffTOF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",2,true,"Psame",0.1,10,1e-6,1,"Protons TOF",24,true);
	PlotTH1FintoGraph(gPad,NaFDB, (TH1F*)FullSet_P->EffNaF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",2,true,"Psame",0.1,10,1e-6,1,"Protons NaF",26,true);
	PlotTH1FintoGraph(gPad,AglDB, (TH1F*)FullSet_P->EffAgl->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",2,true,"Psame",0.1,10,1e-6,1,"Protons Agl",30,true);

	PlotTH1FintoGraph(gPad,ToFDB, (TH1F*)FullSet_D->EffTOF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",4,true,"Psame",0.1,10,1e-6,1,"Deutons TOF",24,true);
	PlotTH1FintoGraph(gPad,NaFDB, (TH1F*)FullSet_D->EffNaF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",4,true,"Psame",0.1,10,1e-6,1,"Deutons NaF",26,true);
	PlotTH1FintoGraph(gPad,AglDB, (TH1F*)FullSet_D->EffAgl->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",4,true,"Psame",0.1,10,1e-6,1,"Deutons Agl",30,true);

		
	Plots.Add(c1);
	Plots.writeObjsInFolder("Efficiencies");

	//trigger.
	TCanvas * c_ = new TCanvas("Trigger Efficiency"); 
        c_->SetCanvasSize(2000,1500);

	DrawEfficiencyComparison(Trigger_P,Trigger_D);

	Plots.Add(c_);
	Plots.writeObjsInFolder("Efficiencies");


	//fragm.
	TCanvas * c = new TCanvas(" Fragmentation Efficiency"); 
        c->SetCanvasSize(2000,1500);

	DrawEfficiencyComparison(Fragmentation_P,Fragmentation_D);

	Plots.Add(c);
	Plots.writeObjsInFolder("Efficiencies");

	//presel.
	TCanvas * d = new TCanvas(" Preselection Efficiency"); 
        d->SetCanvasSize(2000,1500);

	DrawEfficiencyComparison(Preselections_P,Preselections_D);

	Plots.Add(d);
	Plots.writeObjsInFolder("Efficiencies");

	//charge+BDT
	TCanvas * c2 = new TCanvas(" Quality Efficiency"); 
        c2->SetCanvasSize(2000,1500);

	DrawEfficiencyComparison(Quality_P,Quality_D);

	Plots.Add(c2);
	Plots.writeObjsInFolder("Efficiencies");

	//rich rec
	TCanvas * c3 = new TCanvas(" RICH rec. Efficiency"); 
        c3->SetCanvasSize(2000,1500);

	DrawEfficiencyComparison(RICH_P,RICH_D);
	
	Plots.Add(c3);
	Plots.writeObjsInFolder("Efficiencies");


	//rich BDT
	TCanvas * c4 = new TCanvas(" RICH Quality Efficiency"); 
        c4->SetCanvasSize(2000,1500);

	DrawEfficiencyComparison(RICH_PQual,RICH_DQual);
	
	Plots.Add(c4);
	Plots.writeObjsInFolder("Efficiencies");

	//global
	TCanvas * d4 = new TCanvas("Global Efficiency"); 
        d4->SetCanvasSize(2000,1500);

	DrawEfficiencyComparison(FullSetTOT_P,FullSetTOT_D);
	
	Plots.Add(d4);
	Plots.writeObjsInFolder("Efficiencies");



	return 0;
}




void DrawEfficiencyComparison( AllRangesEfficiency * EffP,AllRangesEfficiency * EffD ){

	TPad * c3_up = new TPad("upperPad", "upperPad",0.0,0.3,1.0,1.0);
	c3_up->Draw();

	TPad * c3_do = new TPad("lowerPad", "lowerPad",0.0,0.0,1.0,0.3);
	c3_do->Draw();

	c3_up->cd();
	
	//fix legend
	PlotTH1FintoGraph(gPad,ToFDB, (TH1F*)EffP->EffTOF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",1,true,"Psame",0.1,10,1e-6,1,"TOF range",8);
	PlotTH1FintoGraph(gPad,NaFDB, (TH1F*)EffP->EffNaF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",1,true,"Psame",0.1,10,1e-6,1,"NaF range",22);
	PlotTH1FintoGraph(gPad,AglDB, (TH1F*)EffP->EffAgl->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",1,true,"Psame",0.1,10,1e-6,1,"Agl range",29);

	PlotTH1FintoGraph(gPad,ToFDB, (TH1F*)EffP->EffTOF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",2,true,"Psame",0.1,10,1e-6,1,"Protons",8);
	PlotTH1FintoGraph(gPad,ToFDB, (TH1F*)EffD->EffTOF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",4,true,"Psame",0.1,10,1e-6,1,"Deutons",8);
	//	

	PlotTH1FintoGraph(gPad,ToFDB, (TH1F*)EffP->EffTOF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",2,true,"Psame",0.1,10,1e-6,1,"TOF range",8,true);
	PlotTH1FintoGraph(gPad,NaFDB, (TH1F*)EffP->EffNaF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",2,true,"Psame",0.1,10,1e-6,1,"NaF range",22,true);
	PlotTH1FintoGraph(gPad,AglDB, (TH1F*)EffP->EffAgl->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",2,true,"Psame",0.1,10,1e-6,1,"Agl range",29,true);

	PlotTH1FintoGraph(gPad,ToFDB, (TH1F*)EffD->EffTOF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",4,true,"Psame",0.1,10,1e-6,1,"Deutons TOF",8,true);
	PlotTH1FintoGraph(gPad,NaFDB, (TH1F*)EffD->EffNaF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",4,true,"Psame",0.1,10,1e-6,1,"Deutons NaF",22,true);
	PlotTH1FintoGraph(gPad,AglDB, (TH1F*)EffD->EffAgl->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",4,true,"Psame",0.1,10,1e-6,1,"Deutons Agl",29,true);
	
	c3_do->cd();
	gPad->SetGridy();
	
	TH1F * RatioTOF = (TH1F*)EffP->EffTOF->GetEfficiency()->Clone();
	TH1F * RatioNaF = (TH1F*)EffP->EffNaF->GetEfficiency()->Clone();
	TH1F * RatioAgl = (TH1F*)EffP->EffAgl->GetEfficiency()->Clone();
	RatioTOF ->Divide( (TH1F*)EffD->EffTOF->GetEfficiency()->Clone());
	RatioNaF ->Divide( (TH1F*)EffD->EffNaF->GetEfficiency()->Clone());
	RatioAgl ->Divide( (TH1F*)EffD->EffAgl->GetEfficiency()->Clone());

	PlotTH1FintoGraph(gPad,ToFDB, RatioTOF,"", "P/D Eff. ratio",1,true,"Psame",0.1,10,0.4,1.7,"TOF range",8,true);
	PlotTH1FintoGraph(gPad,NaFDB, RatioNaF,"", "P/D Eff. ratio",1,true,"Psame",0.1,10,0.4,1.7,"NaF range",22,true);
	PlotTH1FintoGraph(gPad,AglDB, RatioAgl,"", "P/D Eff. ratio",1,true,"Psame",0.1,10,0.4,1.7,"Agl range",29,true);



} 


