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
	SetUpTOIBinning();

        cout<<"****************************** PLOTTING EFFICIENCIES ***************************************"<<endl;

	AllRangesEfficiency * FullSetTOT_P = new AllRangesEfficiency(finalHistos,"FullSetTOT_P","FullSetTOT",ToFPB,NaFPB,AglPB);

        AllRangesEfficiency * FullSetTOT_D = new AllRangesEfficiency(finalHistos,"FullSetTOT_D","FullSetTOT",ToFDB,NaFDB,AglDB);

        Efficiency * Cascade1 = new Efficiency(finalHistos,"Cascade1","Cascade1",PRB);
        Efficiency * Cascade2 = new Efficiency(finalHistos,"Cascade2","Cascade2",PRB);
        Efficiency * Cascade3 = new Efficiency(finalHistos,"Cascade3","Cascade3",PRB);
        Efficiency * Cascade4 = new Efficiency(finalHistos,"Cascade4","Cascade4",PRB);
        Efficiency * Cascade5 = new Efficiency(finalHistos,"Cascade5","Cascade5",PRB);
        Efficiency * Cascade6 = new Efficiency(finalHistos,"Cascade6","Cascade6",PRB);
        Efficiency * Cascade7 = new Efficiency(finalHistos,"Cascade7","Cascade7",PRB);
        Efficiency * Cascade8 = new Efficiency(finalHistos,"Cascade8","Cascade8",PRB);

        AllRangesEfficiency * Trigger_P = new AllRangesEfficiency(finalHistos,"Trigger_P","Trigger",ToFPB,NaFPB,AglPB);
        AllRangesEfficiency * Trigger_D = new AllRangesEfficiency(finalHistos,"Trigger_D","Trigger",ToFDB,NaFDB,AglDB);
        AllRangesEfficiency * MinimBias_P = new AllRangesEfficiency(finalHistos,"MinimBias_P","MinimBias",ToFPB,NaFPB,AglPB);
        AllRangesEfficiency * MinimBias_D = new AllRangesEfficiency(finalHistos,"MinimBias_D","MinimBias",ToFDB,NaFDB,AglDB);
        AllRangesEfficiency * Cleaning_P = new AllRangesEfficiency(finalHistos,"Cleaning_P","Cleaning",ToFPB,NaFPB,AglPB);
        AllRangesEfficiency * Cleaning_D = new AllRangesEfficiency(finalHistos,"Cleaning_D","Cleaning",ToFDB,NaFDB,AglDB);
        AllRangesEfficiency * RICH_P = new AllRangesEfficiency(finalHistos,"RICH_P","RICH",ToFPB,NaFPB,AglPB);
        AllRangesEfficiency * RICH_D = new AllRangesEfficiency(finalHistos,"RICH_D","RICH",ToFDB,NaFDB,AglDB);
        AllRangesEfficiency * RICHQual_P = new AllRangesEfficiency(finalHistos,"RICHQual_P","RICHQual",ToFPB,NaFPB,AglPB);
        AllRangesEfficiency * RICHQual_D = new AllRangesEfficiency(finalHistos,"RICHQual_D","RICHQual",ToFDB,NaFDB,AglDB);
	AllRangesEfficiency * GoldenTOF_P = new AllRangesEfficiency(finalHistos,"GoldenTOF_P","GoldenTOF",ToFPB,NaFPB,AglPB);
	AllRangesEfficiency * GoldenTOF_D = new AllRangesEfficiency(finalHistos,"GoldenTOF_D","GoldenTOF",ToFDB,NaFDB,AglDB);

	AllRangesEfficiency * Trigger_P_PID = new AllRangesEfficiency(finalHistos,"Trigger_P_PID","Trigger",ToFPB,NaFPB,AglPB);
        AllRangesEfficiency * Trigger_D_PID = new AllRangesEfficiency(finalHistos,"Trigger_D_PID","Trigger",ToFDB,NaFDB,AglDB);
        AllRangesEfficiency * MinimBias_P_PID = new AllRangesEfficiency(finalHistos,"MinimBias_P_PID","MinimBias",ToFPB,NaFPB,AglPB);
        AllRangesEfficiency * MinimBias_D_PID = new AllRangesEfficiency(finalHistos,"MinimBias_D_PID","MinimBias",ToFDB,NaFDB,AglDB);
        AllRangesEfficiency * Cleaning_P_PID = new AllRangesEfficiency(finalHistos,"Cleaning_P_PID","Cleaning",ToFPB,NaFPB,AglPB);
        AllRangesEfficiency * Cleaning_D_PID = new AllRangesEfficiency(finalHistos,"Cleaning_D_PID","Cleaning",ToFDB,NaFDB,AglDB);
        AllRangesEfficiency * RICH_P_PID = new AllRangesEfficiency(finalHistos,"RICH_P_PID","RICH",ToFPB,NaFPB,AglPB);
        AllRangesEfficiency * RICH_D_PID = new AllRangesEfficiency(finalHistos,"RICH_D_PID","RICH",ToFDB,NaFDB,AglDB);
        AllRangesEfficiency * RICHQual_P_PID = new AllRangesEfficiency(finalHistos,"RICHQual_P_PID","RICHQual",ToFPB,NaFPB,AglPB);
        AllRangesEfficiency * RICHQual_D_PID = new AllRangesEfficiency(finalHistos,"RICHQual_D_PID","RICHQual",ToFDB,NaFDB,AglDB);
	AllRangesEfficiency * GoldenTOF_P_PID = new AllRangesEfficiency(finalHistos,"GoldenTOF_P_PID","GoldenTOF",ToFPB,NaFPB,AglPB);
	AllRangesEfficiency * GoldenTOF_D_PID = new AllRangesEfficiency(finalHistos,"GoldenTOF_D_PID","GoldenTOF",ToFDB,NaFDB,AglDB);
	
	AllRangesEfficiency * Fragmentation_P = new AllRangesEfficiency(finalHistos,"Fragmentation_P","Fragmentation",ToFPB,NaFPB,AglPB);
	AllRangesEfficiency * Fragmentation_D = new AllRangesEfficiency(finalHistos,"Fragmentation_D","Fragmentation",ToFDB,NaFDB,AglDB);



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
	Plots.writeObjsInFolder("Efficiencies");

	//trigger.
	TCanvas * c_ = new TCanvas("Trigger Efficiency"); 
	c_->SetCanvasSize(2000,1500);
	c_->Divide(2,1);

	c_->cd(1);
	DrawEfficiencyComparison(Trigger_P,Trigger_D);
	c_->cd(2);
	DrawEfficiencyComparison(Trigger_P_PID,Trigger_D_PID);


	Plots.Add(c_);
	Plots.writeObjsInFolder("Efficiencies");


	//fragm.
	TCanvas * c = new TCanvas(" Minimum Bias Efficiency"); 
	c->SetCanvasSize(2000,1500);
	c->Divide(2,1);
	
	c->cd(1);
	DrawEfficiencyComparison(MinimBias_P,MinimBias_D);
	c->cd(2);
	DrawEfficiencyComparison(MinimBias_P_PID,MinimBias_D_PID);

	Plots.Add(c);
	Plots.writeObjsInFolder("Efficiencies");

	//basic.
	TCanvas * d_ = new TCanvas(" Cleaning Efficiency"); 
        d_->SetCanvasSize(2000,1500);
	d_->Divide(2,1);

	d_->cd(1);
	DrawEfficiencyComparison(Cleaning_P,Cleaning_D);
	d_->cd(2);
	DrawEfficiencyComparison(Cleaning_P_PID,Cleaning_D_PID);

	Plots.Add(d_);
	Plots.writeObjsInFolder("Efficiencies");

	//rich rec
	TCanvas * c3 = new TCanvas(" RICH rec. Efficiency"); 
        c3->SetCanvasSize(2000,1500);
	c3->Divide(2,1);

	c3->cd(1);
	DrawEfficiencyComparison(RICH_P,RICH_D);
	c3->cd(2);
	DrawEfficiencyComparison(RICH_P_PID,RICH_D_PID);
	
	Plots.Add(c3);
	Plots.writeObjsInFolder("Efficiencies");


	//rich BDT
	TCanvas * c4 = new TCanvas(" RICH Quality Efficiency"); 
        c4->SetCanvasSize(2000,1500);
	c4->Divide(2,1);
	
	c4->cd(1);
	DrawEfficiencyComparison(RICHQual_P,RICHQual_D);
	c4->cd(2);
	DrawEfficiencyComparison(RICHQual_P_PID,RICHQual_D_PID);
	
	Plots.Add(c4);
	Plots.writeObjsInFolder("Efficiencies");

	//goldenTOF	
	TCanvas * c2 = new TCanvas("GoldenTOF Efficiency"); 
        c2->SetCanvasSize(2000,1500);
	c2->Divide(2,1);

	c2->cd(1);
	DrawEfficiencyComparison(GoldenTOF_P,GoldenTOF_D);
	c2->cd(2);
	DrawEfficiencyComparison(GoldenTOF_P_PID,GoldenTOF_D_PID);

	Plots.Add(c2);
	Plots.writeObjsInFolder("Efficiencies");


	//global
	TCanvas * d4 = new TCanvas("Global Efficiency"); 
        d4->SetCanvasSize(2000,1500);

	DrawEfficiencyComparison(FullSetTOT_P,FullSetTOT_D);
	
	Plots.Add(d4);
	Plots.writeObjsInFolder("Efficiencies");

	//fragmentation
	TCanvas * d4_ = new TCanvas("Deuteron Survival");
        d4_->SetCanvasSize(2000,1500);

        DrawEfficiencyComparison(Fragmentation_P,Fragmentation_D);

        Plots.Add(d4_);
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
	PlotTH1FintoGraph(gPad,ToFPB, (TH1F*)EffP->EffTOF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",1,true,"Psame",0.1,10,1e-6,1,"TOF range",8);
	PlotTH1FintoGraph(gPad,NaFPB, (TH1F*)EffP->EffNaF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",1,true,"Psame",0.1,10,1e-6,1,"NaF range",22);
	PlotTH1FintoGraph(gPad,AglPB, (TH1F*)EffP->EffAgl->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",1,true,"Psame",0.1,10,1e-6,1,"Agl range",29);

	PlotTH1FintoGraph(gPad,ToFPB, (TH1F*)EffP->EffTOF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",2,true,"Psame",0.1,10,1e-6,1,"Protons",8);
	PlotTH1FintoGraph(gPad,ToFDB, (TH1F*)EffD->EffTOF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",4,true,"Psame",0.1,10,1e-6,1,"Deutons",8);
	//	

	PlotTH1FintoGraph(gPad,ToFPB, (TH1F*)EffP->EffTOF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",2,true,"Psame",0.1,10,1e-6,1,"TOF range",8,true);
	PlotTH1FintoGraph(gPad,NaFPB, (TH1F*)EffP->EffNaF->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",2,true,"Psame",0.1,10,1e-6,1,"NaF range",22,true);
	PlotTH1FintoGraph(gPad,AglPB, (TH1F*)EffP->EffAgl->GetEfficiency(),"Kinetic Energy [GeV/nucl.]", "Efficiency",2,true,"Psame",0.1,10,1e-6,1,"Agl range",29,true);

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

	PlotTH1FintoGraph(gPad,ToFDB, RatioTOF,"", "P/D Eff. ratio",1,true,"Psame",0.1,10,0.6,1.3,"TOF range",8,true);
	PlotTH1FintoGraph(gPad,NaFDB, RatioNaF,"", "P/D Eff. ratio",1,true,"Psame",0.1,10,0.6,1.3,"NaF range",22,true);
	PlotTH1FintoGraph(gPad,AglDB, RatioAgl,"", "P/D Eff. ratio",1,true,"Psame",0.1,10,0.6,1.3,"Agl range",29,true);



} 


