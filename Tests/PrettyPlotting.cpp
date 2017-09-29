#include <bitset>
#include "TROOT.h"
#include "TNtuple.h"
#include <TSpline.h>
#include "../include/binning.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include <TVector3.h>
#include "TMath.h"
#include <TFile.h>
#include "TFile.h"
#include "TH2.h"
#include "TF2.h"
#include <TVector3.h>
#include "TMath.h"
#include "TGraphErrors.h"
#include "TFractionFitter.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "../include/PlottingFunctions.h"

#include "../include/GlobalBinning.h"

#include "../Ntuple-making/Commonglobals.cpp"
#include "../include/Variables.hpp"
#include "../include/Cuts.h"

#include "../include/filesaver.h"
#include "../include/MCTuning.h"



void DrawTuning(Tuning * MCTuning,FileSaver finalHistos, FileSaver Plots);

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

        cout<<"****************************** PLOTTING RESULTS ***************************************"<<endl;

	Tuning * MCTuning = new Tuning(finalHistos,PRB,10,50,40);
	
	DrawTuning(MCTuning,finalHistos,Plots);	
	
	TH1F * ShiftBest     = (TH1F *)finalHistos.Get("Best Fit Parameters/Shift Best");
	TH1F * SigmaBest     = (TH1F *)finalHistos.Get("Best Fit Parameters/Sigma Best");
	
	TCanvas * c1 = new TCanvas("Parameters");
	c1->SetCanvasSize(2000,1500);

	TPad * c1_up = new TPad("upperPad", "upperPad",0.0,0.5,1.0,1.0);
	c1_up->Draw();

	TPad * c1_do = new TPad("lowerPad", "lowerPad",0.0,0.0,1.0,0.5);
	c1_do->Draw();

	PlotTH1FintoGraph(c1_up,PRB, ShiftBest, "#beta ToF",  "Mean shift [ps]",2,false,"ep",0.45,0.9,-100,100,"Best #chi^{2} Shift");
	PlotTH1FintoGraph(c1_do,PRB, SigmaBest, "#beta ToF",  "Additive #sigma [ps]",4,false,"ep",0.45,0.9,-40,180,"Best #chi^{2} #sigma");
		
	Plots.Add(c1);
        Plots.writeObjsInFolder("Parameters");
		

}


void DrawTuning(Tuning * MCTuning,FileSaver finalHistos,FileSaver Plots){

	TFile * infile = finalHistos.GetFile();	


	for(int i=1; i<MCTuning->GetBinning().size();i++){
		
		TH1F * OrMC = (TH1F *)infile->Get(("Bin "+to_string(i)+"/Results/OriginalMC_"+to_string(i)).c_str());
		TH1F * OrDT = (TH1F *)infile->Get(("Bin "+to_string(i)+"/Results/DistribData_"+to_string(i)).c_str());
		TH1F * OrDTPrimSec = (TH1F *)infile->Get(("Bin "+to_string(i)+"/Results/DistribDataPrimSec_"+to_string(i)).c_str());
		TH1F * MoMC = (TH1F *)infile->Get(("Bin "+to_string(i)+"/Results/Best Fit MC").c_str());
		TH1F * MoMCFiltered = (TH1F *)infile->Get(("Bin "+to_string(i)+"/Results/Best Fit MC (Cutoff filtered)").c_str());
		TH1F * ResidualBest     = (TH1F *)infile->Get(("Bin "+to_string(i)+"/Results/ResidualBest").c_str());
		TH1F * ResidualOriginal = (TH1F *)infile->Get(("Bin "+to_string(i)+"/Results/ResidualOriginal").c_str());
							

		TCanvas * c1 = new TCanvas("Distribution Comparison");
		c1->SetCanvasSize(2000,1500);

		TPad * c1_up = new TPad("upperPad", "upperPad",0.0,0.3,1.0,1.0);
		c1_up->Draw();

		TPad * c1_do = new TPad("lowerPad", "lowerPad",0.0,0.0,1.0,0.3);
		c1_do->Draw();


	        PlotDistribution(c1_up, OrMC,"Reconstructed Mass [GeV/c^2]","Distribution",2,"same",0.001,0.1,2,"Original MC",true);
		PlotDistribution(c1_up, OrDTPrimSec,"Reconstructed Mass [GeV/c^2]","Distribution",1,"e*same",0.001,0.1,2,"ISS data (prim. + sec.)",false,true);
		PlotDistribution(c1_up, OrDT,"Reconstructed Mass [GeV/c^2]","Distribution",1,"ePsame",0.001,0.1,2,"ISS data (prim. only)",false,true);
		PlotDistribution(c1_up, MoMC,"Reconstructed Mass [GeV/c^2]","Distribution",4,"same",0.001,0.1,2,"Best #chi^2 #beta smear");
		if(MoMCFiltered) PlotDistribution(c1_up, MoMCFiltered,"Reconstructed Mass [GeV/c^2]","Distribution",7,"same",0.001,0.1,2,"Best #chi^2 #beta smear + Cutoff Filter");

		PlotDistribution(c1_do, ResidualOriginal,"","Residuals",2,"ePsame",-15,15,2,"Original MC",false,true);
                PlotDistribution(c1_do, ResidualBest    ,"","Residuals",4,"ePsame",-15,15,2,"Best #chi^2 #beta smear",false,true);



		Plots.Add(c1);
		Plots.writeObjsInFolder(("Bin "+to_string(i)).c_str());	

		TH2F * Chi = (TH2F*) infile->Get(("Bin "+to_string(i)+"/ChiSquare/ChiSquare_"+to_string(i)).c_str());
		TCanvas * c2 = new TCanvas("ChiSquare");
		c2->SetCanvasSize(2000,1500);
		Chi->GetZaxis()->SetRangeUser(0.2,4);
		PlotTH2F(gPad, Chi, "Additive #sigma [ps]","Mean shift [ps]", "colz");

		Plots.Add(c2);
                Plots.writeObjsInFolder(("Bin "+to_string(i)).c_str());
		



	}

	return;

}
