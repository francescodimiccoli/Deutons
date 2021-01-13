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
#include "../include/Globals.h"
#include "TKey.h"
#include "TFractionFitter.h"

#include "../include/Variables.hpp"
#include "../include/Cuts.h"
#include "../include/filesaver.h"

#include "../include/Efficiency.h"

#include "../include/PlottingFunctions.h"

#include "../include/AllRangesEfficiency.h"
#include "../include/Flux.h"

#include "TGraphAsymmErrors.h"


#include "RangeMerger.h"

TSpline3 * GetFluxSpline(TGraphAsymmErrors * Graph){
	int AMSpointnr=1;
	double a,b;
	while (Graph->GetPoint(AMSpointnr,a,b)>0) AMSpointnr++;

	double X[AMSpointnr-1];
	double Y[AMSpointnr-1];

	for(int i=0; i<AMSpointnr-1; i++)	int c= Graph->GetPoint(i+1,X[i],Y[i]);
	TSpline3 *AMSFlux = new TSpline3("AMSFlux",X,Y,AMSpointnr-1);
	return AMSFlux;	
}

TSpline3 * GetFluxSpline(TH1F * Graph,int begin=0){
	int AMSpointnr=Graph->GetNbinsX()-begin;

	double X[AMSpointnr-1];
	double Y[AMSpointnr-1];

	for(int i=0; i<AMSpointnr-1; i++) {X[i] = Graph->GetBinCenter(begin+i+1); Y[i]=Graph->GetBinContent(begin+i+1) ;}
	TSpline3 *AMSFlux = new TSpline3("AMSFlux",X,Y,AMSpointnr-1);
	return AMSFlux;	
}



TSpline3 * GetFluxSpline(TH1F * Graph, Binning bins){

	double X[Graph->GetNbinsX()];
	double Y[Graph->GetNbinsX()];

	for(int i=0; i<Graph->GetNbinsX(); i++){	X[i]=bins.EkPerMassBinCent(i); Y[i]=Graph->GetBinContent(i+1);}
	TSpline3 *AMSFlux = new TSpline3("AMSFlux",X,Y,Graph->GetNbinsX());
	return AMSFlux;	
}


int main(int argc, char * argv[]){
	cout<<"****************************** FILES OPENING ***************************************"<<endl;

	string INPUT(argv[1]);
	string OUTPUT(argv[2]);

	FileSaver finalResults;
	FileSaver Plots;
	finalResults.setName(INPUT.c_str());
	Plots.setName(OUTPUT.c_str());


	bool checkfile = finalResults.CheckFile();


	cout<<"****************************** BINS ***************************************"<<endl;

	SetUpTOIBinning();
	cout<<"****************************** READING DATABASES *************************************"<<endl;

	string filename2="./database_P.root";
	TFile * file2 = TFile::Open(filename2.c_str(),"READ");

	string filename3="./database_D.root";
	TFile * file3 = TFile::Open(filename3.c_str(),"READ");

	string filename4="./database_PD.root";
	TFile * file4 = TFile::Open(filename4.c_str(),"READ");

	std::vector<TGraphAsymmErrors *> P_Graphs;
	std::vector<TGraphAsymmErrors *> D_Graphs;
	std::vector<TGraphAsymmErrors *> PD_Graphs;

        TList *ExperimentsP = file2->GetListOfKeys();
        TIter nextP(ExperimentsP);
        TKey * keyP;
	TObject * obj;

        while((keyP = (TKey*)nextP())){
                obj = file2->Get(keyP->GetName());
                if(obj->InheritsFrom("TGraphAsymmErrors")) P_Graphs.push_back((TGraphAsymmErrors *)obj);
        }


	TList *ExperimentsD = file3->GetListOfKeys();
	TIter nextD(ExperimentsD);
	TKey * keyD;

	while((keyD = (TKey*)nextD())){
		obj = file3->Get(keyD->GetName());
		if(obj->InheritsFrom("TGraphAsymmErrors")) D_Graphs.push_back((TGraphAsymmErrors *)obj);
	}

	TList *ExperimentsPD = file4->GetListOfKeys();
	TIter nextPD(ExperimentsPD);
	TKey * keyPD;

	while((keyPD = (TKey*)nextPD())){
		obj = file4->Get(keyPD->GetName());
		if(obj->InheritsFrom("TGraphAsymmErrors")) PD_Graphs.push_back((TGraphAsymmErrors *)obj);
	}


	cout<<"****************************** PLOTTING FLUXES ***************************************"<<endl;

	Flux * HEPFlux = new Flux(finalResults,"PFluxHE","Acceptance_HE"  ,"Acceptance","HEPCounts/HEPCounts/HEPCounts_after","HEExposure"	,PRB);
	Flux * HEPFluxL1 = new Flux(finalResults,"PFluxL1HE","Acceptance_L1HE"  ,"Acceptance","HEPCountsL1/HEPCountsL1/HEPCountsL1_after","HEExposure"	,PRB);
	Flux * HEPFluxQ  = new Flux(finalResults,"PFluxQHE" ,"Acceptance_QualHE","Acceptance","HEPCountsQual/HEPCountsQual/HEPCountsQual_after","HEExposure",PRB);

	Flux * RigPTOF = new Flux(finalResults, "RigPTOF", "Acceptance_RigPTOF","Acceptance","TOFPCounts/TOFPCounts/TOFPCounts_before","ExposureTOF",GlobalRig.GetToFPBins());
	Flux * RigPNaF = new Flux(finalResults, "RigPNaF", "Acceptance_RigPNaF","Acceptance","NaFPCounts/NaFPCounts/NaFPCounts_before","ExposureNaF",GlobalRig.GetNaFPBins());
	Flux * RigPAgl = new Flux(finalResults, "RigPAgl", "Acceptance_RigPAgl","Acceptance","AglPCounts/AglPCounts/AglPCounts_before","ExposureAgl",GlobalRig.GetAglPBins());

	Flux * RigBetaPTOF = new Flux(finalResults, "RigBetaPTOF", "Acceptance_PTOF","Acceptance","TOFPCountsBeta/TOFPCountsBeta/TOFPCountsBeta_before","ExposureTOF",GlobalRig.GetToFPBins());
	Flux * RigBetaPNaF = new Flux(finalResults, "RigBetaPNaF", "Acceptance_PNaF","Acceptance","NaFPCountsBeta/NaFPCountsBeta/NaFPCountsBeta_before","ExposureNaF",GlobalRig.GetNaFPBins());
	Flux * RigBetaPAgl = new Flux(finalResults, "RigBetaPAgl", "Acceptance_PAgl","Acceptance","AglPCountsBeta/AglPCountsBeta/AglPCountsBeta_before","ExposureAgl",GlobalRig.GetAglPBins());

	Flux * FitBetaPTOF = new Flux(finalResults, "FitBetaPTOF", "Acceptance_PTOF","Acceptance","TOFPCountsFit/TOFPCountsFit/TOFPCountsFit_before","ExposureTOF",GlobalRig.GetToFPBins());
	Flux * FitBetaPNaF = new Flux(finalResults, "FitBetaPNaF", "Acceptance_PNaF","Acceptance","NaFPCountsFit/NaFPCountsFit/NaFPCountsFit_before","ExposureNaF",GlobalRig.GetNaFPBins());
	Flux * FitBetaPAgl = new Flux(finalResults, "FitBetaPAgl", "Acceptance_PAgl","Acceptance","AglPCountsFit/AglPCountsFit/AglPCountsFit_before","ExposureAgl",GlobalRig.GetAglPBins());


	SetUpRigTOIBinning();



	TCanvas *c3_ = new TCanvas("Proton Primary Flux");
	c3_->SetCanvasSize(2500,1500);
	TPad * c3_up = new TPad("upperPad", "upperPad",0.0,0.3,1.0,1.0);
	c3_up->Draw();

	TPad * c3_do = new TPad("lowerPad", "lowerPad",0.0,0.0,1.0,0.3);
	c3_do->Draw();

	c3_up->cd();

	gPad->SetLogx();
	gPad->SetLogy();
	gPad->SetGridx();
	gPad->SetGridy();


	PlotTH1FintoGraph(gPad,PRB, HEPFlux->GetFlux_rig(),"R T.o.I. [GV]", "Flux",2,false,"Psame",0.7,50,0.01,10000,"Inner",8);
	PlotTH1FintoGraph(gPad,PRB, HEPFluxL1->GetFlux_rig(),"R T.o.I. [GV]", "Flux",4,false,"Psame",0.7,50,0.01,10000,"Inner+L1",8);
	PlotTH1FintoGraph(gPad,PRB, HEPFluxQ->GetFlux_rig(),"R T.o.I. [GV]", "Flux",3,false,"Psame",0.7,50,0.01,10000,"Interactions",8);

	TFile * Extern = TFile::Open("proton54_20180708V2N_B1200400R2MCY7UFSMOOTH_totalQYAN.root");
	TH1F * AMSflux = (TH1F*) Extern->Get("Z1fluxh_total");
	TH1F * unfold = (TH1F*) Extern->Get("RawVUnfold");

	AMSflux->Multiply(unfold);

	TSpline3 *AMSFlux = GetFluxSpline(AMSflux);
	AMSflux->Draw("same");	

	TH1F * MyFlux = HEPFluxQ ->GetFlux_rig_plot();
	for(int i=0; i< MyFlux->GetNbinsX();i++) if(MyFlux->GetBinContent(i+1)>1000) MyFlux->SetBinContent(i+1,0);
	TSpline3 *MyFlux_spline = GetFluxSpline(MyFlux,2);
	MyFlux_spline->Draw("same");

	
	c3_do->cd();
	gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();


	PlotRatioWithSplineiAvg(gPad,PRB,    HEPFlux->GetFlux_rig() ,AMSFlux, "", "Flux",2,false,false,"Psame",0.7,50,0.5,1.5,"Inner",8);
	PlotRatioWithSplineiAvg(gPad,PRB,    HEPFluxL1->GetFlux_rig() ,AMSFlux, "", "Flux",4,false,false,"Psame",0.7,50,0.5,1.5,"Inner+L1",8);
	PlotRatioWithSplineiAvg(gPad,PRB,    HEPFluxQ ->GetFlux_rig() ,AMSFlux, "", "Flux",3,true,false,"Psame",0.7,50,0.5,1.5,"Interactions",8);


	Plots.Add(c3_);
	Plots.writeObjsInFolder("Fluxes/Fluxes R");


	TCanvas *c4 = new TCanvas("Binning Comparison");



	PlotRatioWithSplineintoGraph(gPad,GlobalRig.GetToFPBins(),RigPTOF->GetFlux_rig(),MyFlux_spline,"","Ratio",1,false,"Psame",0.7,25,0.5,1.5,"",20);	
	PlotRatioWithSplineintoGraph(gPad,GlobalRig.GetNaFPBins(),RigPNaF->GetFlux_rig(),MyFlux_spline,"","Ratio",1,false,"Psame",0.7,25,0.5,1.5,"",22);	
	PlotRatioWithSplineintoGraph(gPad,GlobalRig.GetAglPBins(),RigPAgl->GetFlux_rig(),MyFlux_spline,"","Ratio",1,false,"Psame",0.7,25,0.5,1.5,"",21);	

	PlotRatioWithSplineintoGraph(gPad,GlobalRig.GetToFPBins(),RigBetaPTOF->GetFlux_rig(),MyFlux_spline,"","Ratio",2,false,"Psame",0.7,25,0.5,1.5,"",20);	
	PlotRatioWithSplineintoGraph(gPad,GlobalRig.GetNaFPBins(),RigBetaPNaF->GetFlux_rig(),MyFlux_spline,"","Ratio",2,false,"Psame",0.7,25,0.5,1.5,"",22);	
	PlotRatioWithSplineintoGraph(gPad,GlobalRig.GetAglPBins(),RigBetaPAgl->GetFlux_rig(),MyFlux_spline,"","Ratio",2,false,"Psame",0.7,25,0.5,1.5,"",21);	

	PlotRatioWithSplineintoGraph(gPad,Global.GetToFPBins(),FitBetaPTOF->GetFlux_rig(),MyFlux_spline,"","Ratio",3,false,"Psame",0.7,25,0.5,1.5,"",20);	
	PlotRatioWithSplineintoGraph(gPad,Global.GetNaFPBins(),FitBetaPNaF->GetFlux_rig(),MyFlux_spline,"","Ratio",3,false,"Psame",0.7,25,0.5,1.5,"",22);	
	PlotRatioWithSplineintoGraph(gPad,Global.GetAglPBins(),FitBetaPAgl->GetFlux_rig(),MyFlux_spline,"","Ratio",3,false,"Psame",0.7,25,0.5,1.5,"",21);	


	
	Plots.Add(c4);
        Plots.writeObjsInFolder("Fluxes/Fluxes R");



	return 0;
}

