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

TSpline3 * GetFluxSpline(TH1F * Graph){
	int AMSpointnr=Graph->GetNbinsX();

	double X[AMSpointnr-1];
	double Y[AMSpointnr-1];

	for(int i=0; i<AMSpointnr-1; i++) {X[i] = Graph->GetBinCenter(i+1); Y[i]=Graph->GetBinContent(i+1) ;}
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

	Flux * HEPFluxL1 = new Flux(finalResults,"PFluxL1HE","Acceptance_L1HE"  ,"Acceptance","HEPCountsL1/HEPCountsL1/HEPCountsL1_after","HEExposure"	,PRB);
	Flux * HEPFluxQ  = new Flux(finalResults,"PFluxQHE" ,"Acceptance_QualHE","Acceptance","HEPCountsQual/HEPCountsQual/HEPCountsQual_after","HEExposure",PRB);


	Flux * DFluxTOF  = new Flux(finalResults, "DFluxTOF", "Acceptance_DTOF","Acceptance","TOFDfits/Fit Results/Primary Deuteron Counts","ExposureTOF",Global.GetToFDBins());
	Flux * DFluxNaF  = new Flux(finalResults, "DFluxNaF", "Acceptance_DNaF","Acceptance","NaFDfits/Fit Results/Primary Deuteron Counts","ExposureNaF",Global.GetNaFDBins());
	Flux * DFluxAgl  = new Flux(finalResults, "DFluxAgl", "Acceptance_DAgl","Acceptance","AglDfits/Fit Results/Primary Deuteron Counts","ExposureAgl",Global.GetAglDBins());

	Flux * PFluxTOF  = new Flux(finalResults, "PFluxTOF", "Acceptance_PTOF","Acceptance","TOFPfits/Fit Results/Primary Proton Counts","ExposureTOF",Global.GetToFPBins());
	Flux * PFluxNaF  = new Flux(finalResults, "PFluxNaF", "Acceptance_PNaF","Acceptance","NaFPfits/Fit Results/Primary Proton Counts","ExposureNaF",Global.GetNaFPBins());
	Flux * PFluxAgl  = new Flux(finalResults, "PFluxAgl", "Acceptance_PAgl","Acceptance","AglPfits/Fit Results/Primary Proton Counts","ExposureAgl",Global.GetAglPBins());
	

	TCanvas *c1_ = new TCanvas("Exposure Time");
	c1_->SetCanvasSize(3000,1500);

//	PlotTH1FintoGraph(gPad,PRB, HEPFluxL1->GetExposureTime(),"Kinetic Energy [GeV/nucl.]", "Exposure Time [sec]",4,true,"Psame",0.1,50,1,2*HEPFluxL1->GetExposureTime()->GetBinContent(30),"HE range",8,true);
	
	PlotTH1FintoGraph(gPad,Global.GetToFDBins(), DFluxTOF->GetExposureTime(),"Kinetic Energy [GeV/nucl.]", "Exposure Time [sec]",4,true,"Psame",0.1,50,1,2*HEPFluxL1->GetExposureTime()->GetBinContent(30),"TOF range",8,true);
	PlotTH1FintoGraph(gPad,Global.GetNaFDBins(), DFluxNaF->GetExposureTime(),"Kinetic Energy [GeV/nucl.]", "Exposure Time [sec]",4,true,"Psame",0.1,50,1,2*HEPFluxL1->GetExposureTime()->GetBinContent(30),"NaF range",22,true);
	PlotTH1FintoGraph(gPad,Global.GetAglDBins(), DFluxAgl->GetExposureTime(),"Kinetic Energy [GeV/nucl.]", "Exposure Time [sec]",4,true,"Psame",0.1,50,1,2*HEPFluxL1->GetExposureTime()->GetBinContent(30),"Agl range",29,true);

	PlotTH1FintoGraph(gPad,Global.GetToFPBins(), PFluxTOF->GetExposureTime(),"Kinetic Energy [GeV/nucl.]", "Exposure Time [sec]",2,true,"Psame",0.1,50,1,2*HEPFluxL1->GetExposureTime()->GetBinContent(30),"TOF range",8,true);
	PlotTH1FintoGraph(gPad,Global.GetNaFPBins(), PFluxNaF->GetExposureTime(),"Kinetic Energy [GeV/nucl.]", "Exposure Time [sec]",2,true,"Psame",0.1,50,1,2*HEPFluxL1->GetExposureTime()->GetBinContent(30),"NaF range",22,true);
	PlotTH1FintoGraph(gPad,Global.GetAglPBins(), PFluxAgl->GetExposureTime(),"Kinetic Energy [GeV/nucl.]", "Exposure Time [sec]",2,true,"Psame",0.1,50,1,2*HEPFluxL1->GetExposureTime()->GetBinContent(30),"Agl range",29,true);

	Plots.Add(c1_);
	Plots.writeObjsInFolder("Fluxes");



	float potenza = 0;
	TCanvas *c2 = new TCanvas("Deuteron Primary Flux");
	c2->SetCanvasSize(2000,1500);

	TH2F * Frame = CreateFrame(gPad,0.01,25,0.0001,100,"Kin.En./nucl. [GeV/nucl.]","Flux [(m^2 sr GeV/nucl.)^{-1}]");	
	c2->cd();
	gPad->SetLogx();
	gPad->SetLogy();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraph* galprop3P=new TGraph();
	TGraph* galprop3P2=new TGraph();
	float x,y=0;
	int j=0;
	{
		string filename="./Galprop/Tom/deut_1500.dat";
		cout<<filename<<endl;
		ifstream fp(filename.c_str());
		while (!fp.eof()){
			fp>>x>>y;
			if(x/1e3>0.05&&x/1e3<=100)
				galprop3P->SetPoint(j,x/1e3,y*1e7*pow(x/1e3,potenza));
			j++;
		}
	}

	j=0;
	{
		string filename="./Galprop/Tom/deut_100.dat";
		cout<<filename<<endl;
		ifstream fp(filename.c_str());
		while (!fp.eof()){
			fp>>x>>y;
			if(x/1e3>0.05&&x/1e3<=100)
				galprop3P2->SetPoint(j,x/1e3,y*1e7*pow(x/1e3,potenza));
			j++;
		}
	}
	galprop3P->GetXaxis()->SetRangeUser(0.1,20);
	galprop3P->GetYaxis()->SetRangeUser(1e-3,1e3);

	galprop3P->SetTitle("Deutons Flux: Geo. Zones");
	galprop3P->GetXaxis()->SetTitle("Kin.En./nucl. [GeV/nucl.]");
	galprop3P ->GetYaxis()->SetTitle("Flux [(m^2 sr GeV/nucl.)^-1]");
	galprop3P ->GetXaxis()->SetTitleSize(0.045);
	galprop3P->GetYaxis()->SetTitleSize(0.045);
	galprop3P ->GetYaxis()->SetRangeUser(1e-2,1e4);

	Frame->Draw();
	galprop3P->Draw("sameC");
	galprop3P2->Draw("sameC");

	TLegend * leg =new TLegend(0.8, 0.1,0.95,0.95);
	leg->SetName("leg");
        for(uint n=0;n<D_Graphs.size();n++){
                D_Graphs[n] ->Draw("Psame");
                D_Graphs[n]->SetMarkerSize(2); 
                leg->AddEntry(D_Graphs[n],D_Graphs[n]->GetTitle(),"p");
        }

	leg->SetFillColor(0);
	leg->SetLineWidth(2);
	leg->Draw("same");
       


	PlotTH1FintoGraph(gPad,Global.GetToFDBins(), DFluxTOF->GetFlux(),"Kinetic Energy [GeV/nucl.]", "Flux",1,true,"Psame",0.1,100,0.0001,100,"This Work (TOF)",8,true);
	PlotTH1FintoGraph(gPad,Global.GetNaFDBins(), DFluxNaF->GetFlux(),"Kinetic Energy [GeV/nucl.]", "Flux",1,true,"Psame",0.1,100,0.0001,100,"This Work (NaF)",22,true);
	PlotTH1FintoGraph(gPad,Global.GetAglDBins(), DFluxAgl->GetFlux(),"Kinetic Energy [GeV/nucl.]", "Flux",1,true,"Psame",0.1,100,0.0001,100,"This Work (Agl)",29,true);

	Plots.Add(c2);
	Plots.writeObjsInFolder("Fluxes/Fluxes Ekin");



	TCanvas *c3 = new TCanvas("Deuteron Primary Flux - Merged");
	c3->SetCanvasSize(2000,1500);
	gPad->SetLogx();
	gPad->SetLogy();
	gPad->SetGridx();
	gPad->SetGridy();

	TH1F * MergedRange = Global.MergeSubDResult_D(DFluxTOF->GetFlux(),DFluxNaF->GetFlux(),DFluxAgl->GetFlux());
	TH1F * MergedRangeP = Global.MergeSubDResult_P(PFluxTOF->GetFlux(),PFluxNaF->GetFlux(),PFluxAgl->GetFlux());


	PlotTH1FintoGraph(gPad,Global.GetGlobalDBins(), MergedRange,"Kinetic Energy [GeV/nucl.]", "Flux",2,true,"Psame",0.1,100,0.01,1.2,"This Work",8);


	for(uint n=0;n<D_Graphs.size();n++){
		D_Graphs[n] ->Draw("Psame");
		D_Graphs[n]->SetMarkerSize(2); 
	}

	Plots.Add(c3);
	Plots.writeObjsInFolder("Fluxes/Fluxes Ekin");


	TCanvas *cl = new TCanvas("D/P ratio Ekin - Merged");
        cl->SetCanvasSize(2000,1500);
        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetGridx();
        gPad->SetGridy();


	TH1F * Ratio_Ekin = Global.MergedRatio_Ekin(MergedRange,MergedRangeP);
	PlotTH1FintoGraph(gPad,Global.GetGlobalDBins(), Ratio_Ekin,"Kinetic Energy [GeV/n]", "D/p ratio",2,true,"Psame",0.1,100,0.0001,10*MergedRange->GetBinContent(MergedRange->GetMaximumBin()),"This Work",8);

	Plots.Add(cl);
	Plots.writeObjsInFolder("Fluxes/Fluxes Ekin");







	TCanvas *c3_ = new TCanvas("Proton Primary Flux - Merged");
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

//	TH1F * MergedRange_P = Global.MergeSubDResult_P(PFluxTOF->GetFlux(),PFluxNaF->GetFlux(),PFluxAgl->GetFlux());

//	PlotTH1FintoGraph(gPad,Global.GetGlobalPBins(), MergedRange_P,"Kinetic Energy [GeV/nucl.]", "Flux",2,true,"Psame",0.1,100,0.0001,10*MergedRange_P->GetBinContent(MergedRange_P->GetMaximumBin()),"This Work",8);
	PlotTH1FintoGraph(gPad,HEPFluxL1->GetBins(),HEPFluxL1->GetFlux() ,"Kinetic Energy [GeV/nucl.]", "Flux",1,true,"Psame",0.1,100,0.0001,10*HEPFluxL1->GetFlux()->GetBinContent(HEPFluxL1->GetFlux()->GetMaximumBin()),"HE Baseline + L1",8);
	PlotTH1FintoGraph(gPad,HEPFluxQ->GetBins() ,HEPFluxQ->GetFlux() ,"Kinetic Energy [GeV/nucl.]", "Flux",4,true,"Psame",0.1,100,0.0001,10*HEPFluxL1->GetFlux()->GetBinContent(HEPFluxL1->GetFlux()->GetMaximumBin()),"HE Interactions",8);


	TFile * Extern = TFile::Open("proton54_20180708V2N_B1200400R2MCY7UFSMOOTH_totalQYAN.root");
	TH1F * AMSflux = (TH1F*) Extern->Get("Z1Fluxh_total");
	TH1F * unfold = (TH1F*) Extern->Get("RawVUnfold");

	AMSflux->Multiply(unfold);

	TSpline3 *AMSFlux = GetFluxSpline(AMSflux);
	AMSFlux->Draw();	

	for(uint n=0;n<P_Graphs.size();n++){
		P_Graphs[n] ->Draw("Psame");
		P_Graphs[n]->SetMarkerSize(2); 
	}

	c3_do->cd();
	gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();

	


	TGraphAsymmErrors * ErrorPubl = new TGraphAsymmErrors();
	for(int i=0;i<P_Graphs[1]->GetN();i++){
		int s=0;
		double x,y=1;
		s=P_Graphs[1]->GetPoint(i+1,x,y);
		ErrorPubl->SetPoint(i+1,x,1);
		ErrorPubl->SetPointError(i+1,0,0,1.3*P_Graphs[1]->GetErrorYlow(i+1)/y,1.3*P_Graphs[1]->GetErrorYhigh(i+1)/y);
	}
	ErrorPubl->SetLineColor(1);
	ErrorPubl->SetLineWidth(3);

	PlotRatioWithSplineintoGraph(gPad,PRB,    HEPFluxL1->GetFlux() ,AMSFlux, "Kinetic Energy [GeV/nucl.]", "Flux",1,true,"Psame",0.1,50,0.5,1.5,"(H.E. - Baseline + L1)",8);
	PlotRatioWithSplineintoGraph(gPad,PRB,    HEPFluxQ ->GetFlux() ,AMSFlux, "Kinetic Energy [GeV/nucl.]", "Flux",4,true,"Psame",0.1,50,0.5,1.5,"(H.E. - Interactions)",8);
	ErrorPubl->Draw("Psame");




//	Plots.Add(c3);
	Plots.Add(c3_);
	Plots.writeObjsInFolder("Fluxes/Fluxes Ekin");






	SetUpRigTOIBinning();

	TCanvas *c4 = new TCanvas("Deuteron Primary Flux");
	c4->SetCanvasSize(2000,1500);
	gPad->SetLogx();
	gPad->SetLogy();
	gPad->SetGridx();
	gPad->SetGridy();


	PlotTH1FintoGraph(gPad,Global.GetToFDBins(), DFluxTOF->GetFlux_rig(),"R T.o.I. [GV]", "Flux",1,false,"Psame",0.1,25,0.0001,100,"This Work (TOF)",8);
	PlotTH1FintoGraph(gPad,Global.GetNaFDBins(), DFluxNaF->GetFlux_rig(),"R T.o.I. [GV]", "Flux",1,false,"Psame",0.1,25,0.0001,100,"This Work (NaF)",22);
	PlotTH1FintoGraph(gPad,Global.GetAglDBins(), DFluxAgl->GetFlux_rig(),"R T.o.I. [GV]", "Flux",1,false,"Psame",0.1,25,0.0001,100,"This Work (Agl)",29);


	Plots.Add(c4);
	Plots.writeObjsInFolder("Fluxes/Fluxes R");

	TCanvas *c5 = new TCanvas("Deuteron Primary Flux - Merged");
	c5->SetCanvasSize(2000,1500);
	gPad->SetLogx();
	gPad->SetLogy();
	gPad->SetGridx();
	gPad->SetGridy();

	TH1F * MergedRange_rig = Global.MergeSubDResult_D(DFluxTOF->GetFlux_rig(),DFluxNaF->GetFlux_rig(),DFluxAgl->GetFlux_rig());

	PlotTH1FintoGraph(gPad,Global.GetGlobalDBins(), MergedRange_rig,"R T.o.I. [GV]", "Flux",2,false,"Psame",0.1,25,0.0001,100,"This Work",8);

	Plots.Add(c5);
	Plots.writeObjsInFolder("Fluxes/Fluxes R");


	TCanvas *c4_ = new TCanvas("Proton Primary Flux");
	c4_->SetCanvasSize(2000,1500);
	gPad->SetLogx();
	gPad->SetLogy();
	gPad->SetGridx();
	gPad->SetGridy();


	PlotTH1FintoGraph(gPad,Global.GetToFPBins(), PFluxTOF->GetFlux_rig(),"R T.o.I. [GV]", "Flux",1,false,"Psame",0.1,15,0.0001,2000,"This Work (TOF)",8);
	PlotTH1FintoGraph(gPad,Global.GetNaFPBins(), PFluxNaF->GetFlux_rig(),"R T.o.I. [GV]", "Flux",1,false,"Psame",0.1,15,0.0001,2000,"This Work (NaF)",22);
	PlotTH1FintoGraph(gPad,Global.GetAglPBins(), PFluxAgl->GetFlux_rig(),"R T.o.I. [GV]", "Flux",1,false,"Psame",0.1,15,0.0001,2000,"This Work (Agl)",29);


	Plots.Add(c4_);
	Plots.writeObjsInFolder("Fluxes/Fluxes R");

	TCanvas *c5_ = new TCanvas("Proton Primary Flux - Merged");
	c5_->SetCanvasSize(2000,1500);
	gPad->SetLogx();
	gPad->SetLogy();
	gPad->SetGridx();
	gPad->SetGridy();

	TH1F * MergedRange_rig_P = Global.MergeSubDResult_P(PFluxTOF->GetFlux_rig(),PFluxNaF->GetFlux_rig(),PFluxAgl->GetFlux_rig());

	PlotTH1FintoGraph(gPad,Global.GetGlobalPBins(), MergedRange_rig_P,"R T.o.I. [GV]", "Flux",2,false,"Psame",0.4,35,0.0001,10000,"This Work",8);
	PlotTH1FintoGraph(gPad,PRB, HEPFluxL1->GetFlux_rig(),"R T.o.I. [GV]", "Flux",1,false,"Psame",0.4,35,0.0001,10000,"This Work",8);
	PlotTH1FintoGraph(gPad,PRB, HEPFluxQ->GetFlux_rig(),"R T.o.I. [GV]", "Flux",1,false,"Psame",0.4,35,0.0001,10000,"This Work",8);





	

	Plots.Add(c5_);
	Plots.writeObjsInFolder("Fluxes/Fluxes R");


	TCanvas *c6_ = new TCanvas("D/P ratio");
	c6_->SetCanvasSize(2000,1500);
	gPad->SetLogx();
	gPad->SetLogy();
	gPad->SetGridx();
	gPad->SetGridy();

	TH1F * MergedRatio_rig = Global.MergedRatio(MergedRange_rig,MergedRange_rig_P);

	PlotTH1FintoGraph(gPad,Global.GetGlobalDBins(), MergedRatio_rig,"R T.o.I. [GV]", "Flux",2,false,"Psame",0.1,25,0.0001,100,"This Work",8);
	Plots.Add(c6_);
	Plots.writeObjsInFolder("Fluxes/Fluxes R");


	TCanvas * c7 = new TCanvas("Uncertainty");

	TH1F * Acc_TOF = DFluxTOF->GetEffAcceptance();
	TH1F * Counts_TOF = DFluxTOF->GetCounts();
	TH1F * Acc_TOF_Err    = DFluxTOF->GetEffAcceptance();
	TH1F * Counts_TOF_Err = DFluxTOF->GetCounts();
	Acc_TOF_Err  ->Clear(); 
        Counts_TOF_Err->Clear();
	for(int i=0;i<Acc_TOF->GetNbinsX();i++){
		Acc_TOF_Err  ->SetBinContent(i+1,Acc_TOF->GetBinError(i+1)/Acc_TOF->GetBinContent(i+1));
		Acc_TOF_Err  ->SetBinError(i+1,0);
		Counts_TOF_Err  ->SetBinContent(i+1,Counts_TOF->GetBinError(i+1)/Counts_TOF->GetBinContent(i+1));
		Counts_TOF_Err  ->SetBinError(i+1,0);
	}
	TH1F * Acc_NaF = DFluxNaF->GetEffAcceptance();
	TH1F * Counts_NaF = DFluxNaF->GetCounts();
	TH1F * Acc_NaF_Err    = DFluxNaF->GetEffAcceptance();
	TH1F * Counts_NaF_Err = DFluxNaF->GetCounts();
	Acc_NaF_Err  ->Clear(); 
        Counts_NaF_Err->Clear();
	for(int i=0;i<Acc_NaF->GetNbinsX();i++){
		Acc_NaF_Err  ->SetBinContent(i+1,Acc_NaF->GetBinError(i+1)/Acc_NaF->GetBinContent(i+1));
		Acc_NaF_Err  ->SetBinError(i+1,0);
		Counts_NaF_Err  ->SetBinContent(i+1,Counts_NaF->GetBinError(i+1)/Counts_NaF->GetBinContent(i+1));
		Counts_NaF_Err  ->SetBinError(i+1,0);
	}

	TH1F * Acc_Agl = DFluxAgl->GetEffAcceptance();
	TH1F * Counts_Agl = DFluxAgl->GetCounts();
	TH1F * Acc_Agl_Err    = DFluxAgl->GetEffAcceptance();
	TH1F * Counts_Agl_Err = DFluxAgl->GetCounts();
	Acc_Agl_Err  ->Clear(); 
        Counts_Agl_Err->Clear();
	for(int i=0;i<Acc_Agl->GetNbinsX();i++){
		Acc_Agl_Err  ->SetBinContent(i+1,Acc_Agl->GetBinError(i+1)/Acc_Agl->GetBinContent(i+1));
		Acc_Agl_Err  ->SetBinError(i+1,0);
		Counts_Agl_Err  ->SetBinContent(i+1,Counts_Agl->GetBinError(i+1)/Counts_Agl->GetBinContent(i+1));
		Counts_Agl_Err  ->SetBinError(i+1,0);
	}


	TH1F * MergedRangeAcc = Global.MergeSubDResult_D(Acc_TOF_Err,Acc_NaF_Err,Acc_Agl_Err);
	TH1F * MergedRangeCou = Global.MergeSubDResult_D(Counts_TOF_Err,Counts_NaF_Err,Counts_Agl_Err);
	TH1F * TOT = CreateHisto("Total Error",Global.GetGlobalDBins());

	TH1F * Cou = ConvertBinnedHisto(MergedRangeAcc,"Counts Error",Global.GetGlobalDBins());
	TH1F * Acc = ConvertBinnedHisto(MergedRangeCou,"Acceptance Error",Global.GetGlobalDBins());

	for(int i=0;i<TOT->GetNbinsX();i++){
	TOT->SetBinContent(i+1,pow(pow(Cou->GetBinContent(i+1),2)+  pow(Acc->GetBinContent(i+1),2),0.5));
	}

	TOT->Draw();
	Cou->Draw("same");
	Acc->Draw("same");
	
	Plots.Add(c7);
        Plots.writeObjsInFolder("Fluxes");



	return 0;
}

