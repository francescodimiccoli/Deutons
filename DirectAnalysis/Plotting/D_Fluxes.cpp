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

#include "../include/FitError.h"
#include "../include/Resolution.h"
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

	PlotTH1FintoGraph(gPad,PRB, HEPFluxL1->GetExposureTime(),"Kinetic Energy [GeV/nucl.]", "Exposure Time [sec]",4,true,"Psame",0.1,50,1,2*HEPFluxL1->GetExposureTime()->GetBinContent(30),"HE range",8,true);
	
	PlotTH1FintoGraph(gPad,Global.GetToFDBins(), DFluxTOF->GetExposureTime(),"Kinetic Energy [GeV/nucl.]", "Exposure Time [sec]",4,true,"Psame",0.1,50,1,2*HEPFluxL1->GetExposureTime()->GetBinContent(30),"TOF range",8,true);
	PlotTH1FintoGraph(gPad,Global.GetNaFDBins(), DFluxNaF->GetExposureTime(),"Kinetic Energy [GeV/nucl.]", "Exposure Time [sec]",4,true,"Psame",0.1,50,1,2*HEPFluxL1->GetExposureTime()->GetBinContent(30),"NaF range",22,true);
	PlotTH1FintoGraph(gPad,Global.GetAglDBins(), DFluxAgl->GetExposureTime(),"Kinetic Energy [GeV/nucl.]", "Exposure Time [sec]",4,true,"Psame",0.1,50,1,2*HEPFluxL1->GetExposureTime()->GetBinContent(30),"Agl range",29,true);

	PlotTH1FintoGraph(gPad,Global.GetToFPBins(), PFluxTOF->GetExposureTime(),"Kinetic Energy [GeV/nucl.]", "Exposure Time [sec]",2,true,"Psame",0.1,50,1,2*HEPFluxL1->GetExposureTime()->GetBinContent(30),"TOF range",8,true);
	PlotTH1FintoGraph(gPad,Global.GetNaFPBins(), PFluxNaF->GetExposureTime(),"Kinetic Energy [GeV/nucl.]", "Exposure Time [sec]",2,true,"Psame",0.1,50,1,2*HEPFluxL1->GetExposureTime()->GetBinContent(30),"NaF range",22,true);
	PlotTH1FintoGraph(gPad,Global.GetAglPBins(), PFluxAgl->GetExposureTime(),"Kinetic Energy [GeV/nucl.]", "Exposure Time [sec]",2,true,"Psame",0.1,50,1,2*HEPFluxL1->GetExposureTime()->GetBinContent(30),"Agl range",29,true);

	Plots.Add(c1_);
	Plots.writeObjsInFolder("Fluxes");



	float potenza = 0;
	/*TCanvas *c2 = new TCanvas("Deuteron Primary Flux");
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

	PlotTH1FintoGraph(gPad,Global.GetGlobalDBins(), MergedRange,"Kinetic Energy [GeV/nucl.]", "Flux",2,true,"Psame",0.1,100,0.0001,10*MergedRange->GetBinContent(MergedRange->GetMaximumBin()),"This Work",8);


	for(uint n=0;n<D_Graphs.size();n++){
		D_Graphs[n] ->Draw("Psame");
		D_Graphs[n]->SetMarkerSize(2); 
	}

	Plots.Add(c3);
	Plots.writeObjsInFolder("Fluxes/Fluxes Ekin");

*/

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

	TSpline3 *AMSFlux = GetFluxSpline(P_Graphs[1]);
	AMSFlux->Draw("same");	

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


/*
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

	PlotTH1FintoGraph(gPad,Global.GetGlobalPBins(), MergedRatio_rig,"R T.o.I. [GV]", "Flux",2,false,"Psame",0.1,25,0.0001,100,"This Work",8);
	Plots.Add(c6_);
	Plots.writeObjsInFolder("Fluxes/Fluxes R");


*/
	/*
	   TCanvas *c3 = new TCanvas("D/P ratio");
	   c3->SetCanvasSize(2000,1500);
	   j=0;
	   TGraph* galprop3DP=new TGraph();
	   TGraph* galprop3DP2=new TGraph();
	   {
	   string filename="./Galprop/Tom/dP_1500.dat";
	   cout<<filename<<endl;
	   ifstream fp(filename.c_str());
	   while (!fp.eof()){
	   fp>>x>>y;
	   if(x/1e3>0.05&&x/1e3<=100)
	   galprop3DP->SetPoint(j,x/(0.5*1e3),y);
	   j++;
	   }
	   }

	   j=0;
	   {
	   string filename="./Galprop/Tom/dP_500.dat";
	   cout<<filename<<endl;
	   ifstream fp(filename.c_str());
	   while (!fp.eof()){
	   fp>>x>>y;
	   if(x/1e3>0.05&&x/1e3<=100)
	   galprop3DP2->SetPoint(j,x/(0.5*1e3),y);
	   j++;
	   }
	   }
	   galprop3DP->GetXaxis()->SetRangeUser(0.1,20);

	   galprop3DP->SetTitle("Deutons Flux: Geo. Zones");
	   galprop3DP->GetXaxis()->SetTitle("Kin.En./nucl. [GeV/nucl.]");
	   galprop3DP ->GetYaxis()->SetTitle("Flux [(m^2 sr GeV/nucl.)^-1]");
	   galprop3DP ->GetXaxis()->SetTitleSize(0.045);

	   TH1F * DPRatioTOF = (TH1F *)finalResults.Get("Fluxes/DP ratio TOF");
	   TH1F * DPRatioNaF = (TH1F *)finalResults.Get("Fluxes/DP ratio NaF");
	   TH1F * DPRatioAgl = (TH1F *)finalResults.Get("Fluxes/DP ratio Agl");

	   DPRatioNaF->Smooth();	

	   cout<<DPRatioTOF<<endl;

	   PlotTH1FintoGraph(gPad,Global.GetToFDBins(), DPRatioTOF ,"Kinetic Energy [GeV/nucl.]", "Flux",1,true,"Psame",0.1,10,0.00001,0.12,"This Work (TOF)",8);
	   PlotTH1FintoGraph(gPad,Global.GetNaFDBins(), DPRatioNaF ,"Kinetic Energy [GeV/nucl.]", "Flux",1,true,"Psame",0.1,10,0.00001,0.12,"This Work (NaF)",22);
	   PlotTH1FintoGraph(gPad,Global.GetAglDBins(), DPRatioAgl ,"Kinetic Energy [GeV/nucl.]", "Flux",1,true,"Psame",0.1,10,0.00001,0.12,"This Work (Agl)",29);

	//	PlotMergedRanges(gPad,DPRatioTOF ,DPRatioNaF ,DPRatioAgl ,"Kinetic Energy [GeV/nucl.]", "Flux",1,true,"Psame",0.1,10,0.00001,0.12,"This Work (TOF)",8);	

	galprop3DP->Draw("sameC");
	galprop3DP2->Draw("sameC");

	leg = (TLegend*) gPad->FindObject("leg");

	for(uint n=0;n<PD_Graphs.size();n++){
	PD_Graphs[n] ->Draw("Psame");
	PD_Graphs[n]->SetMarkerSize(2); 
	leg->AddEntry(PD_Graphs[n],PD_Graphs[n]->GetTitle(),"p");
	}


	Plots.Add(c3);
	Plots.writeObjsInFolder("Fluxes");

*/
	return 0;
}
