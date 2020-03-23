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
#include "../include/GlobalPaths.h"
#include "TKey.h"
#include "TFractionFitter.h"
#include "TPaveLabel.h"

#include "../include/Variables.hpp"
#include "../include/Cuts.h"
#include "../include/filesaver.h"
#include "TGraphAsymmErrors.h"
#include "../include/PlottingFunctions.h"
#include <sstream>
#include "TStyle.h"
#include "TMultiGraph.h"
#include "TGaxis.h"
#include "../perl/List.h"
#include "../include/Flux.h"
#include "../include/TemplateFITbetasmear.h"
#include "../include/EffCorr.h"
#include <chrono>
#include <iomanip>
#include <ctime>

int Binselected[6]= {0,1,2,3,4,5};
int Binselected2[6]= {3,7,11,15,27,29};


int start =0;



struct TimeResult{
	TH1F * Mean;
	TH1F * VariabilityUp;
	TH1F * VariabilityDw;

	void Scale(float s){Mean->Scale(s); VariabilityUp->Scale(s); VariabilityDw->Scale(s);};
	void Smooth(){Mean->Smooth(); VariabilityUp->Smooth(); VariabilityDw->Smooth();};

	
};


TimeResult GetTimeMean(std::vector<TH1F*>Ratios){

	TH1F * RatioMean = (TH1F*)Ratios[0] -> Clone();
	TH1F * RatioMin  = (TH1F*)Ratios[0] -> Clone();
	TH1F * RatioMax  = (TH1F*)Ratios[0] -> Clone();

	for(int j=0;j<RatioMean->GetNbinsX();j++){
		for(int i=start+1;i<Ratios.size();i++){
			if(Ratios[i]->GetBinContent(j+1)<RatioMin->GetBinContent(j+1)&&!(Ratios[i]->GetBinContent(j+1)<0.2*Ratios[1] ->GetBinContent(j+1))) 
				RatioMin->SetBinContent(j+1, Ratios[i]->GetBinContent(j+1));
			if(Ratios[i]->GetBinContent(j+1)>RatioMax->GetBinContent(j+1) &&!(Ratios[i]->GetBinContent(j+1)>3*Ratios[1] ->GetBinContent(j+1))) 
				RatioMax->SetBinContent(j+1, Ratios[i]->GetBinContent(j+1));
		}
	}
	
	int n =0;
	for(int i=start+1;i<Ratios.size();i++){
		if(!(Ratios[i]->GetBinContent(Ratios[i]->GetMaximumBin())>3*(Ratios[1]->GetBinContent(Ratios[i]->GetMaximumBin())))){
			if(!(Ratios[i]->GetBinContent(Ratios[i]->GetMinimumBin())<0.2*(Ratios[1]->GetBinContent(Ratios[i]->GetMinimumBin()))))	{
				RatioMean->Add(Ratios[i]);
				n++;
			}	
		}	
	}
	RatioMean->Scale(1./n);

	TimeResult Result;
	Result.Mean = (TH1F*) RatioMean->Clone();
	Result.VariabilityUp = (TH1F*) RatioMax->Clone();
	Result.VariabilityDw = (TH1F*) RatioMean->Clone();

	Result.VariabilityUp->Add(Result.Mean,-1);
	Result.VariabilityDw->Add(RatioMin,-1);
	
	Result.VariabilityUp->Smooth();
	Result.VariabilityDw->Smooth();
		
	return Result;
}


double FitTime(double *x, double *p){
	float value;
	if(x[0]<1423958400) value = p[0];
	else value = p[1]*x[0] + ( p[0] - p[1]*1423958400);
	return value;
}


std::string Convert (float number){
    std::ostringstream buff;
    buff<<number;
    std::string output= buff.str();
    output.erase(4,output.end()-output.begin()-4);	
    return output;   
}


int bartels = (
1307750400,
1310083200,
1312416000,
1314748800,
1317081600,
1319414400,
1321747200,
1324080000,
1326412800,
1328745600,
1331078400,
1333411200,
1335744000,
1338076800,
1340409600,
1342742400,
1345075200,
1347408000,
1349740800,
1352073600,
1354406400,
1356739200,
1359072000,
1361404800,
1363737600,
1366070400,
1368403200,
1370736000,
1373068800,
1375401600,
1377734400,
1380067200,
1382400000,
1384732800,
1387065600,
1389398400,
1391731200,
1394064000,
1396396800,
1398729600,
1401062400,
1403395200,
1405728000,
1408060800,
1410393600,
1412726400,
1415059200,
1417392000,
1419724800,
1422057600,
1424390400,
1426723200,
1429056000,
1431388800,
1433721600,
1436054400,
1438387200,
1440720000,
1443052800,
1445385600,
1447718400,
1450051200,
1452384000,
1454716800,
1457049600,
1459382400,
1461715200,
1464048000,
1466380800,
1468713600,
1471046400,
1473379200,
1475712000,
1478044800,
1480377600,
1482710400,
1485043200,
1487376000,
1489708800,
1492041600,
1494417600,
1496750400,
1499083200,
1501416000,
1503748800,
1506081600,
1508414400,
1510747200,
1513080000,
1515412800,
1517745600,
1520078400,
1522411200,
1524744000,
1527076800,
1529409600,
1531915200,
1534075200,
1536408000,
1538740800,
1541073600,
1543406400,
1545739200,
1548979200
);




void DrawGalpropRatio(TVirtualPad *c){
	c->cd();
	float x,y=0;
        int j=0;

	TGraph* galprop3DP=new TGraph();
        TGraph* galprop3DP2=new TGraph();
        {
                string filename="./Galprop/Tom/dP_1500.dat";
                cout<<filename<<endl;
                ifstream fp(filename.c_str());
                while (!fp.eof()){
                        fp>>x>>y;
                        if(x/1e3>3&&x/1e3<=100)
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
                        if(x/1e3>3&&x/1e3<=100)
                                galprop3DP2->SetPoint(j,x/(0.5*1e3),y);
                        j++;
                }
        }
        galprop3DP->GetXaxis()->SetRangeUser(0.1,20);

        galprop3DP->SetTitle("Deutons Flux: Geo. Zones");
        galprop3DP->GetXaxis()->SetTitle("Kin.En./nucl. [GeV/nucl.]");
        galprop3DP ->GetYaxis()->SetTitle("Flux [(m^2 sr GeV/nucl.)^-1]");
        galprop3DP ->GetXaxis()->SetTitleSize(0.045);

	galprop3DP->Draw("sameC");
        galprop3DP2->Draw("sameC");


}

void DrawDPRatioEkin(FileSaver Plots, std::vector<FileSaver> Files ){
	TStyle* m_gStyle= new TStyle();;
	m_gStyle->SetPalette(55);
	int nColors = m_gStyle->GetNumberOfColors();
	
	string filename4="./database_PD.root";
        TFile * file4 = TFile::Open(filename4.c_str(),"READ");

	std::vector<TGraphAsymmErrors *> PD_Graphs;

	TList *ExperimentsPD = file4->GetListOfKeys();
        TIter nextPD(ExperimentsPD);
        TKey * keyPD;
	TObject * obj;
        
	while((keyPD = (TKey*)nextPD())){
                obj = file4->Get(keyPD->GetName());
                if(obj->InheritsFrom("TGraphAsymmErrors")) PD_Graphs.push_back((TGraphAsymmErrors *)obj);
        }


	TCanvas *c3 = new TCanvas("D/P ratio Ekin");
	c3->cd();	
	std::vector<TH1F *> Ratios;

	for(int i=start;i<Files.size();i++){
		TFile * f = Files[i].GetFile();
		Ratios.push_back( (TH1F*)f->Get("Fluxes/Ratio_Ekin"));
	}
	TimeResult Ratio = GetTimeMean(Ratios);
	Ratio.Mean->Smooth();

	for(int i=start;i<Files.size();i++){
		PlotTH1FintoGraph(gPad,Global.GetGlobalDBins(),Ratios[i],"Kinetic Energy [GeV/nucl.]", "D/p ratio",1,true,"Psame",0.05,12,0.001,0.1,"",8);
	}


//PlotTH1FintoGraph(gPad,Global.GetGlobalDBins(),Ratio.Mean,"Kinetic Energy [GeV/nucl.]", "D/p ratio",2,true,"Psame",0.05,12,0.001,0.1,"",8);
	PlotTH1FintoGraph(gPad,Global.GetGlobalDBins(),Ratio.Mean,"Kinetic Energy [GeV/n]","D/p ratio", 2,true,"e4Psame", 0.05,12,0.001,0.1,"",8,true,false,Ratio.VariabilityDw,Ratio.VariabilityUp);


        TLegend * leg = (TLegend*) gPad->FindObject("leg");

        for(uint n=0;n<PD_Graphs.size();n++){
                PD_Graphs[n] ->Draw("Psame");
                PD_Graphs[n]->SetMarkerSize(2);
                leg->AddEntry(PD_Graphs[n],PD_Graphs[n]->GetTitle(),"p");
        }

	DrawGalpropRatio(c3);

        Plots.Add(c3);
        Plots.writeObjsInFolder("Fluxes");
	


}


void DrawPFluxEkin(FileSaver Plots, std::vector<FileSaver> Files ){
	TStyle* m_gStyle= new TStyle();;
	m_gStyle->SetPalette(55);

	string filename2="./database_P.root";
	TFile * file2 = TFile::Open(filename2.c_str(),"READ");

	std::vector<TGraphAsymmErrors *> P_Graphs;
	 TList *ExperimentsP = file2->GetListOfKeys();
        TIter nextP(ExperimentsP);
        TKey * keyP;
	TObject * obj;

        while((keyP = (TKey*)nextP())){
                obj = file2->Get(keyP->GetName());
                if(obj->InheritsFrom("TGraphAsymmErrors")) P_Graphs.push_back((TGraphAsymmErrors *)obj);
        }



	TCanvas *c3 = new TCanvas("P Flux Ekin");
	c3->cd();	
	gPad->SetLogx();
	gPad->SetLogy();
	gPad->SetGridx();
	gPad->SetGridy();
	
	std::vector<TH1F *> Ratios;

	TH2F * Frame = CreateFrame(gPad,0.01,60,0.01,4000,"Ekin [GeV/nucl.]","Flux [(m^2 sr GeV/nucl.)^{-1}]");	

	Frame->Draw();	
	for(int i=start;i<Files.size();i++){
		TFile * f = Files[i].GetFile();
		Ratios.push_back( (TH1F*)f->Get("Fluxes/PFluxQHE/PFluxQHE_Flux"));
	}
	TimeResult Ratio = GetTimeMean(Ratios);

	TLegend * leg =new TLegend(0.8, 0.1,0.95,0.95);
	leg->SetName("leg");
        for(uint n=0;n<P_Graphs.size();n++){
                P_Graphs[n] ->Draw("Psame");
                P_Graphs[n]->SetMarkerSize(2); 
                leg->AddEntry(P_Graphs[n],P_Graphs[n]->GetTitle(),"p");
        }

	leg->SetFillColor(0);
	leg->SetLineWidth(2);
	leg->Draw("same");
	     
  	Ratio.Scale(0.84);
	PlotTH1FintoGraph(gPad,PRB,Ratio.Mean,"Kinetic Energy [GeV/nucl.]", "D/p ratio",2,true,"Psame",0.1,60,0.001,4000,"",8);
	PlotTH1FintoGraph(gPad,PRB,Ratio.Mean,"Kinetic Energy [GeV/n]","D/p ratio", 2,true,"e4Psame", 0.1,60,0.01,4000,"",8,true,false,Ratio.VariabilityUp,Ratio.VariabilityDw);

        Plots.Add(c3);
        Plots.writeObjsInFolder("Fluxes");
	


}


void DrawDFluxEkin(FileSaver Plots, std::vector<FileSaver> Files ){
	TStyle* m_gStyle= new TStyle();;
	m_gStyle->SetPalette(55);

	string filename2="./database_D.root";
	TFile * file2 = TFile::Open(filename2.c_str(),"READ");

	std::vector<TGraphAsymmErrors *> D_Graphs;
	 TList *ExperimentsD = file2->GetListOfKeys();
        TIter nextD(ExperimentsD);
        TKey * keyD;
	TObject * obj;

        while((keyD = (TKey*)nextD())){
                obj = file2->Get(keyD->GetName());
                if(obj->InheritsFrom("TGraphAsymmErrors")) D_Graphs.push_back((TGraphAsymmErrors *)obj);
        }



	TCanvas *c3 = new TCanvas("D Flux Ekin");
	c3->cd();	
	gPad->SetLogx();
	gPad->SetLogy();
	gPad->SetGridx();
	gPad->SetGridy();
	
	std::vector<TH1F *> Ratios;

	TH2F * Frame = CreateFrame(gPad,0.01,30,0.01,300,"Ekin [GeV/nucl.]","Flux [(m^2 sr GeV/nucl.)^{-1}]");	

	Frame->Draw();	
	for(int i=start;i<Files.size();i++){
		TFile * f = Files[i].GetFile();
		Ratios.push_back( (TH1F*)f->Get("Fluxes/MergedRange_D_Ekin"));
	}
	TimeResult Ratio = GetTimeMean(Ratios);

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
	     
//  	Ratio.Scale(2.38);
	PlotTH1FintoGraph(gPad,Global.GetGlobalDBins(),Ratio.Mean,"Kinetic Energy [GeV/nucl.]", "D/p ratio",4,true,"Psame",0.01,30,0.01,300,"",8);
	PlotTH1FintoGraph(gPad,Global.GetGlobalDBins(),Ratio.Mean,"Kinetic Energy [GeV/n]","D/p ratio", 4,true,"e4Psame", 0.01,30,0.01,300,"",8,true,false,Ratio.VariabilityUp,Ratio.VariabilityDw);

        Plots.Add(c3);
        Plots.writeObjsInFolder("Fluxes");
	


}


void DrawLatWeights(FileSaver Plots, std::vector<FileSaver> Files ){

	TStyle* z_gStyle= new TStyle();
	z_gStyle->SetPalette(55);
	int nColors = z_gStyle->GetNumberOfColors();


	std::vector<TF1*> Ratios;
	for(int i=start;i<Files.size();i++){
		TFile * f = Files[i].GetFile();
		Ratios.push_back( (TF1*)f->Get("Spectra/weightmodel"));
	}

	TCanvas *c3 = new TCanvas("Weights time dep");
	c3->cd();	
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	TH2F * Frame = CreateFrame(gPad,0,150,0.01,20,"R [GV]","Event Weight");	
	Frame->Draw();

	for(int i=start;i<Files.size();i++){
		Ratios[i]->SetLineColor(z_gStyle->GetColorPalette((float)nColors/Files.size()*(i+1) ));
		Ratios[i]->SetLineWidth(3);
		Ratios[i]->Draw("same");
	}
	Plots.Add(c3);
        Plots.writeObjsInFolder("LatWeights");
	
}

	
void PlotTimeDep(TVirtualPad * c, float D[], float D_err[], float P[], float P_err[], std::vector <int> Times, std::string title) {

	TGraphErrors * TimeDepD = new TGraphErrors();
	TGraphErrors * TimeDepP = new TGraphErrors();


	for(int i=start;i<Times.size();i++){
		TimeDepD->SetPoint(i-start,Times[i],D[i]);
		TimeDepD->SetPointError(i-start,0,D_err[i]);
		TimeDepP->SetPoint(i-start,Times[i],P[i]);
		TimeDepP->SetPointError(i-start,0,P_err[i]);
	
	}

	TimeDepD->SetMarkerStyle(8);
	TimeDepD->SetMarkerSize(1.3);
	TimeDepD->SetMarkerColor(4);
	TimeDepD->SetLineColor(4);
	TimeDepD->SetLineWidth(2);

	TimeDepP->SetMarkerStyle(8);
	TimeDepP->SetMarkerSize(1.3);
	TimeDepP->SetMarkerColor(2);
	TimeDepP->SetLineColor(2);
	TimeDepP->SetLineWidth(4);
	TimeDepP->GetYaxis()->SetLabelColor(kRed);
	TimeDepP->GetYaxis()->SetLabelSize(0.035*gPad->GetCanvas()->GetWh()/(float)gPad->YtoPixel(gPad->GetY1()));
	TimeDepP->GetYaxis()->SetLabelFont(21);
	TimeDepP->GetYaxis()->SetAxisColor(2);
	TimeDepP->GetYaxis()->SetRangeUser( 0.2*(TimeDepP->GetHistogram()->GetMinimum() + TimeDepP->GetHistogram()->GetMaximum())/2,2*(TimeDepP->GetHistogram()->GetMinimum() + TimeDepP->GetHistogram()->GetMaximum())/2);


	TimeDepP->GetXaxis()->	SetTimeDisplay(1);
	TimeDepP->GetXaxis()->  SetTimeOffset(0,"gmt");  
	TimeDepP->GetXaxis()->  SetTimeFormat("%b%y");
//	TimeDepP->GetXaxis()->SetLabelSize(0.035*gPad->GetCanvas()->GetWh()/(float)gPad->YtoPixel(gPad->GetY1()));
//	TimeDepP->GetXaxis()->SetLabelFont(21);
	
	TimeDepP->SetTitle(title.c_str());
	gStyle->SetTitleFontSize(0.2);
	TimeDepP->GetXaxis()->  SetLabelSize(0.075);
	

	c->cd();
	TPad *p1 = new TPad("p1", "", 0, 0, 1, 1);
	p1->SetGrid();
	TPad *p2 = new TPad("p2", "", 0, 0, 1, 1);
	p2->SetFillStyle(4000); // will be transparent
	 p2->SetGrid();
	
	p1->Draw();
  	p1->cd();
	TimeDepP->Draw("AP");
  	gPad->Update();


	Double_t xmin = p1->GetUxmin();
	Double_t xmax = p1->GetUxmax();
	Double_t ymin = 0.2*(TimeDepP->GetHistogram()->GetMinimum() + TimeDepP->GetHistogram()->GetMaximum())/2*D[5]/P[5];
	Double_t ymax = 2*(TimeDepP->GetHistogram()->GetMinimum() + TimeDepP->GetHistogram()->GetMaximum())/2*D[5]/P[5];
	Double_t dx = (xmax - xmin) / 0.8; // 10 percent margins left and right
	Double_t dy = (ymax - ymin) / 0.8; // 10 percent margins top and bottom
	
	p2->Range(xmin-0.1*dx, ymin-0.1*dy, xmax+0.1*dx, ymax+0.1*dy);
	p2->Draw();
	p2->cd();
	TimeDepD->GetYaxis()->SetRangeUser(ymin,ymax);
	TimeDepD->Draw("P");
	gPad->Update();

	TGaxis *axis = new TGaxis(xmax, ymin, xmax, ymax, ymin, ymax, 510, "+L");
	axis->SetLineColor(kBlue);
	axis->SetLabelColor(kBlue);
	axis->SetLabelSize(0.05);
	
	axis->Draw("same");
	gPad->Update();

}

void PlotTimeRatio(TVirtualPad * c, float D[], float D_err[], float P[], float P_err[], std::vector <int> Times, int col, std::string title) {
	TGraphErrors * TimeDepP = (TGraphErrors * ) gPad->FindObject("Time");
	TGraphErrors * TimeDepR = new TGraphErrors();
	
	for(int i=start;i<Times.size();i++){
		TimeDepR->SetPoint(i-start,Times[i],D[i]/P[i]);
		TimeDepR->SetPointError(i-start,0,pow(pow(D_err[i]/D[i],2) + pow (P_err[i]/P[i],2),0.5)*D[i]/(P[i]) );
	
	}
	TimeDepR->SetName(title.c_str());
	TimeDepR->SetTitle(title.c_str());
	gStyle->SetTitleFontSize(0.2);

	TimeDepR->SetMarkerStyle(8);
	TimeDepR->SetMarkerSize(1.5);
	TimeDepR->SetMarkerColor(col);
	TimeDepR->SetLineColor(col);
	TimeDepR->SetLineWidth(3);
	TimeDepR->GetYaxis()->SetRangeUser(0,0.035);
	TimeDepR->GetXaxis()->	SetTimeDisplay(1);
	TimeDepR->GetXaxis()->  SetTimeOffset(0,"gmt");  
	TimeDepR->GetXaxis()->  SetTimeFormat("%b%y");
	TimeDepR->GetXaxis()->  SetLabelSize(0.075);
	TimeDepR->GetYaxis()->SetLabelColor(1);
	TimeDepR->GetYaxis()->SetLabelSize(0.015*gPad->GetCanvas()->GetWh()/(float)gPad->YtoPixel(gPad->GetY1()));
	TimeDepR->GetYaxis()->SetLabelFont(11);
	TF1 * TimeModel = new TF1("timemodel",FitTime,0,1723958400,2);
	TimeDepR->Fit("timemodel");	

	c->cd();
	if(!TimeDepP)	TimeDepR->Draw("AP");
	else      TimeDepR->Draw("Psame");
	
	gPad->Update();
}


void DrawPDFluxRatio(FileSaver Plots, std::vector<FileSaver> Files,int binselected[] ){

	TStyle* f_gStyle= new TStyle();
	f_gStyle->SetPalette(55);
	int nColors = f_gStyle->GetNumberOfColors();

	string filename4="./database_PD.root";
        TFile * file4 = TFile::Open(filename4.c_str(),"READ");

	std::vector<TGraphAsymmErrors *> PD_Graphs;

	TList *ExperimentsPD = file4->GetListOfKeys();
        TIter nextPD(ExperimentsPD);
        TKey * keyPD;
	TObject * obj;
        
	while((keyPD = (TKey*)nextPD())){
                obj = file4->Get(keyPD->GetName());
                if(obj->InheritsFrom("TGraphAsymmErrors")) PD_Graphs.push_back((TGraphAsymmErrors *)obj);
        }


	std::vector<int> Times;
	std::vector<std::string> Dates;
	std::vector<Flux *> FluxesPTOF;
	std::vector<Flux *> FluxesDTOF;
	std::vector<Flux *> FluxesPNaF;
	std::vector<Flux *> FluxesDNaF;
	std::vector<Flux *> FluxesPAgl;
	std::vector<Flux *> FluxesDAgl;

	for(int i=0;i<Files.size();i++){

		Flux * DFluxTOF = new Flux(Files[i], "DFluxTOF", "FullsetEff_D_TOF","FullsetEfficiency","TOFfits/Fit Results/Primary Deuteron Counts","ExposureTOF",Global.GetToFDBins());
		Flux * DFluxNaF = new Flux(Files[i], "DFluxNaF", "FullsetEff_D_NaF","FullsetEfficiency","NaFfits/Fit Results/Primary Deuteron Counts","ExposureNaF",Global.GetNaFDBins());
		Flux * DFluxAgl = new Flux(Files[i], "DFluxAgl", "FullsetEff_D_Agl","FullsetEfficiency","Aglfits/Fit Results/Primary Deuteron Counts","ExposureAgl",Global.GetAglDBins());

		Flux * PFluxTOF = new Flux(Files[i], "PFluxTOF", "FullsetEff_P_TOF","FullsetEfficiency","TOFfits/Fit Results/Primary Proton Counts","ExposureTOF",Global.GetToFPBins());
		Flux * PFluxNaF = new Flux(Files[i], "PFluxNaF", "FullsetEff_P_NaF","FullsetEfficiency","NaFfits/Fit Results/Primary Proton Counts","ExposureNaF",Global.GetNaFPBins());
		Flux * PFluxAgl = new Flux(Files[i], "PFluxAgl", "FullsetEff_P_Agl","FullsetEfficiency","Aglfits/Fit Results/Primary Proton Counts","ExposureAgl",Global.GetAglPBins());


		FluxesPTOF.push_back(PFluxTOF);
		FluxesPNaF.push_back(PFluxNaF);
		FluxesPAgl.push_back(PFluxAgl);
	
		FluxesDTOF.push_back(DFluxTOF);
		FluxesDNaF.push_back(DFluxNaF);
		FluxesDAgl.push_back(DFluxAgl);


	}

	for(int i=0;i<Files.size();i++){
		Times.push_back(std::atoi(GroupedFiles[i].substr(GroupedFiles[i].find("-")+1,10).c_str()) );
		std::uint32_t time_date_stamp = Times[i];
		std::time_t temp = time_date_stamp;
		std::tm* t = std::gmtime(&temp);
		char buffer [10];
		sprintf(buffer,"%d/%d", t->tm_mon+1, t->tm_year-100);
		Dates.push_back(buffer);
		std::cout<<Times[i]<<" "<<Dates[i]<<endl;
	}
	
	TCanvas *c3 = new TCanvas("D/P FluxRatio");
	c3->SetCanvasSize(2000,1500);
	gPad->SetLogx();


	TH1F * DErrTOF[Files.size()];
	TH1F * DErrNaF[Files.size()];
	TH1F * DErrAgl[Files.size()];
	TH1F * Merged_DErr[Files.size()];
	
	TFile * file[Files.size()] ;
		
	for(int i=start;i<Files.size();i++){
		file [i]=  Files[i].GetFile();
		DErrTOF[i]=(TH1F*) file[i]->Get("TOFDfits/Fit Results/StatErrorD");
		DErrNaF[i]=(TH1F*) file[i]->Get("NaFDfits/Fit Results/StatErrorD");
		DErrAgl[i]=(TH1F*) file[i]->Get("AglDfits/Fit Results/StatErrorD");
		Merged_DErr[i] = Global.MergeSubDResult_D(DErrTOF[i],DErrNaF[i],DErrAgl[i]);
	}

	TH1F * PErrTOF[Files.size()];
	TH1F * PErrNaF[Files.size()];
	TH1F * PErrAgl[Files.size()];
	TH1F * Merged_PErr[Files.size()];
	
	for(int i=start;i<Files.size();i++){
		PErrTOF[i]=(TH1F*) file[i]->Get("TOFPfits/Fit Results/StatErrorP"); ;
		PErrNaF[i]=(TH1F*) file[i]->Get("NaFPfits/Fit Results/StatErrorP"); ;
		PErrAgl[i]=(TH1F*) file[i]->Get("AglPfits/Fit Results/StatErrorP"); ;
		Merged_PErr[i] = Global.MergeSubDResult_P(PErrTOF[i],PErrNaF[i],PErrAgl[i]);
	}


	TH1F * DCountsTOF[Files.size()];
	TH1F * DCountsNaF[Files.size()];
	TH1F * DCountsAgl[Files.size()];

	for(int i=start;i<Files.size();i++){
		DCountsTOF[i]=(TH1F*) FluxesDTOF[i]->GetFlux();
		DCountsNaF[i]=(TH1F*) FluxesDNaF[i]->GetFlux();
		DCountsAgl[i]=(TH1F*) FluxesDAgl[i]->GetFlux();
	}

	TH1F * PCountsTOF[Files.size()];
	TH1F * PCountsNaF[Files.size()];
	TH1F * PCountsAgl[Files.size()];

	for(int i=start;i<Files.size();i++){
		PCountsTOF[i]=(TH1F*)FluxesPTOF[i]->GetFlux();
		PCountsNaF[i]=(TH1F*)FluxesPNaF[i]->GetFlux();
		PCountsAgl[i]=(TH1F*)FluxesPAgl[i]->GetFlux();
	}


	TH1F * DCountsTOF_rig[Files.size()];
	TH1F * DCountsNaF_rig[Files.size()];
	TH1F * DCountsAgl_rig[Files.size()];
	TH1F * Merged_D[Files.size()];
	
	for(int i=start;i<Files.size();i++){
		DCountsTOF_rig[i]=(TH1F*) FluxesDTOF[i]->GetFlux_rig();
		DCountsNaF_rig[i]=(TH1F*) FluxesDNaF[i]->GetFlux_rig();
		DCountsAgl_rig[i]=(TH1F*) FluxesDAgl[i]->GetFlux_rig();
		cout<<DCountsTOF_rig[i]<<" "<<DCountsNaF_rig[i]<<" "<<DCountsAgl_rig[i]<<endl;
		Merged_D[i] = Global.MergeSubDResult_D(DCountsTOF_rig[i],DCountsNaF_rig[i],DCountsAgl_rig[i]);
	}


	TH1F * PCountsTOF_rig[Files.size()];
	TH1F * PCountsNaF_rig[Files.size()];
	TH1F * PCountsAgl_rig[Files.size()];

	TH1F * Merged_P[Files.size()];
	
	for(int i=start;i<Files.size();i++){
		PCountsTOF_rig[i]=(TH1F*)FluxesPTOF[i]->GetFlux_rig();
		PCountsNaF_rig[i]=(TH1F*)FluxesPNaF[i]->GetFlux_rig();
		PCountsAgl_rig[i]=(TH1F*)FluxesPAgl[i]->GetFlux_rig();
		Merged_P[i] = Global.MergeSubDResult_P(PCountsTOF_rig[i],PCountsNaF_rig[i],PCountsAgl_rig[i]);
	}

	int bin_P[6];
	int bin_D[6];

	for(int j=0;j<6;j++){
		float Rtoicenter = Global.GetGlobalDBins().RigTOIBinsCent()[binselected[j]];
		bin_P[j] = Global.GetGlobalPBins().GetRTOIBin(Rtoicenter);
		bin_D[j] = Global.GetGlobalDBins().GetRTOIBin(Rtoicenter);
	
	}


	float Bin_D[6][Files.size()];
	float Bin_P[6][Files.size()];
	float Bin_D_err[6][Files.size()];
	float Bin_P_err[6][Files.size()];

	TH1F * DERRAverage = (TH1F*)Merged_DErr[2]->Clone();
	TH1F * PERRAverage = (TH1F*)Merged_PErr[2]->Clone();

	for(int i=start;i<Files.size();i++){
		DERRAverage->Add( Merged_DErr[i]);
		PERRAverage->Add( Merged_PErr[i]);
	}
	DERRAverage->Scale(1./Files.size());
	PERRAverage->Scale(1./Files.size());


	for(int j=0;j<6;j++){
		for(int i=start;i<Files.size();i++){
			Bin_D[j][i]= Merged_D[i]->GetBinContent(bin_D[j]+1);
			Bin_P[j][i]= Merged_P[i]->GetBinContent(bin_P[j]+1);
		//	float errD = Merged_DErr[i]->GetBinContent(bin_D[j]+1)/DERRAverage->GetBinContent(bin_D[j]+1);
		//	float errP = Merged_PErr[i]->GetBinContent(bin_P[j]+1)/PERRAverage->GetBinContent(bin_P[j]+1);
		
			float errD = Merged_D[i]->GetBinError(bin_D[j]+1);
			float errP = Merged_P[i]->GetBinError(bin_P[j]+1);
			
			float errsubd=0.01;
			if(Global.GetToFBinD(bin_D[j])>=0) errsubd=0.01;
			else if (Global.GetNaFBinD(bin_D[j])>=0) errsubd=0.05;
			else if (Global.GetAglBinD(bin_D[j])>=0) errsubd=0.022;
			//Bin_D_err[j][i]= errD*errsubd*Bin_D[j][i];
			//Bin_P_err[j][i]= errP*errsubd*Bin_P[j][i];
			Bin_D_err[j][i]= errD;
			Bin_P_err[j][i]= errP;
		
			cout<<"******************ERROR************** "<<i <<" "<<j<<endl;
			cout<<errD<<endl;//Bin_D_err[j][i]/Bin_D[j][i]<<endl;
			cout<<errP<<endl;//Bin_P_err[j][i]/Bin_P[j][i]<<endl;
		}
	}


	TCanvas *c = new TCanvas(("Time Dep." +to_string(binselected[0])+"_"+to_string(binselected[5])).c_str());
	c->SetCanvasSize(1000,700);
	c->Divide(1,6,0,0);
	c->cd(1);	
	PlotTimeDep(gPad,Bin_D[0],Bin_D_err[0],Bin_P[0],Bin_P_err[0], Times, (Convert(Global.GetGlobalPBins().RigTOIBins()[bin_P[0]])  + " < R < " + Convert(Global.GetGlobalPBins().RigTOIBins()[bin_P[0]+1]) + "GV").c_str());
	c->cd(2);                                      
       	PlotTimeDep(gPad,Bin_D[1],Bin_D_err[1],Bin_P[1],Bin_P_err[1], Times,(Convert(Global.GetGlobalPBins().RigTOIBins()[bin_P[1]])  + " < R < " + Convert(Global.GetGlobalPBins().RigTOIBins()[bin_P[1]+1]) + "GV").c_str());
	c->cd(3);                                      
       	PlotTimeDep(gPad,Bin_D[2],Bin_D_err[2],Bin_P[2],Bin_P_err[2], Times,(Convert(Global.GetGlobalPBins().RigTOIBins()[bin_P[2]])  + " < R < " + Convert(Global.GetGlobalPBins().RigTOIBins()[bin_P[2]+1]) + "GV").c_str());
	c->cd(4);	                               
 	PlotTimeDep(gPad,Bin_D[3],Bin_D_err[3],Bin_P[3],Bin_P_err[3], Times,(Convert(Global.GetGlobalPBins().RigTOIBins()[bin_P[3]])  + " < R < " + Convert(Global.GetGlobalPBins().RigTOIBins()[bin_P[3]+1]) + "GV").c_str());
	c->cd(5);                                      
 	PlotTimeDep(gPad,Bin_D[4],Bin_D_err[4],Bin_P[4],Bin_P_err[4], Times,(Convert(Global.GetGlobalPBins().RigTOIBins()[bin_P[4]])  + " < R < " + Convert(Global.GetGlobalPBins().RigTOIBins()[bin_P[4]+1]) + "GV").c_str());
	c->cd(6);                                      
 	PlotTimeDep(gPad,Bin_D[5],Bin_D_err[5],Bin_P[5],Bin_P_err[5], Times,(Convert(Global.GetGlobalPBins().RigTOIBins()[bin_P[5]])  + " < R < " + Convert(Global.GetGlobalPBins().RigTOIBins()[bin_P[5]+1]) + "GV").c_str());
	
	Plots.Add(c);
        Plots.writeObjsInFolder("Flux");

	TCanvas *c1 = new TCanvas(("Time Ratio." +to_string(binselected[0])+"_"+to_string(binselected[5])).c_str());
	c1->SetCanvasSize(1000,700);
	c1->cd();
	c1->Divide(1,6,0,0);
	c1->cd(1);	
	gPad->SetTickx();
	gPad->SetTicky();
	PlotTimeRatio(gPad,Bin_D[0],Bin_D_err[0],Bin_P[0],Bin_P_err[0], Times,2,(Convert(Global.GetGlobalPBins().RigTOIBins()[bin_P[0]])  + " < R < " + Convert(Global.GetGlobalPBins().RigTOIBins()[bin_P[0]+1]) + "GV").c_str());
	c1->cd(2);                                      
       	gPad->SetTickx();
	gPad->SetTicky();
	PlotTimeRatio(gPad,Bin_D[1],Bin_D_err[1],Bin_P[1],Bin_P_err[1], Times,2,(Convert(Global.GetGlobalPBins().RigTOIBins()[bin_P[1]])  + " < R < " + Convert(Global.GetGlobalPBins().RigTOIBins()[bin_P[1]+1]) + "GV").c_str());
	c1->cd(3);                                      
       	gPad->SetTickx();
	gPad->SetTicky();
	PlotTimeRatio(gPad,Bin_D[2],Bin_D_err[2],Bin_P[2],Bin_P_err[2], Times,2,(Convert(Global.GetGlobalPBins().RigTOIBins()[bin_P[2]])  + " < R < " + Convert(Global.GetGlobalPBins().RigTOIBins()[bin_P[2]+1]) + "GV").c_str());
	c1->cd(4);	                               
 	gPad->SetTickx();
	gPad->SetTicky();
	PlotTimeRatio(gPad,Bin_D[3],Bin_D_err[3],Bin_P[3],Bin_P_err[3], Times,2,(Convert(Global.GetGlobalPBins().RigTOIBins()[bin_P[3]])  + " < R < " + Convert(Global.GetGlobalPBins().RigTOIBins()[bin_P[3]+1]) + "GV").c_str());
	c1->cd(5);                                      
 	gPad->SetTickx();
	gPad->SetTicky();
	PlotTimeRatio(gPad,Bin_D[4],Bin_D_err[4],Bin_P[4],Bin_P_err[4], Times,2,(Convert(Global.GetGlobalPBins().RigTOIBins()[bin_P[4]])  + " < R < " + Convert(Global.GetGlobalPBins().RigTOIBins()[bin_P[4]+1]) + "GV").c_str());
	c1->cd(6);                                      
 	gPad->SetTickx();
	gPad->SetTicky();
	PlotTimeRatio(gPad,Bin_D[5],Bin_D_err[5],Bin_P[5],Bin_P_err[5], Times,2,(Convert(Global.GetGlobalPBins().RigTOIBins()[bin_P[5]])  + " < R < " + Convert(Global.GetGlobalPBins().RigTOIBins()[bin_P[5]+1]) + "GV").c_str());
	
	TCanvas * cerr = new TCanvas("cerr");
	cerr->cd();
	for(int i=start;i<Files.size();i++){
		Merged_PErr[i]->Draw("same");	
	}	
	Plots.Add(c1);
        Plots.writeObjsInFolder("Flux");

	SetUpRigTOIBinning();
	TCanvas *c5 = new TCanvas("D/P ratio R");
	c5->cd();	



	std::vector<TH1F *> Ratios;

	for(int i=start;i<Files.size();i++){
		Ratios.push_back(Global.MergedRatio(Merged_D[i],Merged_P[i]) );
	}
	TimeResult Ratio = GetTimeMean(Ratios);

	//Ratio.Mean->Smooth();

       	PlotTH1FintoGraph(gPad,GlobalRig.GetGlobalPBins(),Ratio.Mean,"Kinetic Energy [GeV/nucl.]", "D/p ratio",2,false,"Psame",0.5,20,0.001,0.1,"",8);
	PlotTH1FintoGraph(gPad,GlobalRig.GetGlobalPBins(),Ratio.Mean,"Kinetic Energy [GeV/n]","D/p ratio", 2,false,"e4Psame", 0.5,20,0.001,0.1,"",8,true,false,Ratio.VariabilityUp,Ratio.VariabilityDw);


	Plots.Add(c5);
        Plots.writeObjsInFolder("Fluxes");

	


	return;
}


void DrawParameters(FileSaver Plots, std::vector<FileSaver> Files){

	TH1F * ChiSquareTOF[Files.size()];
	TH1F * ChiSquareNaF[Files.size()];
	TH1F * ChiSquareAgl[Files.size()];
	float bins[5]={2,5,9,11,13};

	for(int i=0;i<Files.size();i++) {
		ChiSquareTOF[i] = (TH1F*) Files[i].Get("TOFDfits/Fit Results/Best ChiSquare");
		ChiSquareNaF[i] = (TH1F*) Files[i].Get("NaFDfits/Fit Results/Best ChiSquare");
		ChiSquareAgl[i] = (TH1F*) Files[i].Get("AglDfits/Fit Results/Best ChiSquare");
	}	
	TH1F * TimeChiSquareTOF[5]; 
	TH1F * TimeChiSquareNaF[5]; 
	TH1F * TimeChiSquareAgl[5]; 

	for(int j=0;j<5;j++){
		TimeChiSquareTOF[j] = new TH1F(("TimeChiSquareTOF_"+to_string(j)).c_str(),("TimeChiSquareTOF_"+to_string(j)).c_str(),Files.size(),0,Files.size());
		TimeChiSquareNaF[j] = new TH1F(("TimeChiSquareNaF_"+to_string(j)).c_str(),("TimeChiSquareNaF_"+to_string(j)).c_str(),Files.size(),0,Files.size());
		TimeChiSquareAgl[j] = new TH1F(("TimeChiSquareAgl_"+to_string(j)).c_str(),("TimeChiSquareAgl_"+to_string(j)).c_str(),Files.size(),0,Files.size());
	}
	for(int j=0;j<5;j++)
		for(int i=0;i<Files.size();i++) {
			TimeChiSquareTOF[j]->SetBinContent(i+1,ChiSquareTOF[i]->GetBinContent(bins[j]+1));
			TimeChiSquareTOF[j]->SetBinError(i+1,0.2);
			TimeChiSquareNaF[j]->SetBinContent(i+1,ChiSquareNaF[i]->GetBinContent(bins[j]+1));
			TimeChiSquareNaF[j]->SetBinError(i+1,0.2);
			TimeChiSquareAgl[j]->SetBinContent(i+1,ChiSquareAgl[i]->GetBinContent(bins[j]+1));
			TimeChiSquareAgl[j]->SetBinError(i+1,0.2);
		}	


		
	TCanvas * c1 = new TCanvas("Chi2 TOF");
	for(int j=4;j>=0;j--)
		TimeChiSquareTOF[j]->Draw("same");
	TCanvas * c2 = new TCanvas("Chi2 NaF");
	for(int j=4;j>=0;j--)
		TimeChiSquareNaF[j]->Draw("same");
	TCanvas * c3 = new TCanvas("Chi2 Agl");
	for(int j=4;j>=0;j--)
		TimeChiSquareAgl[j]->Draw("same");
	
	Plots.Add(c1);
	Plots.Add(c2);
	Plots.Add(c3);
	Plots.writeObjsInFolder("Parameters/FitChi2");
}




void DrawMassDistributions(FileSaver Plots, std::vector<FileSaver> Files,int binselected[]){

	TStyle* l_gStyle= new TStyle();
	l_gStyle->SetPalette(55);
	int nColors = l_gStyle->GetNumberOfColors();

	TH1F * MassDistrD[Files.size()][6];
	TH1F * MassDistrP[Files.size()][6];

	int bin_P[6];
	int bin_D[6];

	for(int j=0;j<6;j++){
		float Rtoicenter = Global.GetGlobalDBins().RigTOIBinsCent()[binselected[j]];
		bin_P[j] = Global.GetGlobalPBins().GetRTOIBin(Rtoicenter);
		bin_D[j] = Global.GetGlobalDBins().GetRTOIBin(Rtoicenter);
	
	}



	for(int i=0;i<Files.size();i++) {
		for(int j=0;j<6;j++)	{

			if(Global.GetToFBinD(bin_D[j])>=0)
				MassDistrD[i][j] = (TH1F*) Files[i].Get(("TOFDfits/Fit Results/Data/Bin"+to_string(Global.GetToFBinD(bin_D[j]))+"/TOFDfits_Data_" + to_string(Global.GetToFBinD(bin_D[j]))+" 0 5").c_str()  );

			else if(Global.GetNaFBinD(bin_D[j])>=0)
				MassDistrD[i][j] = (TH1F*) Files[i].Get(("NaFDfits/Fit Results/Data/Bin"+to_string(Global.GetNaFBinD(bin_D[j]))+"/NaFDfits_Data_" + to_string(Global.GetNaFBinD(bin_D[j]))+" 0 5").c_str()  );

			else if(Global.GetAglBinD(bin_D[j])>=0)
				MassDistrD[i][j] = (TH1F*) Files[i].Get(("AglDfits/Fit Results/Data/Bin"+to_string(Global.GetAglBinD(bin_D[j]))+"/AglDfits_Data_" + to_string(Global.GetAglBinD(bin_D[j]))+" 0 5").c_str()  );

		}
			
	}

	for(int i=0;i<Files.size();i++) {
		for(int j=0;j<6;j++)	{

			if(Global.GetToFBinP(bin_P[j])>=0)
				MassDistrP[i][j] = (TH1F*) Files[i].Get(("TOFPfits/Fit Results/Data/Bin"+to_string(Global.GetToFBinP(bin_P[j]))+"/TOFPfits_Data_" + to_string(Global.GetToFBinP(bin_P[j]))+" 0 5").c_str()  );

			else if(Global.GetNaFBinP(bin_P[j])>=0)
				MassDistrP[i][j] = (TH1F*) Files[i].Get(("NaFPfits/Fit Results/Data/Bin"+to_string(Global.GetNaFBinP(bin_P[j]))+"/NaFPfits_Data_" + to_string(Global.GetNaFBinP(bin_P[j]))+" 0 5").c_str()  );

			else if(Global.GetAglBinP(bin_P[j])>=0)
				MassDistrP[i][j] = (TH1F*) Files[i].Get(("AglPfits/Fit Results/Data/Bin"+to_string(Global.GetAglBinP(bin_P[j]))+"/AglPfits_Data_" + to_string(Global.GetAglBinP(bin_P[j]))+" 0 5").c_str()  );

		}
			
	}




	TCanvas * c[6];
	for(int j=0;j<6;j++) {
		c[j]= new TCanvas(("Bin D" + to_string(j)).c_str());
		float norm = MassDistrD[0][j]->GetBinContent(MassDistrD[0][j]->GetMaximumBin());
		for(int i=0;i<Files.size();i++){
			MassDistrD[i][j]->Scale(norm/MassDistrD[i][j]->GetBinContent(MassDistrD[i][j]->GetMaximumBin()));
			c[j]->cd();
			gPad->SetLogy();
			MassDistrD[i][j]->SetLineColor(l_gStyle->GetColorPalette( ((float)nColors/Files.size()) *(i+1)));
			MassDistrD[i][j]->SetMarkerColor(l_gStyle->GetColorPalette( ((float)nColors/Files.size()) *(i+1)));
			MassDistrD[i][j]->SetLineWidth(3);	
			MassDistrD[i][j]->Draw("same");	
		}
		Plots.Add(c[j]);
	}
	TCanvas * d[6];
	for(int j=0;j<6;j++) {
		d[j]= new TCanvas(("Bin P" + to_string(j)).c_str());
		float norm = MassDistrP[0][j]->GetBinContent(MassDistrP[0][j]->GetMaximumBin());
		for(int i=0;i<Files.size();i++){
			MassDistrP[i][j]->Scale(norm/MassDistrP[i][j]->GetBinContent(MassDistrP[i][j]->GetMaximumBin()));
			d[j]->cd();
			gPad->SetLogy();
			MassDistrP[i][j]->SetLineColor(l_gStyle->GetColorPalette( ((float)nColors/Files.size()) *(i+1)));
			MassDistrP[i][j]->SetMarkerColor(l_gStyle->GetColorPalette( ((float)nColors/Files.size()) *(i+1)));
			MassDistrP[i][j]->SetLineWidth(3);	
			MassDistrP[i][j]->Draw("same");	
		}
		Plots.Add(d[j]);
	}
	

	Plots.writeObjsInFolder(("Mass Distr." +to_string(binselected[0])+"_"+to_string(binselected[5])).c_str());
}



void DrawTimeEffCorr(FileSaver Plots, std::vector<FileSaver> Files ){

	TStyle* r_gStyle= new TStyle();
	r_gStyle->SetPalette(55);
	int nColors = r_gStyle->GetNumberOfColors();

	std::vector<int> Times;
	std::vector<std::string> Dates;
	std::vector<TSpline3 *> TriggerTSpline3_HE  ;
	std::vector<TSpline3 *> L1PickUpTSpline3_HE ;
	std::vector<TSpline3 *> GoodQTrack_HE      ;
	std::vector<TSpline3 *> GoodChi_HE  	     ;
	std::vector<TSpline3 *> TrackerTSpline3_HE  ;
	std::vector<TSpline3 *> StatusL1Check_HE   ;
	std::vector<TSpline3 *> Good1Track_HE      ;
	std::vector<TSpline3 *> GoodLtof_HE  	     ;
	std::vector<TSpline3 *> GoodUtof_HE  	     ;
	std::vector<TSpline3 *> GoodTime_TOF 	     ;
	std::vector<TSpline3 *> Quality_TOF 	     ;
	std::vector<TSpline3 *> RICHTSpline3_NaF    ;
	std::vector<TSpline3 *> RICHTSpline3_Agl    ;
	std::vector<TSpline3 *> RICHQualTSpline3_NaF;
	std::vector<TSpline3 *> RICHQualTSpline3_Agl;


	for(int i=0;i<Files.size();i++){
	TFile * file = Files[i].GetFile();
	TSpline3 *	riggerTSpline3_HE  	= (TSpline3 *)file->Get("Trigger Eff. Corr/TriggerEffCorr_HE/TriggerEffCorr_HE_CorrSpline");
	TSpline3 *	PickUpTSpline3_HE  	= (TSpline3 *)file->Get("L1PickUp Eff. Corr/L1PickUpEffCorr_HE/L1PickUpEffCorr_HE _CorrSpline");
	TSpline3 *	oodQTrack_HE   		= (TSpline3 *)file->Get("GoodQTrack Eff. Corr/GoodQTrackEffCorr_HE/GoodQTrackEffCorr_HE_CorrSpline");
	TSpline3 *	oodChi_HE   		= (TSpline3 *)file->Get("GoodChi Eff. Corr/GoodChiEffCorr_HE/GoodChiEffCorr_HE_CorrSpline");
	TSpline3 *	rackerTSpline3_HE  	= (TSpline3 *)file->Get("Tracker Eff. Corr/TrackerEffCorr_HE/TrackerEffCorr_HE_CorrSpline");
	TSpline3 *	tatusL1Check_HE  	= (TSpline3 *)file->Get("StatusL1Check Eff. Corr/StatusL1Check_HE/StatusL1Check_HE_CorrSpline");
	TSpline3 *	ood1Track_HE  		= (TSpline3 *)file->Get("Good1Track Eff. Corr/Good1TrackEffCorr_HE/Good1TrackEffCorr_HE_CorrSpline");
	TSpline3 *	oodLtof_HE   		= (TSpline3 *)file->Get("GoodLtof Eff. Corr/GoodLTOFEffCorr_HE/GoodLTOFEffCorr_HE_CorrSpline");
	TSpline3 *	oodUtof_HE   		= (TSpline3 *)file->Get("GoodUtof Eff. Corr/GoodUtofEffCorr_HE/GoodUtofEffCorr_HE_CorrSpline");
	TSpline3 *	oodTime_TOF  		= (TSpline3 *)file->Get("GoodTime Eff. Corr/GoodTimeEffCorr_TOF/GoodTimeEffCorr_TOF_CorrSpline");
	TSpline3 *	uality_TOF  		= (TSpline3 *)file->Get("Quality TOF Eff. Corr/QualityEffCorr_TOF/QualityEffCorr_TOF_CorrSpline");

	TSpline3 *	ICHTSpline3_NaF  	= (TSpline3 *)file->Get("RICH Eff. Corr/RICHCorrection_NaF/RICHCorrection_NaF_CorrSpline");
	TSpline3 *	ICHTSpline3_Agl  	= (TSpline3 *)file->Get("RICH Eff. Corr/RICHCorrection_Agl/RICHCorrection_Agl_CorrSpline");
	TSpline3 *	ICHQualTSpline3_NaF 	= (TSpline3 *)file->Get("RICH Qual Eff. Corr/RICHQualCorrection_NaF/RICHQualCorrection_NaF_CorrSpline");
	TSpline3 *	ICHQualTSpline3_Agl 	= (TSpline3 *)file->Get("RICH Qual. Eff. Corr/RICHqualCorrection_Agl/RICHqualCorrection_Agl_CorrSpline");

	TriggerTSpline3_HE .push_back(		riggerTSpline3_HE  ); 
	L1PickUpTSpline3_HE .push_back(          PickUpTSpline3_HE );
	GoodQTrack_HE      	.push_back(     oodQTrack_HE   	  );
	GoodChi_HE  	  .push_back(           oodChi_HE   	  );
	TrackerTSpline3_HE  .push_back(          rackerTSpline3_HE  );
	StatusL1Check_HE   .push_back(          tatusL1Check_HE   );
	Good1Track_HE      .push_back(          ood1Track_HE  	  );
	GoodLtof_HE  	  .push_back(           oodLtof_HE   	  );
	GoodUtof_HE  	  .push_back(           oodUtof_HE   	  );
	GoodTime_TOF 	  .push_back(           oodTime_TOF  	  );
	Quality_TOF       .push_back(           uality_TOF       );
	RICHTSpline3_NaF    .push_back(          ICHTSpline3_NaF    );
	RICHTSpline3_Agl    .push_back(          ICHTSpline3_Agl    );
	RICHQualTSpline3_NaF.push_back(          ICHQualTSpline3_NaF);
	RICHQualTSpline3_Agl.push_back(          ICHQualTSpline3_Agl);

	}


	TCanvas * f = new TCanvas("Track Quality");
	f->Divide(1,3);
	f->cd(1);
	for(int i=0;i<Files.size();i++) {
		GoodQTrack_HE[i]->SetLineWidth(3);
		GoodQTrack_HE[i]->SetLineColor(r_gStyle->GetColorPalette((float)nColors/Files.size()*(i+1) ));
		if(i==0) GoodQTrack_HE[i]->Draw();
		else GoodQTrack_HE[i]->Draw("same");
	}
	f->cd(2);
	for(int i=0;i<Files.size();i++) {
		GoodChi_HE[i]->SetLineWidth(3);
		GoodChi_HE[i]->SetLineColor(r_gStyle->GetColorPalette((float)nColors/Files.size()*(i+1) ));
		if(i==0) GoodChi_HE[i]->Draw();
		else GoodChi_HE[i]->Draw("same");
	}
	f->cd(3);
	for(int i=0;i<Files.size();i++) {
		Good1Track_HE[i]->SetLineWidth(3);
		Good1Track_HE[i]->SetLineColor(r_gStyle->GetColorPalette((float)nColors/Files.size()*(i+1) ));
		if(i==0) Good1Track_HE[i]->Draw();
		else Good1Track_HE[i]->Draw("same");
	}
	
	
	TCanvas * c = new TCanvas("TOF");
	c->Divide(1,2);
	c->cd(1);
	for(int i=0;i<Files.size();i++) {
		GoodUtof_HE[i]->SetLineWidth(3);
		GoodUtof_HE[i]->SetLineColor(r_gStyle->GetColorPalette((float)nColors/Files.size()*(i+1) ));
		if(i==0) GoodUtof_HE[i]->Draw();
		else GoodUtof_HE[i]->Draw("same");
	}
	c->cd(2);
	for(int i=0;i<Files.size();i++) {
		GoodLtof_HE[i]->SetLineWidth(3);
		GoodLtof_HE[i]->SetLineColor(r_gStyle->GetColorPalette((float)nColors/Files.size()*(i+1) ));
		if(i==0) GoodLtof_HE[i]->Draw();
		else GoodLtof_HE[i]->Draw("same");
	}


	TCanvas * c0 = new TCanvas("Time Of Flight");
        c0->Divide(1,2);
        c0->cd(1);
	for(int i=0;i<Files.size();i++) {
		GoodTime_TOF[i]->SetLineWidth(3);
		GoodTime_TOF[i]->SetLineColor(r_gStyle->GetColorPalette((float)nColors/Files.size()*(i+1) ));
		if(i==0) GoodTime_TOF[i]->Draw();
		else GoodTime_TOF[i]->Draw("same");
	}
	c0->cd(2);
	for(int i=0;i<Files.size();i++) {
		Quality_TOF[i]->SetLineWidth(3);
		Quality_TOF[i]->SetLineColor(r_gStyle->GetColorPalette((float)nColors/Files.size()*(i+1) ));
		if(i==0) Quality_TOF[i]->Draw();
		else Quality_TOF[i]->Draw("same");
	}
	
	TCanvas * c1 = new TCanvas("RICH BDT");
	c1->Divide(1,2);
	c1->cd(1);
	for(int i=0;i<Files.size();i++) {
		RICHQualTSpline3_NaF[i]->SetLineWidth(3);
		RICHQualTSpline3_NaF[i]->SetLineColor(r_gStyle->GetColorPalette((float)nColors/Files.size()*(i+1) ));
		if(i==0) RICHQualTSpline3_NaF[i]->Draw();
		else RICHQualTSpline3_NaF[i]->Draw("same");
	}
	c1->cd(2);
	for(int i=0;i<Files.size();i++) {
		RICHQualTSpline3_Agl[i]->SetLineWidth(3);
		RICHQualTSpline3_Agl[i]->SetLineColor(r_gStyle->GetColorPalette((float)nColors/Files.size() *(i+1)));
		if(i==0)RICHQualTSpline3_Agl[i]->Draw();
		else RICHQualTSpline3_Agl[i]->Draw("same");
	}
	
	TCanvas * c2 = new TCanvas("RICH CIEMAT");
	c2->Divide(1,2);
	c2->cd(1);
	for(int i=0;i<Files.size();i++) {
		RICHTSpline3_NaF[i]->SetLineWidth(3);
		RICHTSpline3_NaF[i]->SetLineColor(r_gStyle->GetColorPalette((float)nColors/Files.size()*(i+1) ));
		if(i==0) RICHTSpline3_NaF[i]->Draw();
		else RICHTSpline3_NaF[i]->Draw("same");
	}
	c2->cd(2);
	for(int i=0;i<Files.size();i++) {
		RICHTSpline3_Agl[i]->SetLineWidth(3);
		RICHTSpline3_Agl[i]->SetLineColor(r_gStyle->GetColorPalette((float)nColors/Files.size() *(i+1)));
		if(i==0)RICHTSpline3_Agl[i]->Draw();
		else RICHTSpline3_Agl[i]->Draw("same");
	}

	Plots.Add(f);
	Plots.Add(c);
	Plots.Add(c0);
	Plots.Add(c2);
	Plots.Add(c1);
	Plots.writeObjsInFolder("Eff. Corrections");

}	

int main(int argc, char * argv[]){

	std::vector<TFile *> Files;
	for(int i=0; i< GroupedFiles.size(); i++){
		cout<<GroupedFiles[i].c_str()<<endl;
		TFile * f = TFile::Open(GroupedFiles[i].c_str()); 
		Files.push_back(f);
	}
	cout<<"Files stored: "<<Files.size()<<endl;
	std::vector<FileSaver> files;
	for(int i=0; i< GroupedFiles.size(); i++){
		FileSaver f;
		f.setName(GroupedFiles[i].c_str());
		files.push_back(f);
	}
	cout<<"files stored: "<<files.size()<<endl;
	
	FileSaver Plots;
	Plots.setName("/afs/cern.ch/user/f/fdimicco/Work/Deutons/DirectAnalysis/Plotting/Time.root");
        cout<<"****************************** BINS ***************************************"<<endl;
	SetUpTOIBinning();
	
	cout<<"**************** plotting *****************"<<endl;
	//DrawDPRatio(Plots,Files );
	DrawPDFluxRatio(Plots,files,Binselected2);
	//DrawMassDistributions(Plots,files,Binselected);
	DrawTimeEffCorr(Plots,files);
	//DrawDPRatioEkin(Plots,files);
	//DrawPFluxEkin(Plots,files);
//	DrawDFluxEkin(Plots,files);
//	DrawLatWeights(Plots,files);
	DrawParameters(Plots,files);

	return 0;
}
