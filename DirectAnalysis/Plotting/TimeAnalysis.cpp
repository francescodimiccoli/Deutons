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

#include <chrono>
#include <iomanip>
#include <ctime>
int start =2;
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

void DrawDPRatio(FileSaver Plots, std::vector<TFile *> Files ){
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


	TCanvas *c3 = new TCanvas("D/P ratio");
	c3->SetCanvasSize(2000,1500);
	TH1F * RatiosTOF[Files.size()];
	TH1F * RatiosNaF[Files.size()];
	TH1F * RatiosAgl[Files.size()];

	for(int i=start;i<Files.size();i++){
		RatiosTOF[i]=(TH1F*)Files[i]->Get("Fluxes/DP ratio TOF");
		RatiosNaF[i]=(TH1F*)Files[i]->Get("Fluxes/DP ratio NaF");
		RatiosAgl[i]=(TH1F*)Files[i]->Get("Fluxes/DP ratio Agl");
	}
	
/*	for(int i=start;i<Files.size();i++){
		PlotMergedRanges(gPad,RatiosTOF[i] ,RatiosNaF[i] ,RatiosAgl[i] ,"Kinetic Energy [GeV/nucl.]", "Flux",m_gStyle->GetColorPalette( ((float)nColors/Files.size()) *(i+1)),true,"Psame",0.1,10,0.00001,0.12,"This Work (TOF)",8);	
	}	

        TLegend * leg = (TLegend*) gPad->FindObject("leg");

        for(uint n=0;n<PD_Graphs.size();n++){
                PD_Graphs[n] ->Draw("Psame");
                PD_Graphs[n]->SetMarkerSize(2);
                leg->AddEntry(PD_Graphs[n],PD_Graphs[n]->GetTitle(),"p");
        }

	DrawGalpropRatio(c3);
*/
        Plots.Add(c3);
        Plots.writeObjsInFolder("Fluxes");
	


	TCanvas *c4 = new TCanvas("Double ratio");
	c4->SetCanvasSize(1000,700);
	
	for(int i=start;i<Files.size();i++){
		RatiosTOF[i]->Divide(RatiosTOF[start]);
		RatiosNaF[i]->Divide(RatiosNaF[start]);
		RatiosAgl[i]->Divide(RatiosAgl[start]);

	}
	

	 for(int i=start;i<Files.size();i++){
                PlotMergedRanges(gPad,RatiosTOF[i] ,RatiosNaF[i] ,RatiosAgl[i] ,"Kinetic Energy [GeV/nucl.]", "Flux",m_gStyle->GetColorPalette( ((float)nColors/Files.size()) *(i+1)),true,"Psame",0.1,10,0.1,5,"This Work (TOF)",8);
        }

	Plots.Add(c4);
        Plots.writeObjsInFolder("Fluxes");


}

void PlotTimeDep(TVirtualPad * c, float D[], float D_err[], float P[], float P_err[], std::vector <int> Times) {

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
	TimeDepP->GetYaxis()->SetRangeUser(0.4*(TimeDepP->GetHistogram()->GetMinimum() + TimeDepP->GetHistogram()->GetMaximum())/2,1.6*(TimeDepP->GetHistogram()->GetMinimum() + TimeDepP->GetHistogram()->GetMaximum())/2);


	TimeDepP->GetXaxis()->	SetTimeDisplay(1);
	TimeDepP->GetXaxis()->  SetTimeOffset(0,"gmt");  
	TimeDepP->GetXaxis()->  SetTimeFormat("%b%y");
	TimeDepP->GetXaxis()->SetLabelSize(0.035*gPad->GetCanvas()->GetWh()/(float)gPad->YtoPixel(gPad->GetY1()));
	TimeDepP->GetXaxis()->SetLabelFont(21);
	

	c->cd();
	TPad *p1 = new TPad("p1", "", 0, 0, 1, 1);
	p1->SetGrid();
	TPad *p2 = new TPad("p2", "", 0, 0, 1, 1);
	p2->SetFillStyle(4000); // will be transparent
	// p2->SetGrid();
	
	p1->Draw();
  	p1->cd();
	TimeDepP->Draw("AP");
  	gPad->Update();


	Double_t xmin = p1->GetUxmin();
	Double_t xmax = p1->GetUxmax();
	Double_t ymin = 0.4*(TimeDepD->GetHistogram()->GetMinimum() + TimeDepD->GetHistogram()->GetMaximum())/2;
	Double_t ymax = 1.6*(TimeDepD->GetHistogram()->GetMinimum() + TimeDepD->GetHistogram()->GetMaximum())/2;
	Double_t dx = (xmax - xmin) / 0.8; // 10 percent margins left and right
	Double_t dy = (ymax - ymin) / 0.8; // 10 percent margins top and bottom
	
	p2->Range(xmin-0.1*dx, ymin-0.1*dy, xmax+0.1*dx, ymax+0.1*dy);
	p2->Draw();
	p2->cd();
	TimeDepD->GetYaxis()->SetRangeUser(0.4*ymin,1.4*ymax);
	TimeDepD->Draw("P");
	gPad->Update();

	TGaxis *axis = new TGaxis(xmax, ymin, xmax, ymax, ymin, ymax, 510, "+L");
	axis->SetLineColor(kBlue);
	axis->SetLabelColor(kBlue);
	axis->Draw();
	gPad->Update();

}

void PlotTimeRatio(TVirtualPad * c, float D[], float D_err[], float P[], float P_err[], std::vector <int> Times, int col) {
	TGraphErrors * TimeDepP = (TGraphErrors * ) gPad->FindObject("Time");
	TGraphErrors * TimeDepR = new TGraphErrors();
	
	for(int i=start;i<Times.size();i++){
		TimeDepR->SetPoint(i-start,Times[i],D[i]/P[i]);
		TimeDepR->SetPointError(i-start,0,pow(pow(D_err[i]/D[i],2) + pow (P_err[i]/P[i],2),0.5)*D[i]/(P[i]) );
	
	}
	TimeDepR->SetName("Time");
	 TimeDepR->SetTitle("Time");

	TimeDepR->SetMarkerStyle(8);
	TimeDepR->SetMarkerSize(1.5);
	TimeDepR->SetMarkerColor(col);
	TimeDepR->SetLineColor(col);
	TimeDepR->SetLineWidth(3);
	TimeDepR->GetYaxis()->SetRangeUser(0,0.035);
	TimeDepR->GetXaxis()->	SetTimeDisplay(1);
	TimeDepR->GetXaxis()->  SetTimeOffset(0,"gmt");  
	TimeDepR->GetXaxis()->  SetTimeFormat("%b%y");
	TimeDepR->GetYaxis()->SetLabelColor(1);
	TimeDepR->GetYaxis()->SetLabelSize(0.035*gPad->GetCanvas()->GetWh()/(float)gPad->YtoPixel(gPad->GetY1()));
	TimeDepR->GetYaxis()->SetLabelFont(21);
	

	c->cd();
	if(!TimeDepP)	TimeDepR->Draw("AP");
	else      TimeDepR->Draw("Psame");
	
	gPad->Update();
}

void DrawPDCountsRatio(FileSaver Plots, std::vector<TFile *> Files ){

	TStyle* m_gStyle= new TStyle();;
	m_gStyle->SetPalette(55);
	int nColors = m_gStyle->GetNumberOfColors();
	
	std::vector<int> Times;
	std::vector<std::string> Dates;


	for(int i=0;i<Files.size();i++){
		Times.push_back(std::atoi(TimeFiles[i].substr(TimeFiles[i].find("-")+1,10).c_str()) );
		std::uint32_t time_date_stamp = Times[i];
		std::time_t temp = time_date_stamp;
		std::tm* t = std::gmtime(&temp);
		char buffer [10];
		sprintf(buffer,"%d/%d", t->tm_mon+1, t->tm_year-100);
		Dates.push_back(buffer);
		std::cout<<Times[i]<<" "<<Dates[i]<<endl;
	}
	
	TCanvas *c3 = new TCanvas("D/P Counts");
	c3->SetCanvasSize(2000,1500);
	gPad->SetLogx();

	TH1F * RatioTOF[Files.size()];
	TH1F * RatioNaF[Files.size()];
	TH1F * RatioAgl[Files.size()];

	TH1F * DCountsTOF[Files.size()];
	TH1F * DCountsNaF[Files.size()];
	TH1F * DCountsAgl[Files.size()];

	for(int i=start;i<Files.size();i++){
		DCountsTOF[i]=(TH1F*)Files[i]->Get("TOFfits/Fit Results/Primary Deuteron Counts");
		DCountsNaF[i]=(TH1F*)Files[i]->Get("NaFfits/Fit Results/Primary Deuteron Counts");
		DCountsAgl[i]=(TH1F*)Files[i]->Get("Aglfits/Fit Results/Primary Deuteron Counts");
	}

	TH1F * PCountsTOF[Files.size()];
	TH1F * PCountsNaF[Files.size()];
	TH1F * PCountsAgl[Files.size()];

	for(int i=start;i<Files.size();i++){
		PCountsTOF[i]=(TH1F*)Files[i]->Get("TOFfits/Fit Results/Primary Proton Counts");
		PCountsNaF[i]=(TH1F*)Files[i]->Get("NaFfits/Fit Results/Primary Proton Counts");
		PCountsAgl[i]=(TH1F*)Files[i]->Get("Aglfits/Fit Results/Primary Proton Counts");
	}

	for(int i=start;i<Files.size();i++){
	RatioTOF[i] 	= (TH1F*)DCountsTOF[i]->Clone();
	RatioNaF[i]	= (TH1F*)DCountsNaF[i]->Clone();
	RatioAgl[i]	= (TH1F*)DCountsAgl[i]->Clone();
	RatioTOF[i]->Divide(PCountsTOF[i]	);
	RatioNaF[i]->Divide(PCountsNaF[i]	);
        RatioAgl[i]->Divide(PCountsAgl[i]	);
	}

	for(int i=start;i<Files.size();i++){
		PlotMergedRanges(gPad,RatioTOF[i] ,RatioNaF[i] ,RatioAgl[i] ,"Ekin [GeV/n]", "Primary Counts (D/P)",m_gStyle->GetColorPalette( ((float)nColors/Files.size()) *(i+1)),true,"Psame",0.1,10,0.00001,0.4,Dates[i],8);	
	}	

        Plots.Add(c3);
        Plots.writeObjsInFolder("Counts");


	float DTOF[Files.size()];
	float DTOF_err[Files.size()];
	float PTOF[Files.size()];
	float PTOF_err[Files.size()];

	for(int i=start;i<Files.size();i++){
		DTOF[i]=DCountsTOF[i]->GetBinContent(4);
		DTOF_err[i] = DCountsTOF[i]->GetBinError(4);		
		PTOF[i]=PCountsTOF[i]->GetBinContent(18);
		PTOF_err[i] = PCountsTOF[i]->GetBinError(18);		
	}

	float DNaF[Files.size()];
	float DNaF_err[Files.size()];
	float PNaF[Files.size()];
	float PNaF_err[Files.size()];

	for(int i=start;i<Files.size();i++){
		DNaF[i]=DCountsNaF[i]->GetBinContent(5);
		DNaF_err[i] = DCountsNaF[i]->GetBinError(5);		
		PNaF[i]=PCountsNaF[i]->GetBinContent(15);
		PNaF_err[i] = PCountsNaF[i]->GetBinError(15);		
	}
	float DAgl[Files.size()];
	float DAgl_err[Files.size()];
	float PAgl[Files.size()];
	float PAgl_err[Files.size()];

	for(int i=start;i<Files.size();i++){
		DAgl[i]=DCountsAgl[i]->GetBinContent(4);
		DAgl_err[i] = DCountsAgl[i]->GetBinError(4);		
		PAgl[i]=PCountsAgl[i]->GetBinContent(13);
		PAgl_err[i] = PCountsAgl[i]->GetBinError(13);		
	}


	TCanvas *c4 = new TCanvas("Time Dep. 1.3 GV");
	c4->SetCanvasSize(1000,300);
	PlotTimeDep(c4,DTOF,DTOF_err,PTOF,PTOF_err, Times);
	
	Plots.Add(c4);
	Plots.writeObjsInFolder("Counts");

	TCanvas *c5 = new TCanvas("Time Dep. 3.5 GV");
	c5->SetCanvasSize(1000,300);
	PlotTimeDep(c5,DNaF,DNaF_err,PNaF,PNaF_err, Times);
	
	Plots.Add(c5);
	Plots.writeObjsInFolder("Counts");

	TCanvas *c6 = new TCanvas("Time Dep. 7.8 GV");
	c6->SetCanvasSize(1000,300);
	PlotTimeDep(c6,DAgl,DAgl_err,PAgl,PAgl_err, Times);
	
	Plots.Add(c6);
	Plots.writeObjsInFolder("Counts");


	return;
}


void DrawPDFluxRatio(FileSaver Plots, std::vector<FileSaver> Files ){

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

		Flux * DFluxTOF = new Flux(Files[i], "DFluxTOF", "FullsetEff_D_TOF","FullsetEfficiency","TOFfits/Fit Results/Primary Deuteron Counts","ExposureTOF",ToFDB);
		Flux * DFluxNaF = new Flux(Files[i], "DFluxNaF", "FullsetEff_D_NaF","FullsetEfficiency","NaFfits/Fit Results/Primary Deuteron Counts","ExposureNaF",NaFDB);
		Flux * DFluxAgl = new Flux(Files[i], "DFluxAgl", "FullsetEff_D_Agl","FullsetEfficiency","Aglfits/Fit Results/Primary Deuteron Counts","ExposureAgl",AglDB);

		Flux * PFluxTOF = new Flux(Files[i], "PFluxTOF", "FullsetEff_P_TOF","FullsetEfficiency","TOFfits/Fit Results/Primary Proton Counts","ExposureTOF",ToFPB);
		Flux * PFluxNaF = new Flux(Files[i], "PFluxNaF", "FullsetEff_P_NaF","FullsetEfficiency","NaFfits/Fit Results/Primary Proton Counts","ExposureNaF",NaFPB);
		Flux * PFluxAgl = new Flux(Files[i], "PFluxAgl", "FullsetEff_P_Agl","FullsetEfficiency","Aglfits/Fit Results/Primary Proton Counts","ExposureAgl",AglPB);

		FluxesPTOF.push_back(PFluxTOF);
		FluxesPNaF.push_back(PFluxNaF);
		FluxesPAgl.push_back(PFluxAgl);
	
		FluxesDTOF.push_back(DFluxTOF);
		FluxesDNaF.push_back(DFluxNaF);
		FluxesDAgl.push_back(DFluxAgl);
	}

	for(int i=0;i<Files.size();i++){
		Times.push_back(std::atoi(TimeFiles[i].substr(TimeFiles[i].find("-")+1,10).c_str()) );
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

	TH1F * RatioTOF[Files.size()];
	TH1F * RatioNaF[Files.size()];
	TH1F * RatioAgl[Files.size()];

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

	for(int i=start;i<Files.size();i++){
	RatioTOF[i] 	= (TH1F*)DCountsTOF[i]->Clone();
	RatioNaF[i]	= (TH1F*)DCountsNaF[i]->Clone();
	RatioAgl[i]	= (TH1F*)DCountsAgl[i]->Clone();
	RatioTOF[i]->Divide(PCountsTOF[i]	);
	RatioNaF[i]->Divide(PCountsNaF[i]	);
        RatioAgl[i]->Divide(PCountsAgl[i]	);
	}

	for(int i=start;i<Files.size();i++){
		PlotMergedRanges(gPad,RatioTOF[i] ,RatioNaF[i] ,RatioAgl[i] ,"Ekin [GeV/n]", "Primary Counts (D/P)",f_gStyle->GetColorPalette( ((float)nColors/Dates.size()) *(i+1)),true,"Psame",0.1,10,0.00001,0.15,Dates[i],8);	
	}	

        TLegend * leg2 = new TLegend();

        for(uint n=0;n<PD_Graphs.size();n++){
                PD_Graphs[n] ->Draw("Psame");
                PD_Graphs[n]->SetMarkerSize(2);
                leg2->AddEntry(PD_Graphs[n],PD_Graphs[n]->GetTitle(),"p");
        }
	leg2->Draw("same");

	DrawGalpropRatio(c3);
        Plots.Add(c3);
        Plots.writeObjsInFolder("Flux");

/*	TCanvas *c4_ = new TCanvas("D/P Double Ratio");
	c4_->SetCanvasSize(1000,600);
	gPad->SetLogx();

	TH1F * RatioTOF_start = (TH1F*) RatioTOF[start]->Clone();
	TH1F * RatioNaF_start = (TH1F*) RatioNaF[start]->Clone();
	TH1F * RatioAgl_start = (TH1F*) RatioAgl[start]->Clone();

	for(int i=start;i<Files.size();i++){
		RatioTOF[i] 	->Divide(RatioTOF_start); 
		RatioNaF[i]	->Divide(RatioNaF_start);  
		RatioAgl[i]	->Divide(RatioAgl_start);  
	}


	for(int i=start;i<Files.size();i++){
		PlotMergedRanges(gPad,RatioTOF[i] ,RatioNaF[i] ,RatioAgl[i] ,"Ekin [GeV/n]", "Primary Counts (D/P)",1,true,"Psame",0.1,10,0.1,4,Dates[i],8);	
	}	


	
	Plots.Add(c4_);
        Plots.writeObjsInFolder("Flux");
*/

	TH1F * DCountsTOF_rig[Files.size()];
	TH1F * DCountsNaF_rig[Files.size()];
	TH1F * DCountsAgl_rig[Files.size()];

	for(int i=start;i<Files.size();i++){
		DCountsTOF_rig[i]=(TH1F*) FluxesDTOF[i]->GetFlux_rig();
		DCountsNaF_rig[i]=(TH1F*) FluxesDNaF[i]->GetFlux_rig();
		DCountsAgl_rig[i]=(TH1F*) FluxesDAgl[i]->GetFlux_rig();
	}

	TH1F * PCountsTOF_rig[Files.size()];
	TH1F * PCountsNaF_rig[Files.size()];
	TH1F * PCountsAgl_rig[Files.size()];

	for(int i=start;i<Files.size();i++){
		PCountsTOF_rig[i]=(TH1F*)FluxesPTOF[i]->GetFlux_rig();
		PCountsNaF_rig[i]=(TH1F*)FluxesPNaF[i]->GetFlux_rig();
		PCountsAgl_rig[i]=(TH1F*)FluxesPAgl[i]->GetFlux_rig();
	}



	float DTOF[Files.size()];
	float DTOF_err[Files.size()];
	float PTOF[Files.size()];
	float PTOF_err[Files.size()];

	for(int i=start;i<Files.size();i++){
		DTOF[i]=DCountsTOF_rig[i]->GetBinContent(4);
		DTOF_err[i] = DCountsTOF_rig[i]->GetBinError(4);		
		PTOF[i]=PCountsTOF_rig[i]->GetBinContent(18);
		PTOF_err[i] = PCountsTOF_rig[i]->GetBinError(18);		
	}

	float DNaF[Files.size()];
	float DNaF_err[Files.size()];
	float PNaF[Files.size()];
	float PNaF_err[Files.size()];

	for(int i=start;i<Files.size();i++){
		DNaF[i]=DCountsNaF_rig[i]->GetBinContent(5);
		DNaF_err[i] = DCountsNaF_rig[i]->GetBinError(5);		
		PNaF[i]=PCountsNaF_rig[i]->GetBinContent(15);
		PNaF_err[i] = PCountsNaF_rig[i]->GetBinError(15);		
	}
	float DAgl[Files.size()];
	float DAgl_err[Files.size()];
	float PAgl[Files.size()];
	float PAgl_err[Files.size()];

	for(int i=start;i<Files.size();i++){
		DAgl[i]=DCountsAgl_rig[i]->GetBinContent(4);
		DAgl_err[i] = DCountsAgl_rig[i]->GetBinError(4);		
		PAgl[i]=PCountsAgl_rig[i]->GetBinContent(13);
		PAgl_err[i] = PCountsAgl_rig[i]->GetBinError(13);		
	}

/*
	TCanvas *c4 = new TCanvas("Time Dep. 1.3 GV");
	c4->SetCanvasSize(1000,700);
	c4->Divide(1,2);
	c4->cd(1);	
	PlotTimeDep(gPad,DTOF,DTOF_err,PTOF,PTOF_err, Times);
	c4->cd(2);
	PlotTimeRatio(gPad,DTOF,DTOF_err,PTOF,PTOF_err, Times);

	Plots.Add(c4);
	Plots.writeObjsInFolder("Flux");

	TCanvas *c5 = new TCanvas("Time Dep. 3.5 GV");
	c5->SetCanvasSize(1000,700);
	c5->Divide(1,2);
	c5->cd(1);	
	PlotTimeDep(gPad,DNaF,DNaF_err,PNaF,PNaF_err, Times);
	c5->cd(2);
	PlotTimeRatio(gPad,DNaF,DNaF_err,PNaF,PNaF_err, Times);
	
	Plots.Add(c5);
	Plots.writeObjsInFolder("Flux");


	TCanvas *c6 = new TCanvas("Time Dep. 7.8 GV");
	c6->SetCanvasSize(1000,700);
	c6->Divide(1,2);
	c6->cd(1);	
	PlotTimeDep(gPad,DAgl,DAgl_err,PAgl,PAgl_err, Times);
	c6->cd(2);
	PlotTimeRatio(gPad,DAgl,DAgl_err,PAgl,PAgl_err, Times);


	Plots.Add(c6);
	Plots.writeObjsInFolder("Flux");
*/
	TCanvas *c = new TCanvas("Time Dep. Allbins");
	c->SetCanvasSize(1000,700);
	c->Divide(1,3);
	c->cd(1);	
	PlotTimeDep(gPad,DTOF,DTOF_err,PTOF,PTOF_err, Times);
	c->cd(2);      
        PlotTimeDep(gPad,DNaF,DNaF_err,PNaF,PNaF_err, Times);
	c->cd(3);      
        PlotTimeDep(gPad,DAgl,DAgl_err,PAgl,PAgl_err, Times);		
	
	Plots.Add(c);
        Plots.writeObjsInFolder("Flux");

	TCanvas *c1 = new TCanvas("Time Ratio Allbins");
	c1->SetCanvasSize(1000,700);
	c1->cd();
	PlotTimeRatio(gPad,DTOF,DTOF_err,PTOF,PTOF_err, Times, 1);
	PlotTimeRatio(gPad,DNaF,DNaF_err,PNaF,PNaF_err, Times, 12);
	PlotTimeRatio(gPad,DAgl,DAgl_err,PAgl,PAgl_err, Times, 15);

	Plots.Add(c1);
        Plots.writeObjsInFolder("Flux");


	return;
}



int main(int argc, char * argv[]){

	std::vector<TFile *> Files;
	for(int i=0; i< TimeFiles.size(); i++){
		TFile * f = TFile::Open(TimeFiles[i].c_str()); 
		Files.push_back(f);
	}
	cout<<"Files stored: "<<Files.size()<<endl;
	std::vector<FileSaver> files;
	for(int i=0; i< TimeFiles.size(); i++){
		FileSaver f;
		f.setName(TimeFiles[i].c_str());
		files.push_back(f);
	}
	cout<<"files stored: "<<files.size()<<endl;
	
	FileSaver Plots;
	Plots.setName("/data1/home/fdimicco/Deutons/DirectAnalysis/Plotting/Time.root");
        cout<<"****************************** BINS ***************************************"<<endl;
	SetUpTOIBinning();
	
	cout<<"**************** PLOTTING *****************"<<endl;
	//DrawDPRatio(Plots,Files );
	DrawPDFluxRatio(Plots,files);
	//DrawPDCountsRatio(Plots,Files);
	return 0;
}
