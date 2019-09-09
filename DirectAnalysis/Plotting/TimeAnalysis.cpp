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

int start =0;

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
1548072000
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
//	TimeDepP->GetXaxis()->SetLabelSize(0.035*gPad->GetCanvas()->GetWh()/(float)gPad->YtoPixel(gPad->GetY1()));
//	TimeDepP->GetXaxis()->SetLabelFont(21);
	

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

void PlotTimeRatio(TVirtualPad * c, float D[], float D_err[], float P[], float P_err[], std::vector <int> Times, int col, std::string title) {
	TGraphErrors * TimeDepP = (TGraphErrors * ) gPad->FindObject("Time");
	TGraphErrors * TimeDepR = new TGraphErrors();
	
	for(int i=start;i<Times.size();i++){
		TimeDepR->SetPoint(i-start,Times[i],D[i]/P[i]);
		TimeDepR->SetPointError(i-start,0,pow(pow(D_err[i]/D[i],2) + pow (P_err[i]/P[i],2),0.5)*D[i]/(P[i]) );
	
	}
	TimeDepR->SetName(title.c_str());
	 TimeDepR->SetTitle(title.c_str());

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
	TimeDepR->GetYaxis()->SetLabelSize(0.015*gPad->GetCanvas()->GetWh()/(float)gPad->YtoPixel(gPad->GetY1()));
	TimeDepR->GetYaxis()->SetLabelFont(11);
	

	c->cd();
	if(!TimeDepP)	TimeDepR->Draw("AP");
	else      TimeDepR->Draw("Psame");
	
	gPad->Update();
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

	float Bin_D[6][Files.size()];
	float Bin_P[6][Files.size()];
	float Bin_D_err[6][Files.size()];
	float Bin_P_err[6][Files.size()];
	int binselected[6]= {3,7,11,15,22,29};
	int binselected_P[6]= {3,7,11,15,22,29};



	for(int j=0;j<6;j++){
		for(int i=start;i<Files.size();i++){
			Bin_D[j][i]= Merged_D[i]->GetBinContent(binselected[j]);
			Bin_P[j][i]= Merged_P[i]->GetBinContent(binselected[j]);
			Bin_D_err[j][i]= Merged_D[i]->GetBinError(binselected[j]);
			Bin_P_err[j][i]= Merged_P[i]->GetBinError(binselected[j]);
		}
	}


	TCanvas *c = new TCanvas("Time Dep. Allbins");
	c->SetCanvasSize(1000,700);
	c->Divide(1,6,0,0);
	c->cd(1);	
	PlotTimeDep(gPad,Bin_D[0],Bin_D_err[0],Bin_P[0],Bin_P_err[0], Times);
	c->cd(2);                                      
       	PlotTimeDep(gPad,Bin_D[1],Bin_D_err[1],Bin_P[1],Bin_P_err[1], Times);
	c->cd(3);                                      
       	PlotTimeDep(gPad,Bin_D[2],Bin_D_err[2],Bin_P[2],Bin_P_err[2], Times);
	c->cd(4);	                               
 	PlotTimeDep(gPad,Bin_D[3],Bin_D_err[3],Bin_P[3],Bin_P_err[3], Times);
	c->cd(5);                                      
 	PlotTimeDep(gPad,Bin_D[4],Bin_D_err[4],Bin_P[4],Bin_P_err[4], Times);
	c->cd(6);                                      
 	PlotTimeDep(gPad,Bin_D[5],Bin_D_err[5],Bin_P[5],Bin_P_err[5], Times);
	
	Plots.Add(c);
        Plots.writeObjsInFolder("Flux");

	TCanvas *c1 = new TCanvas("Time Ratio Allbins");
	c1->SetCanvasSize(1000,700);
	c1->cd();
	c1->Divide(1,6,0,0);
	c1->cd(1);	
	gPad->SetTickx();
	gPad->SetTicky();
	PlotTimeRatio(gPad,Bin_D[0],Bin_D_err[0],Bin_P[0],Bin_P_err[0], Times,2,(Convert(Global.GetGlobalPBins().RigTOIBins()[binselected[0]])  + " < R < " + Convert(Global.GetGlobalPBins().RigTOIBins()[binselected[0]+1]) + "GV").c_str());
	c1->cd(2);                                      
       	gPad->SetTickx();
	gPad->SetTicky();
	PlotTimeRatio(gPad,Bin_D[1],Bin_D_err[1],Bin_P[1],Bin_P_err[1], Times,2,(Convert(Global.GetGlobalPBins().RigTOIBins()[binselected[1]])  + " < R < " + Convert(Global.GetGlobalPBins().RigTOIBins()[binselected[1]+1]) + "GV").c_str());
	c1->cd(3);                                      
       	gPad->SetTickx();
	gPad->SetTicky();
	PlotTimeRatio(gPad,Bin_D[2],Bin_D_err[2],Bin_P[2],Bin_P_err[2], Times,2,(Convert(Global.GetGlobalPBins().RigTOIBins()[binselected[2]])  + " < R < " + Convert(Global.GetGlobalPBins().RigTOIBins()[binselected[2]+1]) + "GV").c_str());
	c1->cd(4);	                               
 	gPad->SetTickx();
	gPad->SetTicky();
	PlotTimeRatio(gPad,Bin_D[3],Bin_D_err[3],Bin_P[3],Bin_P_err[3], Times,2,(Convert(Global.GetGlobalPBins().RigTOIBins()[binselected[3]])  + " < R < " + Convert(Global.GetGlobalPBins().RigTOIBins()[binselected[3]+1]) + "GV").c_str());
	c1->cd(5);                                      
 	gPad->SetTickx();
	gPad->SetTicky();
	PlotTimeRatio(gPad,Bin_D[4],Bin_D_err[4],Bin_P[4],Bin_P_err[4], Times,2,(Convert(Global.GetGlobalPBins().RigTOIBins()[binselected[4]])  + " < R < " + Convert(Global.GetGlobalPBins().RigTOIBins()[binselected[4]+1]) + "GV").c_str());
	c1->cd(6);                                      
 	gPad->SetTickx();
	gPad->SetTicky();
	PlotTimeRatio(gPad,Bin_D[5],Bin_D_err[5],Bin_P[5],Bin_P_err[5], Times,2,(Convert(Global.GetGlobalPBins().RigTOIBins()[binselected[5]])  + " < R < " + Convert(Global.GetGlobalPBins().RigTOIBins()[binselected[5]+1]) + "GV").c_str());
	
	Plots.Add(c1);
        Plots.writeObjsInFolder("Flux");


	return;
}



int main(int argc, char * argv[]){

	std::vector<TFile *> Files;
	for(int i=0; i< GroupedFiles.size(); i++){
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
	Plots.setName("/data1/home/fdimicco/Deutons/DirectAnalysis/Plotting/Time.root");
        cout<<"****************************** BINS ***************************************"<<endl;
	SetUpTOIBinning();
	
	cout<<"**************** PLOTTING *****************"<<endl;
	//DrawDPRatio(Plots,Files );
	DrawPDFluxRatio(Plots,files);
	return 0;
}
