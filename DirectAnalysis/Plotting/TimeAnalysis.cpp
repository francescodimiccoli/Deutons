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

#include "../perl/List.h"

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

	for(int i=0;i<Files.size();i++){
		RatiosTOF[i]=(TH1F*)Files[i]->Get("Fluxes/DP ratio TOF");
		RatiosNaF[i]=(TH1F*)Files[i]->Get("Fluxes/DP ratio NaF");
		RatiosAgl[i]=(TH1F*)Files[i]->Get("Fluxes/DP ratio Agl");
	}
	
	for(int i=0;i<Files.size();i++){
		PlotMergedRanges(gPad,RatiosTOF[i] ,RatiosNaF[i] ,RatiosAgl[i] ,"Kinetic Energy [GeV/nucl.]", "Flux",m_gStyle->GetColorPalette( ((float)nColors/Files.size()) *(i+1)),true,"Psame",0.1,10,0.00001,0.12,"This Work (TOF)",8);	
	}	

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


void DrawPDCountsRatio(FileSaver Plots, std::vector<TFile *> Files ){

	TStyle* m_gStyle= new TStyle();;
	m_gStyle->SetPalette(55);
	int nColors = m_gStyle->GetNumberOfColors();

	TCanvas *c3 = new TCanvas("D/P Counts");
	c3->SetCanvasSize(2000,1500);
	TH1F * DCountsTOF[Files.size()];
	TH1F * DCountsNaF[Files.size()];
	TH1F * DCountsAgl[Files.size()];

	for(int i=0;i<Files.size();i++){
		DCountsTOF[i]=(TH1F*)Files[i]->Get("TOFfits/Fit Results/Deuteron Counts");
		DCountsNaF[i]=(TH1F*)Files[i]->Get("NaFfits/Fit Results/Deuteron Counts");
		DCountsAgl[i]=(TH1F*)Files[i]->Get("Aglfits/Fit Results/Deuteron Counts");
	}

	TH1F * PCountsTOF[Files.size()];
	TH1F * PCountsNaF[Files.size()];
	TH1F * PCountsAgl[Files.size()];

	for(int i=0;i<Files.size();i++){
		PCountsTOF[i]=(TH1F*)Files[i]->Get("TOFfits/Fit Results/Proton Counts");
		PCountsNaF[i]=(TH1F*)Files[i]->Get("NaFfits/Fit Results/Proton Counts");
		PCountsAgl[i]=(TH1F*)Files[i]->Get("Aglfits/Fit Results/Proton Counts");
	}

	for(int i=0;i<Files.size();i++){
	DCountsTOF[i]->Divide(PCountsTOF[i]	);
	DCountsNaF[i]->Divide(PCountsNaF[i]	);
        DCountsAgl[i]->Divide(PCountsAgl[i]	);
	}

	for(int i=0;i<Files.size();i++){
		PlotMergedRanges(gPad,DCountsTOF[i] ,DCountsNaF[i] ,DCountsAgl[i] ,"Kinetic Energy [GeV/nucl.]", "Flux",m_gStyle->GetColorPalette( ((float)nColors/Files.size()) *(i+1)),true,"Psame",0.1,10,0.00001,0.12,"This Work (TOF)",8);	
	}	

        Plots.Add(c3);
        Plots.writeObjsInFolder("Counts");

	return;
}

int main(int argc, char * argv[]){

	std::vector<TFile *> Files;
	for(int i=0; i< TimeFiles.size(); i++){
		TFile * f = TFile::Open(TimeFiles[i].c_str()); 
		Files.push_back(f);
	}
	cout<<"Files stored: "<<Files.size()<<endl;
	
	FileSaver Plots;
	Plots.setName("Time.root");
        cout<<"****************************** BINS ***************************************"<<endl;

        SetBins();

	PRB.Print();

        cout<<"**TOF**"<<endl;
        ToFDB.Print();

        cout<<"**NaF**"<<endl;
        NaFDB.Print();

        cout<<"**Agl**"<<endl;
        AglDB.Print();

        ToFDB.UseREdges();
        NaFDB.UseREdges();
        AglDB.UseREdges();
	
	ToFPB.UseREdges();
        NaFPB.UseREdges();
        AglPB.UseREdges();

        PRB.UseREdges();


        cout<<endl;


	cout<<"**************** PLOTTING *****************"<<endl;
	//DrawDPRatio(Plots,Files );
	DrawPDCountsRatio(Plots,Files);
	return 0;
}
