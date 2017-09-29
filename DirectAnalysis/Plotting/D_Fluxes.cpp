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
#include "../include/GlobalBinning.h"
#include "TKey.h"
#include "TFractionFitter.h"

#include "../../Ntuple-making/Commonglobals.cpp"
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

	cout<<"****************************** READING DATABASES *************************************"<<endl;
	
	string filename2="./database_P.root";
        TFile * file2 = TFile::Open(filename2.c_str(),"READ");
        
        string filename3="./database_D.root";
        TFile * file3 = TFile::Open(filename3.c_str(),"READ");

        string filename4="./database_PD.root";
        TFile * file4 = TFile::Open(filename4.c_str(),"READ");

	std::vector<TGraphAsymmErrors *> D_Graphs;
	std::vector<TGraphAsymmErrors *> PD_Graphs;

        TList *ExperimentsD = file3->GetListOfKeys();
        TIter nextD(ExperimentsD);
        TKey * keyD;
	TObject * obj;

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
	
	Flux * DFluxTOF = new Flux(finalHistos, "DFluxTOF", "FullsetEff_D_TOF","FullsetEfficiency","TOFfits/Fit Results/Primary Deuteron Counts","ExposureTOF","Gen. Acceptance",ToFDB);
	Flux * DFluxNaF = new Flux(finalHistos, "DFluxNaF", "FullsetEff_D_NaF","FullsetEfficiency","NaFfits/Fit Results/Primary Deuteron Counts","ExposureNaF","Gen. Acceptance",NaFDB);
	Flux * DFluxAgl = new Flux(finalHistos, "DFluxAgl", "FullsetEff_D_Agl","FullsetEfficiency","Aglfits/Fit Results/Primary Deuteron Counts","ExposureAgl","Gen. Acceptance",AglDB);

	TCanvas *c1 = new TCanvas("Gen. Acceptance");
	c1->SetCanvasSize(2000,1500);

	PlotTH1FintoGraph(gPad,ToFDB, DFluxTOF->GetGenAcceptance(),"Kinetic Energy [GeV/nucl.]", "Gen. Acceptance [m^{2} sr]",1,true,"Psame",0.1,10,0.1,3,"TOF range",8);
	PlotTH1FintoGraph(gPad,NaFDB, DFluxNaF->GetGenAcceptance(),"Kinetic Energy [GeV/nucl.]", "Gen. Acceptance [m^{2} sr]",1,true,"Psame",0.1,10,0.1,3,"NaF range",22);
	PlotTH1FintoGraph(gPad,AglDB, DFluxAgl->GetGenAcceptance(),"Kinetic Energy [GeV/nucl.]", "Gen. Acceptance [m^{2} sr]",1,true,"Psame",0.1,10,0.1,3,"Agl range",29);

	Plots.Add(c1);
	Plots.writeObjsInFolder("Fluxes");


	float potenza = 0;
	TCanvas *c2 = new TCanvas("Deuteron Primary Flux");
	c2->SetCanvasSize(2000,1500);

	TH2F * Frame = CreateFrame(0.01,25,0.01,1000,"Kin.En./nucl. [GeV/nucl.]","Flux [(m^2 sr GeV/nucl.)^{-1}]");	
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

	DFluxNaF->GetFlux()->Smooth();        

	PlotTH1FintoGraph(gPad,ToFDB, DFluxTOF->GetFlux(),"Kinetic Energy [GeV/nucl.]", "Flux",1,true,"Psame",0.1,10,0.01,1000,"This Work (TOF)",8);
	PlotTH1FintoGraph(gPad,NaFDB, DFluxNaF->GetFlux(),"Kinetic Energy [GeV/nucl.]", "Flux",1,true,"Psame",0.1,10,0.01,1000,"This Work (NaF)",22);
	PlotTH1FintoGraph(gPad,AglDB, DFluxAgl->GetFlux(),"Kinetic Energy [GeV/nucl.]", "Flux",1,true,"Psame",0.1,10,0.01,1000,"This Work (Agl)",29);

	
	Plots.Add(c2);
	Plots.writeObjsInFolder("Fluxes");

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

	TH1F * DPRatioTOF = (TH1F *)finalHistos.Get("Fluxes/DP ratio TOF");
	TH1F * DPRatioNaF = (TH1F *)finalHistos.Get("Fluxes/DP ratio NaF");
	TH1F * DPRatioAgl = (TH1F *)finalHistos.Get("Fluxes/DP ratio Agl");

	DPRatioNaF->Smooth();	

	cout<<DPRatioTOF<<endl;

	PlotTH1FintoGraph(gPad,ToFDB, DPRatioTOF ,"Kinetic Energy [GeV/nucl.]", "Flux",1,true,"Psame",0.1,10,0.00001,0.07,"This Work (TOF)",8);
	PlotTH1FintoGraph(gPad,NaFDB, DPRatioNaF ,"Kinetic Energy [GeV/nucl.]", "Flux",1,true,"Psame",0.1,10,0.00001,0.07,"This Work (NaF)",22);
	PlotTH1FintoGraph(gPad,AglDB, DPRatioAgl ,"Kinetic Energy [GeV/nucl.]", "Flux",1,true,"Psame",0.1,10,0.00001,0.07,"This Work (Agl)",29);

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


	return 0;
}
