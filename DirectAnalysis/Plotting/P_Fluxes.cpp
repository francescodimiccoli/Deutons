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

TSpline3 * GetFluxSpline(TGraphAsymmErrors * Graph);

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
        ToFPB.Print();

        cout<<"**NaF**"<<endl;
        NaFPB.Print();

        cout<<"**Agl**"<<endl;
        AglPB.Print();

        ToFPB.UseBetaEdges();
        NaFPB.UseBetaEdges();
        AglPB.UseBetaEdges();

        PRB.UseREdges();


        cout<<endl;

	cout<<"****************************** READING DATABASES *************************************"<<endl;
	
	string filename2="./database_P.root";
        TFile * file2 = TFile::Open(filename2.c_str(),"READ");
        
        string filename3="./database_D.root";
        TFile * file3 = TFile::Open(filename3.c_str(),"READ");

        string filename4="./database_PD.root";
        TFile * file4 = TFile::Open(filename4.c_str(),"READ");

	std::vector<TGraphAsymmErrors *> P_Graphs;

        TList *ExperimentsP = file2->GetListOfKeys();
        TIter nextP(ExperimentsP);
        TKey * keyP;
	TObject * obj;

        while((keyP = (TKey*)nextP())){
                obj = file2->Get(keyP->GetName());
                if(obj->InheritsFrom("TGraphAsymmErrors")) P_Graphs.push_back((TGraphAsymmErrors *)obj);
        }

        cout<<"****************************** PLOTTING FLUXES ***************************************"<<endl;
	
	Flux * HEPFlux  = new Flux(finalHistos,"PFluxHE", "RigBinFullsetEff","RigBinFullsetEff","HEPCounts/HEPCounts/HEPCounts","HEExposure","Gen. Acceptance",PRB);
	
	Flux * PFluxTOF = new Flux(finalHistos, "PFluxTOF", "FullsetEff_P_TOF","FullsetEfficiency","TOFfits/Fit Results/Primary Proton Counts","ExposureTOF","Gen. Acceptance",ToFPB);
	Flux * PFluxNaF = new Flux(finalHistos, "PFluxNaF", "FullsetEff_P_NaF","FullsetEfficiency","NaFfits/Fit Results/Primary Proton Counts","ExposureNaF","Gen. Acceptance",NaFPB);
	Flux * PFluxAgl = new Flux(finalHistos, "PFluxAgl", "FullsetEff_P_Agl","FullsetEfficiency","Aglfits/Fit Results/Primary Proton Counts","ExposureAgl","Gen. Acceptance",AglPB);

	TCanvas *c1 = new TCanvas("Gen. Acceptance (P)");
	c1->SetCanvasSize(2000,1500);

	PlotTH1FintoGraph(gPad,PRB,   HEPFlux ->GetGenAcceptance(),"Kinetic Energy [GeV/nucl.]", "Gen. Acceptance [m^{2} sr]",2,true,"Psame",0.1,100,0.1,3,"H.E. range",24);
	PlotTH1FintoGraph(gPad,ToFDB, PFluxTOF->GetGenAcceptance(),"Kinetic Energy [GeV/nucl.]", "Gen. Acceptance [m^{2} sr]",2,true,"Psame",0.1,10,0.1,3,"TOF range",8);
	PlotTH1FintoGraph(gPad,NaFDB, PFluxNaF->GetGenAcceptance(),"Kinetic Energy [GeV/nucl.]", "Gen. Acceptance [m^{2} sr]",2,true,"Psame",0.1,10,0.1,3,"NaF range",22);
	PlotTH1FintoGraph(gPad,AglDB, PFluxAgl->GetGenAcceptance(),"Kinetic Energy [GeV/nucl.]", "Gen. Acceptance [m^{2} sr]",2,true,"Psame",0.1,10,0.1,3,"Agl range",29);

	Plots.Add(c1);
	Plots.writeObjsInFolder("Fluxes");


	float potenza = 0;
	TCanvas *c2 = new TCanvas("Proton Primary Flux");
	c2->SetCanvasSize(2000,1500);

	TH2F * Frame = CreateFrame(0.01,100,1e-3,1e7,"Kin.En./nucl. [GeV/nucl.]","Flux [(m^2 sr GeV/nucl.)^{-1}]");	
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
                string filename="./Galprop/Tom/prot_1500.dat";
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
                string filename="./Galprop/Tom/prot_100.dat";
                cout<<filename<<endl;
                ifstream fp(filename.c_str());
                while (!fp.eof()){
                        fp>>x>>y;
                        if(x/1e3>0.05&&x/1e3<=100)
                                galprop3P2->SetPoint(j,x/1e3,y*1e7*pow(x/1e3,potenza));
                        j++;
                }
        }
        galprop3P->GetXaxis()->SetRangeUser(0.1,100);
        galprop3P->GetYaxis()->SetRangeUser(1e-3,1e4);

        galprop3P->SetTitle("Protons Flux: Geo. Zones");
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
        for(uint n=0;n<P_Graphs.size();n++){
                P_Graphs[n] ->Draw("Psame");
                P_Graphs[n]->SetMarkerSize(2); 
                leg->AddEntry(P_Graphs[n],P_Graphs[n]->GetTitle(),"p");
        }

	leg->SetFillColor(0);
	leg->SetLineWidth(2);
	leg->Draw("same");
        
	PlotTH1FintoGraph(gPad,PRB,   HEPFlux->GetFlux(), "Kinetic Energy [GeV/nucl.]", "Flux",1,true,"Psame",0.1,100,1e-3,10000,"This Work (H.E.)",24);
	PlotTH1FintoGraph(gPad,ToFPB, PFluxTOF->GetFlux(),"Kinetic Energy [GeV/nucl.]", "Flux",1,true,"Psame",0.1,100,1e-3,10000,"This Work (TOF)",8);
	PlotTH1FintoGraph(gPad,NaFPB, PFluxNaF->GetFlux(),"Kinetic Energy [GeV/nucl.]", "Flux",1,true,"Psame",0.1,100,1e-3,10000,"This Work (NaF)",22);
	PlotTH1FintoGraph(gPad,AglPB, PFluxAgl->GetFlux(),"Kinetic Energy [GeV/nucl.]", "Flux",1,true,"Psame",0.1,100,1e-3,10000,"This Work (Agl)",29);

	
	Plots.Add(c2);
	Plots.writeObjsInFolder("Fluxes");


	TCanvas *c3 = new TCanvas("P Flux Comparison with AMS-02");
	c3->SetCanvasSize(2000,1500);
	c3->cd();
        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetGridx();
        gPad->SetGridy();


	TSpline3 *AMSFlux = GetFluxSpline(P_Graphs[1]);
	AMSFlux->SetLineColor(2);
	AMSFlux->SetLineWidth(2);
	//PlotTH1FintoGraph(gPad,PRB,   HEPFlux->GetFlux(), "Kinetic Energy [GeV/nucl.]", "Flux",1,true,"Psame",0.1,100,1e-3,10000,"This Work (H.E.)",24);
	PlotRatioWithSplineintoGraph(gPad,PRB, HEPFlux->GetFlux(),AMSFlux, "Kinetic Energy [GeV/nucl.]", "Flux",1,true,"Psame",0.1,100,0.1,2,"This Work (H.E.)",24);
	PlotRatioWithSplineintoGraph(gPad,ToFDB, PFluxTOF->GetFlux(),AMSFlux, "Kinetic Energy [GeV/nucl.]", "Flux",1,true,"Psame",0.1,100,0.1,2,"This Work (TOF)",8);
	PlotRatioWithSplineintoGraph(gPad,NaFDB, PFluxNaF->GetFlux(),AMSFlux, "Kinetic Energy [GeV/nucl.]", "Flux",1,true,"Psame",0.1,100,0.1,2,"This Work (NaF)",22);
	PlotRatioWithSplineintoGraph(gPad,AglDB, PFluxAgl->GetFlux(),AMSFlux, "Kinetic Energy [GeV/nucl.]", "Flux",1,true,"Psame",0.1,100,0.1,2,"This Work (Agl)",29);
	
	//AMSFlux->Draw("same");
	 
        Plots.Add(c3);
	Plots.writeObjsInFolder("Fluxes");


	return 0;
}

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
