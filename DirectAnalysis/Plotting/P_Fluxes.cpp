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
#include "TStyle.h"

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
TSpline3 * GetFluxSpline(TH1F * Graph, Binning bins);
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
	SetUpTOIBinning();
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
	
	Flux * FluxTests[10];

        for(int i=0;i<10;i++) FluxTests[i] = new Flux(finalHistos,("PFluxTests_v"+to_string(i)).c_str(),("CountsTestEff_v" + to_string(i)).c_str(),("CountsTestEff_v" + to_string(i)).c_str(),("HEPTests_v"+to_string(i) + "/HEPTests_v"+to_string(i)+"/HEPTests_v"+to_string(i)+"_before").c_str(),"HEExposure",PRB);

	Flux * HEPFlux  = new Flux(finalHistos,"PFluxHE", "RigBinFullsetEff_Trig","RigBinFullsetEff_Trig","HEPCounts/HEPCounts/HEPCounts_before","HEExposure",PRB);
	Flux * HEPFluxQ  = new Flux(finalHistos,"PFluxQHE", "RigBinQualEff","RigBinQualEff","HEPCountsQual/HEPCountsQual/HEPCountsQual_before","HEExposure",PRB);
		Flux * HEPFluxL1  = new Flux(finalHistos,"PFluxL1HE", "RigBinFullsetEffL1_Trig","RigBinFullsetEffL1_Trig","HEPCountsL1/HEPCountsL1/HEPCountsL1_before","HEExposure",PRB);
	
	Flux * PFluxTOF = new Flux(finalHistos, "PFluxTOF", "FullsetEff_P_TOF","FullsetEfficiency","TOFfits/Fit Results/Primary Proton Counts","ExposureTOF",ToFPB);
	Flux * PFluxNaF = new Flux(finalHistos, "PFluxNaF", "FullsetEff_P_NaF","FullsetEfficiency","NaFfits/Fit Results/Primary Proton Counts","ExposureNaF",NaFPB);
	Flux * PFluxAgl = new Flux(finalHistos, "PFluxAgl", "FullsetEff_P_Agl","FullsetEfficiency","Aglfits/Fit Results/Primary Proton Counts","ExposureAgl",AglPB);
	Flux * DummyPTOF = new Flux(finalHistos,"DummyPTOF", "Baseline_P_TOF","Baseline","TOFPCounts/TOFPCounts/TOFPCounts","ExposureTOF",ToFPB);
	Flux * DummyPNaF = new Flux(finalHistos,"DummyPNaF", "Baseline_P_NaF","Baseline","NaFPCounts/NaFPCounts/NaFPCounts","ExposureNaF",NaFPB);
	Flux * DummyPAgl = new Flux(finalHistos,"DummyPAgl", "Baseline_P_Agl","Baseline","AglPCounts/AglPCounts/AglPCounts","ExposureAgl",AglPB);


	TCanvas *c1 = new TCanvas("Effective Acceptance (P)");
	c1->SetCanvasSize(2000,1500);
	c1->cd();
	//cout<<"ACC ENTRIES: "<<HEPFlux  ->GetEffAcceptance()->GetBinContent(HEPFlux  ->GetEffAcceptance()->GetMaximumBin())<<endl;
//	HEPFlux  ->GetEffAcceptance()->Draw();
//	HEPFluxQ ->GetEffAcceptance()->Draw("same");
	PlotTH1FintoGraph(gPad,PRB,   HEPFlux->GetEffAcceptance(), "Kinetic Energy [GeV/nucl.]", "Flux",1,true,"Psame",0.1,50,1e-3,3,"Baseline",8);	
	PlotTH1FintoGraph(gPad,PRB,   HEPFluxL1->GetEffAcceptance(), "Kinetic Energy [GeV/nucl.]", "Flux",4,true,"Psame",0.1,50,1e-3,3,"Layer 1",8);	
	PlotTH1FintoGraph(gPad,PRB,   HEPFluxQ->GetEffAcceptance(), "Kinetic Energy [GeV/nucl.]", "Flux",2,true,"Psame",0.1,50,1e-3,3,"Interactions",8);	


	Plots.Add(c1);
	Plots.writeObjsInFolder("Fluxes");


	float potenza = 0;
	TCanvas *c2 = new TCanvas("Proton Primary Flux");
	c2->SetCanvasSize(2000,1500);

	TH2F * Frame = CreateFrame(gPad,0.01,100,1e-3,104*HEPFlux->GetFlux()->GetBinContent(HEPFlux->GetFlux()->GetMaximumBin()),"Kin.En./nucl. [GeV/nucl.]","Flux [(m^2 sr GeV/nucl.)^{-1}]");	
        c2->cd();
        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetGridx();
        gPad->SetGridy();



	Frame->Draw();
	
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
        galprop3P->GetXaxis()->SetRangeUser(0.1,50);
        galprop3P->GetYaxis()->SetRangeUser(1e-3,1e4);

        galprop3P->SetTitle("Protons Flux: Geo. Zones");
        galprop3P->GetXaxis()->SetTitle("Kin.En./nucl. [GeV/nucl.]");
        galprop3P ->GetYaxis()->SetTitle("Flux [(m^2 sr GeV/nucl.)^-1]");
        galprop3P ->GetXaxis()->SetTitleSize(0.045);
        galprop3P->GetYaxis()->SetTitleSize(0.045);
        galprop3P ->GetYaxis()->SetRangeUser(1e-2,1e4);

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
       
	PlotTH1FintoGraph(gPad,PRB,   HEPFlux->GetFlux(), "Kinetic Energy [GeV/nucl.]", "Flux",1,true,"Psame",0.1,50,1e-3,104*HEPFlux->GetFlux()->GetBinContent(HEPFlux->GetFlux()->GetMaximumBin()),"This Work (H.E.)",24);
	PlotTH1FintoGraph(gPad,PRB,   HEPFluxQ->GetFlux(), "Kinetic Energy [GeV/nucl.]", "Flux",2,true,"Psame",0.1,50,1e-3,104*HEPFlux->GetFlux()->GetBinContent(HEPFlux->GetFlux()->GetMaximumBin()),"This Work (H.E.)",24);
	
	PlotTH1FintoGraph(gPad,ToFPB, DummyPTOF->GetFlux(),"Kinetic Energy [GeV/nucl.]", "Flux",3,true,"Psame",0.1,50,1e-3,104*HEPFlux->GetFlux()->GetBinContent(HEPFlux->GetFlux()->GetMaximumBin()),"This Work (TOF)",8);
	PlotTH1FintoGraph(gPad,NaFPB, DummyPNaF->GetFlux(),"Kinetic Energy [GeV/nucl.]", "Flux",3,true,"Psame",0.1,50,1e-3,104*HEPFlux->GetFlux()->GetBinContent(HEPFlux->GetFlux()->GetMaximumBin()),"This Work (NaF)",22);
	PlotTH1FintoGraph(gPad,AglPB, DummyPAgl->GetFlux(),"Kinetic Energy [GeV/nucl.]", "Flux",3,true,"Psame",0.1,50,1e-3,104*HEPFlux->GetFlux()->GetBinContent(HEPFlux->GetFlux()->GetMaximumBin()),"This Work (Agl)",29);
/*
	PlotTH1FintoGraph(gPad,ToFPB, PFluxTOF->GetFlux(),"Kinetic Energy [GeV/nucl.]", "Flux",4,true,"Psame",0.1,50,1e-3,104*HEPFlux->GetFlux()->GetBinContent(HEPFlux->GetFlux()->GetMaximumBin()),"This Work (TOF)",8);
	PlotTH1FintoGraph(gPad,NaFPB, PFluxNaF->GetFlux(),"Kinetic Energy [GeV/nucl.]", "Flux",4,true,"Psame",0.1,50,1e-3,104*HEPFlux->GetFlux()->GetBinContent(HEPFlux->GetFlux()->GetMaximumBin()),"This Work (NaF)",22);
	PlotTH1FintoGraph(gPad,AglPB, PFluxAgl->GetFlux(),"Kinetic Energy [GeV/nucl.]", "Flux",4,true,"Psame",0.1,50,1e-3,104*HEPFlux->GetFlux()->GetBinContent(HEPFlux->GetFlux()->GetMaximumBin()),"This Work (Agl)",29);
*/
	
	
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

	//PlotRatioWithSplineintoGraph(gPad,PRB,    HEPFluxQ->GetFlux() ,AMSFlux, "Kinetic Energy [GeV/nucl.]", "Flux",2,true,"Psame",0.1,50,1e-5,120,"This Work (H.E. - Interactions)",8);
	//PlotRatioWithSplineintoGraph(gPad,PRB,    HEPFlux->GetFlux()  ,AMSFlux, "Kinetic Energy [GeV/nucl.]", "Flux",1,true,"4Psame",0.1,50,1e-5,120,"This Work (H.E.- Baseline)",8);
//	PlotRatioWithSplineintoGraph(gPad,ToFPB,  DummyPTOF->GetFlux(),AMSFlux, "Kinetic Energy [GeV/nucl.]", "Flux",3,true,"Psame",0.1,50,0.1,2,"This Work (TOF)",24);
//	PlotRatioWithSplineintoGraph(gPad,NaFPB,  DummyPNaF->GetFlux(),AMSFlux, "Kinetic Energy [GeV/nucl.]", "Flux",3,true,"Psame",0.1,50,0.1,2,"This Work (NaF)",22);
//	PlotRatioWithSplineintoGraph(gPad,AglPB,  DummyPAgl->GetFlux(),AMSFlux, "Kinetic Energy [GeV/nucl.]", "Flux",3,true,"Psame",0.1,50,0.1,2,"This Work (Agl)",29);
	

	
	TStyle *st1 = new TStyle("st1","my style");
	st1->cd();
	st1->SetPalette(55);
	int nColors = st1->GetNumberOfColors();
	int nHistos = 10;
	std::string nameplots[10]={"ultra-basic","tracker","+4/4 tof","+beta>0.3","+qU-tof>0","+0.5<qL_tof<2","+tofch1<100","+chitime<10","+1track","+tofclusters"};
        for(int i=0;i<10;i++) {
		int histocolor=(float)nColors / nHistos * i;
		PlotRatioWithSplineintoGraph(gPad,PRB,    FluxTests[i]->GetFlux()  ,AMSFlux, "Kinetic Energy [GeV/n]", "Flux (ratio against published)",st1->GetColorPalette(histocolor),true,"4Psame",0.1,100,0.1,100,nameplots[i],8);
	}
	PlotRatioWithSplineintoGraph(gPad,PRB,    HEPFlux->GetFlux()  ,AMSFlux, "Kinetic Energy [GeV/nucl.]", "Flux",1,true,"4Psame",0.1,50,1e-5,120,"This Work (H.E.- Baseline)",8);
	ErrorPubl->Draw("Psame"); 

        Plots.Add(c3);
	Plots.writeObjsInFolder("Fluxes");

	TCanvas *c3_ = new TCanvas("MY P Flux Comparison");
	c3_->SetCanvasSize(2000,1500);
	c3_->cd();
	gPad->SetLogy();
        PlotRatioWithSplineintoGraph(gPad,PRB,    HEPFlux->GetFlux()  ,AMSFlux, "Kinetic Energy [GeV/nucl.]", "Flux",1,true,"4Psame",0.1,50,1e-5,120,"This Work (H.E.- Baseline)",8);
 	PlotRatioWithSplineintoGraph(gPad,PRB,    HEPFluxL1->GetFlux()  ,AMSFlux, "Kinetic Energy [GeV/nucl.]", "Flux",4,true,"4Psame",0.1,50,1e-5,120,"This Work (H.E.- Layer 1)",8);
        PlotRatioWithSplineintoGraph(gPad,PRB,    HEPFluxQ->GetFlux()  ,AMSFlux, "Kinetic Energy [GeV/nucl.]", "Flux",2,true,"4Psame",0.1,50,1e-5,120,"This Work (H.E.- Interactions)",8);
	ErrorPubl->Draw("Psame"); 


	Plots.Add(c3_);
        Plots.writeObjsInFolder("Fluxes");
	
	TCanvas * c6 = new TCanvas("Uncertainty Break-down");
        c6->SetCanvasSize(5000,1000);

	TH1F * CountErr = (TH1F *) HEPFlux->GetCounts()->Clone();
	TH1F * TotErr = (TH1F *) HEPFlux->GetFlux()->Clone();
	TH1F * AccErr = (TH1F *) HEPFlux->GetEffAcceptance()->Clone();

	if(HEPFlux->GetAcc_StatErr()) {
	TH1F * AccErr_Stat = (TH1F *) HEPFlux->GetAcc_StatErr()->Clone();
	TH1F * AccErr_Syst = (TH1F *) HEPFlux->GetAcc_SystErr()->Clone();
	for(int i =0; i< TotErr->GetNbinsX();i++) {
		if(TotErr->GetBinContent(i+1)>0&&TotErr->GetBinError(i+1)>0){
		TotErr->SetBinContent(i+1, TotErr->GetBinError(i+1)/TotErr->GetBinContent(i+1));
		TotErr->SetBinError(i+1,0);
		AccErr->SetBinContent(i+1, AccErr->GetBinError(i+1)/AccErr->GetBinContent(i+1));
		AccErr->SetBinError(i+1,0);
		CountErr->SetBinContent(i+1, CountErr->GetBinError(i+1)/CountErr->GetBinContent(i+1));
		CountErr->SetBinError(i+1,0);
	
		}
	}
	c6->cd();
	PlotDistribution(gPad,TotErr,"Bin nr.","Relative error",2,"same",1e-4,1.1,10,"Total Error");
	PlotDistribution(gPad,AccErr_Stat,"Bin nr.","Relative error",4,"same",1e-4,1.1,4,"Acc. Stat.",true);
	PlotDistribution(gPad,AccErr_Syst,"Bin nr.","Relative error",3,"same",1e-4,1.1,4,"Acc. Syst.",true);
	PlotDistribution(gPad,CountErr,"Bin nr.","Relative error",1,"same",1e-4,1.1,4,"Count stat.");
	}

       Plots.Add(c6);
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


TSpline3 * GetFluxSpline(TH1F * Graph, Binning bins){

	double X[Graph->GetNbinsX()];
	double Y[Graph->GetNbinsX()];

	for(int i=0; i<Graph->GetNbinsX(); i++){	X[i]=bins.EkPerMassBinCent(i); Y[i]=Graph->GetBinContent(i+1);}
	TSpline3 *AMSFlux = new TSpline3("AMSFlux",X,Y,Graph->GetNbinsX());
	return AMSFlux;	
}
