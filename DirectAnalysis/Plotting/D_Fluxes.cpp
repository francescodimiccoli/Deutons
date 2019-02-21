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

        Flux * HEPFlux  = new Flux(finalHistos,"PFluxHE", "RigBinFullsetEff","RigBinFullsetEff","HEPCounts/HEPCounts/HEPCounts","HEExposure",PRB);
	
	Flux * DFluxTOF = new Flux(finalHistos, "DFluxTOF", "FullsetEff_D_TOF","FullsetEfficiency","TOFfits/Fit Results/Primary Deuteron Counts","ExposureTOF",ToFDB);
	Flux * DFluxNaF = new Flux(finalHistos, "DFluxNaF", "FullsetEff_D_NaF","FullsetEfficiency","NaFfits/Fit Results/Primary Deuteron Counts","ExposureNaF",NaFDB);
	Flux * DFluxAgl = new Flux(finalHistos, "DFluxAgl", "FullsetEff_D_Agl","FullsetEfficiency","Aglfits/Fit Results/Primary Deuteron Counts","ExposureAgl",AglDB);

	Flux * PFluxTOF = new Flux(finalHistos, "PFluxTOF", "FullsetEff_P_TOF","FullsetEfficiency","TOFfits/Fit Results/Primary Proton Counts","ExposureTOF",ToFPB);
	Flux * PFluxNaF = new Flux(finalHistos, "PFluxNaF", "FullsetEff_P_NaF","FullsetEfficiency","NaFfits/Fit Results/Primary Proton Counts","ExposureNaF",NaFPB);
	Flux * PFluxAgl = new Flux(finalHistos, "PFluxAgl", "FullsetEff_P_Agl","FullsetEfficiency","Aglfits/Fit Results/Primary Proton Counts","ExposureAgl",AglPB);

	Flux * DummyDTOF = new Flux(finalHistos, "DummyDTOF", "Baseline_D_TOF","Baseline","TOFfits/Fit Results/Primary Deuteron Counts","ExposureTOF",ToFDB);
	Flux * DummyDNaF = new Flux(finalHistos, "DummyDNaF", "Baseline_D_NaF","Baseline","NaFfits/Fit Results/Primary Deuteron Counts","ExposureNaF",NaFDB);
	Flux * DummyDAgl = new Flux(finalHistos, "DummyDAgl", "Baseline_D_Agl","Baseline","Aglfits/Fit Results/Primary Deuteron Counts","ExposureAgl",AglDB);
                                                                        
	Flux * DummyPTOF = new Flux(finalHistos,"DummyPTOF", "Baseline_P_TOF","Baseline","TOFfits/Fit Results/Primary Proton Counts","ExposureTOF",ToFPB);
	Flux * DummyPNaF = new Flux(finalHistos,"DummyPNaF", "Baseline_P_NaF","Baseline","NaFfits/Fit Results/Primary Proton Counts","ExposureNaF",NaFPB);
	Flux * DummyPAgl = new Flux(finalHistos,"DummyPAgl", "Baseline_P_Agl","Baseline","Aglfits/Fit Results/Primary Proton Counts","ExposureAgl",AglPB);



	TCanvas *c2_ = new TCanvas("Eff. Acceptance");
	c2_->SetCanvasSize(2000,1500);

	PlotTH1FintoGraph(gPad,ToFDB, DFluxTOF->GetEffAcceptance(),"Kinetic Energy [GeV/nucl.]", "Eff. Acceptance [m^{2} sr]",4,true,"Lsame",0.1,10,0.0001,0.3,"Deuterons",8);
	PlotTH1FintoGraph(gPad,ToFPB, PFluxTOF->GetEffAcceptance(),"Kinetic Energy [GeV/nucl.]", "Eff. Acceptance [m^{2} sr]",2,true,"Lsame",0.1,10,0.0001,0.3,"Protons",8);

	PlotTH1FintoGraph(gPad,ToFDB, DFluxTOF->GetEffAcceptance(),"Kinetic Energy [GeV/nucl.]", "Eff. Acceptance [m^{2} sr]",4,true,"Lsame",0.1,10,0.0001,0.6,"TOF range",8,true);
	PlotTH1FintoGraph(gPad,NaFDB, DFluxNaF->GetEffAcceptance(),"Kinetic Energy [GeV/nucl.]", "Eff. Acceptance [m^{2} sr]",4,true,"Lsame",0.1,10,0.0001,0.6,"NaF range",22,true);
	PlotTH1FintoGraph(gPad,AglDB, DFluxAgl->GetEffAcceptance(),"Kinetic Energy [GeV/nucl.]", "Eff. Acceptance [m^{2} sr]",4,true,"Lsame",0.1,10,0.0001,0.6,"Agl range",29,true);

	PlotTH1FintoGraph(gPad,ToFPB, PFluxTOF->GetEffAcceptance(),"Kinetic Energy [GeV/nucl.]", "Eff. Acceptance [m^{2} sr]",2,true,"Lsame",0.1,10,0.0001,0.6,"TOF range",8,true);
        PlotTH1FintoGraph(gPad,NaFPB, PFluxNaF->GetEffAcceptance(),"Kinetic Energy [GeV/nucl.]", "Eff. Acceptance [m^{2} sr]",2,true,"Lsame",0.1,10,0.0001,0.6,"NaF range",22,true);
        PlotTH1FintoGraph(gPad,AglPB, PFluxAgl->GetEffAcceptance(),"Kinetic Energy [GeV/nucl.]", "Eff. Acceptance [m^{2} sr]",2,true,"Lsame",0.1,10,0.0001,0.6,"Agl range",29,true);

/*	PlotTH1FintoGraph(gPad,ToFDB, DummyDTOF->GetEffAcceptance(),"Kinetic Energy [GeV/nucl.]", "Eff. Acceptance [m^{2} sr]",4,true,"Lsame",0.1,10,0.0001,0.6,"TOF range",8,true);
	PlotTH1FintoGraph(gPad,NaFDB, DummyDNaF->GetEffAcceptance(),"Kinetic Energy [GeV/nucl.]", "Eff. Acceptance [m^{2} sr]",4,true,"Lsame",0.1,10,0.0001,0.6,"NaF range",22,true);
	PlotTH1FintoGraph(gPad,AglDB, DummyDAgl->GetEffAcceptance(),"Kinetic Energy [GeV/nucl.]", "Eff. Acceptance [m^{2} sr]",4,true,"Lsame",0.1,10,0.0001,0.6,"Agl range",29,true);

	PlotTH1FintoGraph(gPad,ToFPB, DummyPTOF->GetEffAcceptance(),"Kinetic Energy [GeV/nucl.]", "Eff. Acceptance [m^{2} sr]",2,true,"Lsame",0.1,10,0.0001,0.6,"TOF range",8,true);
        PlotTH1FintoGraph(gPad,NaFPB, DummyPNaF->GetEffAcceptance(),"Kinetic Energy [GeV/nucl.]", "Eff. Acceptance [m^{2} sr]",2,true,"Lsame",0.1,10,0.0001,0.6,"NaF range",22,true);
        PlotTH1FintoGraph(gPad,AglPB, DummyPAgl->GetEffAcceptance(),"Kinetic Energy [GeV/nucl.]", "Eff. Acceptance [m^{2} sr]",2,true,"Lsame",0.1,10,0.0001,0.6,"Agl range",29,true);
*/

	Plots.Add(c2_);
	Plots.writeObjsInFolder("Fluxes");

/*
	TCanvas * c6 = new TCanvas("Acceptance Uncertainty Break-down");
        c6->SetCanvasSize(5000,1000);
        c6->Divide(3,1);

	TH1F * TotErrTOF = (TH1F *) DFluxTOF->GetEffAcceptance()->Clone();
	TH1F * TotErrNaF = (TH1F *) DFluxNaF->GetEffAcceptance()->Clone();
	TH1F * TotErrAgl = (TH1F *) DFluxAgl->GetEffAcceptance()->Clone();
	
	for(int i =0; i< TotErrTOF->GetNbinsX();i++) {
		if(TotErrTOF->GetBinContent(i+1)>0&&TotErrTOF->GetBinError(i+1)>0){
		TotErrTOF->SetBinContent(i+1, TotErrTOF->GetBinError(i+1)/TotErrTOF->GetBinContent(i+1));
		TotErrTOF->SetBinError(i+1,0);
		}
	}
	for(int i =0; i< TotErrNaF->GetNbinsX();i++) {
		if(TotErrNaF->GetBinContent(i+1)>0&&TotErrNaF->GetBinError(i+1)>0){
		TotErrNaF->SetBinContent(i+1, TotErrNaF->GetBinError(i+1)/TotErrNaF->GetBinContent(i+1));
		TotErrNaF->SetBinError(i+1,0);
		}
	}
	for(int i =0; i< TotErrAgl->GetNbinsX();i++) {
		if(TotErrAgl->GetBinContent(i+1)>0&&TotErrAgl->GetBinError(i+1)>0){
		TotErrAgl->SetBinContent(i+1, TotErrAgl->GetBinError(i+1)/TotErrAgl->GetBinContent(i+1));
		TotErrAgl->SetBinError(i+1,0);
		}
	}

	c6->cd(1);
	PlotDistribution(gPad,TotErrTOF,"TOF Range Bin","Relative error",2,"same",1e-4,1.1,10,"Total Error");
	PlotDistribution(gPad,DFluxTOF->GetAcc_StatErr(),"TOF Range Bin","Relative error",4,"same",1e-4,1.1,4,"Stat. Error",true);
	PlotDistribution(gPad,DFluxTOF->GetAcc_SystErr(),"TOF Range Bin","Relative error",1,"same",1e-4,1.1,4,"Syst. D/P",true);
	c6->cd(2);
	PlotDistribution(gPad,TotErrNaF,"NaF Range Bin","Relative error",2,"same",1e-4,1.1,10,"Total Error");
	PlotDistribution(gPad,DFluxNaF->GetAcc_StatErr(),"NaF Range Bin","Relative error",4,"same",1e-4,1.1,4,"Stat. Error",true);
	PlotDistribution(gPad,DFluxNaF->GetAcc_SystErr(),"NaF Range Bin","Relative error",1,"same",1e-4,1.1,4,"Syst. D/P",true);
	c6->cd(3);
	PlotDistribution(gPad,TotErrAgl,"Agl Range Bin","Relative error",2,"same",1e-4,1.1,10,"Total Error");
	PlotDistribution(gPad,DFluxAgl->GetAcc_StatErr(),"Agl Range Bin","Relative error",4,"same",1e-4,1.1,4,"Stat. Error",true);
	PlotDistribution(gPad,DFluxAgl->GetAcc_SystErr(),"Agl Range Bin","Relative error",1,"same",1e-4,1.1,4,"Syst. D/P",true);

	
	Plots.Add(c6);
        Plots.writeObjsInFolder("Fluxes");	

*/
	
	TCanvas *c1_ = new TCanvas("Exposure Time");
	c1_->SetCanvasSize(3000,1500);



	PlotTH1FintoGraph(gPad,PRB, HEPFlux->GetExposureTime(),"Kinetic Energy [GeV/nucl.]", "Exposure Time [sec]",1,true,"Psame",0.1,50,1,2*   HEPFlux->GetExposureTime()->GetBinContent(30),"HE range",8);
/*	
	PlotTH1FintoGraph(gPad,ToFDB, DFluxTOF->GetExposureTime(),"Kinetic Energy [GeV/nucl.]", "Exposure Time [sec]",1,true,"Psame",0.1,50,1,2*HEPFlux->GetExposureTime()->GetBinContent(30),"TOF range",8);
	PlotTH1FintoGraph(gPad,NaFDB, DFluxNaF->GetExposureTime(),"Kinetic Energy [GeV/nucl.]", "Exposure Time [sec]",1,true,"Psame",0.1,50,1,2*HEPFlux->GetExposureTime()->GetBinContent(30),"NaF range",22);
	PlotTH1FintoGraph(gPad,AglDB, DFluxAgl->GetExposureTime(),"Kinetic Energy [GeV/nucl.]", "Exposure Time [sec]",1,true,"Psame",0.1,50,1,2*HEPFlux->GetExposureTime()->GetBinContent(30),"Agl range",29);

	PlotTH1FintoGraph(gPad,ToFDB, DFluxTOF->GetExposureTime(),"Kinetic Energy [GeV/nucl.]", "Exposure Time [sec]",4,true,"Psame",0.1,50,1,2*HEPFlux->GetExposureTime()->GetBinContent(30),"Deuterons",8);
	PlotTH1FintoGraph(gPad,ToFPB, PFluxTOF->GetExposureTime(),"Kinetic Energy [GeV/nucl.]", "Exposure Time [sec]",2,true,"Psame",0.1,50,1,2*HEPFlux->GetExposureTime()->GetBinContent(30),"Protons",8);

	PlotTH1FintoGraph(gPad,ToFDB, DFluxTOF->GetExposureTime(),"Kinetic Energy [GeV/nucl.]", "Exposure Time [sec]",4,true,"Psame",0.1,50,1,2*HEPFlux->GetExposureTime()->GetBinContent(30),"TOF range",8,true);
	PlotTH1FintoGraph(gPad,NaFDB, DFluxNaF->GetExposureTime(),"Kinetic Energy [GeV/nucl.]", "Exposure Time [sec]",4,true,"Psame",0.1,50,1,2*HEPFlux->GetExposureTime()->GetBinContent(30),"NaF range",22,true);
	PlotTH1FintoGraph(gPad,AglDB, DFluxAgl->GetExposureTime(),"Kinetic Energy [GeV/nucl.]", "Exposure Time [sec]",4,true,"Psame",0.1,50,1,2*HEPFlux->GetExposureTime()->GetBinContent(30),"Agl range",29,true);

	PlotTH1FintoGraph(gPad,ToFPB, PFluxTOF->GetExposureTime(),"Kinetic Energy [GeV/nucl.]", "Exposure Time [sec]",2,true,"Psame",0.1,50,1,2*HEPFlux->GetExposureTime()->GetBinContent(30),"TOF range",8,true);
	PlotTH1FintoGraph(gPad,NaFPB, PFluxNaF->GetExposureTime(),"Kinetic Energy [GeV/nucl.]", "Exposure Time [sec]",2,true,"Psame",0.1,50,1,2*HEPFlux->GetExposureTime()->GetBinContent(30),"NaF range",22,true);
	PlotTH1FintoGraph(gPad,AglPB, PFluxAgl->GetExposureTime(),"Kinetic Energy [GeV/nucl.]", "Exposure Time [sec]",2,true,"Psame",0.1,50,1,2*HEPFlux->GetExposureTime()->GetBinContent(30),"Agl range",29,true);
*/
	Plots.Add(c1_);
	Plots.writeObjsInFolder("Fluxes");



	float potenza = 0;
	TCanvas *c2 = new TCanvas("Deuteron Primary Flux");
	c2->SetCanvasSize(2000,1500);

	TH2F * Frame = CreateFrame(gPad,0.01,25,0.01,1000,"Kin.En./nucl. [GeV/nucl.]","Flux [(m^2 sr GeV/nucl.)^{-1}]");	
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

	PlotTH1FintoGraph(gPad,ToFDB, DPRatioTOF ,"Kinetic Energy [GeV/nucl.]", "Flux",1,true,"Psame",0.1,10,0.00001,0.12,"This Work (TOF)",8);
	PlotTH1FintoGraph(gPad,NaFDB, DPRatioNaF ,"Kinetic Energy [GeV/nucl.]", "Flux",1,true,"Psame",0.1,10,0.00001,0.12,"This Work (NaF)",22);
	PlotTH1FintoGraph(gPad,AglDB, DPRatioAgl ,"Kinetic Energy [GeV/nucl.]", "Flux",1,true,"Psame",0.1,10,0.00001,0.12,"This Work (Agl)",29);

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


	return 0;
}
