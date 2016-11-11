#include "TFile.h"
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"
#include "TStyle.h"
#include <TVector3.h>
#include <fstream>
#include <sstream>
#include <vector>
#include "TCanvas.h"
#include "TFrame.h"
#include "TLegend.h"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <cstdlib>
#include <stdio.h>
#include <stdarg.h>
#include <TSpline.h>
#include "TFractionFitter.h"
#include "THStack.h"
#include "TNtuple.h"
#include "TKey.h"
#include "TObject.h"
#include "TGraphAsymmErrors.h"
#include <TGraph.h>
#include "TGraphErrors.h"
#include <cstring>
#include "Mesi.h"
#include "TimePlots.h"

using namespace std;

int colorbase=55;

void CheckFileIntegrity(int i,std::string mesi[]);

std::string ConvertString(std::string mese);

void DrawCalibration(TVirtualPad * c, TSpline3 * Corr[],const string &name,const string &xaxisname,const string &yaxisname){
	c -> cd();
	gPad->SetGridx();
	gPad->SetGridy();
	gStyle->SetPalette(1);
	TLegend* leg =new TLegend(0.91,0.1,1.0,0.9);
	Corr[0]->SetLineWidth(2);
	Corr[0]->SetLineColor(1);
	leg->AddEntry(Corr[1],mesi[1].c_str(), "ep");
	Corr[0]->Draw("");
	TH1F *frame=c->DrawFrame(0.4,0.85,1,1.1,name.c_str());
	frame->GetXaxis()->SetTitle(xaxisname.c_str());
	frame->GetYaxis()->SetTitle(yaxisname.c_str());
	for(int i=1;i<num_mesi;i++) {
		Corr[i]->SetLineWidth(2);
		Corr[i]->SetLineColor(colorbase+2*i);
		Corr[i]->SetMarkerColor(colorbase+2*i);
		Corr[i]->SetMarkerStyle(8);
		leg->AddEntry(Corr[i],ConvertString(mesi[i]).c_str(), "ep");
		Corr[i]->Draw("same");
	}
	leg->Draw("same");

}


void DrawCorrection(TVirtualPad * c, TGraphErrors * Corr[],const string &name,const string &xaxisname,const string &yaxisname){

	gPad->SetGridx();
	gPad->SetGridy();
	gPad->SetLogx();
	gStyle->SetPalette(1);
	TLegend* leg =new TLegend(0.91,0.1,1.0,0.9);
	Corr[0]->SetLineWidth(2);
	Corr[0]->SetLineColor(1);
	Corr[0]->SetMarkerStyle(8);
	Corr[0]->SetMarkerColor(1);
	leg->AddEntry(Corr[0],ConvertString(mesi[0]).c_str(), "ep");
	Corr[0]->Draw("AC");
	Corr[0]->GetXaxis()->SetTitle(xaxisname.c_str());
	Corr[0]->GetYaxis()->SetTitle(yaxisname.c_str());
	Corr[0]->GetYaxis()->SetRangeUser(0.7,1.3);
	for(int i=1;i<num_mesi;i++) {
		Corr[i]->SetLineWidth(2);
		Corr[i]->SetLineColor(colorbase+2*i);
		Corr[i]->SetMarkerStyle(8);
		Corr[i]->SetMarkerColor(colorbase+2*i);
		leg->AddEntry(Corr[i],ConvertString(mesi[i]).c_str(), "ep");
		Corr[i]->Draw("Csame");
	}
	leg->Draw("same");
}

TGraphErrors * FluxesMean(TGraphErrors * Fluxes[],int bins){
	TGraphErrors * FluxMean = new TGraphErrors();	
	double x,y,ey=0;
        double x0,y0=0;
        int u,v=0;

	double flux[bins];
	double err[bins];
	
	for(int p=1;p<bins;p++){
		flux[p]=0;
		err[p]=0;
		for(int i=0;i<num_mesi;i++) {
			u=Fluxes[i]->GetPoint(p,x,y);
			ey=Fluxes[i]->GetErrorY(p);
			flux[p]+=y/pow(ey,2);
			err[p]+=1/pow(ey,2);	
		}
	flux[p] = flux[p]/err[p];
	err[p] = pow(1/err[p],0.5); 
	}
	
	for(int p=1;p<bins;p++){
		v=Fluxes[0]->GetPoint(p,x,y);
		FluxMean->SetPoint(p,x,flux[p]);
		FluxMean->SetPointError(p,0,err[p]);
	}
	return FluxMean;
}

void DrawFluxRatio(TVirtualPad * c, TGraphErrors * Fluxes[], TGraphErrors * FluxMean, const string &name,const string &xaxisname,const string &yaxisname, bool same = false){

        double x,y,ey=0;
        double x0,y0=0;
	int u,v=0;
	TGraphErrors *Fluxesratio[num_mesi];
        {       gPad->SetGridx();
                gPad->SetGridy();
                gPad->SetLogx();
                TLegend* leg =new TLegend(0.91,0.1,1.0,0.9);
                for(int i=0;i<num_mesi;i++) {
                Fluxesratio[i]=new TGraphErrors();
                for(int p=1;p<43;p++){
			u=Fluxes[i]->GetPoint(p,x,y);
                        v=FluxMean->GetPoint(p,x0,y0);
                        ey=Fluxes[i]->GetErrorY(p);
                        Fluxesratio[i]->SetPoint(p-1,x,y/y0);
                        Fluxesratio[i]->SetPointError(p-1,0,ey/y0);
                        }
                }
                Fluxesratio[0]->SetLineWidth(2);
                Fluxesratio[0]->SetLineColor(1);
                Fluxesratio[0]->SetMarkerStyle(8);
		Fluxesratio[0]->SetMarkerSize(1.5);
                Fluxesratio[0]->SetMarkerColor(1);
                Fluxesratio[0]->SetFillStyle(3002);
                leg->AddEntry(Fluxesratio[0],ConvertString(mesi[0]).c_str(), "ep");
                Fluxesratio[0]->GetYaxis()->SetRangeUser(0.2,1.3);
                if(!same) Fluxesratio[0]->Draw("AP");
		else Fluxesratio[0]->Draw("Psame");
                Fluxesratio[0]->GetXaxis()->SetTitle(xaxisname.c_str());
                Fluxesratio[0]->GetYaxis()->SetTitle(yaxisname.c_str());
                for(int i=0;i<num_mesi;i++) {
                        Fluxesratio[i]->SetLineWidth(2);
                        Fluxesratio[i]->SetLineColor(colorbase+2*i);
                        Fluxesratio[i]->SetMarkerStyle(8);
                        Fluxesratio[i]->SetMarkerSize(1.5);
			Fluxesratio[i]->SetMarkerColor(colorbase+2*i);
                        Fluxesratio[i]->SetFillStyle(3002);
                        leg->AddEntry(Fluxesratio[i],ConvertString(mesi[i]).c_str(), "ep");
                        Fluxesratio[i]->Draw("Psame");
                }
                leg->Draw("same");
        }

}

int main()
{
	int colorbase = 55;
	gStyle->SetPalette(1);
	cout<<"Analyzed Months: "<<endl;
	for(int i=0;i<num_mesi;i++){
		cout<<mesi[i].c_str()<<endl;
	}
	TFile * calib[40];
	TFile * result[40];
	string nomefile;
	cout<<" *************** LETTURA FILES *************** "<<endl;
	for(int i=0;i<num_mesi;i++){
		nomefile="./CALIBRAZIONI/"+mesi[i]+".root";
		calib[i]=TFile::Open(nomefile.c_str());
		nomefile="./Final_plots/"+mesi[i]+".root";
		result[i]=TFile::Open(nomefile.c_str());
		cout<<"calib: "<<calib[i]<<" results: "<<result[i]<<endl;
	}
	cout<<endl;

	string filename2="../database_P.root";
        TFile * file2 = TFile::Open(filename2.c_str(),"READ");
	
	string filename3="../database_D.root";
        TFile * file3 = TFile::Open(filename3.c_str(),"READ");

	cout<<"*** Protons ***"<<endl;
	TList *Experiments = file2->GetListOfKeys();
	TIter next(Experiments);
	TKey * key;
	TObject * obj;	

	std::vector<TGraphAsymmErrors *> P_Graphs;	
	
	while((key = (TKey*)next())){
		obj = file2->Get(key->GetName());
		if(obj->InheritsFrom("TGraphAsymmErrors")) P_Graphs.push_back((TGraphAsymmErrors *)obj); 
	}
	

	cout<<"*** Deutons ***"<<endl;

	std::vector<TGraphAsymmErrors *> D_Graphs;
	
	TList *ExperimentsD = file3->GetListOfKeys();
        TIter nextD(ExperimentsD);
        TKey * keyD;
	
        while((keyD = (TKey*)nextD())){
                obj = file3->Get(keyD->GetName());
		if(obj->InheritsFrom("TGraphAsymmErrors")) D_Graphs.push_back((TGraphAsymmErrors *)obj);
        }



	cout<<" ************** RISULTATI  ******************** "<<endl;

	for(int i=0;i<num_mesi;i++){
		Corr_L1[i] 		=  (TSpline3 *) 	calib[i]->Get("Fit Results/Splines/Corr_L1");
		Corr_TOFU[i] 		=  (TSpline3 *) 	calib[i]->Get("Fit Results/Splines/Corr_TOFU");
		Corr_Track[i] 		=  (TSpline3 *) 	calib[i]->Get("Fit Results/Splines/Corr_Track");
		Corr_TOFD[i] 		=  (TSpline3 *) 	calib[i]->Get("Fit Results/Splines/Corr_TOFD");	
		
		CorrLATpre_Spl[0][i] 	=  (TGraphErrors *) 	result[i]->Get("Export/Matching TOF");
		CorrLATpre_Spl[1][i] 	=  (TGraphErrors *) 	result[i]->Get("Export/Chi^2 R");
		CorrLATpre_Spl[2][i] 	=  (TGraphErrors *) 	result[i]->Get("Export/1 Tr. Track");
		CorrLAT_LikTOF_Spl[i] 	=    (TGraphErrors *) 	result[i]->Get("Export/CorrLAT_Lik_Spl");
		CorrLAT_LikNaF_Spl[i] 	=    (TGraphErrors *) 	result[i]->Get("Export/CorrLAT_LikNaF_Spl");
		CorrLAT_LikAgl_Spl[i] 	=    (TGraphErrors *) 	result[i]->Get("Export/CorrLAT_LikAgl_Spl");
		CorrLAT_DistTOF_Spl[i] 	=   (TGraphErrors *) 	result[i]->Get("Export/CorrLAT_Dist_Spl");
		CorrLAT_DistNaF_Spl[i] 	=   (TGraphErrors *) 	result[i]->Get("Export/CorrLAT_DistNaF_Spl");
		CorrLAT_DistAgl_Spl[i] 	=   (TGraphErrors *) 	result[i]->Get("Export/CorrLAT_DistAgl_Spl");

		preDVSMC_P[0][i] 	=  (TGraphErrors *) 	result[i]->Get("Export/DvsMC/DvsMC: Matching TOF_R");
                preDVSMC_P[1][i] 	=  (TGraphErrors *) 	result[i]->Get("Export/DvsMC/DvsMC: Chi^2 R_R");
                preDVSMC_P[2][i] 	=  (TGraphErrors *) 	result[i]->Get("Export/DvsMC/DvsMC: 1 Tr. Track_R");
		LikDVSMC_P[i]    	=  (TGraphErrors *) 	result[i]->Get("Export/DvsMC/LikDVSMCFit_P_Graph");
                DistDVSMC_P[i]   	=  (TGraphErrors *) 	result[i]->Get("Export/DvsMC/DistDVSMCFit_P_Graph");
		
		preDVSMC_PTOF[0][i] 	=  (TGraphErrors *) 	result[i]->Get("Export/DvsMC/DvsMC: Matching TOF_TOF");
                preDVSMC_PTOF[1][i] 	=  (TGraphErrors *) 	result[i]->Get("Export/DvsMC/DvsMC: Chi^2 R_TOF");
                preDVSMC_PTOF[2][i] 	=  (TGraphErrors *) 	result[i]->Get("Export/DvsMC/DvsMC: 1 Tr. Track_TOF");
		LikDVSMC_PTOF[i]    	=  (TGraphErrors *) 	result[i]->Get("Export/DvsMC/LikDVSMCFit_P_GraphTOF");
                DistDVSMC_PTOF[i]   	=  (TGraphErrors *) 	result[i]->Get("Export/DvsMC/DistDVSMCFit_P_GraphTOF");
		
		preDVSMC_PNaF[0][i] 	=  (TGraphErrors *) 	result[i]->Get("Export/DvsMC/DvsMC: Matching TOF_NaF");
                preDVSMC_PNaF[1][i] 	=  (TGraphErrors *) 	result[i]->Get("Export/DvsMC/DvsMC: Chi^2 R_NaF");
                preDVSMC_PNaF[2][i] 	=  (TGraphErrors *) 	result[i]->Get("Export/DvsMC/DvsMC: 1 Tr. Track_NaF");
		LikDVSMC_PNaF[i]    	=  (TGraphErrors *) 	result[i]->Get("Export/DvsMC/LikDVSMCFit_P_GraphNaF");
                DistDVSMC_PNaF[i]   	=  (TGraphErrors *) 	result[i]->Get("Export/DvsMC/DistDVSMCFit_P_GraphNaF");
		
		preDVSMC_PAgl[0][i] 	=  (TGraphErrors *) 	result[i]->Get("Export/DvsMC/DvsMC: Matching TOF_Agl");
                preDVSMC_PAgl[1][i] 	=  (TGraphErrors *) 	result[i]->Get("Export/DvsMC/DvsMC: Chi^2 R_Agl");
                preDVSMC_PAgl[2][i] 	=  (TGraphErrors *) 	result[i]->Get("Export/DvsMC/DvsMC: 1 Tr. Track_Agl");
		LikDVSMC_PAgl[i]    	=  (TGraphErrors *) 	result[i]->Get("Export/DvsMC/LikDVSMCFit_P_GraphAgl");
                DistDVSMC_PAgl[i]   	=  (TGraphErrors *) 	result[i]->Get("Export/DvsMC/DistDVSMCFit_P_GraphAgl");
		

		TrackerEff[i] 		 =  (TH1F *) 	   	result[i]->Get("Export/TrackerEfficiencyData");	
		TriggerEff[i]  		=  (TH1F *) 	   	result[i]->Get("Export/TriggerGlobalFactor");
		P_Fluxes[i]    		=  (TGraphErrors *) 	result[i]->Get("Export/PFluxes/Protons Primary Flux");
		D_FluxesTOF[i] 		=  (TGraphErrors *) 	result[i]->Get("Export/DFluxes/Deutons Flux: Primaries TOF");
		D_FluxesNaF[i] 		=  (TGraphErrors *) 	result[i]->Get("Export/DFluxes/Deutons Flux: Primaries NaF");
		D_FluxesAgl[i] 		=  (TGraphErrors *) 	result[i]->Get("Export/DFluxes/Deutons Flux: Primaries Agl");


	}
	cout<<" *************** FILES CHECK  ************************ "<<endl;	
	for(int i=0;i<num_mesi;i++) CheckFileIntegrity(i,mesi);
	cout<<endl;
        cout<<" *************** DRAW  ************************ "<<endl;
	TCanvas *c1=new TCanvas("L1 Energy Deposition");
	TCanvas *c2=new TCanvas("U. TOF Energy Deposition");
	TCanvas *c3=new TCanvas("Tracker Energy Deposition");
	TCanvas *c4=new TCanvas("L. TOF Energy Deposition");
	TCanvas *c5=new TCanvas("Matching TOF");
	TCanvas *c6=new TCanvas("Chi^2");
	TCanvas *c7=new TCanvas("1 Tr. Track");
	TCanvas *c8=new TCanvas("Likelihood");
	TCanvas *c9=new TCanvas("Distance");
	TCanvas *c10=new TCanvas("DvsMC: Matching TOF");
        TCanvas *c11=new TCanvas("DvsMC: Chi^2");
        TCanvas *c12=new TCanvas("DvsMC: 1 Tr. Track");
	TCanvas *c13=new TCanvas("DvsMC: Likelihood");
        TCanvas *c14=new TCanvas("DvsMC: Distance");
	TCanvas *c  =new TCanvas("Tracker Efficiency");
	TCanvas *c15=new TCanvas("Protons Fluxes Ratio");
	TCanvas *c16=new TCanvas("Deutons Fluxes Ratio");
	TCanvas *c15_bis=new TCanvas("Protons Fluxes");
        TCanvas *c16_bis=new TCanvas("Deutons Fluxes");
	TCanvas *c17=new TCanvas("Time dependence");

	c1->cd();
	DrawCalibration(c1,Corr_L1,"L1 E.dep. Corr. factors","Beta","Corr. Factor");
	
	c2->cd();
	DrawCalibration(c2,Corr_TOFU,"Upper TOF E.dep. Corr. factors","Beta","Corr. Factor");
	
	c3->cd();
	DrawCalibration(c3,Corr_Track,"Tracker E.dep. Corr. factors","Beta","Corr. Factor");

	c4->cd();
	DrawCalibration(c4,Corr_TOFD,"Lower E.dep. Corr. factors","Beta","Corr. Factor");
	c5->cd();
	DrawCorrection(c5,CorrLATpre_Spl[0],"","Latitude","Corr. Factor");

	c6->cd();
	DrawCorrection(c6,CorrLATpre_Spl[1],"","Latitude","Corr. Factor");

	c7->cd();
	DrawCorrection(c6,CorrLATpre_Spl[2],"","Latitude","Corr. Factor");
	c8->Divide(1,3);
	c8->cd(1);
	DrawCorrection(c8,CorrLAT_LikTOF_Spl,"","latitude","Corr. Factor");

	c8->cd(2);
	DrawCorrection(c8,CorrLAT_LikNaF_Spl,"","Latitude","Corr- Factor");

	c8->cd(3);
	DrawCorrection(c8,CorrLAT_LikAgl_Spl,"nome","Latitude","Corr. Factor");
	c9->Divide(1,3);
	c9->cd(1);
	DrawCorrection(c9,CorrLAT_DistTOF_Spl,"","Latitude","Corr. Factor");

	c9->cd(2);
	DrawCorrection(c9,CorrLAT_DistNaF_Spl,"","Latitude","Corr. Factor");

	c9->cd(3);
	DrawCorrection(c9,CorrLAT_DistNaF_Spl,"","Latitude","Corr. Factor");

	c10->Divide(2,2);
	c11->Divide(2,2);
	c12->Divide(2,2);
	c13->Divide(2,2);
	c14->Divide(2,2);

	c10->cd(1);
	DrawCorrection(c10,preDVSMC_P[0],"","R [GV]","Corr. Factor");
	c10->cd(2);
	DrawCorrection(c10,preDVSMC_PTOF[0],"","Kin. En. / nucl. [GeV/nucl.]","Corr. Factor");
	c10->cd(3);
	DrawCorrection(c10,preDVSMC_PNaF[0],"","Kin. En. / nucl. [GeV/nucl.]","Corr. Factor");
	c10->cd(4);
	DrawCorrection(c10,preDVSMC_PAgl[0],"","Kin. En. / nucl. [GeV/nucl.]","Corr. Factor");
			
	c11->cd(1);
	DrawCorrection(c11,preDVSMC_P[1],"","R [GV]","Corr. Factor");
	c11->cd(2);
	DrawCorrection(c10,preDVSMC_PTOF[1],"","Kin. En. / nucl. [GeV/nucl.]","Corr. Factor");
	c11->cd(3);
	DrawCorrection(c10,preDVSMC_PNaF[1],"","Kin. En. / nucl. [GeV/nucl.]","Corr. Factor");
	c11->cd(4);
	DrawCorrection(c10,preDVSMC_PAgl[1],"","Kin. En. / nucl. [GeV/nucl.]","Corr. Factor");
		
	c12->cd(1);
	DrawCorrection(c12,preDVSMC_P[2],"","R [GV]","Corr. Factor");
	c12->cd(2);
	DrawCorrection(c10,preDVSMC_PTOF[2],"","Kin. En. / nucl. [GeV/nucl.]","Corr. Factor");
	c12->cd(3);
	DrawCorrection(c10,preDVSMC_PNaF[2],"","Kin. En. / nucl. [GeV/nucl.]","Corr. Factor");
	c12->cd(4);
	DrawCorrection(c10,preDVSMC_PAgl[2],"","Kin. En. / nucl. [GeV/nucl.]","Corr. Factor");
	
	c13->cd(1);
	DrawCorrection(c13,LikDVSMC_P,"","R [GV]","Corr. Factor");
	c13->cd(2);
	DrawCorrection(c10,LikDVSMC_PTOF,"","Kin. En. / nucl. [GeV/nucl.]","Corr. Factor");
	c13->cd(3);
	DrawCorrection(c10,LikDVSMC_PNaF,"","Kin. En. / nucl. [GeV/nucl.]","Corr. Factor");
	c13->cd(4);
	DrawCorrection(c10,LikDVSMC_PAgl,"","Kin. En. / nucl. [GeV/nucl.]","Corr. Factor");
	
	c14->cd(1);
	DrawCorrection(c14,DistDVSMC_P,"","R [GV]","Corr. Factor");
	c14->cd(2);
	DrawCorrection(c10,DistDVSMC_PTOF,"","Kin. En. / nucl. [GeV/nucl.]","Corr. Factor");
	c14->cd(3);
	DrawCorrection(c10,DistDVSMC_PNaF,"","Kin. En. / nucl. [GeV/nucl.]","Corr. Factor");
	c14->cd(4);
	DrawCorrection(c10,DistDVSMC_PAgl,"","Kin. En. / nucl. [GeV/nucl.]","Corr. Factor");
	
	c->cd();
	gPad->SetGridx();
	gPad->SetGridy();
	double effy,effey=0;
	TH1F * Trackeff_time = new TH1F("","",num_mesi,0,num_mesi);
        TH1F * Triggeff_time = new TH1F("","",num_mesi,0,num_mesi);
	for(int i=0;i<num_mesi;i++) {
		effy = TrackerEff[i]-> GetBinContent(33);
		Trackeff_time -> SetBinContent(i+1,effy);
                effey = TrackerEff[i]-> GetBinError(33);
                Trackeff_time -> SetBinError(i+1,effey);
                Trackeff_time -> GetXaxis() -> SetBinLabel(i+1,ConvertString(mesi[i]).c_str());
        	effy = TriggerEff[i]-> GetBinContent(1);                
                Triggeff_time -> SetBinContent(i+1,effy);
                effey = TriggerEff[i]-> GetBinError(1);
                Triggeff_time -> SetBinError(i+1,effey);
		Trackeff_time -> GetXaxis() -> SetBinLabel(i+1,ConvertString(mesi[i]).c_str());
	}
	Trackeff_time -> SetMarkerStyle(8);
	Trackeff_time -> LabelsOption("v");
	Triggeff_time -> SetMarkerStyle(8);
	Trackeff_time -> LabelsOption();
	Trackeff_time -> GetXaxis() -> SetLabelSize(0.085);
	Trackeff_time -> SetMarkerColor(1);
	Trackeff_time -> SetMarkerColor(4);
	Triggeff_time -> SetMarkerColor(2);
	Trackeff_time ->GetYaxis()->SetRangeUser(0,1); 
	Trackeff_time -> Draw();
	Triggeff_time -> Draw("same");
	Trackeff_time -> Draw();


	c15_bis-> Divide(2,1);
	c15_bis->cd(1);
	gPad->SetGridx();
	gPad->SetGridy();
	gPad->SetLogx();
	gPad->SetLogy();
	P_Fluxes[0]->Draw("AP");
	P_Fluxes[0]->SetMarkerColor(colorbase);
	P_Fluxes[0]->SetLineColor(colorbase);
	for(int i=1;i<num_mesi;i++){
		P_Fluxes[i]->SetMarkerColor(colorbase+2*i);
		P_Fluxes[i]->SetLineColor(colorbase+2*i);
		P_Fluxes[i]->SetMarkerStyle(8);
		P_Fluxes[i]->SetMarkerSize(1);
		P_Fluxes[i]->Draw("Psame");
	}
	{
		TLegend* leg =new TLegend(0.91,0.1,1.0,0.9);
		for(int i=0;i<num_mesi;i++) leg->AddEntry(P_Fluxes[i],ConvertString(mesi[i]).c_str(),"ep");
		leg->Draw("same");
	}

	c15_bis->cd(2);
        gPad->SetGridx();
        gPad->SetGridy();
        gPad->SetLogx();
        gPad->SetLogy();
	TGraphErrors * P_Mean = FluxesMean(P_Fluxes,43);
        P_Mean->SetMarkerColor(colorbase);
        P_Mean->SetLineColor(colorbase);
	P_Mean->SetTitle("Protons Primary Flux");
        P_Mean->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
        P_Mean->GetYaxis()->SetTitle("Flux [(m^2 sec sr GeV/nucl.)^-1]");
	P_Mean->GetXaxis()->SetTitleSize(0.045);
        P_Mean->GetYaxis()->SetTitleSize(0.045);
	P_Mean->SetLineWidth(3);
	P_Mean->SetMarkerStyle(8);
        P_Mean->SetMarkerSize(2);
	P_Mean->Draw("AP");
	{
	TLegend* leg =new TLegend(0.91,0.1,1.0,0.9);
                leg->AddEntry(P_Mean,"This Work", "ep");

        for(uint n=0;n<P_Graphs.size();n++){
                P_Graphs[n] -> SetMarkerSize(1.5);
		P_Graphs[n] -> SetLineWidth(1.5);
		P_Graphs[n] ->Draw("Psame");
                 leg->AddEntry(P_Graphs[n],P_Graphs[n]->GetTitle(),"ep");
        }
	leg->Draw("same");
	}

	c15-> Divide(1,2);
        c15->cd(1);
        DrawFluxRatio(c15,P_Fluxes,P_Mean,"R [GV]","Proton Flux","Proton Flux (norm. to mean flux)");

	c16_bis-> Divide(1,2);
	c16_bis->cd(1);
        gPad->SetGridx();
        gPad->SetGridy();
        gPad->SetLogx();
        gPad->SetLogy();

	TH2F * Frame = new TH2F("Frame","Frame",1000,0,30,1e4,0.01,1000);
	Frame->SetTitle("Deuterons Primary Flux");
	Frame->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
	Frame->GetYaxis()->SetTitle("Flux [(m^2 sec sr GeV/nucl.)^-1]");
	Frame->GetXaxis()->SetTitleSize(0.045);
	Frame->GetYaxis()->SetTitleSize(0.045);
	Frame->Draw();
	D_FluxesTOF[0]->Draw("Psame");
        D_FluxesTOF[0]->SetMarkerColor(colorbase);
        D_FluxesTOF[0]->SetLineColor(colorbase);
        for(int i=1;i<num_mesi;i++){
                D_FluxesTOF[i]->SetMarkerColor(colorbase+2*i);
                D_FluxesTOF[i]->SetLineColor(colorbase+2*i);
                D_FluxesTOF[i]->Draw("Psame");
        }
        {
                TLegend* leg =new TLegend(0.91,0.1,1.0,0.9);
                for(int i=0;i<num_mesi;i++) leg->AddEntry(D_FluxesTOF[i],ConvertString(mesi[i]).c_str(),"ep");
                leg->Draw("same");
        }

	D_FluxesNaF[0]->Draw("Psame");
        D_FluxesNaF[0]->SetMarkerColor(colorbase);
        D_FluxesNaF[0]->SetLineColor(colorbase);
        for(int i=1;i<num_mesi;i++){
                D_FluxesNaF[i]->SetMarkerColor(colorbase+2*i);
                D_FluxesNaF[i]->SetLineColor(colorbase+2*i);
                D_FluxesNaF[i]->Draw("Psame");
        }
        D_FluxesAgl[0]->Draw("Psame");
        D_FluxesAgl[0]->SetMarkerColor(colorbase);
        D_FluxesAgl[0]->SetLineColor(colorbase);
        for(int i=1;i<num_mesi;i++){
                D_FluxesAgl[i]->SetMarkerColor(colorbase+2*i);
                D_FluxesAgl[i]->SetLineColor(colorbase+2*i);
                D_FluxesAgl[i]->Draw("Psame");
        }	
		

	c16_bis->cd(2);
        gPad->SetGridx();
        gPad->SetGridy();
        gPad->SetLogx();
        gPad->SetLogy();

        TH2F * Frame2 = new TH2F("Frame2","Frame2",1000,0,30,1e4,0.01,1000);
        Frame2->SetTitle("Deuterons Primary Flux");
        Frame2->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
        Frame2->GetYaxis()->SetTitle("Flux [(m^2 sec sr GeV/nucl.)^-1]");
        Frame2->GetXaxis()->SetTitleSize(0.045);
        Frame2->GetYaxis()->SetTitleSize(0.045);
        Frame2->Draw();
	TGraphErrors * D_MeanTOF = FluxesMean(D_FluxesTOF,18);
        D_MeanTOF->SetMarkerColor(colorbase);
        D_MeanTOF->SetLineColor(colorbase);
        D_MeanTOF->SetMarkerStyle(8);
	D_MeanTOF->SetLineWidth(3);
        D_MeanTOF->SetMarkerSize(2);
        D_MeanTOF->Draw("Psame");
	TGraphErrors * D_MeanNaF = FluxesMean(D_FluxesNaF,18);
        D_MeanNaF->SetMarkerColor(colorbase);
        D_MeanNaF->SetLineColor(colorbase);
        D_MeanNaF->SetLineWidth(3);
	D_MeanNaF->SetMarkerStyle(8);
        D_MeanNaF->SetMarkerSize(2);
        D_MeanNaF->Draw("Psame");
	TGraphErrors * D_MeanAgl = FluxesMean(D_FluxesAgl,18);
        D_MeanAgl->SetMarkerColor(colorbase);
        D_MeanAgl->SetLineColor(colorbase);
        D_MeanAgl->SetLineWidth(3);
	D_MeanAgl->SetMarkerStyle(8);
        D_MeanAgl->SetMarkerSize(2);
        D_MeanAgl->Draw("Psame");
	{
        TLegend* leg =new TLegend(0.91,0.1,1.0,0.9);
                leg->AddEntry(D_MeanTOF,"This Work", "ep");

        for(uint n=0;n<D_Graphs.size();n++){
                D_Graphs[n] ->Draw("Psame");
                D_Graphs[n] -> SetMarkerSize(1.5);
                D_Graphs[n] -> SetLineWidth(1.5);
		leg->AddEntry(D_Graphs[n],D_Graphs[n]->GetTitle(),"ep");
        }
        leg->Draw("same");
        }

	c16->cd();
        bool same = true;
	TH2F * Frame3 = new TH2F("Frame3","Frame3",1000,0,30,1000,0.4,1.4);
	Frame3->Draw();
	DrawFluxRatio(c16,D_FluxesAgl,D_MeanAgl,"Kin. En./nucl.","Deuton Flux","Deutons Flux (norm. to mean flux)",same);
	DrawFluxRatio(c16,D_FluxesTOF,D_MeanTOF,"Kin. En./nucl.","Deuton Flux","Deutons Flux (norm. to mean flux)",same);	
	DrawFluxRatio(c16,D_FluxesNaF,D_MeanNaF,"Kin. En./nucl.","Deuton Flux","Deutons Flux (norm. to mean flux)",same);



	
	c17->Divide(1,3);

	TH1F * Time_depD1 = new TH1F("","",num_mesi,0,num_mesi);
	TH1F * Time_depP1 = new TH1F("","",num_mesi,0,num_mesi);
	TH1F * Time_depD2 = new TH1F("","",num_mesi,0,num_mesi);
        TH1F * Time_depP2 = new TH1F("","",num_mesi,0,num_mesi);
	TH1F * Time_depD3 = new TH1F("","",num_mesi,0,num_mesi);
        TH1F * Time_depP3 = new TH1F("","",num_mesi,0,num_mesi);
	
	TH1F * Time_depP4 = new TH1F("","",num_mesi,0,num_mesi);
	TH1F * errorP4    = new TH1F("","",num_mesi,0,num_mesi);
	

	double x01,x02,x03,x04;
	double y01,y02,y03,y04;
	double x1,x2,x3,x4;
	double y1,y2,y3,y4;
	double ey1,ey2,ey3,ey4;

	for(int i=0;i<num_mesi;i++) {
		D_FluxesTOF[i]->GetPoint(6,x1,y1);
		D_MeanTOF->GetPoint(6,x01,y01);
		
		D_FluxesTOF[i]->GetPoint(12,x2,y2);
                D_MeanTOF->GetPoint(12,x02,y02);
		
		D_FluxesAgl[i]->GetPoint(9,x3,y3);
                D_MeanAgl->GetPoint(9,x03,y03);

		Time_depD1 -> SetBinContent(i+1,y1/y01);
		Time_depD2 -> SetBinContent(i+1,y2/y02);		
		Time_depD3 -> SetBinContent(i+1,y3/y03);		

		ey1 = D_FluxesTOF[i]->GetErrorY(4);
		ey2 = D_FluxesTOF[i]->GetErrorY(12);
		ey3 = D_FluxesAgl[i]->GetErrorY(9);
		
		Time_depD1 -> SetBinError(i+1,ey1/y01);
		Time_depD2 -> SetBinError(i+1,ey2/y02);
		Time_depD3 -> SetBinError(i+1,ey3/y03);

		Time_depD1 -> GetXaxis() -> SetBinLabel(i+1,ConvertString(mesi[i]).c_str());
		Time_depD2 -> GetXaxis() -> SetBinLabel(i+1,ConvertString(mesi[i]).c_str());
		Time_depD3 -> GetXaxis() -> SetBinLabel(i+1,ConvertString(mesi[i]).c_str());
	}
	for(int i=0;i<num_mesi;i++) {
                P_Fluxes[i]->GetPoint(4,x1,y1);
                P_Mean->GetPoint(4,x01,y01);
                
		P_Fluxes[i]->GetPoint(11,x2,y2);
                P_Mean->GetPoint(11,x02,y02);
       
		P_Fluxes[i]->GetPoint(25,x3,y3);
                P_Mean->GetPoint(25,x03,y03);
                  
		P_Fluxes[i]->GetPoint(40,x4,y4);
                P_Mean->GetPoint(40,x04,y04);
                
		Time_depP1 -> SetBinContent(i+1,y1/y01);
		Time_depP2 -> SetBinContent(i+1,y2/y02);         
		Time_depP3 -> SetBinContent(i+1,y3/y03);         
		Time_depP4 -> SetBinContent(i+1,y4/y04);
		errorP4 -> SetBinContent(i+1,y04/y04);


	        ey1 = P_Fluxes[i]->GetErrorY(6);
		ey2 = P_Fluxes[i]->GetErrorY(11);	
		ey3 = P_Fluxes[i]->GetErrorY(25);	
		ey4 = P_Fluxes[i]->GetErrorY(40);		

		Time_depP1 -> SetBinError(i+1,ey1/y01);
        	Time_depP2 -> SetBinError(i+1,ey2/y02);
		Time_depP3 -> SetBinError(i+1,ey3/y03);
		errorP4    -> SetBinError(i+1,ey4/y04);
		Time_depP4 -> GetXaxis() -> SetBinLabel(i+1,ConvertString(mesi[i]).c_str());

		
	}
	Time_depD1 -> SetMarkerStyle(8);
	Time_depP1 -> SetMarkerStyle(8);
	Time_depD1 -> SetMarkerColor(4);
        Time_depP1 -> SetMarkerColor(2);
	Time_depD1 -> SetLineColor(4);
        Time_depP1 -> SetLineColor(2);

	Time_depD2 -> SetMarkerStyle(8);
        Time_depP2 -> SetMarkerStyle(8);
        Time_depD2 -> SetMarkerColor(4);
        Time_depP2 -> SetMarkerColor(2);
	Time_depD2 -> SetLineColor(4);
        Time_depP2 -> SetLineColor(2);

	Time_depD3 -> SetMarkerStyle(8);
        Time_depP3 -> SetMarkerStyle(8);
        Time_depD3 -> SetMarkerColor(4);
        Time_depP3 -> SetMarkerColor(2);
	Time_depD3 -> SetLineColor(4);
        Time_depP3 -> SetLineColor(2);

	Time_depP4 -> SetMarkerStyle(8);
        Time_depP4 -> SetMarkerColor(2);
        errorP4	   -> SetFillColor(1);
	

	gStyle -> SetTitleFontSize(0.3);
	gStyle -> SetLegendFillColor(0);	
	c17->cd(1);
        gPad->SetGridx();
        gPad->SetGridy();
	Time_depD1 -> GetYaxis() -> SetRangeUser(0.6,1.3);
	Time_depD1 -> LabelsOption("h");
	Time_depD1 -> SetTitle("R = 1 GV");
	Time_depD1 -> GetYaxis() -> SetLabelSize(0.085);
	Time_depD1 -> GetXaxis() -> SetLabelSize(0.085);
	Time_depD1 -> GetYaxis() -> SetTitle("Flux / Mean");
	Time_depD1 -> GetYaxis() -> SetTitleSize(0.085);
	Time_depD1 -> SetMarkerSize(2);
	Time_depP1 -> SetMarkerSize(2);
	Time_depD1 -> Draw();
	Time_depP1 -> Draw("same");	
	 TLegend* leg =new TLegend(0.7,0.1,0.9,0.3);
         leg->AddEntry(Time_depD1,"Deuterons", "ep");
	leg->AddEntry(Time_depP1,"Protons", "ep");
	leg ->Draw("same");

	c17->cd(2);
        gPad->SetGridx();
        gPad->SetGridy();
        Time_depD2 -> GetYaxis() -> SetRangeUser(0.6,1.3);
	Time_depD2 -> SetTitle("R = 2 GV");
	Time_depD2 -> LabelsOption("h");
	Time_depD2 -> GetYaxis() -> SetLabelSize(0.085);
	Time_depD2 -> GetXaxis() -> SetLabelSize(0.085);
	Time_depD2 -> GetYaxis() -> SetTitle("Flux / Mean");
	Time_depD2 -> GetYaxis() -> SetTitleSize(0.085);
	Time_depD2 -> SetMarkerSize(2);
	Time_depP2 -> SetMarkerSize(2); 
	Time_depD2-> Draw();
        Time_depP2 -> Draw("same");

	c17->cd(3);
        gPad->SetGridx();
        gPad->SetGridy();
        Time_depD3 -> GetYaxis() -> SetRangeUser(0.6,1.3);
	Time_depD3 -> SetTitle("R = 10 GV");
	Time_depD3 -> LabelsOption("h");
	Time_depD3 -> GetYaxis() -> SetLabelSize(0.085);
	Time_depD3 -> GetXaxis() -> SetLabelSize(0.085);
	Time_depD3 -> GetYaxis() -> SetTitle("Flux / Mean");
	Time_depD3 -> GetYaxis() -> SetTitleSize(0.085);
	Time_depD3 -> SetMarkerSize(2);
	Time_depP3 -> SetMarkerSize(2);
	Time_depD3 -> Draw();
        Time_depP3 -> Draw("same");

	c15->cd(2);
        gPad->SetGridx();
        gPad->SetGridy();
        Time_depP4 -> GetYaxis() -> SetRangeUser(0.8,1.2);
	Time_depP4 -> SetTitle("High Energy Flux stability (R=80 GV)");
	Time_depP4 -> GetYaxis() -> SetLabelSize(0.085);
	Time_depP4 -> GetXaxis() -> SetLabelSize(0.085);
	Time_depP4 -> GetYaxis() -> SetTitle("Flux / Mean");
	Time_depP4 -> GetYaxis() -> SetTitleSize(0.085);
	Time_depP4 -> SetMarkerSize(2);
	Time_depP4 -> Draw("P");
	errorP4 -> SetFillStyle(3001);
	errorP4 -> Draw("E4same");








	cout<<"*************** OUTPUT **************************"<<endl;
	nomefile="./Final_plots/TimeAnalysis.root";
	TFile *f_out=new TFile(nomefile.c_str(), "RECREATE");
	f_out->mkdir("Calibrations");
	f_out->mkdir("Lat. Dependence");
	f_out->mkdir("Tracker Eff.");
	f_out->mkdir("Data vs MC");
	f_out->mkdir("Fluxes");

	f_out->cd("Calibrations");
	c1->Write();
	c2->Write();
	c3->Write();
        c4->Write();
	f_out->cd("Lat. Dependence");
	c5->Write();
	c6->Write();
	c7->Write();
	c8->Write();
	c9->Write();
	f_out->cd("Data vs MC");
	c10->Write();
	c11->Write();
        c12->Write();
	c13->Write();
	c14->Write();
	f_out->cd("Tracker Eff.");
	c -> Write();
	f_out->cd("Fluxes");
	c15->Write();
	c15_bis->Write();
	c16->Write();	
	c16_bis->Write();
	c17->Write();
	f_out->Write();
	f_out->Close();
	
	return 0;
}


void CheckFileIntegrity(int i,std::string mesi[]){
	if(Corr_L1[i] 			&&	
           Corr_TOFU[i] 		&&
           Corr_Track[i] 		&&
           Corr_TOFD[i] 		&&
          
           CorrLATpre_Spl[0][i] 	&& 
           CorrLATpre_Spl[1][i] 	&&
           CorrLATpre_Spl[2][i] 	&&
           CorrLAT_LikTOF_Spl[i] 	&&
           CorrLAT_LikNaF_Spl[i] 	&& 
           CorrLAT_LikAgl_Spl[i] 	&&
           CorrLAT_DistTOF_Spl[i] 	&&
           CorrLAT_DistNaF_Spl[i] 	&&
           CorrLAT_DistAgl_Spl[i] 	&&	
                                   
           preDVSMC_P[0][i] 	&& 
           preDVSMC_P[1][i] 	&&
           preDVSMC_P[2][i] 	&&
           LikDVSMC_P[i]    	&&
           DistDVSMC_P[i]   	&&	
           
           preDVSMC_PTOF[0][i] 	&& 
           preDVSMC_PTOF[1][i] 	&&
           preDVSMC_PTOF[2][i] 	&&
           LikDVSMC_PTOF[i]    	&&
           DistDVSMC_PTOF[i]   &&	
           
           preDVSMC_PNaF[0][i] 	&& 
           preDVSMC_PNaF[1][i] 	&&
           preDVSMC_PNaF[2][i] 	&&
           LikDVSMC_PNaF[i]    	&&
           DistDVSMC_PNaF[i]   &&	
           
           preDVSMC_PAgl[0][i] 	&& 
           preDVSMC_PAgl[1][i] 	&&
           preDVSMC_PAgl[2][i] 	&&
           LikDVSMC_PAgl[i]    	&&
           DistDVSMC_PAgl[i]    &&	
           
                                   
           TrackerEff[i] 	&& 	
           TriggerEff[i]  	&& 	
           P_Fluxes[i]    	&& 	
           D_FluxesTOF[i] 	&& 	
           D_FluxesNaF[i] 	&&	
           D_FluxesAgl[i] 	)  cout<<"file "<<ConvertString(mesi[i]).c_str()<<"( "<<i<<") is complete"<<endl;
		else cout<<"file "<<ConvertString(mesi[i]).c_str()<<"( "<<i<<") is broken"<<endl;	
}

std::string ConvertString(std::string mese){
	string Mese;
	if(mese=="2011_01") Mese="Jan '11";
	if(mese=="2011_02") Mese="Feb '11";
	if(mese=="2011_03") Mese="Mar '11";
	if(mese=="2011_04") Mese="Apr '11";
	if(mese=="2011_05") Mese="May '11";
	if(mese=="2011_06") Mese="Jun '11";
	if(mese=="2011_07") Mese="Jul '11";
	if(mese=="2011_08") Mese="Aug '11";
	if(mese=="2011_09") Mese="Sep '11";
	if(mese=="2011_10") Mese="Oct '11";
	if(mese=="2011_11") Mese="Nov '11";
	if(mese=="2011_12") Mese="Dec '11";
	if(mese=="2012_01") Mese="Jan '12";
	if(mese=="2012_02") Mese="Feb '12";
	if(mese=="2012_03") Mese="Mar '12";
	if(mese=="2012_04") Mese="Apr '12";
	if(mese=="2012_05") Mese="May '12";
	if(mese=="2012_06") Mese="Jun '12";
	if(mese=="2012_07") Mese="Jul '12";
	if(mese=="2012_08") Mese="Aug '12";
	if(mese=="2012_09") Mese="Sep '12";
	if(mese=="2012_10") Mese="Oct '12";
	if(mese=="2012_11") Mese="Nov '12";
	if(mese=="2012_12") Mese="Dec '12";
	if(mese=="2013_01") Mese="Jan '13";
	if(mese=="2013_02") Mese="Feb '13";
	if(mese=="2013_03") Mese="Mar '13";
	if(mese=="2013_04") Mese="Apr '13";
	if(mese=="2013_05") Mese="May '13";
	if(mese=="2013_06") Mese="Jun '13";
	if(mese=="2013_07") Mese="Jul '13";
	if(mese=="2013_08") Mese="Aug '13";
	if(mese=="2013_09") Mese="Sep '13";
	if(mese=="2013_10") Mese="Oct '13";
	if(mese=="2013_11") Mese="Nov '13";
	if(mese=="2013_12") Mese="Dec '13";
	return Mese;
}












