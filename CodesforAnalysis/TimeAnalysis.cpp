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
#include "TObject.h"
#include "TGraphAsymmErrors.h"
#include <TGraph.h>
#include "TGraphErrors.h"
#include <cstring>
#include "Mesi.h"

using namespace std;

int main()
{
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
	cout<<" ************** RISULTATI  ******************** "<<endl;

	TSpline3 *Corr_L1[40];
	TSpline3 *Corr_TOFU[40];
	TSpline3 *Corr_Track[40];
	TSpline3 *Corr_TOFD[40];
	TGraphErrors *CorrLATpre_Spl[3][40];
	TGraphErrors *CorrLAT_DistTOF_Spl[40];
	TGraphErrors *CorrLAT_LikTOF_Spl[40];	
	TGraphErrors *CorrLAT_DistNaF_Spl[40];
        TGraphErrors *CorrLAT_LikNaF_Spl[40];
	TGraphErrors *CorrLAT_DistAgl_Spl[40];
        TGraphErrors *CorrLAT_LikAgl_Spl[40];
	TGraphErrors *preDVSMC_P[3][40];
	TGraphErrors *LikDVSMC_P[40];
        TGraphErrors *DistDVSMC_P[40];
	TGraphErrors *P_Fluxes[40];
	TGraphErrors *P_Fluxesratio[40];
	TGraphErrors *D_FluxesTOF[40];
        TGraphErrors *D_FluxesratioTOF[40];
	TGraphErrors *D_FluxesNaF[40];
        TGraphErrors *D_FluxesratioNaF[40];
	TGraphErrors *D_FluxesAgl[40];
        TGraphErrors *D_FluxesratioAgl[40];	
	for(int i=0;i<num_mesi;i++){
		Corr_L1[i] =  (TSpline3 *) calib[i]->Get("Fit Results/Splines/Corr_L1");
		Corr_TOFU[i] =  (TSpline3 *) calib[i]->Get("Fit Results/Splines/Corr_TOFU");
		Corr_Track[i] =  (TSpline3 *) calib[i]->Get("Fit Results/Splines/Corr_Track");
		Corr_TOFD[i] =  (TSpline3 *) calib[i]->Get("Fit Results/Splines/Corr_TOFD");	
		
		CorrLATpre_Spl[0][i] =  (TGraphErrors *) result[i]->Get("Export/Matching TOF");
		CorrLATpre_Spl[1][i] =  (TGraphErrors *) result[i]->Get("Export/Chi^2 R");
		CorrLATpre_Spl[2][i] =  (TGraphErrors *) result[i]->Get("Export/1 Tr. Track");
		CorrLAT_LikTOF_Spl[i] =    (TGraphErrors *) result[i]->Get("Export/CorrLAT_LikTOF_Spl");
		CorrLAT_LikNaF_Spl[i] =    (TGraphErrors *) result[i]->Get("Export/CorrLAT_LikNaF_Spl");
		CorrLAT_LikAgl_Spl[i] =    (TGraphErrors *) result[i]->Get("Export/CorrLAT_LikAgl_Spl");
		CorrLAT_DistTOF_Spl[i] =   (TGraphErrors *) result[i]->Get("Export/CorrLAT_LikTOF_Spl");
		CorrLAT_DistNaF_Spl[i] =   (TGraphErrors *) result[i]->Get("Export/CorrLAT_LikNaF_Spl");
		CorrLAT_DistAgl_Spl[i] =   (TGraphErrors *) result[i]->Get("Export/CorrLAT_LikAgl_Spl");
		/*preDVSMC_P[0][i] =  (TGraphErrors *) result[i]->Get("Export/DvsMC: Matching TOF_Graph");
                preDVSMC_P[1][i] =  (TGraphErrors *) result[i]->Get("Export/DvsMC: Chi^2 R_Graph");
                preDVSMC_P[2][i] =  (TGraphErrors *) result[i]->Get("Export/DvsMC: 1 Tr. Track_Graph");
		LikDVSMC_P[i] =  (TGraphErrors *) result[i]->Get("Export/LikDVSMC_P_Graph");
                DistDVSMC_P[i] =  (TGraphErrors *) result[i]->Get("Export/DistDVSMC_P_Graph");*/
		P_Fluxes[i]    =  (TGraphErrors *) result[i]->Get("Export/Protons Primary Flux");
		D_FluxesTOF[i] =  (TGraphErrors *) result[i]->Get("Export/Deutons Primary Flux: TOF");
		D_FluxesNaF[i] =  (TGraphErrors *) result[i]->Get("Export/Deutons Primary Flux: NaF");
		D_FluxesAgl[i] =  (TGraphErrors *) result[i]->Get("Export/Deutons Primary Flux: Agl");
	}	
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
	/*TCanvas *c10=new TCanvas("DvsMC: Matching TOF");
        TCanvas *c11=new TCanvas("DvsMC: Chi^2");
        TCanvas *c12=new TCanvas("DvsMC: 1 Tr. Track");
	TCanvas *c13=new TCanvas("DvsMC: Likelihood");
        TCanvas *c14=new TCanvas("DvsMC: Distance");*/
	TCanvas *c15=new TCanvas("Protons Fluxes");
	TCanvas *c16=new TCanvas("Deutons Fluxes");
	c1->cd();
	{gPad->SetGridx();
		gPad->SetGridy();
		gStyle->SetPalette(0);
		TLegend* leg =new TLegend(0.91,0.1,1.0,0.9);
		Corr_L1[0]->SetLineWidth(2); 
		Corr_L1[0]->SetLineColor(1);
		leg->AddEntry(Corr_L1[0],mesi[0].c_str(), "ep");
		Corr_L1[0]->Draw("");
		TH1F *frame=c1->DrawFrame(0.4,0.85,1,1.1,"L1 E.dep. Corr. factors");
		frame->GetXaxis()->SetTitle("Beta");
		frame->GetYaxis()->SetTitle("Corr. Factor");
		for(int i=1;i<num_mesi;i++) {
			Corr_L1[i]->SetLineWidth(2);    
			Corr_L1[i]->SetLineColor(i+1);
			leg->AddEntry(Corr_L1[i],mesi[i].c_str(), "ep");
			Corr_L1[i]->Draw("same");
		}
		leg->Draw("same");		
	}
	cout<<"1"<<endl;	
	c2->cd();
        {gPad->SetGridx();
                gPad->SetGridy();
                gStyle->SetPalette(0);
                TLegend* leg =new TLegend(0.91,0.1,1.0,0.9);
                Corr_TOFU[0]->SetLineWidth(2);
                Corr_TOFU[0]->SetLineColor(1);
                leg->AddEntry(Corr_TOFU[0],mesi[0].c_str(), "ep");
                Corr_TOFU[0]->Draw("");
                TH1F *frame=c2->DrawFrame(0.4,0.85,1,1.1,"Upper TOF E.dep. Corr. factors");
                frame->GetXaxis()->SetTitle("Beta");
                frame->GetYaxis()->SetTitle("Corr. Factor");
                for(int i=1;i<num_mesi;i++) {
                        Corr_TOFU[i]->SetLineWidth(2);
                        Corr_TOFU[i]->SetLineColor(i+1);
                        leg->AddEntry(Corr_L1[i],mesi[i].c_str(), "ep");
                        Corr_TOFU[i]->Draw("same");
                }
                leg->Draw("same");       
        }
	cout<<"2"<<endl;
	c3->cd();
        {gPad->SetGridx();
                gPad->SetGridy();
                gStyle->SetPalette(0);
                TLegend* leg =new TLegend(0.91,0.1,1.0,0.9);
                Corr_Track[0]->SetLineWidth(2);
                Corr_Track[0]->SetLineColor(1);
                leg->AddEntry(Corr_Track[0],mesi[0].c_str(), "ep");
                Corr_Track[0]->Draw("");
                TH1F *frame=c3->DrawFrame(0.4,0.85,1,1.1,"Upper Track E.dep. Corr. factors");
                frame->GetXaxis()->SetTitle("Beta");
                frame->GetYaxis()->SetTitle("Corr. Factor");
                for(int i=1;i<num_mesi;i++) {
                        Corr_Track[i]->SetLineWidth(2);
                        Corr_Track[i]->SetLineColor(i+1);
                        leg->AddEntry(Corr_Track[i],mesi[i].c_str(), "ep");
                        Corr_Track[i]->Draw("same");
                }
                leg->Draw("same");
        }
	cout<<"3"<<endl;
	c4->cd();
        {gPad->SetGridx();
                gPad->SetGridy();
                gStyle->SetPalette(0);
                TLegend* leg =new TLegend(0.91,0.1,1.0,0.9);
                Corr_TOFD[0]->SetLineWidth(2);
                Corr_TOFD[0]->SetLineColor(1);
		leg->AddEntry(Corr_TOFD[0],mesi[0].c_str(), "ep");
                Corr_TOFD[0]->Draw("AL");
		TH1F *frame=c4->DrawFrame(0.4,0.85,1,1.1,"Lower TOF E.dep. Corr. factors");
                frame->GetXaxis()->SetTitle("Beta");
                frame->GetYaxis()->SetTitle("Corr. Factor");
                for(int i=1;i<num_mesi;i++) {
                        Corr_TOFD[i]->SetLineWidth(2);
                        Corr_TOFD[i]->SetLineColor(i+1);
                        leg->AddEntry(Corr_TOFD[i],mesi[i].c_str(), "ep");
                        Corr_TOFD[i]->Draw("Lsame");
                }
                leg->Draw("same");
        }
	cout<<"4"<<endl;
	c5->cd();
        {gPad->SetGridx();
                gPad->SetGridy();
                gStyle->SetPalette(0);
                TLegend* leg =new TLegend(0.91,0.1,1.0,0.9);
               	cout<<CorrLATpre_Spl[0][0]<<endl;
		CorrLATpre_Spl[0][0]->SetLineWidth(2);
                CorrLATpre_Spl[0][0]->SetLineColor(1);
		CorrLATpre_Spl[0][0]->SetMarkerStyle(8);
		CorrLATpre_Spl[0][0]->SetMarkerColor(1);
		leg->AddEntry(CorrLATpre_Spl[0][0],mesi[0].c_str(), "ep");
                CorrLATpre_Spl[0][0]->Draw("AC");
                CorrLATpre_Spl[0][0]->GetXaxis()->SetTitle("Latitude");
                CorrLATpre_Spl[0][0]->GetYaxis()->SetTitle("Corr. Factor");
                for(int i=1;i<num_mesi;i++) {
                        CorrLATpre_Spl[0][i]->SetLineWidth(2);
                        CorrLATpre_Spl[0][i]->SetLineColor(i+1);
			CorrLATpre_Spl[0][i]->SetMarkerStyle(8);
	                CorrLATpre_Spl[0][i]->SetMarkerColor(i+1);
                        leg->AddEntry(CorrLATpre_Spl[0][i],mesi[i].c_str(), "ep");
                        CorrLATpre_Spl[0][i]->Draw("Csame");
                }
                leg->Draw("same");
        }
	cout<<"5"<<endl;
	c6->cd();
        {gPad->SetGridx();
                gPad->SetGridy();
                gStyle->SetPalette(0);
                TLegend* leg =new TLegend(0.91,0.1,1.0,0.9);
                CorrLATpre_Spl[1][0]->SetLineWidth(2);
                CorrLATpre_Spl[1][0]->SetLineColor(1);
                CorrLATpre_Spl[1][0]->SetMarkerStyle(8);
                CorrLATpre_Spl[1][0]->SetMarkerColor(1);
		leg->AddEntry(CorrLATpre_Spl[1][0],mesi[0].c_str(), "ep");
                CorrLATpre_Spl[1][0]->Draw("AC");
                CorrLATpre_Spl[1][0]->GetXaxis()->SetTitle("Latitude");
                CorrLATpre_Spl[1][0]->GetYaxis()->SetTitle("Corr. Factor");
                for(int i=1;i<num_mesi;i++) {
                        CorrLATpre_Spl[1][i]->SetLineWidth(2);
                        CorrLATpre_Spl[1][i]->SetLineColor(i+1);
			CorrLATpre_Spl[1][i]->SetMarkerStyle(8);
	                CorrLATpre_Spl[1][i]->SetMarkerColor(i+1);
                        leg->AddEntry(CorrLATpre_Spl[1][i],mesi[i].c_str(), "ep");
                        CorrLATpre_Spl[1][i]->Draw("Csame");
                }
                leg->Draw("same");
        }
	cout<<"6"<<endl;
	c7->cd();
        {gPad->SetGridx();
                gPad->SetGridy();
                gStyle->SetPalette(0);
                TLegend* leg =new TLegend(0.91,0.1,1.0,0.9);
                CorrLATpre_Spl[2][0]->SetLineWidth(2);
                CorrLATpre_Spl[2][0]->SetLineColor(1);
		CorrLATpre_Spl[2][0]->SetMarkerStyle(8);
                CorrLATpre_Spl[2][0]->SetMarkerColor(1);
                leg->AddEntry(CorrLATpre_Spl[2][0],mesi[0].c_str(), "ep");
                CorrLATpre_Spl[2][0]->Draw("ACP");
                CorrLATpre_Spl[2][0]->GetXaxis()->SetTitle("Latitude");
                CorrLATpre_Spl[2][0]->GetYaxis()->SetTitle("Corr. Factor");
                for(int i=1;i<num_mesi;i++) {
                        CorrLATpre_Spl[2][i]->SetLineWidth(2);
                        CorrLATpre_Spl[2][i]->SetLineColor(i+1);
			CorrLATpre_Spl[2][i]->SetMarkerStyle(8);
	                CorrLATpre_Spl[2][i]->SetMarkerColor(i+1);
                        leg->AddEntry(CorrLATpre_Spl[2][i],mesi[i].c_str(), "ep");
                        CorrLATpre_Spl[2][i]->Draw("CPsame");
                }
                leg->Draw("same");
        }
	cout<<"7"<<endl;
	c8->Divide(1,3);
	c8->cd(1);
        {gPad->SetGridx();
                gPad->SetGridy();
                gStyle->SetPalette(0);
                TLegend* leg =new TLegend(0.91,0.1,1.0,0.9);
                CorrLAT_LikTOF_Spl[0]->SetLineWidth(2);
                CorrLAT_LikTOF_Spl[0]->SetLineColor(1);
                CorrLAT_LikTOF_Spl[0]->SetMarkerStyle(8);
		CorrLAT_LikTOF_Spl[0]->SetMarkerColor(1);
		leg->AddEntry(CorrLAT_LikTOF_Spl[0],mesi[0].c_str(), "ep");
                CorrLAT_LikTOF_Spl[0]->Draw("ACP");
                CorrLAT_LikTOF_Spl[0]->GetXaxis()->SetTitle("Latitude");
                CorrLAT_LikTOF_Spl[0]->GetYaxis()->SetTitle("Corr. Factor");
                for(int i=1;i<num_mesi;i++) {
                        CorrLAT_LikTOF_Spl[i]->SetLineWidth(2);
                        CorrLAT_LikTOF_Spl[i]->SetLineColor(i+1);
			CorrLAT_LikTOF_Spl[i]->SetMarkerStyle(8);
                	CorrLAT_LikTOF_Spl[i]->SetMarkerColor(i+1);
                        leg->AddEntry(CorrLAT_LikTOF_Spl[i],mesi[i].c_str(), "ep");
                        CorrLAT_LikTOF_Spl[i]->Draw("CPsame");
                }
                leg->Draw("same");
        }
	c8->cd(2);
        {gPad->SetGridx();
                gPad->SetGridy();
                gStyle->SetPalette(0);
                TLegend* leg =new TLegend(0.91,0.1,1.0,0.9);
                CorrLAT_LikNaF_Spl[0]->SetLineWidth(2);
                CorrLAT_LikNaF_Spl[0]->SetLineColor(1);
                CorrLAT_LikNaF_Spl[0]->SetMarkerStyle(8);
                CorrLAT_LikNaF_Spl[0]->SetMarkerColor(1);
                leg->AddEntry(CorrLAT_LikNaF_Spl[0],mesi[0].c_str(), "ep");
                CorrLAT_LikNaF_Spl[0]->Draw("ACP");
                CorrLAT_LikNaF_Spl[0]->GetXaxis()->SetTitle("Latitude");
                CorrLAT_LikNaF_Spl[0]->GetYaxis()->SetTitle("Corr. Factor");
                for(int i=1;i<num_mesi;i++) {
                        CorrLAT_LikNaF_Spl[i]->SetLineWidth(2);
                        CorrLAT_LikNaF_Spl[i]->SetLineColor(i+1);
                        CorrLAT_LikNaF_Spl[i]->SetMarkerStyle(8);
                        CorrLAT_LikNaF_Spl[i]->SetMarkerColor(i+1);
                        leg->AddEntry(CorrLAT_LikNaF_Spl[i],mesi[i].c_str(), "ep");
                        CorrLAT_LikNaF_Spl[i]->Draw("CPsame");
                }
                leg->Draw("same");
        }
	c8->cd(3);
        {gPad->SetGridx();
                gPad->SetGridy();
                gStyle->SetPalette(0);
                TLegend* leg =new TLegend(0.91,0.1,1.0,0.9);
                CorrLAT_LikAgl_Spl[0]->SetLineWidth(2);
                CorrLAT_LikAgl_Spl[0]->SetLineColor(1);
                CorrLAT_LikAgl_Spl[0]->SetMarkerStyle(8);
                CorrLAT_LikAgl_Spl[0]->SetMarkerColor(1);
                leg->AddEntry(CorrLAT_LikAgl_Spl[0],mesi[0].c_str(), "ep");
                CorrLAT_LikAgl_Spl[0]->Draw("ACP");
                CorrLAT_LikAgl_Spl[0]->GetXaxis()->SetTitle("Latitude");
                CorrLAT_LikAgl_Spl[0]->GetYaxis()->SetTitle("Corr. Factor");
                for(int i=1;i<num_mesi;i++) {
                        CorrLAT_LikAgl_Spl[i]->SetLineWidth(2);
                        CorrLAT_LikAgl_Spl[i]->SetLineColor(i+1);
                        CorrLAT_LikAgl_Spl[i]->SetMarkerStyle(8);
                        CorrLAT_LikAgl_Spl[i]->SetMarkerColor(i+1);
                        leg->AddEntry(CorrLAT_LikAgl_Spl[i],mesi[i].c_str(), "ep");
                        CorrLAT_LikAgl_Spl[i]->Draw("CPsame");
                }
                leg->Draw("same");
        }

	cout<<"8"<<endl;
	c9->Divide(1,3);
	c9->cd(1);
	{gPad->SetGridx();
                gPad->SetGridy();
                gStyle->SetPalette(0);
                TLegend* leg =new TLegend(0.91,0.1,1.0,0.9);
                CorrLAT_DistTOF_Spl[0]->SetLineWidth(2);
                CorrLAT_DistTOF_Spl[0]->SetLineColor(1);
		CorrLAT_DistTOF_Spl[0]->SetMarkerStyle(8);
		CorrLAT_DistTOF_Spl[0]->SetMarkerColor(1);
                leg->AddEntry(CorrLAT_DistTOF_Spl[0],mesi[0].c_str(), "ep");
                CorrLAT_DistTOF_Spl[0]->Draw("ACP");
                CorrLAT_DistTOF_Spl[0]->GetXaxis()->SetTitle("Latitude");
                CorrLAT_DistTOF_Spl[0]->GetYaxis()->SetTitle("Corr. Factor");
                for(int i=1;i<num_mesi;i++) {
                        CorrLAT_DistTOF_Spl[i]->SetLineWidth(2);
                        CorrLAT_DistTOF_Spl[i]->SetLineColor(i+1);
                        CorrLAT_DistTOF_Spl[i]->SetMarkerStyle(8);
	                CorrLAT_DistTOF_Spl[i]->SetMarkerColor(i+1);
			leg->AddEntry(CorrLAT_DistTOF_Spl[i],mesi[i].c_str(), "ep");
                        CorrLAT_DistTOF_Spl[i]->Draw("CPsame");
                }
                leg->Draw("same");
        }
	c9->cd(2);
        {gPad->SetGridx();
                gPad->SetGridy();
                gStyle->SetPalette(0);
                TLegend* leg =new TLegend(0.91,0.1,1.0,0.9);
                CorrLAT_DistNaF_Spl[0]->SetLineWidth(2);
                CorrLAT_DistNaF_Spl[0]->SetLineColor(1);
                CorrLAT_DistNaF_Spl[0]->SetMarkerStyle(8);
                CorrLAT_DistNaF_Spl[0]->SetMarkerColor(1);
                leg->AddEntry(CorrLAT_DistNaF_Spl[0],mesi[0].c_str(), "ep");
                CorrLAT_DistNaF_Spl[0]->Draw("ACP");
                CorrLAT_DistNaF_Spl[0]->GetXaxis()->SetTitle("Latitude");
                CorrLAT_DistNaF_Spl[0]->GetYaxis()->SetTitle("Corr. Factor");
                for(int i=1;i<num_mesi;i++) {
                        CorrLAT_DistNaF_Spl[i]->SetLineWidth(2);
                        CorrLAT_DistNaF_Spl[i]->SetLineColor(i+1);
                        CorrLAT_DistNaF_Spl[i]->SetMarkerStyle(8);
                        CorrLAT_DistNaF_Spl[i]->SetMarkerColor(i+1);
                        leg->AddEntry(CorrLAT_DistNaF_Spl[i],mesi[i].c_str(), "ep");
                        CorrLAT_DistNaF_Spl[i]->Draw("CPsame");
                }
                leg->Draw("same");
        }
	c9->cd(3);
        {gPad->SetGridx();
                gPad->SetGridy();
                gStyle->SetPalette(0);
                TLegend* leg =new TLegend(0.91,0.1,1.0,0.9);
                CorrLAT_DistAgl_Spl[0]->SetLineWidth(2);
                CorrLAT_DistAgl_Spl[0]->SetLineColor(1);
                CorrLAT_DistAgl_Spl[0]->SetMarkerStyle(8);
                CorrLAT_DistAgl_Spl[0]->SetMarkerColor(1);
                leg->AddEntry(CorrLAT_DistAgl_Spl[0],mesi[0].c_str(), "ep");
                CorrLAT_DistAgl_Spl[0]->Draw("ACP");
                CorrLAT_DistAgl_Spl[0]->GetXaxis()->SetTitle("Latitude");
                CorrLAT_DistAgl_Spl[0]->GetYaxis()->SetTitle("Corr. Factor");
                for(int i=1;i<num_mesi;i++) {
                        CorrLAT_DistAgl_Spl[i]->SetLineWidth(2);
                        CorrLAT_DistAgl_Spl[i]->SetLineColor(i+1);
                        CorrLAT_DistAgl_Spl[i]->SetMarkerStyle(8);
                        CorrLAT_DistAgl_Spl[i]->SetMarkerColor(i+1);
                        leg->AddEntry(CorrLAT_DistAgl_Spl[i],mesi[i].c_str(), "ep");
                        CorrLAT_DistAgl_Spl[i]->Draw("CPsame");
                }
                leg->Draw("same");
        }

	cout<<"9"<<endl;
	/*c10->cd();
        {gPad->SetGridx();
                gPad->SetGridy();
                gStyle->SetPalette(0);
		gPad->SetLogx();
                TLegend* leg =new TLegend(0.91,0.1,1.0,0.9);
                preDVSMC_P[0][0]->SetLineWidth(2);
                preDVSMC_P[0][0]->SetLineColor(1);
                preDVSMC_P[0][0]->SetMarkerStyle(8);
                preDVSMC_P[0][0]->SetMarkerColor(1);
		preDVSMC_P[0][0]->SetFillStyle(3002);
                leg->AddEntry(preDVSMC_P[0][0],mesi[0].c_str(), "ep");
                preDVSMC_P[0][0]->GetYaxis()->SetRangeUser(0.8,1.25);
		preDVSMC_P[0][0]->Draw("AC4");
                preDVSMC_P[0][0]->GetXaxis()->SetTitle("R [GV]");
                preDVSMC_P[0][0]->GetYaxis()->SetTitle("Corr. Factor");
                for(int i=1;i<num_mesi;i++) {
                        preDVSMC_P[0][i]->SetLineWidth(2);
                        preDVSMC_P[0][i]->SetLineColor(i+1);
                        preDVSMC_P[0][i]->SetMarkerStyle(8);
                        preDVSMC_P[0][i]->SetMarkerColor(i+1);
			preDVSMC_P[0][i]->SetFillStyle(3002);
                        leg->AddEntry(CorrLATpre_Spl[0][i],mesi[i].c_str(), "ep");
                        preDVSMC_P[0][i]->Draw("C4same");
                }
                leg->Draw("same");
        }
	
	c11->cd();
        {gPad->SetGridx();
                gPad->SetGridy();
                gStyle->SetPalette(0);
		gPad->SetLogx();
                TLegend* leg =new TLegend(0.91,0.1,1.0,0.9);
                preDVSMC_P[1][0]->SetLineWidth(2);
                preDVSMC_P[1][0]->SetLineColor(1);
                preDVSMC_P[1][0]->SetMarkerStyle(8);
                preDVSMC_P[1][0]->SetMarkerColor(1);
                preDVSMC_P[1][0]->SetFillStyle(3002);
		leg->AddEntry(preDVSMC_P[1][0],mesi[0].c_str(), "ep");
                preDVSMC_P[1][0]->GetYaxis()->SetRangeUser(0.7,1.15);
		preDVSMC_P[1][0]->Draw("AC4");
                preDVSMC_P[1][0]->GetXaxis()->SetTitle("R [GV]");
                preDVSMC_P[1][0]->GetYaxis()->SetTitle("Corr. Factor");
                for(int i=1;i<num_mesi;i++) {
                        preDVSMC_P[1][i]->SetLineWidth(2);
                        preDVSMC_P[1][i]->SetLineColor(i+1);
                        preDVSMC_P[1][i]->SetMarkerStyle(8);
                        preDVSMC_P[1][i]->SetMarkerColor(i+1);
                        preDVSMC_P[1][i]->SetFillStyle(3002);
			leg->AddEntry(CorrLATpre_Spl[1][i],mesi[i].c_str(), "ep");
                        preDVSMC_P[1][i]->Draw("C4same");
                }
                leg->Draw("same");
        }
	
	c12->cd();
        {gPad->SetGridx();
                gPad->SetGridy();
                gStyle->SetPalette(0);
                gPad->SetLogx();
		TLegend* leg =new TLegend(0.91,0.1,1.0,0.9);
                preDVSMC_P[2][0]->SetLineWidth(2);
                preDVSMC_P[2][0]->SetLineColor(1);
                preDVSMC_P[2][0]->SetMarkerStyle(8);
                preDVSMC_P[2][0]->SetMarkerColor(1);
                preDVSMC_P[2][0]->SetFillStyle(3002);
		leg->AddEntry(preDVSMC_P[2][0],mesi[0].c_str(), "ep");
		preDVSMC_P[2][0]->GetYaxis()->SetRangeUser(0.8,1.15);
                preDVSMC_P[2][0]->Draw("AC4");
                preDVSMC_P[2][0]->GetXaxis()->SetTitle("R [GV]");
                preDVSMC_P[2][0]->GetYaxis()->SetTitle("Corr. Factor");
                for(int i=1;i<num_mesi;i++) {
                        preDVSMC_P[2][i]->SetLineWidth(2);
                        preDVSMC_P[2][i]->SetLineColor(i+1);
                        preDVSMC_P[2][i]->SetMarkerStyle(8);
                        preDVSMC_P[2][i]->SetMarkerColor(i+1);
                        preDVSMC_P[2][i]->SetFillStyle(3002);
			leg->AddEntry(CorrLATpre_Spl[2][i],mesi[i].c_str(), "ep");
                        preDVSMC_P[2][i]->Draw("C4same");
                }
                leg->Draw("same");
        }

	c13->cd();
        {gPad->SetGridx();
                gPad->SetGridy();
                gStyle->SetPalette(0);
                gPad->SetLogx();
		TLegend* leg =new TLegend(0.91,0.1,1.0,0.9);
                LikDVSMC_P[0]->SetLineWidth(2);
                LikDVSMC_P[0]->SetLineColor(1);
		LikDVSMC_P[0]->SetMarkerStyle(8);
                LikDVSMC_P[0]->SetMarkerColor(1);
                LikDVSMC_P[0]->SetFillStyle(3002);

                leg->AddEntry(LikDVSMC_P[0],mesi[0].c_str(), "ep");
                LikDVSMC_P[0]->GetYaxis()->SetRangeUser(0.8,1.4);
		LikDVSMC_P[0]->Draw("AC4");
                LikDVSMC_P[0]->GetXaxis()->SetTitle("R [GV]");
                LikDVSMC_P[0]->GetYaxis()->SetTitle("Corr. Factor");
                for(int i=1;i<num_mesi;i++) {
                        LikDVSMC_P[i]->SetLineWidth(2);
                        LikDVSMC_P[i]->SetLineColor(i+1);
                        LikDVSMC_P[i]->SetMarkerStyle(8);
                	LikDVSMC_P[i]->SetMarkerColor(i+1);
                	LikDVSMC_P[i]->SetFillStyle(3002);
			leg->AddEntry(LikDVSMC_P[i],mesi[i].c_str(), "ep");
                        LikDVSMC_P[i]->Draw("C4same");
                }
                leg->Draw("same");
        }

	c14->cd();
        {gPad->SetGridx();
                gPad->SetGridy();
                gStyle->SetPalette(0);
                gPad->SetLogx();
                TLegend* leg =new TLegend(0.91,0.1,1.0,0.9);
                DistDVSMC_P[0]->SetLineWidth(2);
                DistDVSMC_P[0]->SetLineColor(1);
                DistDVSMC_P[0]->SetMarkerStyle(8);
                DistDVSMC_P[0]->SetMarkerColor(1);
                DistDVSMC_P[0]->SetFillStyle(3002);

                leg->AddEntry(DistDVSMC_P[0],mesi[0].c_str(), "ep");
                DistDVSMC_P[0]->GetYaxis()->SetRangeUser(0.7,1.3);
                DistDVSMC_P[0]->Draw("AC4");
                DistDVSMC_P[0]->GetXaxis()->SetTitle("R [GV]");
                DistDVSMC_P[0]->GetYaxis()->SetTitle("Corr. Factor");
                for(int i=1;i<num_mesi;i++) {
                        DistDVSMC_P[i]->SetLineWidth(2);
                        DistDVSMC_P[i]->SetLineColor(i+1);
                        DistDVSMC_P[i]->SetMarkerStyle(8);
                        DistDVSMC_P[i]->SetMarkerColor(i+1);
                        DistDVSMC_P[i]->SetFillStyle(3002);
                        leg->AddEntry(DistDVSMC_P[i],mesi[i].c_str(), "ep");
                        DistDVSMC_P[i]->Draw("C4same");
                }
                leg->Draw("same");
        }*/
	cout<<"10"<<endl;
	c15->cd();
	double x,y,ey=0;
	double x0,y0=0;
	int u,v=0;
        {       gPad->SetGridx();
                gPad->SetGridy();
                gStyle->SetPalette(0);
                gPad->SetLogx();
                gPad->SetLogy();
		TLegend* leg =new TLegend(0.91,0.1,1.0,0.9);
                for(int i=0;i<num_mesi;i++) {
                P_Fluxesratio[i]=new TGraphErrors();
		for(int p=1;p<42;p++){
                        u=P_Fluxes[i]->GetPoint(p,x,y);
                        v=P_Fluxes[0]->GetPoint(p,x0,y0);
			P_Fluxesratio[i]->SetPoint(p,x,y/y0);
                	}
		}
		P_Fluxesratio[0]->SetLineWidth(2);
                P_Fluxesratio[0]->SetLineColor(1);
                P_Fluxesratio[0]->SetMarkerStyle(8);
                P_Fluxesratio[0]->SetMarkerColor(1);
                P_Fluxesratio[0]->SetFillStyle(3002);
                leg->AddEntry(P_Fluxesratio[0],mesi[0].c_str(), "ep");
                P_Fluxesratio[0]->GetYaxis()->SetRangeUser(0.4,2.3);
                P_Fluxesratio[0]->Draw("AP");
                P_Fluxesratio[0]->GetXaxis()->SetTitle("R [GV]");
                P_Fluxesratio[0]->GetYaxis()->SetTitle("Flux Ratio (vs first month)");
                for(int i=1;i<num_mesi;i++) {
                        P_Fluxesratio[i]->SetLineWidth(2);
                        P_Fluxesratio[i]->SetLineColor(i+1);
                        P_Fluxesratio[i]->SetMarkerStyle(8);
                        P_Fluxesratio[i]->SetMarkerColor(i+1);
                        P_Fluxesratio[i]->SetFillStyle(3002);
                        leg->AddEntry(P_Fluxesratio[i],mesi[i].c_str(), "ep");
                        P_Fluxesratio[i]->Draw("Psame");
                }
                leg->Draw("same");
        }

	c16->cd();
        x,y,ey=0;
        x0,y0=0;
        u,v=0;
        {       gPad->SetGridx();
                gPad->SetGridy();
                gStyle->SetPalette(0);
                gPad->SetLogx();
                gPad->SetLogy();
                TLegend* leg =new TLegend(0.91,0.1,1.0,0.9);
                for(int i=0;i<num_mesi;i++) {
                D_FluxesratioTOF[i]=new TGraphErrors();
                for(int p=1;p<18;p++){
                        u=D_FluxesTOF[i]->GetPoint(p,x,y);
                        v=D_FluxesTOF[0]->GetPoint(p,x0,y0);
                        D_FluxesratioTOF[i]->SetPoint(p,x,y/y0);
                        }
                }
		D_FluxesratioTOF[0]->SetPoint(18,20,1e-3);
                D_FluxesratioTOF[0]->SetLineWidth(2);
                D_FluxesratioTOF[0]->SetMarkerStyle(8);
		D_FluxesratioTOF[0]->SetLineColor(1);
                D_FluxesratioTOF[0]->SetMarkerColor(1);
                D_FluxesratioTOF[0]->SetFillStyle(3002);
                leg->AddEntry(D_FluxesratioTOF[0],mesi[0].c_str(), "ep");
                D_FluxesratioTOF[0]->GetYaxis()->SetRangeUser(0.4,2.3);
                D_FluxesratioTOF[0]->Draw("AP");
                D_FluxesratioTOF[0]->GetXaxis()->SetTitle("R [GV]");
                D_FluxesratioTOF[0]->GetYaxis()->SetTitle("Flux Ratio (vs first month)");
                for(int i=1;i<num_mesi;i++) {
                        D_FluxesratioTOF[i]->SetLineWidth(2);
                        D_FluxesratioTOF[i]->SetLineColor(i+1);
                        D_FluxesratioTOF[i]->SetMarkerColor(i+1);
                        D_FluxesratioTOF[i]->SetMarkerStyle(8);
			D_FluxesratioTOF[i]->SetFillStyle(3002);
                        leg->AddEntry(D_FluxesratioTOF[i],mesi[i].c_str(), "ep");
                        D_FluxesratioTOF[i]->Draw("PLsame");
                }
                leg->Draw("same");
        	for(int i=0;i<num_mesi;i++) {
                D_FluxesratioNaF[i]=new TGraphErrors();
                for(int p=1;p<18;p++){
                        u=D_FluxesNaF[i]->GetPoint(p,x,y);
                        v=D_FluxesNaF[0]->GetPoint(p,x0,y0);
                        D_FluxesratioNaF[i]->SetPoint(p,x,y/y0);
			}
                }
                D_FluxesratioNaF[0]->SetLineWidth(2);
                D_FluxesratioNaF[0]->SetMarkerStyle(4);
                D_FluxesratioNaF[0]->SetLineColor(1);
                D_FluxesratioNaF[0]->SetMarkerColor(1);
                D_FluxesratioNaF[0]->SetFillStyle(3002);
		D_FluxesratioNaF[0]->Draw("Psame");
		for(int i=1;i<num_mesi;i++) {
                        D_FluxesratioNaF[i]->SetLineWidth(2);
                        D_FluxesratioNaF[i]->SetLineColor(i+1);
                        D_FluxesratioNaF[i]->SetMarkerColor(i+1);
                        D_FluxesratioNaF[i]->SetMarkerStyle(4);
                        D_FluxesratioNaF[i]->SetFillStyle(3002);
                        D_FluxesratioNaF[i]->Draw("PLsame");
                }
		 for(int i=0;i<num_mesi;i++) {
                D_FluxesratioAgl[i]=new TGraphErrors();
                for(int p=1;p<18;p++){
                        u=D_FluxesAgl[i]->GetPoint(p,x,y);
                        v=D_FluxesAgl[0]->GetPoint(p,x0,y0);
                        D_FluxesratioAgl[i]->SetPoint(p,x,y/y0);
                        }
                }
                D_FluxesratioAgl[0]->SetLineWidth(2);
                D_FluxesratioAgl[0]->SetMarkerStyle(3);
                D_FluxesratioAgl[0]->SetLineColor(1);
                D_FluxesratioAgl[0]->SetMarkerColor(1);
                D_FluxesratioAgl[0]->SetFillStyle(3002);
                D_FluxesratioAgl[0]->Draw("Psame");
                for(int i=1;i<num_mesi;i++) {
                        D_FluxesratioAgl[i]->SetLineWidth(2);
                        D_FluxesratioAgl[i]->SetLineColor(i+1);
                        D_FluxesratioAgl[i]->SetMarkerColor(i+1);
                        D_FluxesratioAgl[i]->SetMarkerStyle(3);
                        D_FluxesratioAgl[i]->SetFillStyle(3002);
                        D_FluxesratioAgl[i]->Draw("PLsame");
                }

	
	}

	


	cout<<"*************** OUTPUT **************************"<<endl;
	nomefile="./TimeAnalysis.root";
	TFile *f_out=new TFile(nomefile.c_str(), "RECREATE");
	f_out->mkdir("Calibrations");
	f_out->mkdir("Lat. Dependence");
	f_out->mkdir("Data vs MC");
	f_out->mkdir("P Fluxes");

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
	/*f_out->cd("Data vs MC");
	c10->Write();
	c11->Write();
        c12->Write();
	c13->Write();
	c14->Write();*/
	f_out->cd("P Fluxes");
	c15->Write();
	c16->Write();	
	f_out->Write();
	f_out->Close();
	
	return 0;
}



