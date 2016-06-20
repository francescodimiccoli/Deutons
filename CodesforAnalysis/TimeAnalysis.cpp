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

int colorbase=30;

void DrawCalibration(TVirtualPad * c, TSpline3 * Corr[],const string &name,const string &xaxisname,const string &yaxisname){
	c -> cd();
	gPad->SetGridx();
	gPad->SetGridy();
	gStyle->SetPalette(1);
	TLegend* leg =new TLegend(0.91,0.1,1.0,0.9);
	Corr[0]->SetLineWidth(2);
	Corr[0]->SetLineColor(1);
	leg->AddEntry(Corr[0],mesi[0].c_str(), "ep");
	Corr[0]->Draw("");
	TH1F *frame=c->DrawFrame(0.4,0.85,1,1.1,name.c_str());
	frame->GetXaxis()->SetTitle(xaxisname.c_str());
	frame->GetYaxis()->SetTitle(yaxisname.c_str());
	for(int i=1;i<num_mesi;i++) {
		Corr[i]->SetLineWidth(2);
		Corr[i]->SetLineColor(colorbase+i);
		Corr[i]->SetMarkerColor(colorbase+i);
		Corr[i]->SetMarkerStyle(8);
		leg->AddEntry(Corr[i],mesi[i].c_str(), "ep");
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
	leg->AddEntry(Corr[0],mesi[0].c_str(), "ep");
	Corr[0]->Draw("AC");
	Corr[0]->GetXaxis()->SetTitle(xaxisname.c_str());
	Corr[0]->GetYaxis()->SetTitle(yaxisname.c_str());
	Corr[0]->GetYaxis()->SetRangeUser(0.7,1.3);
	for(int i=1;i<num_mesi;i++) {
		Corr[i]->SetLineWidth(2);
		Corr[i]->SetLineColor(colorbase+i);
		Corr[i]->SetMarkerStyle(8);
		Corr[i]->SetMarkerColor(colorbase+i);
		leg->AddEntry(Corr[i],mesi[i].c_str(), "ep");
		Corr[i]->Draw("Csame");
	}
	leg->Draw("same");
}


void DrawFluxRatio(TVirtualPad * c, TGraphErrors * Fluxes[],const string &name,const string &xaxisname,const string &yaxisname){

	gStyle->SetPalette(0);
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
                        v=Fluxes[0]->GetPoint(p,x0,y0);
                        ey=Fluxes[i]->GetErrorY(p);
                        Fluxesratio[i]->SetPoint(p,x,y/y0);
                        Fluxesratio[i]->SetPointError(p,0,ey/y0);
                        }
                }
                Fluxesratio[0]->SetLineWidth(2);
                Fluxesratio[0]->SetLineColor(1);
                Fluxesratio[0]->SetMarkerStyle(8);
                Fluxesratio[0]->SetMarkerColor(1);
                Fluxesratio[0]->SetFillStyle(3002);
                leg->AddEntry(Fluxesratio[0],mesi[0].c_str(), "ep");
                Fluxesratio[0]->GetYaxis()->SetRangeUser(0.2,1.3);
                Fluxesratio[0]->Draw("AP");
                Fluxesratio[0]->GetXaxis()->SetTitle(xaxisname.c_str());
                Fluxesratio[0]->GetYaxis()->SetTitle(yaxisname.c_str());
                for(int i=1;i<num_mesi;i++) {
                        Fluxesratio[i]->SetLineWidth(2);
                        Fluxesratio[i]->SetLineColor(38+i);
                        Fluxesratio[i]->SetMarkerStyle(8);
                        Fluxesratio[i]->SetMarkerColor(colorbase+i);
                        Fluxesratio[i]->SetFillStyle(3002);
                        leg->AddEntry(Fluxesratio[i],mesi[i].c_str(), "ep");
                        Fluxesratio[i]->Draw("Psame");
                }
                leg->Draw("same");
        }

}

int main()
{
	int colorbase = 30;

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
	TH1F * TrackerEff[40];
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
		preDVSMC_P[0][i] =  (TGraphErrors *) result[i]->Get("Export/DvsMC/DvsMC: Matching TOF_R");
                preDVSMC_P[1][i] =  (TGraphErrors *) result[i]->Get("Export/DvsMC/DvsMC: Chi^2 R_R");
                preDVSMC_P[2][i] =  (TGraphErrors *) result[i]->Get("Export/DvsMC/DvsMC: 1 Tr. Track_R");
		LikDVSMC_P[i] =  (TGraphErrors *) result[i]->Get("Export/DvsMC/LikDVSMC_P_Graph");
                DistDVSMC_P[i] =  (TGraphErrors *) result[i]->Get("Export/DvsMC/DistDVSMC_P_Graph");
		TrackerEff[i]  =  (TH1F *) 	   result[i]->Get("Export/TrakerEfficiencyData");	
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
	TCanvas *c10=new TCanvas("DvsMC: Matching TOF");
        TCanvas *c11=new TCanvas("DvsMC: Chi^2");
        TCanvas *c12=new TCanvas("DvsMC: 1 Tr. Track");
	TCanvas *c13=new TCanvas("DvsMC: Likelihood");
        TCanvas *c14=new TCanvas("DvsMC: Distance");
	TCanvas *c  =new TCanvas("Tracker Efficiency");
	TCanvas *c15=new TCanvas("Protons Fluxes");
	TCanvas *c16=new TCanvas("Deutons Fluxes");
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

	c10->cd();
	DrawCorrection(c10,preDVSMC_P[0],"","R [GV]","Corr. Factor");
		
	c11->cd();
	DrawCorrection(c11,preDVSMC_P[1],"","R [GV]","Corr. Factor");
	
	c12->cd();
	DrawCorrection(c12,preDVSMC_P[2],"","R [GV]","Corr. Factor");
	
	c13->cd();
	DrawCorrection(c12,LikDVSMC_P,"","R [GV]","Corr. Factor");
	
	c14->cd();
	DrawCorrection(c12,DistDVSMC_P,"","R [GV]","Corr. Factor");
	
	c->cd();
	gPad->SetGridx();
	gPad->SetGridy();
	double effy,effey=0;
	TH1F * Trackeff_time = new TH1F("","",num_mesi,0,num_mesi);
        for(int i=0;i<num_mesi;i++) {
                cout<<TrackerEff[i]<<endl;
		effy = TrackerEff[i]-> GetBinContent(30);
		Trackeff_time -> SetBinContent(i+1,effy);
                effey = TrackerEff[i]-> GetBinError(30);
                Trackeff_time -> SetBinError(i+1,effey);
                Trackeff_time -> GetXaxis() -> SetBinLabel(i+1,mesi[i].c_str());
        }
	Trackeff_time -> SetMarkerStyle(8);
	Trackeff_time -> LabelsOption("v");
	Trackeff_time -> GetXaxis() -> SetLabelSize(0.085);
	Trackeff_time -> SetMarkerColor(1);
	Trackeff_time -> Draw();


	c15-> Divide(1,2);
	c15->cd(1);
	DrawFluxRatio(c15,P_Fluxes,"R [GV]","Proton Flux","Proton Flux (norm. to first month)");

	c16->cd();
        DrawFluxRatio(c16,D_FluxesTOF,"Kin. En./nucl.","Deuton Flux","Deutons Flux (norm. to first month)");	
	DrawFluxRatio(c16,D_FluxesNaF,"Kin. En./nucl.","Deuton Flux","Deutons Flux (norm. to first month)");
	DrawFluxRatio(c16,D_FluxesAgl,"Kin. En./nucl.","Deuton Flux","Deutons Flux (norm. to first month)");
	
	c17->Divide(1,3);

	TH1F * Time_depD1 = new TH1F("","",num_mesi,0,num_mesi);
	TH1F * Time_depP1 = new TH1F("","",num_mesi,0,num_mesi);
	TH1F * Time_depD2 = new TH1F("","",num_mesi,0,num_mesi);
        TH1F * Time_depP2 = new TH1F("","",num_mesi,0,num_mesi);
	TH1F * Time_depD3 = new TH1F("","",num_mesi,0,num_mesi);
        TH1F * Time_depP3 = new TH1F("","",num_mesi,0,num_mesi);
	
	TH1F * Time_depP4 = new TH1F("","",num_mesi,0,num_mesi);
	TH1F * errorP4 = new TH1F("","",num_mesi,0,num_mesi);
	

	double x01,x02,x03,x04;
	double y01,y02,y03,y04;
	double x1,x2,x3,x4;
	double y1,y2,y3,y4;
	double ey1,ey2,ey3,ey4;

	for(int i=0;i<num_mesi;i++) {
		D_FluxesTOF[i]->GetPoint(4,x1,y1);
		D_FluxesTOF[0]->GetPoint(4,x01,y01);
		
		D_FluxesTOF[i]->GetPoint(12,x2,y2);
                D_FluxesTOF[0]->GetPoint(12,x02,y02);
		
		D_FluxesAgl[i]->GetPoint(9,x3,y3);
                D_FluxesAgl[0]->GetPoint(9,x03,y03);

		Time_depD1 -> SetBinContent(i+1,y1/y01);
		Time_depD2 -> SetBinContent(i+1,y2/y02);		
		Time_depD3 -> SetBinContent(i+1,y3/y03);		

		ey1 = D_FluxesTOF[i]->GetErrorY(4);
		ey2 = D_FluxesTOF[i]->GetErrorY(12);
		ey3 = D_FluxesAgl[i]->GetErrorY(9);
		
		Time_depD1 -> SetBinError(i+1,ey1/y01);
		Time_depD2 -> SetBinError(i+1,ey2/y02);
		Time_depD3 -> SetBinError(i+1,ey3/y03);

		Time_depD1 -> GetXaxis() -> SetBinLabel(i+1,mesi[i].c_str());
		Time_depD2 -> GetXaxis() -> SetBinLabel(i+1,mesi[i].c_str());
		Time_depD3 -> GetXaxis() -> SetBinLabel(i+1,mesi[i].c_str());
	}
	for(int i=0;i<num_mesi;i++) {
                P_Fluxes[i]->GetPoint(6,x1,y1);
                P_Fluxes[0]->GetPoint(6,x01,y01);
                
		P_Fluxes[i]->GetPoint(11,x2,y2);
                P_Fluxes[0]->GetPoint(11,x02,y02);
       
		P_Fluxes[i]->GetPoint(25,x3,y3);
                P_Fluxes[0]->GetPoint(25,x03,y03);
                  
		P_Fluxes[i]->GetPoint(40,x4,y4);
                P_Fluxes[0]->GetPoint(40,x04,y04);
                
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
		errorP4    -> SetBinError(i+1,2*ey4/y04);
		Time_depP4 -> GetXaxis() -> SetBinLabel(i+1,mesi[i].c_str());

		
	}
	Time_depD1 -> SetMarkerStyle(8);
	Time_depP1 -> SetMarkerStyle(8);
	Time_depD1 -> SetMarkerColor(4);
        Time_depP1 -> SetMarkerColor(2);

	Time_depD2 -> SetMarkerStyle(8);
        Time_depP2 -> SetMarkerStyle(8);
        Time_depD2 -> SetMarkerColor(4);
        Time_depP2 -> SetMarkerColor(2);

	Time_depD3 -> SetMarkerStyle(8);
        Time_depP3 -> SetMarkerStyle(8);
        Time_depD3 -> SetMarkerColor(4);
        Time_depP3 -> SetMarkerColor(2);

	Time_depP4 -> SetMarkerStyle(8);
        Time_depP4 -> SetMarkerColor(2);
        errorP4	   -> SetFillColor(1);
	

	gStyle -> SetTitleFontSize(0.3);
	gStyle -> SetLegendFillColor(0);	
	c17->cd(1);
        gPad->SetGridx();
        gPad->SetGridy();
	Time_depD1 -> GetYaxis() -> SetRangeUser(0.65,1.2);
	Time_depD1 -> LabelsOption("h");
	Time_depD1 -> SetTitle("R = 1 GV");
	Time_depD1 -> GetYaxis() -> SetLabelSize(0.085);
	Time_depD1 -> GetXaxis() -> SetLabelSize(0.085);
	Time_depD1 -> Draw();
	Time_depP1 -> Draw("same");	
	 TLegend* leg =new TLegend(0.7,0.1,0.9,0.3);
         leg->AddEntry(Time_depD1,"Deuterons", "ep");
	leg->AddEntry(Time_depP1,"Proons", "ep");
	leg ->Draw("same");

	c17->cd(2);
        gPad->SetGridx();
        gPad->SetGridy();
        Time_depD2 -> GetYaxis() -> SetRangeUser(0.65,1.2);
	Time_depD2 -> SetTitle("R = 2 GV");
	Time_depD2 -> LabelsOption("h");
	Time_depD2 -> GetYaxis() -> SetLabelSize(0.085);
	Time_depD2 -> GetXaxis() -> SetLabelSize(0.085);
	Time_depD2 -> Draw();
        Time_depP2 -> Draw("same");

	c17->cd(3);
        gPad->SetGridx();
        gPad->SetGridy();
        Time_depD3 -> GetYaxis() -> SetRangeUser(0.65,1.2);
	Time_depD3 -> SetTitle("R = 10 GV");
	Time_depD3 -> LabelsOption("h");
	Time_depD3 -> GetYaxis() -> SetLabelSize(0.085);
	Time_depD3 -> GetXaxis() -> SetLabelSize(0.085);
	Time_depD3 -> Draw();
        Time_depP3 -> Draw("same");

	c15->cd(2);
        gPad->SetGridx();
        gPad->SetGridy();
        Time_depP4 -> GetYaxis() -> SetRangeUser(0.9,1.1);
	Time_depP4 -> SetTitle("High Energy Flux stability (R=80 GV)");
	Time_depP4 -> GetYaxis() -> SetLabelSize(0.085);
	Time_depP4 -> GetXaxis() -> SetLabelSize(0.085);
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
	c16->Write();	
	c17->Write();
	f_out->Write();
	f_out->Close();
	
	return 0;
}



