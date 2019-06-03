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

#include "../include/EffCorr.h"


void DrawCorrection(EffCorr * Correction, FileSaver Plots, Binning Bins, std::string CorrName,std::string RangeName,TVirtualPad * c3 = 0x0, TVirtualPad * c3_up = 0x0,  TVirtualPad * c3_do = 0x0, TVirtualPad * c4 =0x0,TVirtualPad * c4_up =0x0,TVirtualPad * c4_do =0x0, float rangemin=-1, float rangemax=-1,bool skipleg=false) {

	std::string latitudes[10]={"0.0<|#theta_{M}|<0.2","0.2<|#theta_{M}|<0.3","0.3<|#theta_{M}|<0.4","0.4<|#theta_{M}|<0.5","0.5<|#theta_{M}|<0.6","0.6<|#theta_{M}|<0.7","0.7<|#theta_{M}|<0.8","0.8<|#theta_{M}|<0.9","0.9<|#theta_{M}|<1","|#theta_{M}|>1"};


	if(!c3||!c3_up) {
		c3 = new TCanvas((CorrName + " Latitude").c_str());
	        c3->SetCanvasSize(2500,1920);
       		c3->cd();
	 	c3_up = new TPad("upperPad", "upperPad",0.0,0.3,1.0,1.0);
        	c3_up->Draw();
        	c3_do = new TPad("lowerPad", "lowerPad",0.0,0.0,1.0,0.3);
        	c3_do->Draw();
	}

	TH1F * MCEff = (TH1F*)Correction->GetMCEfficiency();
	TH1F * MCEff2 = (TH1F*)Correction->GetMCEfficiency2();
	TH1F * MCEffPID = (TH1F*)Correction->GetMCEfficiency_noPID();
	cout<<MCEff<<endl;
	if(rangemin==-1&&rangemax==-1) {
		rangemin =0.8*Bins.EkPerMasBins()[0];
		rangemax =1.2*Bins.EkPerMasBins()[Bins.EkPerMasBins().size()-1];
	}

	c3_up->cd();
	gPad->SetLogx();
	gPad->SetGridx();
	for(int lat=0;lat<10;lat++){
		TH1F * LatEff = (TH1F*)Correction->GetEfficiencyLat(lat);
		cout<<LatEff<<endl;
		PlotTH1FintoGraph(gPad,Bins, LatEff,"R [GV]", "Efficiency",55+5*lat,false,"Psame",rangemin,rangemax,0.5*MCEff->GetBinContent(2),1.3*MCEff->GetBinContent(13),latitudes[lat].c_str(),8,skipleg,true);
	}

	c3_do->cd();
	gPad->SetLogx();
	gPad->SetGridx();
	for(int lat=0;lat<10;lat++){ 
	PlotTH1FintoGraph(gPad,Bins, (TH1F*)Correction->GetCorrectionLat(lat),"", "Eff. (Data/MC)",55+5*lat,false,"Psame",rangemin,rangemax,0.65*(Correction->GetGlobCorrection_noPID()->GetBinContent(10)),1.2*(Correction->GetGlobCorrection_noPID()->GetBinContent(10)),latitudes[lat].c_str(),8,true,false);
	}

	if(!c4||!c4_up) {
		c4 = new TCanvas((CorrName + " Global").c_str());
	        c4->SetCanvasSize(2500,1920);
       		c4->cd();
	 	c4_up = new TPad("upperPad", "upperPad",0.0,0.3,1.0,1.0);
        	c4_up->Draw();
        	c4_do = new TPad("lowerPad", "lowerPad",0.0,0.0,1.0,0.3);
        	c4_do->Draw();
	}

	TH1F * DataEff = (TH1F*)Correction->GetGlobEfficiency();
	cout<<DataEff<<endl;
	c4_up->cd();
	gPad->SetLogx();
	gPad->SetGridx();
	PlotTH1FintoGraph(gPad,Bins, MCEff,"R [GV]", "Efficiency",2,false,"e4Psame",rangemin,rangemax,0.5*DataEff->GetBinContent(2),1.3*MCEff->GetBinContent(13),"MC (P)",8,skipleg);
	PlotTH1FintoGraph(gPad,Bins, MCEff2,"R [GV]", "Efficiency",4,false,"e4Psame",rangemin,rangemax,0.5*DataEff->GetBinContent(2),1.3*MCEff->GetBinContent(13),"MC (D)",8,skipleg);
	PlotTH1FintoGraph(gPad,Bins, MCEffPID,"R [GV]", "Efficiency",2,false,"e4Psame",rangemin,rangemax,0.5*DataEff->GetBinContent(2),1.3*MCEff->GetBinContent(13),"MC (P - PID)",3,skipleg);
	
	PlotTH1FintoGraph(gPad,Bins, DataEff,"R [GV]", "Efficiency",1,false,"Psame",rangemin,rangemax,0.5*DataEff->GetBinContent(2),1.3*MCEff->GetBinContent(13),"Data",8,skipleg);



	c4_do->cd();
	gPad->SetLogx();
	gPad->SetGridx();
	PlotTH1FintoGraph(gPad,Bins, (TH1F*) Correction->GetGlobCorrection(),"", "Eff. (Data/MC)",1,false,"Psame",rangemin,rangemax,0.65*(Correction->GetGlobCorrection_noPID()->GetBinContent(10)),1.2*(Correction->GetGlobCorrection_noPID()->GetBinContent(10)),"Data",8,true,false);	
	
	Plots.Add(c3);
	Plots.Add(c4);
	Plots.writeObjsInFolder(("Eff. Correction/" +CorrName+"/"+RangeName).c_str());

}


void DrawCorrection(EffCorr * Correction, FileSaver Plots, Binning Bins, std::string CorrName,std::string RangeName,float rangemin=-1, float rangemax=-1,bool skipleg=false) {
	DrawCorrection(Correction,Plots,Bins,CorrName,RangeName,0x0,0x0,0x0,0x0,0x0,0x0,rangemin,rangemax,skipleg);
}


void DrawAllRangeCorrection( EffCorr * CorrectionTOF,EffCorr * CorrectionNaF,EffCorr * CorrectionAgl,FileSaver Plots,std::string CorrName){

	TCanvas * c3 = new TCanvas((CorrName + " Latitude").c_str()); 
        c3->SetCanvasSize(2000,1500);
	TPad * c3_up = new TPad("upperPad", "upperPad",0.0,0.3,1.0,1.0);
        c3_up->Draw();
        TPad * c3_do = new TPad("lowerPad", "lowerPad",0.0,0.0,1.0,0.3);
        c3_do->Draw();

	TCanvas * c4 = new TCanvas((CorrName + " Global").c_str()); 
        c4->SetCanvasSize(2000,1500);
	TPad * c4_up = new TPad("upperPad", "upperPad",0.0,0.3,1.0,1.0);
        c4_up->Draw();
        TPad * c4_do = new TPad("lowerPad", "lowerPad",0.0,0.0,1.0,0.3);
        c4_do->Draw();


	DrawCorrection(CorrectionTOF,Plots,ForEffCorr,CorrName,"TOF",c3,c3_up,c3_do,c4,c4_up,c4_do,0.4,12);
	DrawCorrection(CorrectionNaF,Plots,ForEffCorr,CorrName,"NaF",c3,c3_up,c3_do,c4,c4_up,c4_do,0.4,12,true);
	DrawCorrection(CorrectionAgl,Plots,ForEffCorr,CorrName,"Agl",c3,c3_up,c3_do,c4,c4_up,c4_do,0.4,12,true);
	
	Plots.Add(c3);
	Plots.Add(c4);
	Plots.writeObjsInFolder(("Eff. Correction/" +CorrName).c_str());
}


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
	cout<<"**************************** PLOTTING ***************************************"<<endl;

	EffCorr * HEPPresEffCorr = new EffCorr(finalHistos,"HEPPresEffCorr","HEPPresEffCorr",ForEffCorr,"IsPositive&IsMinimumBias&IsLooseCharge1","IsPositive&IsMinimumBias&IsLooseCharge1&IsGolden","","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * HEPQualEffCorr = new EffCorr(finalHistos,"HEPQualEffCorr","HEPQualEffCorr",ForEffCorr,"IsPositive&IsMinimumBias&IsLooseCharge1","IsPositive&IsMinimumBias&IsLooseCharge1&IsGolden","","IsPurePMC","IsPureDMC","IsDeutonMC");
	
	std::string before;
        std::string after;
        before = "";
        after  = ""; 
	EffCorr * TriggerEffCorr_HE = new EffCorr(finalHistos,"TriggerEffCorr_HE","Trigger Eff. Corr",ForEffCorr,before,after,"IsPrimary",    "IsProtonMC","IsPureDMC","IsDeutonMC"); 
//	EffCorr * Trigger2EffCorr_HE = new EffCorr(finalHistos,"Trigger2EffCorr_HE","Trigger2 Eff. Corr",ForEffCorr,before,after,"IsPrimary",    "IsProtonMC","IsPureDMC","IsDeutonMC"); 

	EffCorr * L1PickUpEffCorr_HE = new EffCorr(finalHistos,"L1PickUpEffCorr_HE","L1PickUp Eff. Corr",ForEffCorr,before,after,"IsPrimary",    "IsPurePMC","IsPureDMC","IsDeutonMC"); 
	EffCorr * TrackerEffCorr_HE = new EffCorr(finalHistos,"TrackerEffCorr_HE","Tracker Eff. Corr",ForEffCorr,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodChi_HE  = new EffCorr(finalHistos,"GoodChiEffCorr_HE","GoodChi Eff. Corr",ForEffCorr,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");

	EffCorr * GoodUtof_HE  = new EffCorr(finalHistos,"GoodUtofEffCorr_HE","GoodUtof Eff. Corr",ForEffCorr,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");

	EffCorr * GoodLtof_HE  = new EffCorr(finalHistos,"GoodLTOFEffCorr_HE","GoodLtof Eff. Corr",ForEffCorr,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");

	EffCorr * Good1Track_HE  = new EffCorr(finalHistos,"Good1TrackEffCorr_HE","Good1Track Eff. Corr",ForEffCorr,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");

	EffCorr * GoodQTrack_HE  = new EffCorr(finalHistos,"GoodQTrackEffCorr_HE","GoodQTrack Eff. Corr",ForEffCorr,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");

	EffCorr * GoodTime_TOF = new EffCorr(finalHistos,"GoodTimeEffCorr_TOF","GoodTime Eff. Corr",ForEffCorr,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");

	EffCorr * RICHEffCorr_NaF = new EffCorr(finalHistos,"RICHCorrection_NaF","RICH Eff. Corr",ForEffCorr,before,(after+"&IsFromNaF").c_str(),"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * RICHEffCorr_Agl = new EffCorr(finalHistos,"RICHCorrection_Agl","RICH Eff. Corr",ForEffCorr,before,(after+"&IsFromAgl").c_str(),"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * RICHQualEffCorr_NaF = new EffCorr(finalHistos,"RICHQualCorrection_NaF","RICH Qual Eff. Corr",ForEffCorr,(before+"&IsFromNaF").c_str(),(after+"&IsFromNaF&RICHBDTCut").c_str(),"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * RICHQualEffCorr_Agl = new EffCorr(finalHistos,"RICHqualCorrection_Agl","RICH Qual. Eff. Corr",ForEffCorr,(before+"&IsFromNaF").c_str(),(after+"&IsFromAgl&RICHBDTCut").c_str(),"","IsPurePMC","IsPureDMC","IsDeutonMC");

	DrawCorrection(TriggerEffCorr_HE,Plots,ForEffCorr,"TriggerEffCorr","HE",0.5,70);
	DrawCorrection(L1PickUpEffCorr_HE,Plots,ForEffCorr,"L1PickUpEffCorr","HE",0.5,70);
	DrawCorrection(TrackerEffCorr_HE,Plots,ForEffCorr,"TrackerEffCorr","HE",0.5,70);
	DrawCorrection(GoodChi_HE,Plots,ForEffCorr,"GoodChi_HE_EffCorr","HE",0.5,70);
	DrawCorrection(GoodQTrack_HE,Plots,ForEffCorr,"GoodQTrack_HE_EffCorr","HE",0.5,70);

	DrawCorrection(GoodUtof_HE,Plots,ForEffCorr,"GoodUtof_HE_EffCorr","HE",0.5,70);
	DrawCorrection(GoodLtof_HE,Plots,ForEffCorr,"GoodLtof_HE_EffCorr","HE",0.5,70);
	DrawCorrection(Good1Track_HE,Plots,ForEffCorr,"Good1Track_HE_EffCorr","HE",0.5,70);
	
	DrawCorrection(GoodTime_TOF,Plots, ForEffCorr, "GoodTime","TOF",0.5,70);
	DrawCorrection(RICHEffCorr_NaF,Plots, ForEffCorr, "RICHEffCorr","NaF",0.5,70);
	DrawCorrection(RICHEffCorr_Agl,Plots, ForEffCorr, "RICHEffCorr","Agl",0.5,70);
	DrawCorrection(RICHQualEffCorr_NaF,Plots, ForEffCorr, "RICHQualEffCorr","NaF",0.5,70);
	DrawCorrection(RICHQualEffCorr_Agl,Plots, ForEffCorr, "RICHQualEffCorr","Agl",0.5,70);
	

	return 0;
}

