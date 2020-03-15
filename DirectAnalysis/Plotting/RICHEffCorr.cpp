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
		if(Correction->IsEkin()){
			rangemin =0.8*Bins.EkPerMasBins()[0];
			rangemax =1.2*Bins.EkPerMasBins()[Bins.EkPerMasBins().size()-1];
		}
		else{
			rangemin =0.8*Bins.GetBinCenter(0);
			rangemax =1.2*Bins.GetBinCenter(Bins.EkPerMasBins().size()-1);
		}	
	}

	c3_up->cd();
	gPad->SetLogx();
	gPad->SetGridx();
	for(int lat=0;lat<10;lat++){
		TH1F * LatEff = (TH1F*)Correction->GetEfficiencyLat(lat);
		cout<<LatEff<<endl;
		if(Correction->IsEkin()) PlotTH1FintoGraph(gPad,Bins, LatEff,"Ekin [GeV/n]", "Efficiency",55+5*lat,true,"Psame",rangemin,rangemax,0.5*MCEff->GetBinContent(2),1.3*MCEff->GetBinContent(13),latitudes[lat].c_str(),8,skipleg,true);
		else PlotTH1FintoGraph(gPad,Bins, LatEff,"R [GV]", "Efficiency",55+5*lat,false,"Psame",rangemin,rangemax,0.5*MCEff->GetBinContent(2),1.3*MCEff->GetBinContent(13),latitudes[lat].c_str(),8,skipleg,true);
	}

	c3_do->cd();
	gPad->SetLogx();
	gPad->SetGridx();
	for(int lat=0;lat<10;lat++){ 
	if(Correction->IsEkin()) PlotTH1FintoGraph(gPad,Bins, (TH1F*)Correction->GetCorrectionLat(lat),"", "Eff. (Data/MC)",55+5*lat,true,"Psame",rangemin,rangemax,0.65*(Correction->GetGlobCorrection_noPID()->GetBinContent(10)),1.2*(Correction->GetGlobCorrection_noPID()->GetBinContent(10)),latitudes[lat].c_str(),8,true,false);
	else PlotTH1FintoGraph(gPad,Bins, (TH1F*)Correction->GetCorrectionLat(lat),"", "Eff. (Data/MC)",55+5*lat,false,"Psame",rangemin,rangemax,0.65*(Correction->GetGlobCorrection_noPID()->GetBinContent(10)),1.2*(Correction->GetGlobCorrection_noPID()->GetBinContent(10)),latitudes[lat].c_str(),8,true,false);
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
	if(Correction->IsEkin()){
		PlotTH1FintoGraph(gPad,Bins, MCEff,"Ekin [GeV/n]", "Efficiency",2,true,"Psame",rangemin,rangemax,0.5*DataEff->GetBinContent(2),1.3*MCEff->GetBinContent(MCEff->GetMaximumBin()),"MC (P)",8,skipleg);
		PlotTH1FintoGraph(gPad,Bins, MCEff2,"Ekin [GeV/n]", "Efficiency",4,true,"Psame",rangemin,rangemax,0.5*DataEff->GetBinContent(2),1.3*MCEff->GetBinContent(MCEff->GetMaximumBin()),"MC (D)",8,skipleg);
	//	PlotTH1FintoGraph(gPad,Bins, MCEffPID,"Ekin [GeV/n]", "Efficiency",2,true,"Psame",rangemin,rangemax,0.5*DataEff->GetBinContent(2),1.3*MCEff->GetBinContent(MCEff->GetMaximumBin()),"MC (P - PID)",3,skipleg);
		PlotTH1FintoGraph(gPad,Bins, DataEff,"Ekin [GeV/n]", "Efficiency",1,true,"Psame",rangemin,rangemax,0.5*DataEff->GetBinContent(2),1.3*MCEff->GetBinContent(MCEff->GetMaximumBin()),"Data",8,skipleg);
	}
	else{
		PlotTH1FintoGraph(gPad,Bins, MCEff,"R [GV]", "Efficiency",2,false,"Psame",rangemin,rangemax,0.5*DataEff->GetBinContent(2),1.3*MCEff->GetBinContent(MCEff->GetMaximumBin()),"MC (P)",8,skipleg);
		PlotTH1FintoGraph(gPad,Bins, MCEff2,"R [GV]", "Efficiency",4,false,"Psame",rangemin,rangemax,0.5*DataEff->GetBinContent(2),1.3*MCEff->GetBinContent(MCEff->GetMaximumBin()),"MC (D)",8,skipleg);
	//	PlotTH1FintoGraph(gPad,Bins, MCEffPID,"R [GV]", "Efficiency",2,false,"Psame",rangemin,rangemax,0.5*DataEff->GetBinContent(2),1.3*MCEff->GetBinContent(MCEff->GetMaximumBin()),"MC (P - PID)",3,skipleg);
		PlotTH1FintoGraph(gPad,Bins, DataEff,"R [GV]", "Efficiency",1,false,"Psame",rangemin,rangemax,0.5*DataEff->GetBinContent(2),1.3*MCEff->GetBinContent(MCEff->GetMaximumBin()),"Data",8,skipleg);
	}


	c4_do->cd();
	gPad->SetLogx();
	gPad->SetGridx();
	if(Correction->IsEkin()){
	float ymin = Correction->GetGlobCorrection_noPID()->GetBinContent(Correction->GetGlobCorrection_noPID()->GetMinimumBin());
	float ymax = Correction->GetGlobCorrection_noPID()->GetBinContent(Correction->GetGlobCorrection_noPID()->GetMaximumBin());
	PlotTH1FintoGraph(gPad,Bins, (TH1F*) Correction->GetGlobCorrection(),"", "Eff. (Data/MC)",1,true,"Psame",rangemin,rangemax,0.65*ymin,1.2*ymax,"Data",8,true,false);
	}	
	else{
	float ymin = Correction->GetGlobCorrection_noPID()->GetBinContent(Correction->GetGlobCorrection_noPID()->GetMinimumBin());
	float ymax = Correction->GetGlobCorrection_noPID()->GetBinContent(Correction->GetGlobCorrection_noPID()->GetMaximumBin());
	PlotTH1FintoGraph(gPad,Bins, (TH1F*) Correction->GetGlobCorrection(),"", "Eff. (Data/MC)",1,false,"Psame",rangemin,rangemax,0.65*ymin,1.2*ymax,"Data",8,true,false); 
	}
	TSpline3 * model = (TSpline3 *) Correction->GetCorrectionModel();
	model->SetLineWidth(3);
	model->SetLineColor(2);
	model->Draw("same");


	Plots.Add(c3);
	Plots.Add(c4);
	Plots.writeObjsInFolder(("Eff. Correction/" +CorrName+"/"+RangeName).c_str());

}


void DrawCorrection(EffCorr * Correction, FileSaver Plots, std::string CorrName,std::string RangeName,float rangemin=-1, float rangemax=-1,bool skipleg=false) {
	DrawCorrection(Correction,Plots,Correction->GetBins(),CorrName,RangeName,0x0,0x0,0x0,0x0,0x0,0x0,rangemin,rangemax,skipleg);
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

	std::string before;
        std::string after;
        before = "";
        after  = ""; 

	EffCorr * TriggerEffCorr_HE  = new EffCorr(finalHistos,"TriggerEffCorr_HE" ,"Trigger Eff. Corr",true  ,before,after,"IsPrimary","IsProtonMC","IsPureDMC","IsPurePMC");
	EffCorr * TriggerFullSpan_HE  = new EffCorr(finalHistos,"TriggerFullSpan_HE" ,"Trigger Full Span",true  ,before,after,"IsPrimary","IsProtonMC","IsPureDMC","IsPurePMC");
	EffCorr * L1PickUpEffCorr_HE = new EffCorr(finalHistos,"L1PickUpEffCorr_HE","L1PickUp Eff. Corr",false,before,after,"IsPrimary",    "IsProtonMC","IsPureDMC","IsPurePMC");
	EffCorr * GoodQTrack_HE  = new EffCorr(finalHistos,"GoodQTrackEffCorr_HE","GoodQTrack Eff. Corr",true,before,after,"IsPrimary","IsProtonMC","IsPureDMC","IsPurePMC");
	EffCorr * GoodChi_HE  = new EffCorr(finalHistos,"GoodChiEffCorr_HE","GoodChi Eff. Corr",false,before,after,"IsPrimary","IsProtonMC","IsPureDMC","IsPurePMC");
	EffCorr * TrackerEffCorr_HE = new EffCorr(finalHistos,"TrackerEffCorr_HE","Tracker Eff. Corr",false,before,after,"IsPrimary","IsProtonMC","IsPureDMC","IsProtonPMC");
	EffCorr * StatusL1Check_HE = new EffCorr(finalHistos,"StatusL1Check_HE","StatusL1Check Eff. Corr",true,before,after,"IsPrimary","IsProtonMC","IsPureDMC","IsPurePMC");
	EffCorr * Good1Track_HE  = new EffCorr(finalHistos,"Good1TrackEffCorr_HE","Good1Track Eff. Corr",true,before,after,"IsPrimary","IsProtonMC","IsPureDMC","IsPurePMC");
	EffCorr * GoodLtof_HE  = new EffCorr(finalHistos,"GoodLTOFEffCorr_HE" ,"GoodLtof Eff. Corr",true,before,after,"IsPrimary","IsProtonMC","IsPureDMC","IsPurePMC");
	EffCorr * GoodUtof_HE  = new EffCorr(finalHistos,"GoodUtofEffCorr_HE" ,"GoodUtof Eff. Corr",true,before,after,"IsPrimary","IsProtonMC","IsPureDMC","IsPurePMC");
	EffCorr * GoodTime_TOF = new EffCorr(finalHistos,"GoodTimeEffCorr_TOF","GoodTime Eff. Corr",true,before,after,"IsPrimary","IsProtonMC","IsPureDMC","IsPurePMC");
	EffCorr * Quality_TOF = new EffCorr(finalHistos,"QualityEffCorr_TOF","Quality TOF Eff. Corr",true,before,after,"IsPrimary","IsProtonMC","IsPureDMC","IsPurePMC");
	EffCorr * RICHEffCorr_NaF = new EffCorr(finalHistos,"RICHCorrection_NaF","RICH Eff. Corr",true,before,(after+"&IsFromNaF").c_str(),"IsPrimary","IsProtonMC","IsPureDMC","IsPurePMC");
	EffCorr * RICHEffCorr_Agl = new EffCorr(finalHistos,"RICHCorrection_Agl","RICH Eff. Corr",true,before,(after+"&IsFromAgl").c_str(),"IsPrimary","IsProtonMC","IsPureDMC","IsPurePMC");
	EffCorr * RICHQualEffCorr_NaF = new EffCorr(finalHistos,"RICHQualCorrection_NaF","RICH Qual Eff. Corr",true,(before+"&IsFromNaF").c_str(),(after+"&IsFromNaF&RICHBDTCut").c_str(),"IsPrimary","IsProtonMC","IsPureDMC","IsPurePMC");
	EffCorr * RICHQualEffCorr_Agl = new EffCorr(finalHistos,"RICHqualCorrection_Agl","RICH Qual. Eff. Corr",true,(before+"&IsFromAgl").c_str(),(after+"&IsFromAgl&RICHBDTCut").c_str(),"IsPrimary","IsProtonMC","IsPureDMC","IsPurePMC");

	DrawCorrection(TriggerEffCorr_HE,Plots,"TriggerEffCorr","HE",0.1,40);
	DrawCorrection(TriggerFullSpan_HE,Plots,"TriggerFullSpan","HE",0.1,40);
	DrawCorrection(L1PickUpEffCorr_HE,Plots,"L1PickUpEffCorr","HE",0.1,40);
	DrawCorrection(TrackerEffCorr_HE,Plots,"TrackerEffCorr","HE",0.1,40);
	DrawCorrection(GoodChi_HE,Plots,"GoodChi_HE_EffCorr","HE",0.1,40);
	DrawCorrection(GoodQTrack_HE,Plots,"GoodQTrack_HE_EffCorr","HE",0.1,40);

	DrawCorrection(GoodUtof_HE,Plots,"GoodUtof_HE_EffCorr","HE",0.1,40);
	DrawCorrection(GoodLtof_HE,Plots,"GoodLtof_HE_EffCorr","HE",0.1,40);
	DrawCorrection(Good1Track_HE,Plots,"Good1Track_HE_EffCorr","HE",0.1,40);
	
	DrawCorrection(GoodTime_TOF,Plots, "GoodTime","TOF",0.1,40);
	DrawCorrection(Quality_TOF,Plots, "QualityTOF","TOF",0.1,40);
	DrawCorrection(RICHEffCorr_NaF,Plots,"RICHEffCorr","NaF",0.1,40);
	DrawCorrection(RICHEffCorr_Agl,Plots,"RICHEffCorr","Agl",0.1,40);
	DrawCorrection(RICHQualEffCorr_NaF,Plots,"RICHQualEffCorr","NaF",0.1,40);
	DrawCorrection(RICHQualEffCorr_Agl,Plots,"RICHQualEffCorr","Agl",0.1,40);
	

	return 0;
}

