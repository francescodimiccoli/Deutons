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

	if(rangemin==-1&&rangemax==-1) {
		rangemin =0.8*Bins.EkPerMasBins()[0];
		rangemax =1.2*Bins.EkPerMasBins()[Bins.EkPerMasBins().size()-1];
	}

	c3_up->cd();
	gPad->SetLogx();
	gPad->SetGridx();
	for(int lat=0;lat<10;lat++){
		TH1F * LatEff = (TH1F*)Correction->GetEfficiencyLat(lat);
		PlotTH1FintoGraph(gPad,Bins, LatEff,"R [GV]", "Efficiency",55+5*lat,false,"Psame",rangemin,rangemax,0.5*MCEff->GetBinContent(2),1.3*MCEff->GetBinContent(13),latitudes[lat].c_str(),8,skipleg,true);
	}

	c3_do->cd();
	gPad->SetLogx();
	gPad->SetGridx();
	for(int lat=0;lat<10;lat++){ 
	PlotTH1FintoGraph(gPad,Bins, (TH1F*)Correction->GetCorrectionLat(lat),"", "Eff. (Data/MC)",55+5*lat,false,"Psame",rangemin,rangemax,0.65*(Correction->GetGlobCorrection()->GetBinContent(10)),1.2*(Correction->GetGlobCorrection()->GetBinContent(10)),latitudes[lat].c_str(),8,true,true);
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
	c4_up->cd();
	gPad->SetLogx();
	gPad->SetGridx();
	PlotTH1FintoGraph(gPad,Bins, MCEff,"R [GV]", "Efficiency",2,false,"e4Psame",rangemin,rangemax,0.5*DataEff->GetBinContent(2),1.3*MCEff->GetBinContent(13),"MC",8,skipleg,true);
	PlotTH1FintoGraph(gPad,Bins, DataEff,"R [GV]", "Efficiency",1,false,"Psame",rangemin,rangemax,0.5*DataEff->GetBinContent(2),1.3*MCEff->GetBinContent(13),"Data",8,skipleg,true);

	c4_do->cd();
	gPad->SetLogx();
	gPad->SetGridx();
	PlotTH1FintoGraph(gPad,Bins, (TH1F*) Correction->GetGlobCorrection(),"", "Eff. (Data/MC)",1,false,"Psame",rangemin,rangemax,0.65*(Correction->GetGlobCorrection()->GetBinContent(10)),1.2*(Correction->GetGlobCorrection()->GetBinContent(10)),"Data",8,true,true);	
	
	Plots.Add(c3);
	Plots.Add(c4);
	Plots.writeObjsInFolder(("Eff. Correction/" +CorrName+"/"+RangeName).c_str());

}


void DrawCorrection(EffCorr * Correction, FileSaver Plots, Binning Bins, std::string CorrName,std::string RangeName,float rangemin=-1, float rangemax=-1,bool skipleg=false) {
	DrawCorrection(Correction,Plots,Bins,CorrName,RangeName,0x0,0x0,0x0,0x0,0x0,0x0,rangemin,rangemax,skipleg);
}


void DrawAllRangeCorrection(EffCorr * CorrectionTOF,EffCorr * CorrectionNaF,EffCorr * CorrectionAgl,FileSaver Plots,std::string CorrName){

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


	DrawCorrection(CorrectionTOF,Plots,ToFPB,CorrName,"TOF",c3,c3_up,c3_do,c4,c4_up,c4_do,0.4,12);
	DrawCorrection(CorrectionNaF,Plots,NaFPB,CorrName,"NaF",c3,c3_up,c3_do,c4,c4_up,c4_do,0.4,12,true);
	DrawCorrection(CorrectionAgl,Plots,AglPB,CorrName,"Agl",c3,c3_up,c3_do,c4,c4_up,c4_do,0.4,12,true);
	
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

	cout<<"**************************** PLOTTING ***************************************"<<endl;

	EffCorr * HEPPresEffCorr = new EffCorr(finalHistos,"HEPPresEffCorr","HEPPresEffCorr",PRB,"IsPositive&IsMinimumBias&IsLooseCharge1","IsPositive&IsMinimumBias&IsLooseCharge1&IsGolden","","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * HEPQualEffCorr = new EffCorr(finalHistos,"HEPQualEffCorr","HEPQualEffCorr",PRB,"IsPositive&IsMinimumBias&IsLooseCharge1","IsPositive&IsMinimumBias&IsLooseCharge1&IsGolden","","IsPurePMC","IsPureDMC","IsDeutonMC");
	
	std::string before;
        std::string after;
        before = "";
        after  = ""; 
	EffCorr * TriggerEffCorr_HE = new EffCorr(finalHistos,"TriggerEffCorr_HE","Trigger Eff. Corr",PRB,before,after,"IsPrimary",    "IsProtonMC","IsPureDMC","IsDeutonMC"); 
	EffCorr * TriggerEffCorr_TOF = new EffCorr(finalHistos,"TriggerEffCorr_TOF","Trigger Eff. Corr",ToFPB,before,after,"IsPrimary","IsProtonMC","IsPureDMC","IsDeutonMC");
	EffCorr * TriggerEffCorr_NaF = new EffCorr(finalHistos,"TriggerEffCorr_NaF","Trigger Eff. Corr",NaFPB,before,after,"IsPrimary","IsProtonMC","IsPureDMC","IsDeutonMC");
	EffCorr * TriggerEffCorr_Agl = new EffCorr(finalHistos,"TriggerEffCorr_Agl","Trigger Eff. Corr",AglPB,before,after,"IsPrimary","IsProtonMC","IsPureDMC","IsDeutonMC");

	EffCorr * Trigger2EffCorr_HE = new EffCorr(finalHistos,"Trigger2EffCorr_HE","Trigger2 Eff. Corr",PRB,before,after,"IsPrimary",    "IsProtonMC","IsPureDMC","IsDeutonMC"); 
	

	EffCorr * L1PickUpEffCorr_HE = new EffCorr(finalHistos,"L1PickUpEffCorr_HE","L1PickUp Eff. Corr",PRB,before,after,"IsPrimary",    "IsPurePMC","IsPureDMC","IsDeutonMC"); 
	EffCorr * L1PickUpEffCorr_TOF = new EffCorr(finalHistos,"L1PickUpEffCorr_TOF","L1PickUp Eff. Corr",ToFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * L1PickUpEffCorr_NaF = new EffCorr(finalHistos,"L1PickUpEffCorr_NaF","L1PickUp Eff. Corr",NaFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * L1PickUpEffCorr_Agl = new EffCorr(finalHistos,"L1PickUpEffCorr_Agl","L1PickUp Eff. Corr",AglPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");

	EffCorr * MinTOFEffCorr_HE  = new EffCorr(finalHistos,"MinTOFEffCorr_HE" ,"Min TOF Eff. Corr",PRB,before,after,  "IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC"); 
	EffCorr * MinTOFEffCorr_TOF = new EffCorr(finalHistos,"MinTOFEffCorr_TOF","Min TOF Eff. Corr",ToFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * MinTOFEffCorr_NaF = new EffCorr(finalHistos,"MinTOFEffCorr_NaF","Min TOF Eff. Corr",NaFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * MinTOFEffCorr_Agl = new EffCorr(finalHistos,"MinTOFEffCorr_Agl","Min TOF Eff. Corr",AglPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");

	EffCorr * TrackerEffCorr_HE = new EffCorr(finalHistos,"TrackerEffCorr_HE","Tracker Eff. Corr",ToFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	
	EffCorr * GoodChi_TOF = new EffCorr(finalHistos,"GoodChiEffCorr_TOF","GoodChi Eff. Corr",ToFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodChi_NaF = new EffCorr(finalHistos,"GoodChiEffCorr_NaF","GoodChi Eff. Corr",NaFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodChi_Agl = new EffCorr(finalHistos,"GoodChiEffCorr_Agl","GoodChi Eff. Corr",AglPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodUtof_TOF = new EffCorr(finalHistos,"GoodUtofEffCorr_TOF","GoodUtof Eff. Corr",ToFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodUtof_NaF = new EffCorr(finalHistos,"GoodUtofEffCorr_NaF","GoodUtof Eff. Corr",NaFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodUtof_Agl = new EffCorr(finalHistos,"GoodUtofEffCorr_Agl","GoodUtof Eff. Corr",AglPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
        EffCorr * GoodLtof_TOF = new EffCorr(finalHistos,"GoodLTOFEffCorr_TOF","GoodQTrack Eff. Corr",ToFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
        EffCorr * GoodLtof_NaF = new EffCorr(finalHistos,"GoodLTOFEffCorr_NaF","GoodQTrack Eff. Corr",NaFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
        EffCorr * GoodLtof_Agl = new EffCorr(finalHistos,"GoodLTOFEffCorr_Agl","GoodQTrack Eff. Corr",AglPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * Good1Track_TOF = new EffCorr(finalHistos,"Good1TrackEffCorr_TOF","Good1Track Eff. Corr",ToFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * Good1Track_NaF = new EffCorr(finalHistos,"Good1TackEffCorr_NaF", "Good1Track Eff. Corr",NaFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * Good1Track_Agl = new EffCorr(finalHistos,"Good1TrackEffCorr_Agl","Good1Track Eff. Corr",AglPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodQTrack_TOF = new EffCorr(finalHistos,"GoodQTrackEffCorr_TOF","GoodQTrack Eff. Corr",ToFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodQTrack_NaF = new EffCorr(finalHistos,"GoodQTrackEffCorr_NaF","GoodQTrack Eff. Corr",NaFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodQTrack_Agl = new EffCorr(finalHistos,"GoodQTrackEffCorr_Agl","GoodQTrack Eff. Corr",AglPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodTime_TOF = new EffCorr(finalHistos,"GoodTimeEffCorr_TOF","GoodTime Eff. Corr",ToFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * RICHEffCorr_NaF = new EffCorr(finalHistos,"RICHCorrection_NaF","RICH Eff. Corr",NaFPB,before,(after+"&IsFromNaF").c_str(),"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * RICHEffCorr_Agl = new EffCorr(finalHistos,"RICHCorrection_Agl","RICH Eff. Corr",AglPB,before,(after+"&IsFromAgl").c_str(),"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * RICHQualEffCorr_NaF = new EffCorr(finalHistos,"RICHQualCorrection_NaF","RICH Qual Eff. Corr",NaFPB,(before+"&IsFromNaF").c_str(),(after+"&IsFromNaF&RICHBDTCut").c_str(),"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * RICHQualEffCorr_Agl = new EffCorr(finalHistos,"RICHqualCorrection_Agl","RICH Qual. Eff. Corr",AglPB,(before+"&IsFromNaF").c_str(),(after+"&IsFromAgl&RICHBDTCut").c_str(),"","IsPurePMC","IsPureDMC","IsDeutonMC");

/*DrawCorrection(TrackerEffCorr_TOF,Plots, ToFPB,"TrackerEffCorr" ,"TOF");
	DrawCorrection(TrackerEffCorr_NaF,Plots, NaFPB,"TrackerEffCorr" ,"NaF");
	DrawCorrection(TrackerEffCorr_Agl,Plots, AglPB,"TrackerEffCorr" ,"Agl");
	DrawCorrection(GoodChi_TOF,Plots, ToFPB, "GoodChi","TOF");
	DrawCorrection(GoodChi_NaF,Plots, NaFPB, "GoodChi","NaF");
	DrawCorrection(GoodChi_Agl,Plots, AglPB, "GoodChi","Agl");
	DrawCorrection(GoodUtof_TOF,Plots, ToFPB, "GoodUtof","TOF");
	DrawCorrection(GoodUtof_NaF,Plots, NaFPB, "GoodUtof","NaF");
	DrawCorrection(GoodUtof_Agl,Plots, AglPB, "GoodUtof","Agl");
	DrawCorrection(GoodLtof_TOF,Plots, ToFPB, "GoodLtof","TOF");
	DrawCorrection(GoodLtof_NaF,Plots, NaFPB, "GoodLtof","NaF");
	DrawCorrection(GoodLtof_Agl,Plots, AglPB, "GoodLtof","Agl");
	DrawCorrection(Good1Track_TOF,Plots, ToFPB, "Good1Track","TOF");
	DrawCorrection(Good1Track_NaF,Plots, NaFPB, "Good1Track","NaF");
	DrawCorrection(Good1Track_Agl,Plots, AglPB, "Good1Track","Agl");
	DrawCorrection(GoodQTrack_TOF,Plots, ToFPB, "GoodQTrack","TOF");
	DrawCorrection(GoodQTrack_NaF,Plots, NaFPB, "GoodQTrack","NaF");
	DrawCorrection(GoodQTrack_Agl,Plots, AglPB, "GoodQTrack","Agl");
	*/
	DrawAllRangeCorrection(GoodChi_TOF,GoodChi_NaF,GoodChi_Agl,Plots,"GoodChiEffCorr");
	DrawAllRangeCorrection(GoodUtof_TOF,GoodUtof_NaF,GoodUtof_Agl,Plots,"GoodUtofEffCorr");
	DrawAllRangeCorrection(GoodLtof_TOF,GoodLtof_NaF,GoodLtof_Agl,Plots,"GoodLtofEffCorr");
	DrawAllRangeCorrection(Good1Track_TOF,Good1Track_NaF,Good1Track_Agl,Plots,"Good1TrackEffCorr");
	DrawAllRangeCorrection(GoodQTrack_TOF,GoodQTrack_NaF,GoodQTrack_Agl,Plots,"GoodQTrackEffCorr");

	DrawCorrection(TriggerEffCorr_HE,Plots,PRB,"TriggerEffCorr","HE",0,110);
	DrawCorrection(L1PickUpEffCorr_HE,Plots,PRB,"L1PickUpEffCorr","HE",0,110);
	DrawCorrection(TrackerEffCorr_HE,Plots,PRB,"TrackerEffCorr","HE",0,110);
	DrawCorrection(MinTOFEffCorr_HE,Plots,PRB,"MinTOFEffCorr","HE",0,110);
	
	DrawCorrection(GoodTime_TOF,Plots, ToFPB, "GoodTime","TOF",0.5,3);
	DrawCorrection(RICHEffCorr_NaF,Plots, NaFPB, "RICHEffCorr","NaF",1,7);
	DrawCorrection(RICHEffCorr_Agl,Plots, AglPB, "RICHEffCorr","Agl",3,10);
	DrawCorrection(RICHQualEffCorr_NaF,Plots, NaFPB, "RICHQualEffCorr","NaF",1,7);
	DrawCorrection(RICHQualEffCorr_Agl,Plots, AglPB, "RICHQualEffCorr","Agl",3,10);
	

	return 0;
}

