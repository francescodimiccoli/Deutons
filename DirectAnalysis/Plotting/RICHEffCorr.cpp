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

#include "../include/EffCorr.h"


void DrawCorrection(EffCorr * Correction, FileSaver Plots, Binning Bins, std::string CorrName,std::string RangeName, float rangemin=-1, float rangemax=-1,bool skipleg=false) {



	TCanvas * c4 = new TCanvas((CorrName + " Global").c_str());
	c4->cd();
	TVirtualPad * c4_up = new TPad("upperPad", "upperPad",0.0,0.3,1.0,1.0);
	c4_up->Draw();
	TVirtualPad * c4_do = new TPad("lowerPad", "lowerPad",0.0,0.0,1.0,0.3);
	c4_do->Draw();


	TH1F * MCEff = (TH1F*)Correction->GetMCEfficiency();
	TH1F * DataEff = (TH1F*)Correction->GetGlobEfficiency();

	cout<<MCEff<<" "<<DataEff<<endl;

	TH1F*	MCEff_plot   = ConvertBinnedHisto(MCEff,   "Protons MC",Correction->GetBins(),Correction->IsEkin());
	TH1F*	DataEff_plot = ConvertBinnedHisto(DataEff, "ISS data"  ,Correction->GetBins(),Correction->IsEkin());

	cout<<MCEff_plot<<" "<<DataEff_plot<<endl;

	c4_up->cd();
	gPad->SetLogx();
	gPad->SetGridx();
	MCEff_plot->SetLineColor(2);
	DataEff_plot->SetLineColor(1);
	MCEff_plot->SetLineWidth(2);
	DataEff_plot->SetLineWidth(2);
	
	MCEff_plot->GetYaxis()->SetRangeUser(0.2,1.1);
	MCEff_plot->Draw();
	DataEff_plot->Draw("same");


	c4_do->cd();
	gPad->SetLogx();
	gPad->SetGridx();
	TH1F * Correct = (TH1F*) Correction->GetGlobCorrection();
	TH1F * Correct_plot = ConvertBinnedHisto(Correct,   "Eff. ratio",Correction->GetBins(),Correction->IsEkin());
	
	TSpline3 * model = (TSpline3 *) Correction->GetCorrectionModel();
	model->SetLineWidth(3);
	model->SetLineColor(2);
	Correct_plot->SetLineWidth(2);
	Correct_plot->SetLineColor(1);
	Correct_plot->GetYaxis()->SetRangeUser(0.95,1.05);
	Correct_plot->GetYaxis()->SetLabelSize(0.045);
	Correct_plot->Draw();
	model->Draw("same");

	Plots.Add(c4);
	Plots.writeObjsInFolder(("Eff. Correction/" +CorrName+"/"+RangeName).c_str());

}


void DrawCorrection(EffCorr * Correction, FileSaver Plots, std::string CorrName,std::string RangeName,float rangemin=-1, float rangemax=-1,bool skipleg=false) {
	DrawCorrection(Correction,Plots,Correction->GetBins(),CorrName,RangeName,rangemin,rangemax,skipleg);
}

void DrawBlockMCcomparison(std::vector<EffCorr *> block, std::string Blockname,FileSaver Plots){

	TH1F * Efficiency_P   = block[0]->GetMCEfficiency();
	TH1F * Efficiency_D   = block[0]->GetMCEfficiency2();
	TH1F * Efficiency_Ppid= block[0]->GetMCEfficiency_noPID();
	TH1F * Efficiency_Dpid= block[0]->GetMCEfficiency_noPID2();

	for(int i=1;i<block.size();i++){
		Efficiency_P->Multiply(block[i]->GetMCEfficiency());	
		Efficiency_D->Multiply(block[i]->GetMCEfficiency());	
		Efficiency_Ppid->Multiply(block[i]->GetMCEfficiency_noPID());	
		Efficiency_Dpid->Multiply(block[i]->GetMCEfficiency_noPID2());	
	}	

	TH1F * Efficiency_P_plot   = ConvertBinnedHisto(Efficiency_P,   "Protons MC",block[0]->GetBins(),block[0]->IsEkin());
	TH1F * Efficiency_D_plot   = ConvertBinnedHisto(Efficiency_D,   "Deutons MC",block[0]->GetBins(),block[0]->IsEkin());
	TH1F * Efficiency_Ppid_plot= ConvertBinnedHisto(Efficiency_Ppid,"Protons MC (PID)",block[0]->GetBins(),block[0]->IsEkin());
	TH1F * Efficiency_Dpid_plot= ConvertBinnedHisto(Efficiency_Dpid,"Deutons MC (PID)",block[0]->GetBins(),block[0]->IsEkin());

	Efficiency_P_plot   ->SetLineColor(1);
	Efficiency_D_plot   ->SetLineColor(2);
        Efficiency_Ppid_plot->SetLineColor(1);
	Efficiency_Dpid_plot->SetLineColor(2);

	Efficiency_P_plot   ->SetLineWidth(2);
	Efficiency_D_plot   ->SetLineWidth(2);
        Efficiency_Ppid_plot->SetLineWidth(2);
	Efficiency_Dpid_plot->SetLineWidth(2);

	Efficiency_P_plot   ->SetLineStyle(1);
	Efficiency_D_plot   ->SetLineStyle(1);
        Efficiency_Ppid_plot->SetLineStyle(2);
	Efficiency_Dpid_plot->SetLineStyle(2);

	TCanvas * c1 = new TCanvas("MC comparison");
	Efficiency_P_plot   ->Draw("Psame");
       	Efficiency_D_plot   ->Draw("Psame");
        Efficiency_Ppid_plot->Draw("Psame");
       	Efficiency_Dpid_plot->Draw("Psame");

	Plots.Add(c1);
	Plots.writeObjsInFolder(("Blocks/MC comparison/"+Blockname).c_str());	

}	

void DrawBlockCascade(std::vector<EffCorr *> block, std::string Blockname,FileSaver Plots,float rangemin, float rangemax){

	TH1F * Efficiencies[block.size()];
	for(int i=0;i<block.size();i++){
		Efficiencies[i] = block[i]->GetGlobEfficiency();
		if(i!=0) Efficiencies[i]->Multiply(Efficiencies[i-1]);
	}	
	
	TCanvas * c1 = new TCanvas("Cascade");
	TH1F * Efficiencies_plot[block.size()];
	for(int i=0;i<block.size();i++){
		Efficiencies_plot[i] = CreateHisto(Efficiencies[i]->GetTitle(),block[0]->GetBins(),block[0]->IsEkin()); 
		for(int bin=0;bin<Efficiencies[i]->GetNbinsX();bin++){
			Efficiencies_plot[i]->SetBinContent(bin+1,Efficiencies[i]->GetBinContent(bin+1));
			Efficiencies_plot[i]->SetBinError(bin+1,Efficiencies[i]->GetBinError(bin+1));
		}
		Efficiencies_plot[i]->SetLineColor(i+2);
		Efficiencies_plot[i]->SetLineWidth(2);
		Efficiencies_plot[i]->Draw("same");
	}


	Plots.Add(c1);
	Plots.writeObjsInFolder(("Blocks/Cascade/"+Blockname).c_str());	
}


void DrawCorrectionBlock(std::vector<EffCorr *> block, std::string Blockname,FileSaver Plots){
	TStyle* m_gStyle= new TStyle();;
        m_gStyle->SetPalette(55);
	int nColors = m_gStyle->GetNumberOfColors();


	TGraph * TotalCorrection = new TGraph();
	float y=1;
	for(int bin=0;bin<block[0]->GetBins().size();bin++){
		y=1;
		for(int i=0;i<block.size();i++){
			if(block[i]->IsEkin()) y*=block[i]->GetCorrectionModel()->Eval(block[i]->GetBins().EkPerMassBinCent(bin));
			else y*=block[i]->GetCorrectionModel()->Eval(block[i]->GetBins().RigBinCent(bin));
		}
		if(block[0]->IsEkin()) TotalCorrection->SetPoint(bin,block[0]->GetBins().EkPerMassBinCent(bin),y);
		else TotalCorrection->SetPoint(bin,block[0]->GetBins().RigBinCent(bin),y); 	
	}
	TotalCorrection->SetLineColor(1);
	TotalCorrection->SetLineWidth(2);	
	TotalCorrection->SetLineStyle(2);	
	TotalCorrection->SetTitle((Blockname+ " Corrections").c_str());
	TotalCorrection->GetYaxis()->SetRangeUser(0.8*y,1.1*block[0]->GetCorrectionModel()->Eval(5));
	TotalCorrection->GetYaxis()->SetTitle("Correction");
	if(block[0]->IsEkin())	TotalCorrection->GetXaxis()->SetTitle("Ekin [GeV/n]");
	else TotalCorrection->GetXaxis()->SetTitle("R [GV]");
	
	TCanvas * c1 = new TCanvas("Correction Block");
	TotalCorrection->Draw("AC");
	for(int i=0;i<block.size();i++) block[i]->GetCorrectionModel()->SetLineColor(i+2);	
	for(int i=0;i<block.size();i++) block[i]->GetCorrectionModel()->Draw("same");

	Plots.Add(c1);
	Plots.writeObjsInFolder(("Blocks/"+Blockname).c_str());	


//	TH1F * Total_StatError = (TH1F *) block[0]->GetStat_Err()->Clone("Total_StatError");
//	TH1F * Total_SystError = (TH1F *) block[0]->GetSyst_Err()->Clone("Total_SystError");

	
	TH1F * Total_StatError = CreateHisto("Total_StatError",block[0]->GetBins(),block[0]->IsEkin()); 
	TH1F * Total_SystError = CreateHisto("Total_SystError",block[0]->GetBins(),block[0]->IsEkin());
	TH1F * Total_SystStat  = CreateHisto("Total_SystStat",block[0]->GetBins(),block[0]->IsEkin());



	for(int bin=0;bin<block[0]->GetBins().size();bin++){
		float errstat=0;
		float errsyst=0;
		float errsyststat=0;
		
		for(int i=0;i<block.size();i++){
			errstat+=pow(block[i]->GetStat_Err()->GetBinContent(bin+1),2);			
			errsyst+=pow(block[i]->GetSyst_Err()->GetBinContent(bin+1),2);			
			errsyststat+=pow(block[i]->GetSyst_Stat()->GetBinContent(bin+1),2);			
		}	
		cout<<"systatat: "<<errsyststat<<endl;
		Total_StatError->SetBinContent(bin+1,pow(errstat,0.5));
		Total_SystError->SetBinContent(bin+1,pow(errsyst,0.5));
		Total_SystStat->SetBinContent(bin+1,pow(errsyststat,0.5));
	}

	TCanvas * c2 = new TCanvas("Block Error");
	
	Total_SystError->SetLineColor(2);
        Total_StatError->SetLineColor(1);
	Total_SystError->SetLineWidth(2);
        Total_StatError->SetLineWidth(1);
	Total_SystStat->SetLineColor(4);
        Total_SystStat->SetLineWidth(2);
	Total_SystError->Draw("hist");
	Total_StatError->Draw("hist,same");	
	Total_SystStat->Draw("hist,same");	


	Plots.Add(c2);
	Plots.writeObjsInFolder(("Blocks/"+Blockname).c_str());	


}


void DrawTotalUncertainty(std::vector<std::vector<EffCorr *>> Total,FileSaver Plots, std::string name){
	TH1F * Total_Uncertainty[Total.size()];

	for(int i =0;i<Total.size();i++) 
		Total_Uncertainty[i]=CreateHisto(("Total_Uncertainty"+to_string(i)).c_str(),Total[0][0]->GetBins(),Total[0][0]->IsEkin());

			
	for(int bin=0;bin<Total[0][0]->GetBins().size();bin++){
		for(int i =0;i<Total.size();i++){
				float errcorrection=0;
				for(int j=0;j<Total[i].size();j++){
					errcorrection+=pow(Total[i][j]->GetStat_Err()->GetBinContent(bin+1),2);
					errcorrection+=pow(Total[i][j]->GetSyst_Err()->GetBinContent(bin+1),2);
					errcorrection+=pow(Total[i][j]->GetSyst_Stat()->GetBinContent(bin+1),2);
				}
				Total_Uncertainty[i]->SetBinContent(bin+1,pow(errcorrection,0.5));
			}			
	  }

	TH1F * TOTAL = CreateHisto("Total Uncertainty",Total[0][0]->GetBins(),Total[0][0]->IsEkin());
	for(int bin=0;bin<Total_Uncertainty[0]->GetNbinsX();bin++){
		float content=0;
		for(int i =0;i<Total.size();i++){
			content += pow(Total_Uncertainty[i]->GetBinContent(bin+1),2);
			}
		TOTAL->SetBinContent(bin+1,pow(content,0.5) );
	}


	TCanvas * c2 = new TCanvas(("Global Error "+name).c_str());
	TOTAL->SetLineColor(1);
	TOTAL->SetLineWidth(3);
	TOTAL->Draw("hist");
	
	for(int i =0;i<Total.size();i++){

		Total_Uncertainty[i]->SetLineColor(i+2);
		Total_Uncertainty[i]->SetLineWidth(2);
		Total_Uncertainty[i]->Draw("hist,same");
	}		
	Plots.Add(c2);
        Plots.writeObjsInFolder("Global Errors");
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

	EffCorr * TriggerEffCorr_HE  = new EffCorr(finalHistos,"TriggerEffCorr_HE" ,"Trigger Eff. Corr",true  ,before,after,"IsPrimary"); 
	EffCorr * TriggerFullSpan_HE  = new EffCorr(finalHistos,"TriggerFullSpan_HE" ,"Trigger Full Span",true  ,before,after,"IsPrimary"); 
	EffCorr * L1PickUpGeom_HE = new EffCorr(finalHistos,"L1PickUpGeom_HE","L1PickUp Geom.",false,before,after,"IsPrimaryInner",true); 
	EffCorr * L1PickUpEffCorr_HE = new EffCorr(finalHistos,"L1PickUpEffCorr_HE","L1PickUp Eff. Corr",false,before,after,"IsPrimaryInner",true); 
	EffCorr * TrackerEffCorr_HE = new EffCorr(finalHistos,"TrackerEffCorr_HE","Tracker Eff. Corr",false,before,after,"",true);
	EffCorr * GoodQTrack_HE  = new EffCorr(finalHistos,"GoodQTrackEffCorr_HE","GoodQTrack Eff. Corr",false,before,after,"",true);
	EffCorr * GoodChi_HE  = new EffCorr(finalHistos,"GoodChiEffCorr_HE","GoodChi Eff. Corr",false,before,after,"",true);
	EffCorr * KalmanEffCorr_HE = new EffCorr(finalHistos,"KalmanEffCorr_HE","Kalman Eff. Corr",false,before,after,"",true);
	EffCorr * Good1Track_HE  = new EffCorr(finalHistos,"Good1TrackEffCorr_HE","Good1Track Eff. Corr",true,before,after,"");
        EffCorr * GoodLtof_HE  = new EffCorr(finalHistos,"GoodLTOFEffCorr_HE" ,"GoodLtof Eff. Corr",true,before,after,"");
	EffCorr * GoodUtof_HE  = new EffCorr(finalHistos,"GoodUtofEffCorr_HE" ,"GoodUtof Eff. Corr",true,before,after,"");
	EffCorr * GoodTime_TOF = new EffCorr(finalHistos,"GoodTimeEffCorr_TOF","GoodTime Eff. Corr",true,before,after,"");
	EffCorr * Quality_TOF = new EffCorr(finalHistos,"QualityEffCorr_TOF","Quality TOF Eff. Corr",true,before,after,"");
	EffCorr * RICHEffCorr_NaF = new EffCorr(finalHistos,"RICHCorrection_NaF","RICH Eff. Corr",true,before,(after+"&IsFromNaF").c_str(),"IsPrimary");
	EffCorr * RICHEffCorr_Agl = new EffCorr(finalHistos,"RICHCorrection_Agl","RICH Eff. Corr",true,before,(after+"&IsFromAgl").c_str(),"IsPrimary");
	EffCorr * RICHQualEffCorr_NaF = new EffCorr(finalHistos,"RICHQualCorrection_NaF","RICH Qual Eff. Corr",true,(before+"&IsFromNaF").c_str(),(after+"&IsFromNaF&RICHBDTCut").c_str(),"IsPrimary");
	EffCorr * RICHQualEffCorr_Agl = new EffCorr(finalHistos,"RICHqualCorrection_Agl","RICH Qual. Eff. Corr",true,(before+"&IsFromAgl").c_str(),(after+"&IsFromAgl&RICHBDTCut").c_str(),"IsPrimary");



	//single corrections
	DrawCorrection(TriggerEffCorr_HE,Plots,"TriggerEffCorr","HE",0.1,40);
	DrawCorrection(TriggerFullSpan_HE,Plots,"TriggerFullSpan","HE",0.1,40);
	DrawCorrection(L1PickUpGeom_HE,Plots,"L1PickUpGeom","HE",0.1,40);
	DrawCorrection(L1PickUpEffCorr_HE,Plots,"L1PickUpEffCorr","HE",0.1,40);
	DrawCorrection(TrackerEffCorr_HE,Plots,"TrackerEffCorr","HE",0.1,40);
	DrawCorrection(GoodChi_HE,Plots,"GoodChi_HE_EffCorr","HE",0.1,40);
	DrawCorrection(GoodQTrack_HE,Plots,"GoodQTrack_HE_EffCorr","HE",0.1,40);
	DrawCorrection(KalmanEffCorr_HE,Plots,"KalmanEffCorr","HE",0.1,40);
	
	DrawCorrection(GoodUtof_HE,Plots,"GoodUtof_HE_EffCorr","HE",0.1,40);
	DrawCorrection(GoodLtof_HE,Plots,"GoodLtof_HE_EffCorr","HE",0.1,40);
	DrawCorrection(Good1Track_HE,Plots,"Good1Track_HE_EffCorr","HE",0.1,40);
	
	DrawCorrection(GoodTime_TOF,Plots, "GoodTime","TOF",0.1,40);
	DrawCorrection(Quality_TOF,Plots, "QualityTOF","TOF",0.1,40);
	DrawCorrection(RICHEffCorr_NaF,Plots,"RICHEffCorr","NaF",0.1,40);
	DrawCorrection(RICHEffCorr_Agl,Plots,"RICHEffCorr","Agl",0.1,40);
	DrawCorrection(RICHQualEffCorr_NaF,Plots,"RICHQualEffCorr","NaF",0.1,40);
	DrawCorrection(RICHQualEffCorr_Agl,Plots,"RICHQualEffCorr","Agl",0.1,40);
	

	//correction blocks
	std::vector<EffCorr *> Trigg;
        Trigg.push_back(TriggerEffCorr_HE);
        DrawBlockCascade(Trigg,"Trigger",Plots,0.1,40);
        DrawCorrectionBlock(Trigg,"Trigger",Plots);
        DrawBlockMCcomparison(Trigg,"Trigger",Plots);  

	std::vector<EffCorr *> TriggFS;
        TriggFS.push_back(TriggerFullSpan_HE);
        DrawBlockCascade(TriggFS,"TriggerFS",Plots,0.1,40);
        DrawCorrectionBlock(TriggFS,"TriggerFS",Plots);
        DrawBlockMCcomparison(TriggFS,"TriggerFS",Plots);  

	std::vector<EffCorr *> Track;
	Track.push_back(TrackerEffCorr_HE);
	Track.push_back(GoodQTrack_HE);
	Track.push_back(GoodChi_HE);
	Track.push_back(KalmanEffCorr_HE);
	DrawBlockCascade(Track	   ,"Tracking",Plots,0.1,40);
        DrawCorrectionBlock(Track  ,"Tracking",Plots);
        DrawBlockMCcomparison(Track,"Tracking",Plots);  

        std::vector<EffCorr *> L1PickUpGeom;
        L1PickUpGeom.push_back(L1PickUpGeom_HE);
        DrawBlockCascade(L1PickUpGeom,"L1PickUpGeom",Plots,0.1,40);       
        DrawCorrectionBlock(L1PickUpGeom,"L1PickUpGeom",Plots); 

        std::vector<EffCorr *> L1PickUpAss;
        L1PickUpAss.push_back(L1PickUpEffCorr_HE);
        DrawBlockCascade(L1PickUpAss,"L1PickUpAss",Plots,0.1,40);       
        DrawCorrectionBlock(L1PickUpAss,"L1PickUpAss",Plots); 

	std::vector<EffCorr *> NoInteractions;
	NoInteractions.push_back(Good1Track_HE);
	NoInteractions.push_back(GoodUtof_HE);
	NoInteractions.push_back(GoodLtof_HE);
	DrawBlockCascade(NoInteractions,"NoInteractions",Plots,0.1,40);	
	DrawCorrectionBlock(NoInteractions,"NoInteractions",Plots);	
	DrawBlockMCcomparison(NoInteractions,"NoInteractions",Plots);	

	std::vector<EffCorr *> GoodBetaTOF;
	GoodBetaTOF.push_back(GoodTime_TOF);
	GoodBetaTOF.push_back(Quality_TOF);
	DrawBlockCascade(GoodBetaTOF,"GoodBetaTOF",Plots,0.1,40);	
	DrawCorrectionBlock(GoodBetaTOF,"GoodBetaTOF",Plots);	
	DrawBlockMCcomparison(GoodBetaTOF,"GoodBetaTOF",Plots);	


	std::vector<EffCorr *> GoodBetaNaF;
	GoodBetaNaF.push_back(RICHEffCorr_NaF);
	GoodBetaNaF.push_back(RICHQualEffCorr_NaF);
	DrawBlockCascade(GoodBetaNaF,"GoodBetaNaF",Plots,0.1,40);	
	DrawCorrectionBlock(GoodBetaNaF,"GoodBetaNaF",Plots);	
	DrawBlockMCcomparison(GoodBetaNaF,"GoodBetaNaF",Plots);	


	std::vector<EffCorr *> GoodBetaAgl;
	GoodBetaAgl.push_back(RICHEffCorr_Agl);
	GoodBetaAgl.push_back(RICHQualEffCorr_Agl);
	DrawBlockCascade(GoodBetaAgl,"GoodBetaAgl",Plots,0.1,40);	
	DrawCorrectionBlock(GoodBetaAgl,"GoodBetaAgl",Plots);	
	DrawBlockMCcomparison(GoodBetaAgl,"GoodBetaAgl",Plots);	


	std::vector<std::vector<EffCorr *>> TotalTOF;
	TotalTOF.push_back(Trigg);
	TotalTOF.push_back(L1PickUpGeom);
	TotalTOF.push_back(L1PickUpAss);
	TotalTOF.push_back(NoInteractions);
	TotalTOF.push_back(GoodBetaTOF);
	DrawTotalUncertainty(TotalTOF,Plots,"TOF");

	std::vector<std::vector<EffCorr *>> TotalNaF;
	TotalNaF.push_back(Trigg);
	TotalNaF.push_back(L1PickUpGeom);
	TotalNaF.push_back(L1PickUpAss);
	TotalNaF.push_back(Track);
	TotalNaF.push_back(NoInteractions);
	TotalNaF.push_back(GoodBetaNaF);
	DrawTotalUncertainty(TotalNaF,Plots,"NaF");

	std::vector<std::vector<EffCorr *>> TotalAgl;
	TotalAgl.push_back(Trigg);
	TotalAgl.push_back(L1PickUpGeom);
	TotalAgl.push_back(L1PickUpAss);
	TotalAgl.push_back(Track);
	TotalAgl.push_back(NoInteractions);
	TotalAgl.push_back(GoodBetaAgl);
	DrawTotalUncertainty(TotalAgl,Plots,"Agl");


	//Fragmentation
	TH1F * FragmP_TOF = Quality_TOF->GetMCAfter();
	FragmP_TOF -> Divide(Quality_TOF->GetMCAfter_noPID());

	Plots.Add(FragmP_TOF);

	Plots.writeObjsInFolder("Fragmentation");

	return 0;
}

