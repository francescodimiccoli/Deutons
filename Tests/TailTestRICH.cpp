#include <iostream>

#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TH3.h"
#include "TF2.h"
#include <vector>
#include <string>
#include <sstream>

#include "TVector3.h"
#include "TMath.h"
#include "TKey.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TNtuple.h"
#include "TObject.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TObjArray.h"

#include "../include/filesaver.h"

#include "../include/GlobalBinning.h"

#include "../Ntuple-making/Commonglobals.cpp"
#include "../include/Variables.hpp"
#include "../include/Cuts.h"

int frac =5;
using namespace std;

float SmearBeta(float Beta, float sigma, float shift){

        float time = 1.2/(Beta*3e-4);
        time = time + shift + Rand->Gaus(0,sigma);
        return 1.2/(time*3e-4);

}

float SimulateBadEvents(int EvNum,float Beta,int Level){
	if(EvNum%Level==0) { float beta= Rand->Uniform(0.96,1); return beta;}
	return Beta; 
	
};

int TailTestRICH(){
	cout<<"************************ READING DATA ***************************"<<endl;

	string inputfileMC = "../Ntuple-making/Ntuples/MC/NtupleMC.root";
	TFile * inputMC = TFile::Open(inputfileMC.c_str());
	string inputfileDT = "../Ntuple-making/Ntuples/1314835200/NtupleData.root";
	TFile * inputDT = TFile::Open(inputfileDT.c_str());


	TNtuple * treeMC = (TNtuple *)inputMC->Get("Q");
	TNtuple * treeDT = (TNtuple *)inputDT->Get("Q");
	
	cout<<inputMC<<" "<<inputDT<<endl;

	FileSaver finalHistos;
        finalHistos.setName("TailTest_Plots.root");

	Variables * vars= new Variables();

	cout<<"************************ CUTS & VARIABLES DEFINITIONS ***************************"<<endl;

	float rangemin=0;
	float rangemax=6;

	TH1F * MassMC = new TH1F("MassMC","MassMC",100,rangemin,rangemax);
	TH1F * MassDT = new TH1F("MassDT","MassDT",100,rangemin,rangemax);
	TH1F * MassMCSmear1 = new TH1F("MassMCSmear1","MassMCSmear1",100,rangemin,rangemax);
	TH1F * MassMCSmear2 = new TH1F("MassMCSmear2","MassMCSmear2",100,rangemin,rangemax);
	TH1F * MassMCSmear3 = new TH1F("MassMCSmear3","MassMCSmear3",100,rangemin,rangemax);
			

	vars->ReadAnalysisBranches(treeDT);

        for(int i=0;i<treeDT->GetEntries()/frac;i++){
                vars->AnalysisVariablseReset();
                UpdateProgressBar(i, treeDT->GetEntries()/frac);
                treeDT->GetEvent(i);
		if(ControlSample(vars)&&IsFromAgl(vars)&&vars->BetaRICH_new<0.985&&LikelihoodCut(vars)&&DistanceCut(vars))
			MassDT->Fill(vars->R/vars->BetaRICH_new*pow(1-pow(vars->BetaRICH_new,2),0.5));
		
	}

	vars->ReadAnalysisBranches(treeMC);
	int passed=0;
        for(int i=0;i<treeMC->GetEntries()/frac;i++){
                vars->AnalysisVariablseReset();
                UpdateProgressBar(i, treeMC->GetEntries()/frac);
                treeMC->GetEvent(i);
		if(IsFromAgl(vars)&&IsProtonMC(vars))passed++;
		if(ControlSample(vars)&&IsFromAgl(vars)&&vars->BetaRICH_new<0.985&&IsProtonMC(vars)&&LikelihoodCut(vars)&&DistanceCut(vars)){
			MassMC->Fill(vars->R/vars->BetaRICH_new*pow(1-pow(vars->BetaRICH_new,2),0.5),vars->mcweight);
		}
		float betasmear1,betasmear2,betasmear3;
		if(IsFromAgl(vars)&&IsProtonMC(vars))  { betasmear1=SimulateBadEvents(passed,vars->BetaRICH_new,200);
							 betasmear2=SimulateBadEvents(passed,vars->BetaRICH_new,350);
							 betasmear3=SimulateBadEvents(passed,vars->BetaRICH_new,150);
							}			
		if(ControlSample(vars)&&IsFromAgl(vars)&&betasmear1<0.985&&IsProtonMC(vars)&&LikelihoodCut(vars)&&DistanceCut(vars)){
			MassMCSmear1->Fill(vars->R/betasmear1*pow(1-pow(betasmear1,2),0.5),vars->mcweight);
		}	
		if(ControlSample(vars)&&IsFromAgl(vars)&&betasmear2<0.985&&IsProtonMC(vars)&&LikelihoodCut(vars)&&DistanceCut(vars)){
			MassMCSmear2->Fill(vars->R/betasmear2*pow(1-pow(betasmear2,2),0.5),vars->mcweight);
		}	
		if(ControlSample(vars)&&IsFromAgl(vars)&&betasmear3<0.985&&IsProtonMC(vars)&&LikelihoodCut(vars)&&DistanceCut(vars)){
			MassMCSmear3->Fill(vars->R/betasmear3*pow(1-pow(betasmear3,2),0.5),vars->mcweight);
		}	
			

	}

	MassDT->Scale(1/(float)MassDT->GetBinContent(MassDT->GetMaximumBin()));
	MassMC->Scale(1/(float)MassMC->GetBinContent(MassDT->GetMaximumBin()));	
	MassMCSmear1->Scale(1/(float)MassMCSmear1->GetBinContent(MassDT->GetMaximumBin()));	
	MassMCSmear2->Scale(1/(float)MassMCSmear2->GetBinContent(MassDT->GetMaximumBin()));	
	MassMCSmear3->Scale(1/(float)MassMCSmear3->GetBinContent(MassDT->GetMaximumBin()));	



	MassDT->SetLineColor(1);
	MassMC->SetLineColor(2);
	MassMCSmear1->SetLineColor(3);
	MassMCSmear2->SetLineColor(4);
	MassMCSmear3->SetLineColor(5);



	MassDT->SetLineWidth(6);
	MassMC->SetLineWidth(6);
	MassMCSmear1->SetLineWidth(6);
	MassMCSmear2->SetLineWidth(6);
	MassMCSmear3->SetLineWidth(6);




	MassDT->Draw();
	MassMC->Draw("same");
	MassMCSmear1->Draw("same");
	MassMCSmear2->Draw("same");
	MassMCSmear3->Draw("same");





	return 0;
}
