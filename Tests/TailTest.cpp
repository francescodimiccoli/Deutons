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

#include "../DirectAnalysis/include/filesaver.h"

#include "../DirectAnalysis/include/GlobalBinning.h"

#include "../DirectAnalysis/include/Commonglobals.cpp"
#include "../DirectAnalysis/include/Variables.hpp"

int frac =20;
using namespace std;

float SmearBeta(float Beta, float sigma, float shift){

        float time = 1.2/(Beta*3e-4);
        time = time + shift + Rand->Gaus(0,sigma);
        return 1.2/(time*3e-4);

}


int TailTest(){
	cout<<"************************ READING DATA ***************************"<<endl;

	string inputfileMC = "~/fdimicco/MAIN/sommaMC/temp/sommaMC123.root";
	TFile * inputMC = TFile::Open(inputfileMC.c_str());
	string inputfileDT = "~/fdimicco/MAIN/sommadati/sommadati123.root";
	TFile * inputDT = TFile::Open(inputfileDT.c_str());


	TTree * treeMC = (TTree *)inputMC->Get("parametri_geo");
	TTree * treeDT = (TTree *)inputDT->Get("parametri_geo");
	
	cout<<inputMC<<" "<<inputDT<<endl;

	FileSaver finalHistos;
        finalHistos.setName("TailTest_Plots.root");

	Variables * vars= new Variables();

	cout<<"************************ CUTS & VARIABLES DEFINITIONS ***************************"<<endl;


	TH1F * MassMC = new TH1F("MassMC","MassMC",100,0,5);
	TH1F * MassDT = new TH1F("MassDT","MassDT",100,0,5);
	TH1F * MassMCSmear1 = new TH1F("MassMCSmear1","MassMCSmear1",100,0,5);
	TH1F * MassMCSmear2 = new TH1F("MassMCSmear2","MassMCSmear2",100,0,5);
	TH1F * MassMCSmear3 = new TH1F("MassMCSmear3","MassMCSmear3",100,0,5);
			

	vars->ReadBranches(treeDT);

        for(int i=0;i<treeDT->GetEntries()/frac;i++){
                UpdateProgressBar(i, treeDT->GetEntries()/frac);
                treeDT->GetEvent(i);
		vars->Update();
		if(vars->qL1>0&&vars->qL1<1.7&&vars->qUtof>0.8&&vars->qUtof<1.3&&vars->qLtof>0.8&&vars->qLtof<1.3&&vars->qInner>0.8&&vars->qInner<1.3&&vars->Beta>0.83&&vars->Beta<0.866)
			MassDT->Fill(vars->R/vars->Beta*pow(1-pow(vars->Beta,2),0.5));
		
	}

	vars->ReadBranches(treeMC);

        for(int i=0;i<treeMC->GetEntries()/frac;i++){
                UpdateProgressBar(i, treeMC->GetEntries()/frac);
                treeMC->GetEvent(i);
		vars->Update();
		if(vars->qL1>0&&vars->qUtof>0.8&&vars->qUtof<1.3&&vars->qLtof>0.8&&vars->qLtof<1.3&&vars->qInner>0.8&&vars->qInner<1.3&&vars->Beta>0.83&&vars->Beta<0.866&&vars->Massa_gen<1)
			MassMC->Fill(vars->R/vars->Beta*pow(1-pow(vars->Beta,2),0.5),vars->mcweight);
		

		float betasmear=SmearBeta(vars->Beta,90,0);
		

		if(vars->R>2.7) betasmear = SmearBeta(vars->Beta,90,0);
		if(vars->qL1>0&&vars->qL1<1.7&&vars->qUtof>0.8&&vars->qUtof<1.3&&vars->qLtof>0.8&&vars->qLtof<1.3&&vars->qInner>0.8&&vars->qInner<1.3&&betasmear>0.83&&betasmear<0.866&&vars->Massa_gen<1)
                        MassMCSmear1->Fill(vars->R/betasmear*pow(1-pow(betasmear,2),0.5),vars->mcweight);

		if(vars->R>2.7) betasmear = SmearBeta(vars->Beta,110,0);
		if(vars->qL1>0&&vars->qL1<1.7&&vars->qUtof>0.8&&vars->qUtof<1.3&&vars->qLtof>0.8&&vars->qLtof<1.3&&vars->qInner>0.8&&vars->qInner<1.3&&betasmear>0.83&&betasmear<0.866&&vars->Massa_gen<1)
                        MassMCSmear2->Fill(vars->R/betasmear*pow(1-pow(betasmear,2),0.5),vars->mcweight);

		if(vars->R>2.7) betasmear = SmearBeta(vars->Beta,150,0);
		if(vars->qL1>0&&vars->qL1<1.7&&vars->qUtof>0.8&&vars->qUtof<1.3&&vars->qLtof>0.8&&vars->qLtof<1.3&&vars->qInner>0.8&&vars->qInner<1.3&&betasmear>0.83&&betasmear<0.866&&vars->Massa_gen<1)
                        MassMCSmear3->Fill(vars->R/betasmear*pow(1-pow(betasmear,2),0.5),vars->mcweight);



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
