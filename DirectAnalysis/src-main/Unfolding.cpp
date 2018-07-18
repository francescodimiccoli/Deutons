#include <bitset>
#include "TROOT.h"
#include "TNtuple.h"
#include <TSpline.h>
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include <TVector3.h>
#include "TMath.h"
#include <TFile.h>
#include <string>
#include "TFile.h"
#include "TH2.h"
#include "TF2.h"
#include <TVector3.h>
#include "TChain.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TFractionFitter.h"
#include "TRandom3.h"

#include "InputFileReader.h"
#include "DBarReader.h"

#include "filesaver.h"
#include "binning.h"
#include "Globals.h"
#include "Variables.hpp"
#include "ParallelFiller.h"
#include "Cuts.h"

#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char * argv[])
{


	cout<<"****************************** FILES OPENING ***************************************"<<endl;

	TH1::SetDefaultSumw2();     	
	cout<<"****************************** FILES OPENING ***************************************"<<endl;

	string INPUT1(argv[1]);
	string INPUT2(argv[2]);
	string INPUT3(argv[3]);

	TChain * chainDT = InputFileReader(INPUT1.c_str(),"parametri_geo");
	TChain * chainMC = InputFileReader(INPUT2.c_str(),"parametri_MC");

       FileSaver finalResults;
       finalResults.setName((INPUT3+"_Results").c_str());



	cout<<"****************************** BINS ***************************************"<<endl;
	SetUpUsualBinning();

	
	Variables * vars = new Variables();

	cout<<"****************************** ANALYIS ******************************************"<<endl;
	
	TH1F * GenBeta = new TH1F ("GenBeta","GenBeta",100,0,1);
	TH1F * MeasBeta = new TH1F ("MeasBeta","MeasBeta",100,0,1);
	
	TH1F * GenP = new TH1F("GenP","GenP",ToFDB.size(),0,ToFDB.size());
	TH1F * GenD = new TH1F("GenD","GenD",ToFDB.size(),0,ToFDB.size());
	TH1F * MeasP = new TH1F("MeasP","MeasP",ToFDB.size(),0,ToFDB.size());
	TH1F * MeasD = new TH1F("MeasD","MeasD",ToFDB.size(),0,ToFDB.size());


	vars->ReadBranches(chainMC);
	
	for(int i=0;i<chainMC->GetEntries();i++){
		UpdateProgressBar(i, chainMC->GetEntries());
		if(i%(int)FRAC!=0) continue;
		vars->ResetVariables();
            	chainMC->GetEvent(i);
		vars->Update();
		if(ApplyCuts("IsPositive&IsMinimumBias&IsLooseCharge1",vars)){
				int kbin;
				kbin =ToFDB.GetBin(GetBetaGen(vars));
				if(kbin>0){
			 		if(ApplyCuts("IsPurePMC",vars)) GenP->Fill(kbin);
					if(ApplyCuts("IsPureDMC",vars)) GenD->Fill(kbin);
				}
				kbin =ToFDB.GetBin(GetBetaTOF(vars));
				if(kbin>0){
			 		if(ApplyCuts("IsPurePMC",vars)) MeasP->Fill(kbin);
					if(ApplyCuts("IsPureDMC",vars)) MeasD->Fill(kbin);
				}
	                        if(ApplyCuts("IsPurePMC",vars)) GenBeta->Fill(GetBetaGen(vars));
		 		 if(ApplyCuts("IsPurePMC",vars)) MeasBeta->Fill(GetBetaTOF(vars));


		}
	}

	TH1F *GenRatio = (TH1F *) GenD->Clone();
	GenRatio->Divide(GenP);
	GenRatio->SetName("GenRatio");
	GenRatio->SetTitle("GenRatio");

	TH1F *MeasRatio = (TH1F *) MeasD->Clone();
	MeasRatio->Divide(MeasP);
	MeasRatio->SetName("MeasRatio");
	MeasRatio->SetTitle("MeasRatio");


	finalResults.Add(GenP);
	finalResults.Add(GenD);
	finalResults.Add(GenRatio);
	finalResults.Add(MeasP);
	finalResults.Add(MeasD);
	finalResults.Add(MeasRatio);
		finalResults.Add(GenBeta);
		finalResults.Add(MeasBeta);
	
	finalResults.writeObjsInFolder("");

	return 0;
}

