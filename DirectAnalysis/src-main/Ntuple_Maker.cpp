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

int main(int argc, char * argv[])
{


	cout<<"****************************** FILES OPENING ***************************************"<<endl;

	TH1::SetDefaultSumw2();     	
	cout<<"****************************** FILES OPENING ***************************************"<<endl;

	string INPUT1(argv[1]);
	string INPUT2(argv[2]);
	string INPUT3(argv[3]);

	TChain * chain_RTI = InputFileReader(INPUT1.c_str(),"RTI");

	TChain * chainDT = InputFileReader(INPUT1.c_str(),"Event");
	TChain * chainMC = InputFileReader(INPUT2.c_str(),"Event");



	cout<<"****************************** BINS ***************************************"<<endl;
	SetUpUsualBinning();

	
	Variables * vars = new Variables();

	cout<<"****************************** ANALYIS ******************************************"<<endl;
	
	TFile * File = new TFile(INPUT3.c_str(), "RECREATE");

        TTree * measure_stuff= new TTree("parametri_geo","parametri_geo");
	TTree * template_stuff= new TTree("template_stuff","template_stuff");
	
	vars->RegisterBranches(measure_stuff);
	vars->RegisterTemplatesBranches(template_stuff);

	DBarReader readerDT(chainDT, false,chain_RTI);

	for(int i=0;i<readerDT.GetTreeEntries()/FRAC;i++){
		UpdateProgressBar(i, readerDT.GetTreeEntries()/FRAC);
		//if(i%(int)FRAC!=0) continue;
		readerDT.FillVariables(i,vars);
		vars->Update();
		if(ApplyCuts("IsPositive&IsMinimumBias",vars)) measure_stuff->Fill();
		if(
			(ApplyCuts("IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsGoodTime",vars))||
			(ApplyCuts("IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut",vars))||	
			(ApplyCuts("IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut",vars))	
		)		
		template_stuff->Fill();	
	}
	File->Write();
	File->Close();
		
	/*
	TFile * File2 = new TFile((INPUT3+"_MC").c_str(), "RECREATE");

	TTree * measure_stuffiMC= new TTree("parametri_MC","parametri_MC");
	TTree * template_stuffMC= new TTree("template_stuffMC","template_stuffMC");
		
	vars->RegisterBranches(measure_stuffMC);
	vars->RegisterTemplatesBranches(template_stuffMC);

	DBarReader readerMC(chainMC, true);

	for(int i=0;i<readerMC.GetTreeEntries();i++){
		UpdateProgressBar(i, readerMC.GetTreeEntries());
		if(i%(int)FRAC!=0) continue;
		readerMC.FillVariables(i,vars);
		vars->Update();
		if(ApplyCuts("IsPositive&IsMinimumBias",vars)) measure_stuffMC->Fill();
		if(
			(ApplyCuts("IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsGoodTime",vars))||
			(ApplyCuts("IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut",vars))||	
			(ApplyCuts("IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut",vars))	
		)		
			   template_stuffMC->Fill();	
	}
	File2->Write();
	File2->Close();
	*/

	return 0;
}

