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

	TChain * chainDT = InputFileReader(INPUT1.c_str(),"Event");
	TChain * chainMC = InputFileReader(INPUT2.c_str(),"Event");



	cout<<"****************************** BINS ***************************************"<<endl;
	SetUpUsualBinning();

	
	Variables * vars = new Variables();

	cout<<"****************************** ANALYIS ******************************************"<<endl;
	
	ofstream myfile;
	myfile.open(INPUT3.c_str());

	DBarReader readerMC(chainMC, true);

	for(int i=0;i<readerMC.GetTreeEntries();i++){
		UpdateProgressBar(i, readerMC.GetTreeEntries());
		if(i%(int)FRAC!=0) continue;
		readerMC.FillVariables(i,vars);
		vars->Update();
			if((ApplyCuts("IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF&NafBetaSafetyCut&RICHBDTCut",vars))||	
			(ApplyCuts("IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl&AglBetaSafetyCut&RICHBDTCut",vars)))	
			if(ApplyCuts("IsProtonMC",vars)&&GetRecMassRICH(vars)>3.5) { cout<<endl; myfile<<"Run: "<<vars->Run<<"; Event: "<<vars->Event<<endl;	}
	}
	
	myfile.close();

	return 0;
}

