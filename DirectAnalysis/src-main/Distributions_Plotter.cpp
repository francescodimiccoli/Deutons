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

#include "HistoBooker.h"

#include "Plotter.h"

int main(int argc, char * argv[])
{


	cout<<"****************************** FILES OPENING ***************************************"<<endl;

	TH1::SetDefaultSumw2();     	
	cout<<"****************************** FILES OPENING ***************************************"<<endl;

	string INPUT1(argv[1]);
	string INPUT2(argv[2]);
	string OUTPUT(argv[3]);

	string refill="";
	if(argc > 4 ) 	refill = argv[4];	

	bool Refill = false;
	if(refill!="") Refill=true;

	TChain * chainDT = InputFileReader(INPUT1.c_str(),"Event");
	TChain * chainMC = InputFileReader(INPUT2.c_str(),"Event");

	FileSaver finalHistos;
	finalHistos.setName(OUTPUT.c_str());

	FileSaver finalResults;
        finalResults.setName((OUTPUT+"_Results").c_str());


	bool checkfile = finalHistos.CheckFile();

	cout<<"****************************** BINS ***************************************"<<endl;
	SetUpUsualBinning();

	cout<<"****************************** VARIABLES ***************************************"<<endl;
	Variables * vars = new Variables();

	cout<<"****************************** ANALYIS ******************************************"<<endl;

	Plotter Plots(finalHistos,finalResults);

/*	Plots.BookMassAnalysis();
        Plots.BookCleaningCutsAnalysis();
        Plots.BookRichBDTAnalysis();
        Plots.BookBetaResMatrixAnalysis();
*/
	Plots.BookGenAcceptanceAnalysis();	

       if(Refill)	Plots.FillAllAnalyses(chainDT,chainMC);
	Plots.DoAllAnalyses();

	return 0;
}




