#include <bitset>
#include "TROOT.h"
#include "TNtuple.h"
#include <TSpline.h>
#include "../include/binning.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include <TVector3.h>
#include "TMath.h"
#include <TFile.h>
#include "TFile.h"
#include "TH2.h"
#include "TF2.h"
#include <TVector3.h>
#include "TMath.h"
#include "TGraphErrors.h"
#include "TFractionFitter.h"
#include "TRandom3.h"
#include "TChain.h"
#include "Globals.h"

#include "filesaver.h"
#include "Analyzer.h"

int main(int argc, char * argv[])
{
    TH1::SetDefaultSumw2();     	
     cout<<"********************** FILES OPENING ***********************************"<<endl;

     string INPUT1 = "";
     string INPUT2 = "";
     string OUTPUT = "";

     if(argc<=2) { 
	     OUTPUT = argv[1];
     }	

     else {
	     INPUT1 = argv[1];
	     INPUT2 = argv[2];
	     OUTPUT = argv[3];
     }
     string refill="";
     if(argc > 4 ) 	refill = argv[4];	

     bool Refill = false;
     if(refill!="") Refill=true;
	
    FileSaver finalHistosCounts;
    finalHistosCounts.setName((OUTPUT + "_Counts").c_str());
    
    FileSaver finalResults;
    finalHistosCounts.setName((OUTPUT + "_Results").c_str());


    cout<<"****************************** BINS ***************************************"<<endl;
    SetUpEffCorrBinning();
   
    cout<<"****************************** Analysis*************************************"<<endl;
    	 
    Analyzer analyzer(INPUT1,INPUT2);	 

    analyzer.BookCountsAnalysis(finalHistosCounts,finalResults,Refill);	
   
    analyzer.FillAll();

} 
