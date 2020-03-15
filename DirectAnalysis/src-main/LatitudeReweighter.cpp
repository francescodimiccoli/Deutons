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
#include "../include/InputFileReader.h"

#include "../include/Variables.hpp"
#include "../include/Cuts.h"
#include "../include/ParallelFiller.h"

#include "../include/filesaver.h"
#include "../include/LatReweighter.h"

int main(int argc, char * argv[])
{

    TH1::SetDefaultSumw2();     	
    cout<<"****************************** FILES OPENING ***************************************"<<endl;
     string INPUT1 = "";
     string INPUT2 = "";
     string OUTPUT = "";

  
    if(argc<=2) { 
	    OUTPUT = argv[1];
    }	

    else{
	    INPUT1=argv[1];
	    INPUT2=argv[2];
	    OUTPUT=argv[3];
    }

    string refill="";
    if(argc > 4 ) 	refill = argv[4];	

    bool Refill = false;
    if(refill!="") Refill=true;
    
   TChain * chainRTI = InputFileReader(INPUT1.c_str(),"RTI");
   TChain * chainDT = InputFileReader(INPUT1.c_str(),"Event");


    FileSaver finalHistos;
    finalHistos.setName(OUTPUT.c_str());

    FileSaver finalResults;
    finalResults.setName((OUTPUT+"_Results").c_str());


    bool checkfile = finalHistos.CheckFile();

    TTree *TreeDT = NULL;

    cout<<"****************************** BINS ***************************************"<<endl;
    SetUpUsualBinning();
    
    cout<<"****************************** VARIABLES ***************************************"<<endl;
//    Variables * vars = new Variables(1);


    cout<<"****************************** ANALYSIS ***************************************"<<endl;
    LatReweighter * weighter = new LatReweighter("LatWeights","IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",500,0,150);

    if(Refill){	
	//bisogna fare un loop senza variables
//    	weighter->LoopOnRTI(DBarReader(chainRTI, false ),vars,Refill);
//     	weighter->LoopOnData(DBarReader(chainDT, false ),vars,Refill);
    }
    else weighter = new LatReweighter(finalHistos,"LatWeights");	
	
    weighter->CalculateWeights();	
    if(Refill) weighter->Save(finalHistos);
    else  weighter->SaveResults(finalResults);
	
    return 0;
}



