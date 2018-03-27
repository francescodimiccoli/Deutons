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
    
    string INPUT1(argv[1]);
    string INPUT2(argv[2]);
    string OUTPUT(argv[3]);

    string refill="";
    if(argc > 4 ) 	refill = argv[4];	
    
    bool Refill = false;
    if(refill!="") Refill=true;
    
   //TChain * chainDT = InputFileReader(INPUT1.c_str(),"Event");
  TChain * chainDT = InputFileReader(INPUT1.c_str(),"template_stuff");
  //TChain * chainDT = InputFileReader(INPUT1.c_str(),"parametri_geo");


    FileSaver finalHistos;
    finalHistos.setName(OUTPUT.c_str());

    FileSaver finalResults;
    finalResults.setName((OUTPUT+"_Results").c_str());


    bool checkfile = finalHistos.CheckFile();

    TTree *TreeDT = NULL;

    cout<<"****************************** BINS ***************************************"<<endl;
    SetUpUsualBinning();
    
    cout<<"****************************** VARIABLES ***************************************"<<endl;
    Variables * vars = new Variables();


    cout<<"****************************** ANALYSIS ***************************************"<<endl;
    LatReweighter * weighter = new LatReweighter("LatWeights","IsPositive&IsPreselected&LikelihoodCut&DistanceCut&IsOnlyFromToF",500,0,150);

    if(Refill){	
    	weighter->LoopOnData(chainDT,vars,Refill);
    }
    else weighter = new LatReweighter(finalHistos,"LatWeights");	
	
    weighter->CalculateWeights();	
    weighter->Save(finalHistos);

    return 0;
}



