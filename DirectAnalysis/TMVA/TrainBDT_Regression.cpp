//
// Compile it with:
//
//  g++ -o regression  TrainBDT_Regression.cpp `root-config --cflags --libs` -lTMVA


#include <TFile.h>
#include <TMVA/Factory.h>
#include "TChain.h"
#include "string.h"
#include "include/InputFileReader.h"


int main(int argc, char * argv[])
{
    //Processing input options
    int c;
    std::string outFname;
    outFname = std::string("QualityRegr.root");

    // Open  input files, get the trees
    TChain *mc = InputFileReader("FileListNtuples_red.txt","parametri_MC");
    // Preparing options for the TMVA::Factory
    std::string options( 
        "!V:" 
        "!Silent:"
        "Color:"
        "DrawProgressBar:"
        "Transformations=I;D;P;G,D:"
        "AnalysisType=Regression"
    );

    //Creating the factory
    TFile *   ldFile = new TFile(outFname.c_str(),"RECREATE");
    TMVA::Factory * factory = new TMVA::Factory("QualityRegr", ldFile, options.c_str());

    //Preparing variables 
    //general
    factory->AddVariable("Chisquare", 'F');
    factory->AddVariable("EdepTOFU", 'F');
    factory->AddVariable("EdepTOFD", 'F');
    factory->AddVariable("EdepTrack", 'F');
//  factory->AddSpectator("Beta", 'F');
//  factory->AddSpectator("R", 'F');

    factory->AddTarget("Beta", 'F');


    //Preselection cuts
    std::string PreSelection    = "qL1>0&&(joinCutmask&187)==187&&qL1<1.75&&R>0";
    std::string ChargeCut 	= "qUtof>0.8&&qUtof<1.3&&qLtof>0.8&&qLtof<1.3";
  std::string VelocityCut 	= "Beta<0.85";
    std::string signalCut 	= "GenMass<1&&GenMass<2";	

    //factory->AddTree(mc,"Signal"    ,1,(PreSelection +"&&"+ ChargeCut + "&&" + VelocityCut + "&&"+ signalCut).c_str());
    factory->AddRegressionTree(mc,1);//,(PreSelection +"&&"+ ChargeCut + "&&" + VelocityCut + "&&"+ signalCut).c_str());
    factory->SetCut(static_cast<TString>((PreSelection +"&&"+ ChargeCut + "&&" + VelocityCut + "&&"+ signalCut).c_str()));
    // Preparing
    std::string preselection = "";
    std::string inputparams(
        "SplitMode=Random:"
        "NormMode=NumEvents:"
        "!V"
    );
    factory->PrepareTrainingAndTestTree(preselection.c_str(),inputparams.c_str());

    // Training
    std::string trainparams ="!H:!V:MaxDepth=3";
    factory->BookMethod(TMVA::Types::kBDT, "BDT", trainparams.c_str());

    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();
}
