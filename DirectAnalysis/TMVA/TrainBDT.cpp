//
// Compile it with:
//
//  g++ -o likelihood  TrainBDT.cpp `root-config --cflags --libs` -lTMVA


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
    outFname = std::string("QualityAgl.root");

    // Open  input files, get the trees
    TChain *mc = InputFileReader("FileListNtuples_ext.txt","parametri_geo");
    // Preparing options for the TMVA::Factory
    std::string options( 
        "!V:" 
        "!Silent:"
        "Color:"
        "DrawProgressBar:"
        "Transformations=I;D;P;G,D:"
        "AnalysisType=Classification"
    );

    //Creating the factory
    TFile *   ldFile = new TFile(outFname.c_str(),"RECREATE");
    TMVA::Factory * factory = new TMVA::Factory("QualityAgl", ldFile, options.c_str());

    //Preparing variables 
    //general
    /*factory->AddVariable("Chisquare", 'F');
    factory->AddVariable("Layernonusati", 'I');
    factory->AddVariable("NTofUsed", 'I');
    factory->AddVariable("diffR", 'F');
    factory->AddVariable("TOF_Up_Down", 'F');*/
    //Tof	
    //factory->AddVariable("TOFchisq_s", 'F');
    //factory->AddVariable("TOFchisq_t", 'F');

    //RICH	
    factory->AddVariable("Richtotused", 'F');	
    factory->AddVariable("RichPhEl", 'F');
    factory->AddVariable("RICHprob", 'F');
    factory->AddVariable("RICHcollovertotal");
    factory->AddVariable("RICHLipBetaConsistency");  
    factory->AddVariable("RICHTOFBetaConsistency");  
    factory->AddVariable("RICHChargeConsistency");
    
    factory->AddVariable("RICHPmts");
    factory->AddVariable("RICHgetExpected");		
    factory->AddVariable("tot_hyp_p_uncorr");
    factory->AddVariable("Bad_ClusteringRICH");
    factory->AddVariable("NSecondariesRICHrich");


    //Preselection cuts
    std::string PreSelection    = "qL1>0&&(joinCutmask&187)==187&&qL1<1.75&&R>0";
    std::string ChargeCut 	= "qUtof>0.8&&qUtof<1.3&&qLtof>0.8&&qLtof<1.3";
    std::string VelocityCut 	= /*"Beta<0.8";*/"((joinCutmask>>11))==0&&BetaRICH_new>0&&BetaRICH_new<0.99";
    std::string signalCut 	= /*"(R/Beta)*(1-Beta^2)^0.5>1.65&&GenMass>1&&GenMass<2";*/"(R/BetaRICH_new)*(1-BetaRICH_new^2)^0.5>0.5&&(R/BetaRICH_new)*(1-BetaRICH_new^2)^0.5<1.5";	
    std::string bkgndCut 	= /*"(R/Beta)*(1-Beta^2)^0.5>1.65&&GenMass>0&&GenMass<1";*/"(R/BetaRICH_new)*(1-BetaRICH_new^2)^0.5>3";		 

    factory->AddTree(mc,"Signal"    ,1,(PreSelection +"&&"+ ChargeCut + "&&" + VelocityCut + "&&"+ signalCut).c_str());
    factory->AddTree(mc,"Background",1,(PreSelection +"&&"+ ChargeCut + "&&" + VelocityCut + "&&"+ bkgndCut).c_str());

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

    trainparams ="!H:!V";
    factory->BookMethod(TMVA::Types::kLikelihood, "Likelihood", trainparams.c_str());

    trainparams ="!H:!V:VarTransform=Decorrelate";
    //factory->BookMethod(TMVA::Types::kLikelihood, "LikelihoodD", trainparams.c_str());

    trainparams ="!H:!V";
    //factory->BookMethod(TMVA::Types::kCuts, "Cuts", trainparams.c_str());



    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();
}
