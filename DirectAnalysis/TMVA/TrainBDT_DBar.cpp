//
// Compile it with:
//
//  g++ -o likelihood  TrainBDT_DBar.cpp `root-config --cflags --libs` -lTMVA


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
    if(outFname.empty()) outFname = std::string("QualityTOF.root");

    // Open  input files, get the trees

    TChain *mc = InputFileReader("FileListMC.txt","Event");

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
    factory->AddVariable("nanti", 'I');
    factory->AddVariable("chisqn[1][0]", 'F');
    factory->AddVariable("layernonusati := 7 - (pattxy&1) - (pattxy&2)*1/2 - (pattxy&4)*1/4 - (pattxy&8)*1/8 - (pattxy&16)*1/16 - (pattxy&32)*1/32 - (pattxy&64)*1/64", 'I');
    factory->AddVariable("NTofUsed := ntofh - beta_ncl", 'I');
    factory->AddVariable("diffR := TMath::Abs(rig[2]-rig[3])/rig[1]", 'F');
    factory->AddVariable("TOF_Up_Down := TMath::Abs(edep[2][0]+edep[3][0]-edep[0][0]-edep[1][0])", 'F');
   // factory->AddVariable("Richtotused := nhit_used", 'I');	
   // factory->AddVariable("RichPhEl:= np_exp_uncorr/np", 'F');
	
   std::string IsPreselected = "rig[0]>0&&(trigpatt&0x2)!=0&&chisqcn<10&&chisqtn<10&&rig[4]!=0.0&&chisqn[1][0] < 10&&chisqn[1][1] < 10&&nparticle==1&&flagp[1]==0&&flagp[2]==0&&flagp[3]==0&&flag==0";
	
    //Preselection cuts
    std::string signalCut    = "q_lay[1][0]>0&&beta>0.4&&beta<0.77&&(rig[1]/beta)*(1-beta^2)^0.5>1.7&&mass>1&&mass<2";
    std::string backgnCut    = "q_lay[1][0]>0&&beta>0.4&&beta<0.77&&(rig[1]/beta)*(1-beta^2)^0.5>1.7&&mass>0&&mass<1";
    std::string signalCutNaF = "q_lay[1][0]>0&&beta_refit>0&&beta_refit<0.96&&(rig[1]/beta_refit)*(1-beta_refit^2)^0.5>1.7&&mass>1&&mass<2&&selection==1";
    std::string backgnCutNaF = "q_lay[1][0]>0&&beta_refit>0&&beta_refit<0.96&&(rig[1]/beta_refit)*(1-beta_refit^2)^0.5>1.7&&mass>0&&mass<1&&selection==1";
    std::string signalCutAgl = "q_lay[1][0]>0&&beta_refit>0&&beta_refit<0.96&&(rig[1]/beta_refit)*(1-beta_refit^2)^0.5>1.7&&mass>1&&mass<2&&selection==0";
    std::string backgnCutAgl = "q_lay[1][0]>0&&beta_refit>0&&beta_refit<0.96&&(rig[1]/beta_refit)*(1-beta_refit^2)^0.5>1.7&&mass>0&&mass<1&&selection==0";
		 

    factory->AddTree(mc,"Signal"    ,1,(signalCut+"&&"+IsPreselected).c_str());
    factory->AddTree(mc,"Background",1,(backgnCut+"&&"+IsPreselected).c_str());

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


    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();
}
