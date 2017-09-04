//
// Compile it with:
//
//  g++ -o likelihood  TrainBDT.cpp `root-config --cflags --libs` -lTMVA


#include <TFile.h>
#include <TMVA/Factory.h>


int main(int argc, char * argv[])
{
    //Processing input options
    int c;
    std::string outFname;
    std::string  inFname;
    /*while((c = getopt(argc, argv, "o:")) != -1) {
        if(c == 'o') outFname = std::string(optarg);
    }
    if (optind < argc) inFname = std::string(argv[optind++]); else return 1;
    */	
    if(inFname.empty())  inFname = std::string("/home/AMS/fdimicco/fdimicco/MAIN/sommaMC/temp/sommaMC.root");
    if(outFname.empty()) outFname = std::string("QualityAgl.root");

    // Open  input files, get the trees
    TFile * mcFile = new TFile(inFname.c_str());
    TTree * mc = (TTree*)mcFile->Get("parametri_geo");

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
    factory->AddVariable("NAnticluster", 'I');
    factory->AddVariable("Chisquare", 'F');
    factory->AddVariable("layernonusati := 7 - (hitbits&1) - (hitbits&2)*1/2 - (hitbits&4)*1/4 - (hitbits&8)*1/8 - (hitbits&16)*1/16 - (hitbits&32)*1/32 - (hitbits&64)*1/64", 'I');
    factory->AddVariable("NTofUsed := NTofClusters - NTofClustersusati", 'I');
    factory->AddVariable("diffR := TMath::Abs(Rup-Rdown)/R", 'F');
    factory->AddVariable("TOF_Up_Down := TMath::Abs(TOFEndep[2]+TOFEndep[3]-TOFEndep[0]-TOFEndep[1])", 'F');
    factory->AddVariable("Richtotused");	
    factory->AddVariable("RichPhEl");
	
	
    //Preselection cuts
    std::string signalCut = "qL1>0&&BetaHR>0.4&&BetaHR<0.77&&(R/BetaHR)*(1-BetaHR^2)^0.5>1.7&&GenMass>1&&GenMass<2&&(CUTMASK&187)==187";
    std::string backgnCut = "qL1>0&&BetaHR>0.4&&BetaHR<0.77&&(R/BetaHR)*(1-BetaHR^2)^0.5>1.7&&GenMass<1&&GenMass<2&&(CUTMASK&187)==187";
    std::string signalCutNaF = "qL1>0&&BetaRICH>0&&BetaRICH<0.96&&(R/BetaRICH)*(1-BetaRICH^2)^0.5>1.7&&GenMass>1&&GenMass<2&&(CUTMASK&187)==187&&(RICHmask&1023)==512";
    std::string backgnCutNaF = "qL1>0&&BetaRICH>0&&BetaRICH<0.96&&(R/BetaRICH)*(1-BetaRICH^2)^0.5>1.7&&GenMass<1&&GenMass<2&&(CUTMASK&187)==187&&(RICHmask&1023)==512";
    std::string signalCutAgl = "qL1>0&&BetaRICH>0&&BetaRICH<0.985&&(R/BetaRICH)*(1-BetaRICH^2)^0.5>1.7&&GenMass>1&&GenMass<2&&(CUTMASK&187)==187&&(RICHmask&1023)==0";
    std::string backgnCutAgl = "qL1>0&&BetaRICH>0&&BetaRICH<0.985&&(R/BetaRICH)*(1-BetaRICH^2)^0.5>1.7&&GenMass<1&&GenMass<2&&(CUTMASK&187)==187&&(RICHmask&1023)==0";
		 

    factory->AddTree(mc,"Signal"    ,1,signalCutAgl.c_str());
    factory->AddTree(mc,"Background",1,backgnCutAgl.c_str());

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
