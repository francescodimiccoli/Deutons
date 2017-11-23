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
//#include "../include/GlobalBinning.h"
#include "../include/InputFileReader.h"
#include "../include/Commonglobals.cpp"
#include "../include/Resolution.h"

#include "../include/Variables.hpp"
#include "../include/Cuts.h"
#include "../include/ParallelFiller.h"

#include "../include/filesaver.h"
#include "../include/TemplateFITbetasmear.h"

void ExtractSimpleCountNr(FileSaver finalhistos, FileSaver finalResults, TTree* tree,Binning Bins,float (*discr_var) (Variables * vars),std::string name,std::string cut, bool refill);

int main(int argc, char * argv[])
{


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

    TTree *TreeDT = NULL;//(TTree *)fileDT->Get("parametri_geo");

    bool TRDCalibfound = ReadCalibration();

    cout<<"****************************** BINS ***************************************"<<endl;

        SetBins();

        PRB.Print();

        cout<<"**TOF**"<<endl;
        ToFDB.Print();

        cout<<"**NaF**"<<endl;
        NaFPB.Print();

        cout<<"**Agl**"<<endl;
        AglDB.Print();

        ToFDB.UseBetaEdges();
        NaFDB.UseBetaEdges();
        AglDB.UseBetaEdges();

        PRB.UseREdges();


        cout<<endl;


    cout<<"****************************** VARIABLES ***************************************"<<endl;

    Variables * vars = new Variables();
    TF1 * HeContTOF=0x0;
    TF1 * HeContNaF=0x0;
    TF1 * HeContAgl=0x0;

    cout<<"****************************** ANALYIS ******************************************"<<endl;
    if(finalResults.GetFile()){
        HeContTOF = (TF1 *) finalResults.Get("HeContTemplate/Model");
        HeContNaF = (TF1 *) finalResults.Get("HeContTemplate/Model");
        HeContAgl = (TF1 *) finalResults.Get("HeContTemplate/Model");
    }

    BadEventSimulator * NaFBadEvSimulator= new BadEventSimulator("IsFromNaF",22,0.72,1); 
    BadEventSimulator * AglBadEvSimulator= new BadEventSimulator("IsFromAgl",250,0.95,1); 

    TemplateFIT * TOFfits= new TemplateFIT("TOFfits",ToFDB,"IsPreselected&LikelihoodCut&DistanceCut",100,0.1,4);
    TemplateFIT * NaFfits= new TemplateFIT("NaFfits",NaFDB,"IsPreselected&LikelihoodCut&DistanceCut&IsFromNaF",100,0.1,4,true,11,2000,1000);
    TemplateFIT * Aglfits= new TemplateFIT("Aglfits",AglDB,"IsPreselected&LikelihoodCut&DistanceCut&IsFromAgl",100,0.1,4,true,11,600,500);	

    NaFfits->SetUpBadEventSimulator(NaFBadEvSimulator);
    Aglfits->SetUpBadEventSimulator(AglBadEvSimulator);
    NaFfits->SetFitWithNoiseMode();
    Aglfits->SetFitWithNoiseMode();


    if((!checkfile)||Refill){

        ParallelFiller<TemplateFIT *> Filler;

        Filler.AddObject2beFilled(TOFfits,GetRecMassTOF ,GetBetaTOF);
        Filler.AddObject2beFilled(NaFfits,GetRecMassRICH,GetBetaRICH);
        Filler.AddObject2beFilled(Aglfits,GetRecMassRICH,GetBetaRICH);

        //main loops
        Filler.LoopOnMC  (DBarReader(chainMC, true ),vars);
        Filler.LoopOnData(DBarReader(chainDT, false),vars);
        //

        TOFfits->DisableFit();
        TOFfits->Save(finalHistos);

        NaFfits->DisableFit();
        NaFfits->Save(finalHistos);

        Aglfits->DisableFit();
        Aglfits->Save(finalHistos);
    }

    else { 

        TOFfits= new TemplateFIT(finalHistos,"TOFfits",ToFDB);
        //NaFfits= new TemplateFIT(finalHistos,"NaFfits",NaFDB,true,11,400,200);
        //Aglfits= new TemplateFIT(finalHistos,"Aglfits",AglDB,true,11,110,80);

        NaFfits= new TemplateFIT(finalHistos,"NaFfits",NaFDB,true,11,2000,1000);
        Aglfits= new TemplateFIT(finalHistos,"Aglfits",AglDB,true,11,600,500);


        //TOFfits->DisableFit();
        TOFfits->SetHeliumContamination(HeContTOF);				
        TOFfits->ExtractCounts(finalHistos);	
        TOFfits->SaveFitResults(finalResults);

        NaFfits->SetFitRange(0.6,4);
        //       NaFfits->DisableFit();
        NaFfits->SetHeliumContamination(HeContNaF);
        NaFfits->ExtractCounts(finalHistos);
        NaFfits->SaveFitResults(finalResults);

        Aglfits->SetFitRange(0.6,4);
        // Aglfits->DisableFit();
        Aglfits->SetHeliumContamination(HeContAgl);
        Aglfits->ExtractCounts(finalHistos);
        Aglfits->SaveFitResults(finalResults);	
    }

    ExtractSimpleCountNr(finalHistos,finalResults,TreeDT,PRB,GetRigidity,"HEPCounts","IsPreselected&LikelihoodCut&DistanceCut&IsPrimary",false);

    return 0;
}


void ExtractSimpleCountNr(FileSaver finalhistos, FileSaver finalResults, TTree* tree,Binning Bins,float (*discr_var) (Variables * vars),std::string name,std::string cut, bool refill){

	TH1F * Counts;
	if(refill){
		Counts = new TH1F(name.c_str(),name.c_str(),Bins.size(),0,Bins.size());
		cout<<name.c_str()<<" Filling ... "<< endl;
		Variables * vars = new Variables;
		vars->ReadBranches(tree);
	
		for(int i=0;i<tree->GetEntries()/FRAC;i++){
			UpdateProgressBar(i, tree->GetEntries()/FRAC);
			tree->GetEvent(i);
			vars->Update();
			int kbin;
			kbin = 	Bins.GetBin(discr_var(vars));
			if(ApplyCuts(cut.c_str(),vars)&&kbin>0)
				Counts->Fill(kbin);
		}
	}
	else Counts = (TH1F*) finalhistos.Get((name + "/" + name+ "/" + name).c_str());
	
      if(Counts){
		finalhistos.Add(Counts);
      		finalhistos.writeObjsInFolder((name + "/" + name).c_str());

      		finalResults.Add(Counts);
      		finalResults.writeObjsInFolder((name + "/" + name).c_str());
	}
      return;	
}
	


