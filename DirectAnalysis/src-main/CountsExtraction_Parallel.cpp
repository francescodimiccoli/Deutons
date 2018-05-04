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
#include "../include/LatReweighter.h"

#include "../include/Variables.hpp"
#include "../include/Cuts.h"
#include "../include/ParallelFiller.h"

#include "../include/filesaver.h"
#include "../include/TemplateFITbetasmear.h"

void ExtractSimpleCountNr(FileSaver finalhistos, FileSaver finalResults, DBarReader readerDT, Binning Bins,float (*discr_var) (Variables * vars),std::string name,std::string cut, bool refill);

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
    
  TChain * chainDT = InputFileReader(INPUT1.c_str(),"Event");
  TChain * chainMC = InputFileReader(INPUT2.c_str(),"Event");
  //TChain * chainDT = InputFileReader(INPUT1.c_str(),"template_stuff");
  //TChain * chainMC = InputFileReader(INPUT2.c_str(),"template_stuffMC");
  //TChain * chainDT = InputFileReader(INPUT1.c_str(),"parametri_geo");
  //TChain * chainMC = InputFileReader(INPUT2.c_str(),"parametri_MC");

    FileSaver LatWeights;
    LatWeights.setName("/afs/cern.ch/work/f/fdimicco/private/Deutons/DirectAnalysis/LatWeights/Weights.root");
    LatReweighter * weighter = new LatReweighter(LatWeights,"LatWeights");	
    

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

  //  TemplateFIT * SmearingCheck = new TemplateFIT("SmearingCheck",PRB,"IsPositive&IsPreselected&LikelihoodCut&DistanceCut&IsOnlyFromToF",60,0.3,1.6);	
    TemplateFIT * TOFfits= new TemplateFIT("TOFfits",ToFDB,"IsPositive&IsPreselected&LikelihoodCut&DistanceCut&IsOnlyFromToF",150,0.4,7.5);
    TemplateFIT * NaFfits= new TemplateFIT("NaFfits",NaFDB,"IsPositive&IsPreselected&LikelihoodCut&DistanceCut&IsFromNaF&RICHBDTCut",60,0.4,5,true,11,400,200);
    TemplateFIT * Aglfits= new TemplateFIT("Aglfits",AglDB,"IsPositive&IsPreselected&LikelihoodCut&DistanceCut&IsFromAgl&RICHBDTCut",60,0.4,5,true,11,110,80);	

  //  SmearingCheck->SetLatitudeReweighter(weighter);
    TOFfits->SetLatitudeReweighter(weighter);	
    NaFfits->SetLatitudeReweighter(weighter);	
    Aglfits->SetLatitudeReweighter(weighter);	

    NaFfits->SetUpBadEventSimulator(NaFBadEvSimulator);
    Aglfits->SetUpBadEventSimulator(AglBadEvSimulator);
    NaFfits->SetFitWithNoiseMode();
    Aglfits->SetFitWithNoiseMode();


    if((!checkfile)||Refill){

        ParallelFiller<TemplateFIT *> Filler;

//	Filler.AddObject2beFilled(SmearingCheck,GetBetaTOF,GetRigidity);
        Filler.AddObject2beFilled(TOFfits,GetRecMassTOF ,GetBetaTOF);
        Filler.AddObject2beFilled(NaFfits,GetRecMassRICH,GetBetaRICH);
        Filler.AddObject2beFilled(Aglfits,GetRecMassRICH,GetBetaRICH);

        //main loops
        Filler.LoopOnMC  (DBarReader(chainMC, true ),vars);
        Filler.LoopOnData(DBarReader(chainDT, false),vars);
        //

//	SmearingCheck->DisableFit();
  //      SmearingCheck->Save(finalHistos);
	
        TOFfits->DisableFit();
        TOFfits->Save(finalHistos);

        NaFfits->DisableFit();
        NaFfits->Save(finalHistos);

        Aglfits->DisableFit();
        Aglfits->Save(finalHistos);
    }

    else { 
	//TemplateFIT * SmearingCheck = new TemplateFIT(finalHistos,"SmearingCheck",PRB);
        TOFfits= new TemplateFIT(finalHistos,"TOFfits",ToFDB);
        NaFfits= new TemplateFIT(finalHistos,"NaFfits",NaFDB,true,11,400,200);
        Aglfits= new TemplateFIT(finalHistos,"Aglfits",AglDB,true,11,110,80);
//	NaFfits->SetFitWithNoiseMode();
    //	Aglfits->SetFitWithNoiseMode();


	//SmearingCheck->SetFitRangeByQuantiles(0.05,0.95);
	//SmearingCheck->SetFitConstraints(0.99,1,0.0001,0.001,0.01,0.001);
	//SmearingCheck->DisableFit();
    	//SmearingCheck->ExtractCounts(finalHistos);	
        //SmearingCheck->SaveFitResults(finalResults);

        //TOFfits->DisableFit();
        TOFfits->SetFitRange(0.6,4);
        TOFfits->SetFitConstraints(0.9,1,0.015,0.06,0.005,0.015,true);
	TOFfits->SetHeliumContamination(HeContTOF);
	TOFfits->ExtractCounts(finalHistos);	
        TOFfits->SaveFitResults(finalResults);

        NaFfits->SetFitRange(0.6,5);
        //NaFfits->DisableFit();
        NaFfits->SetFitConstraints(0.9,1,0.001,0.1,0.0001,0.0005,true);
	NaFfits->SetHeliumContamination(HeContNaF);
       NaFfits->ExtractCounts(finalHistos);
	NaFfits->SaveFitResults(finalResults);

        Aglfits->SetFitRange(0.6,5);
        // Aglfits->DisableFit();
	 Aglfits->SetFitConstraints(0.9,1,0.001,0.1,0.0001,0.0005);

        Aglfits->SetHeliumContamination(HeContAgl);
  	Aglfits->ExtractCounts(finalHistos);
        Aglfits->SaveFitResults(finalResults);	
 	
    }

    ExtractSimpleCountNr(finalHistos,finalResults,DBarReader(chainDT, false),PRB,GetRigidity,"HEPCounts","IsPositive&IsPreselected&LikelihoodCut&DistanceCut&IsPrimary",Refill);

    return 0;
}


void ExtractSimpleCountNr(FileSaver finalhistos, FileSaver finalResults, DBarReader readerDT, Binning Bins,float (*discr_var) (Variables * vars),std::string name,std::string cut, bool refill){

	TH1F * Counts;
	if(refill){
		Counts = new TH1F(name.c_str(),name.c_str(),Bins.size(),0,Bins.size());
		cout<<name.c_str()<<" Filling ... "<< endl;
		Variables * vars = new Variables;
		vars->ReadBranches(readerDT.GetTree());	
		for(int i=0;i<readerDT.GetTreeEntries()/FRAC;i++){
			UpdateProgressBar(i, readerDT.GetTreeEntries()/FRAC);
			vars->ResetVariables();	
			
			readerDT.GetTree()->GetEvent(i);
			int kbin;
			kbin = 	Bins.GetBin(discr_var(vars));
			if(ApplyCuts(cut.c_str(),vars)&&kbin>0)
				Counts->Fill(kbin,vars->PrescaleFactor);
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
	


