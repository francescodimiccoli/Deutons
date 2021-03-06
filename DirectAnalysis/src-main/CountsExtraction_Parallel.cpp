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
#include "TProof.h"
#include "../include/InputFileReader.h"
#include "../include/LatReweighter.h"

#include "../include/Variables.hpp"
#include "../include/Cuts.h"
#include "../include/ParallelFiller.h"
#include "../include/Efficiency.h"

#include "../include/filesaver.h"
#include "../include/TemplateFITbetasmear.h"

void ExtractSimpleCountNr(FileSaver finalhistos, FileSaver finalResults, DBarReader readerDT, bool refill);

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
	
	else {
	INPUT1 = argv[1];
	INPUT2 = argv[2];
	OUTPUT = argv[3];
	}
	string refill="";
	if(argc > 4 ) 	refill = argv[4];	
	
	bool Refill = false;
	if(refill!="") Refill=true;

 
  TChain * chain_RTI = InputFileReader(INPUT1.c_str(),"RTI");
 
  //TChain * chainDT = InputFileReader(INPUT1.c_str(),"template_stuff");
  //TChain * chainMC = InputFileReader(INPUT2.c_str(),"template_stuffMC");
  TChain * chainDT = InputFileReader(INPUT1.c_str(),"Event");
  TChain * chainMC = InputFileReader(INPUT2.c_str(),"Event");
  //TChain * chainDT = InputFileReader(INPUT1.c_str(),"parametri_geo");
  //TChain * chainMC = InputFileReader(INPUT2.c_str(),"parametri_MC");


    FileSaver LatWeights;
    LatWeights.setName("./LatWeights/Weights.root");
    LatReweighter * weighter = new LatReweighter(LatWeights,"LatWeights");	
    

    FileSaver finalHistos;
    finalHistos.setName(OUTPUT.c_str());

    FileSaver finalResults;
    finalResults.setName((OUTPUT+"_Results").c_str());


    bool checkfile = finalHistos.CheckFile();

    TTree *TreeDT = NULL;

    cout<<"****************************** BINS ***************************************"<<endl;
	SetUpEffCorrBinning();
    
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

    ExtractSimpleCountNr(finalHistos,finalResults,DBarReader(chainDT, false,chain_RTI),Refill);



  //  TemplateFIT * SmearingCheck = new TemplateFIT("SmearingCheck",PRB,"IsPositive&IsPreselected&LikelihoodCut&DistanceCut&IsOnlyFromToF",60,0.3,1.6);	
    TemplateFIT * TOFfits= new TemplateFIT("TOFfits",ToFDB,"IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsGoodTime&IsOnlyFromToF"       ,150,0.4,7.5);
    TemplateFIT * NaFfits= new TemplateFIT("NaFfits",NaFDB,"IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut"           ,60,0.4,5,true,11,400,200);
    TemplateFIT * Aglfits= new TemplateFIT("Aglfits",AglDB,"IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut"           ,60,0.4,5,true,11,110,80);	

  // SmearingCheck->SetLatitudeReweighter(weighter);
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
        Filler.LoopOnData(DBarReader(chainDT, false,chain_RTI),vars);
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

        TOFfits->DisableFit();
        //TOFfits->SetFitRange(0.6,4);
        TOFfits->SetFitConstraints(0.9,1,0.015,0.06,0.005,0.015,true);
	TOFfits->SetHeliumContamination(HeContTOF);
	TOFfits->ExtractCounts(finalHistos);	
        TOFfits->SaveFitResults(finalResults);

        NaFfits->SetFitRange(0.6,5);
        //NaFfits->DisableFit();
        NaFfits->SetFitConstraints(0.9,1,0.001,0.1,0.0001,0.0005);
	NaFfits->SetHeliumContamination(HeContNaF);
        NaFfits->ExtractCounts(finalHistos);
	NaFfits->SaveFitResults(finalResults);

        Aglfits->SetFitRange(0.6,5);
        //Aglfits->DisableFit();
	Aglfits->SetFitConstraints(0.9,1,0.001,0.1,0.0001,0.0005);
        Aglfits->SetHeliumContamination(HeContAgl);
  	Aglfits->ExtractCounts(finalHistos);
        Aglfits->SaveFitResults(finalResults);	
 	
    }

      return 0;
}


void ExtractSimpleCountNr(FileSaver finalhistos, FileSaver finalResults, DBarReader readerDT, bool refill){

	std::string nameHE      = "HEPCounts" ;
	std::string namequalHE  = "HEPCountsQual" ;
	std::string nameTOF 	= "TOFPCounts";
	std::string nameNaF 	= "NaFPCounts";
	std::string nameAgl 	= "AglPCounts";


	Variables * vars = new Variables();
	Efficiency * CountsHE = new Efficiency(finalhistos,nameHE ,nameHE ,PRB,"IsPositive&IsPrimary&IsMinimumBias&IsLooseCharge1","IsPositive&IsPrimary&IsMinimumBias&IsLooseCharge1");
	Efficiency * CountsQualHE = new Efficiency(finalhistos,namequalHE ,namequalHE ,PRB,"IsPositive&IsPrimary&IsMinimumBias&IsLooseCharge1&IsCleaning","IsPositive&IsPrimary&IsMinimumBias&IsLooseCharge1&IsCleaning");
	Efficiency * CountsTOF= new Efficiency(finalhistos,nameTOF,nameTOF,ToFPB,"IsPositive&IsPrimary&IsMinimumBias&IsLooseCharge1&IsCleaning&IsGoodTime&IsOnlyFromToF","IsPositive&IsPrimary&IsMinimumBias&IsLooseCharge1&IsCleaning&IsGoodTime&IsOnlyFromToF");
	Efficiency * CountsNaF= new Efficiency(finalhistos,nameNaF,nameNaF,NaFPB,"IsPositive&IsPrimary&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut"    ,"IsPositive&IsPrimary&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut"    );
	Efficiency * CountsAgl= new Efficiency(finalhistos,nameAgl,nameAgl,AglPB,"IsPositive&IsPrimary&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut"    ,"IsPositive&IsPrimary&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut"    );

	ParallelFiller<Efficiency *> Filler;
	Filler.AddObject2beFilled(CountsHE,GetRigidity,GetRigidity); 
	Filler.AddObject2beFilled(CountsQualHE,GetRigidity,GetRigidity); 
	Filler.AddObject2beFilled(CountsTOF,GetRigidity,GetRigidity); 
	Filler.AddObject2beFilled(CountsNaF,GetRigidity,GetRigidity); 
	Filler.AddObject2beFilled(CountsAgl,GetRigidity,GetRigidity); 
	Filler.ReinitializeAll(refill);

	Filler.LoopOnData(readerDT,vars);
		
	CountsHE->Save(finalhistos);
      	CountsQualHE->Save(finalhistos);
     	CountsTOF->Save(finalhistos);
      	CountsNaF->Save(finalhistos);
      	CountsAgl->Save(finalhistos);
      
	TH1F * Counts_P = (TH1F*) CountsHE->GetBefore();
	Counts_P->SetName(nameHE.c_str());
	Counts_P->SetTitle(nameHE.c_str());

	if(Counts_P){
      		finalResults.Add(Counts_P);
      		finalResults.writeObjsInFolder((nameHE + "/" + nameHE).c_str());
	}

	TH1F * Counts_QualP = (TH1F*) CountsQualHE->GetBefore();
	Counts_QualP->SetName(nameHE.c_str());
	Counts_QualP->SetTitle(nameHE.c_str());

	if(Counts_QualP){
      		finalResults.Add(Counts_QualP);
      		finalResults.writeObjsInFolder((namequalHE + "/" + nameHE).c_str());
	}

	TH1F * Counts_TOFP = (TH1F*) CountsTOF->GetBefore();
	Counts_TOFP->SetName(nameTOF.c_str());
	Counts_TOFP->SetTitle(nameTOF.c_str());

	if(Counts_TOFP){
      		finalResults.Add(Counts_TOFP);
      		finalResults.writeObjsInFolder((nameTOF + "/" + nameTOF).c_str());
	}

	TH1F * Counts_NaFP = (TH1F*) CountsNaF->GetBefore();
	Counts_NaFP->SetName(nameNaF.c_str());
	Counts_NaFP->SetTitle(nameNaF.c_str());

	if(Counts_NaFP){
      		finalResults.Add(Counts_NaFP);
      		finalResults.writeObjsInFolder((nameNaF + "/" + nameNaF).c_str());
	}
	TH1F * Counts_AglP = (TH1F*) CountsAgl->GetBefore();
	Counts_AglP->SetName(nameAgl.c_str());
	Counts_AglP->SetTitle(nameAgl.c_str());

	if(Counts_AglP){
      		finalResults.Add(Counts_AglP);
      		finalResults.writeObjsInFolder((nameAgl + "/" + nameAgl).c_str());
	}

 
	return;	
}
	


