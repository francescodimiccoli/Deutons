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
#include "Efficiency.h"
#include "../include/TemplateFITbetasmear.h"
#include "ChargeFitter.h"

void ExtractContaminationWeight(TemplateFIT * HeContTemplate, FileSaver finalResults);

int main(int argc, char * argv[])
{


	cout<<"****************************** FILES OPENING ***************************************"<<endl;

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

	TChain *		chain_RTI  	= InputFileReader(INPUT1.c_str(),"RTI");
	TChain *		chainDT    	= InputFileReader(INPUT1.c_str(),"Event");
	TChain *		chainMC    	= InputFileReader(INPUT2.c_str(),"Event");
	TChain *		chainDT_Cpct    = InputFileReader(INPUT1.c_str(),"Compact");
	TChain *		chainMC_Cpct    = InputFileReader(INPUT2.c_str(),"Compact");

	FileSaver finalHistos;
	finalHistos.setName(OUTPUT.c_str());

	FileSaver finalResults;
    	finalResults.setName((OUTPUT+"_Results").c_str());
	
	bool checkfile = finalHistos.CheckFile();
	if(!checkfile) Refill=true;

	cout<<"****************************** VARIABLES ***************************************"<<endl;

        Variables * vars = new Variables;

	cout<<"****************************** BINS ***************************************"<<endl;
	SetUpUsualBinning();
	
	cout<<"****************************** ANALYIS ******************************************"<<endl;
	
	BadEventSimulator * NaFBadEvSimulator= new BadEventSimulator("IsFromNaF",22,0.72,1); 
    	BadEventSimulator * AglBadEvSimulator= new BadEventSimulator("IsFromAgl",250,0.95,1); 


	/*
	Efficiency * TestReweigthingP  = new Efficiency(finalHistos,"TestReweigthingP","TestReweigthingP",PRB,"IsProtonMC" ,"IsProtonMC");
	Efficiency * TestReweigthingHe = new Efficiency(finalHistos,"TestReweigthingHe","TestReweigthingHe",PRB,"IsHeliumMC" ,"IsHeliumMC");

	ParallelFiller<Efficiency *> Filler1;
	Filler1.AddObject2beFilled(TestReweigthingP,GetRigidity,GetRigidity);
	Filler1.AddObject2beFilled(TestReweigthingHe,GetRigidity,GetRigidity);
	Filler1.ReinitializeAll(false);
	//main loop
	Filler1.LoopOnMC(DBarReader(chainMC, true),vars);
	
	TestReweigthingP ->Save(finalHistos);
	TestReweigthingHe->Save(finalHistos);
	
	TestReweigthingP ->Save(finalResults);
	TestReweigthingHe->Save(finalResults);
	*/

	ChargeFitter * HeContTOF = new ChargeFitter(finalHistos,"HeContTOF",Global.GetToFDBins(),"IsPositive&IsPreselectedInner&TofBetaSafetyCut","IsPositive&IsPreselected&IsClearQ1ExceptL2&TofBetaSafetyCut","IsPositive&IsPreselectedHe&IsClearQ2ExceptL2&TofBetaSafetyCut","IsPositive&IsPreselectedHe&IsCleaning&TofBetaSafetyCut","IsPositive&IsPreselected&IsCleaning&TofBetaSafetyCut",300,0,4);

	ChargeFitter * HeContNaF = new ChargeFitter(finalHistos,"HeContNaF",Global.GetNaFDBins(),"IsPositive&IsPreselectedInner&IsFromNaF&RICHBDTCut&NafBetaSafetyCut","IsPositive&IsPreselected&IsClearQ1ExceptL2&IsFromNaF&RICHBDTCut&NafBetaSafetyCut","IsPositive&IsPreselectedHe&IsClearQ2ExceptL2&IsFromNaF&RICHBDTCut&NafBetaSafetyCut","IsPositive&IsPreselectedHe&IsCleaning&IsFromNaF&RICHBDTCut&NafBetaSafetyCut","IsPositive&IsPreselected&IsCleaning&IsFromNaF&RICHBDTCut&NafBetaSafetyCut",300,0,4);

	ChargeFitter * HeContAgl = new ChargeFitter(finalHistos,"HeContAgl",Global.GetAglDBins(),"IsPositive&IsPreselectedInner&IsFromAgl&RICHBDTCut&AglBetaSafetyCut","IsPositive&IsPreselected&IsClearQ1ExceptL2&IsFromAgl&RICHBDTCut&AglBetaSafetyCut","IsPositive&IsPreselectedHe&IsClearQ2ExceptL2&IsFromAgl&RICHBDTCut&AglBetaSafetyCut","IsPositive&IsPreselectedHe&IsCleaning&IsFromAgl&RICHBDTCut&AglBetaSafetyCut","IsPositive&IsPreselected&IsCleaning&IsFromAgl&RICHBDTCut&AglBetaSafetyCut",300,0,4);

	ChargeFitter * HeContTOFCheck[10];
	for(int i=0;i<10;i++){
		 std::string basename = "HeContCheck" + to_string(i);
		 std::string FragmentCut = "IsPositive&IsCleaning&IsPreselectedHe" + to_string(i);
		 HeContTOFCheck[i] = new ChargeFitter(finalHistos,basename,Global.GetToFDBins(),"IsPositive&IsPreselectedInner","IsPositive&IsPreselected&IsClearQ1ExceptL2","IsPositive&IsPreselectedHe&IsClearQ2ExceptL2",FragmentCut,"IsPositive&IsPreselected&IsCleaning",300,0,4);
	}


	HeContNaF->SetUpBadEventSimulator(NaFBadEvSimulator);
	HeContAgl->SetUpBadEventSimulator(AglBadEvSimulator);

	ParallelFiller<ChargeFitter *> Filler2;
	Filler2.AddObject2beFilled(HeContTOF,GetBetaTOF,GetLoweredBetaTOF);
	Filler2.AddObject2beFilled(HeContNaF,GetBetaRICH,GetLoweredBetaNaF);
	Filler2.AddObject2beFilled(HeContAgl,GetBetaRICH,GetLoweredBetaAgl);
	for(int i=0;i<10;i++) Filler2.AddObject2beFilled(HeContTOFCheck[i],GetBetaTOF,GetLoweredBetaTOF);
	Filler2.ReinitializeAll(Refill);
	//main loops
	Filler2.LoopOnData(DBarReader(chainDT, false ,chain_RTI,chainDT_Cpct),vars);
	Filler2.LoopOnMC(DBarReader(chainMC, true ,chain_RTI,chainMC_Cpct),vars);

	HeContTOF->Save(finalHistos);
	HeContNaF->Save(finalHistos);
	HeContAgl->Save(finalHistos);
	for(int i=0;i<10;i++) HeContTOFCheck[i]->Save(finalHistos);

	if(!Refill) {
	HeContTOF->FitDistributions(0.8, 2.3);
	HeContNaF->FitDistributions(0.8, 2.3);
	HeContAgl->FitDistributions(0.8, 2.3);
	for(int i=0;i<10;i++) HeContTOFCheck[i]->FitDistributions(0.8, 2.3);

	HeContTOF->CreateFragmentsMass(1.85);
	HeContNaF->CreateFragmentsMass(1.85);
	HeContAgl->CreateFragmentsMass(1.85);
	for(int i=0;i<10;i++) HeContTOFCheck[i]->CreateFragmentsMass(1.65+0.05*i);

	HeContTOF->SaveResults(finalResults);
	HeContNaF->SaveResults(finalResults);
	HeContAgl->SaveResults(finalResults);
	for(int i=0;i<10;i++) HeContTOFCheck[i]->SaveResults(finalResults);
	}
	return 0;
}


