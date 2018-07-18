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

	string INPUT1(argv[1]);
	string INPUT2(argv[2]);
	string OUTPUT(argv[3]);

	string refill="";
	if(argc > 4 ) 	refill = argv[4];	

	bool Refill = false;
	if(refill!="") Refill=true;

	TChain * chain_RTI = InputFileReader(INPUT1.c_str(),"RTI");

	TChain * chainDT = InputFileReader(INPUT1.c_str(),"parametri_geo");
	TChain * chainMC = InputFileReader(INPUT2.c_str(),"parametri_MC");
	//TChain * chainDT = InputFileReader(INPUT1.c_str(),"Event");
	//TChain * chainMC = InputFileReader(INPUT2.c_str(),"Event");


	FileSaver finalHistos;
	finalHistos.setName(OUTPUT.c_str());

	FileSaver finalResults;
    	finalResults.setName((OUTPUT+"_Results").c_str());
	
	bool checkfile = finalHistos.CheckFile();

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

	ChargeFitter * HeContTOF = new ChargeFitter(finalHistos,"HeContTOF",ToFDB,"IsPositive&IsPreselectedInner","IsPositive&IsPreselected&IsClearQ1ExceptL2","IsPositive&IsPreselectedHe&IsClearQ2ExceptL2","IsPositive&IsPreselectedHe&LikelihoodCut&DistanceCut","IsPositive&IsPreselected&LikelihoodCut&DistanceCut",300,0,4);
	ChargeFitter * HeContNaF = new ChargeFitter(finalHistos,"HeContNaF",NaFDB,"IsPositive&IsPreselectedInner&IsFromNaF&RICHBDTCut","IsPositive&IsPreselected&IsClearQ1ExceptL2&IsFromNaF&RICHBDTCut","IsPositive&IsPreselectedHe&IsClearQ2ExceptL2&IsFromNaF&RICHBDTCut","IsPositive&IsPreselectedHe&LikelihoodCut&DistanceCut&IsFromNaF&RICHBDTCut","IsPositive&IsPreselected&LikelihoodCut&DistanceCut&IsFromNaF&RICHBDTCut",300,0,4);
	ChargeFitter * HeContAgl = new ChargeFitter(finalHistos,"HeContAgl",AglDB,"IsPositive&IsPreselectedInner&IsFromAgl&RICHBDTCut","IsPositive&IsPreselected&IsClearQ1ExceptL2&IsFromAgl&RICHBDTCut","IsPositive&IsPreselectedHe&IsClearQ2ExceptL2&IsFromAgl&RICHBDTCut","IsPositive&IsPreselectedHe&LikelihoodCut&DistanceCut&IsFromAgl&RICHBDTCut","IsPositive&IsPreselected&LikelihoodCut&DistanceCut&IsFromAgl&RICHBDTCut",300,0,4);


	ChargeFitter * HeContTOFCheck[10];
	for(int i=0;i<10;i++){
		 std::string basename = "HeContCheck" + to_string(i);
		 std::string FragmentCut = "IsPositive&LikelihoodCut&DistanceCut&IsPreselectedHe" + to_string(i);
		 HeContTOFCheck[i] = new ChargeFitter(finalHistos,basename,ToFDB,"IsPositive&IsPreselectedInner","IsPositive&IsPreselected&IsClearQ1ExceptL2","IsPositive&IsPreselectedHe&IsClearQ2ExceptL2",FragmentCut,"IsPositive&IsPreselected&LikelihoodCut&DistanceCut",300,0,4);
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
	Filler2.LoopOnData(DBarReader(chainDT, false,chain_RTI),vars);
	Filler2.LoopOnMC(DBarReader(chainMC, true),vars);

	HeContTOF->Save(finalHistos);
	HeContNaF->Save(finalHistos);
	HeContAgl->Save(finalHistos);
	for(int i=0;i<10;i++) HeContTOFCheck[i]->Save(finalHistos);

	if(!Refill) {
	HeContTOF->FitDistributions(0.8, 2.3);
	HeContNaF->FitDistributions(0.8, 2.3);
	HeContAgl->FitDistributions(0.8, 2.3);
	for(int i=0;i<10;i++) HeContTOFCheck[i]->FitDistributions(0.8, 2.3);

	HeContTOF->CreateFragmentsMass(1.75);
	HeContNaF->CreateFragmentsMass(1.75);
	HeContAgl->CreateFragmentsMass(1.75);
	for(int i=0;i<10;i++) HeContTOFCheck[i]->CreateFragmentsMass(1.65+0.05*i);

	HeContTOF->SaveResults(finalResults);
	HeContNaF->SaveResults(finalResults);
	HeContAgl->SaveResults(finalResults);
	for(int i=0;i<10;i++) HeContTOFCheck[i]->SaveResults(finalResults);
	}
	return 0;
}


