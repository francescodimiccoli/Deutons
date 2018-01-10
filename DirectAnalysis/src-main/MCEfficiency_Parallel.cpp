#include <bitset>
#include "TROOT.h"
#include "TNtuple.h"
#include <TSpline.h>
#include "../../include/binning.h"
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
#include "../include/InputFileReader.h"
#include "../include/Globals.h"

#include "../include/Variables.hpp"
#include "../include/Cuts.h"
#include "../include/ParallelFiller.h"

#include "../include/filesaver.h"

#include "../include/Efficiency.h"
#include "../include/AllRangesEfficiency.h"

int main(int argc, char * argv[])
{


        cout<<"****************************** FILES OPENING ***************************************"<<endl;

        string INPUT1(argv[1]);
        string INPUT2(argv[2]);
        string OUTPUT(argv[3]);

	string refill="";
        if(argc > 4 )   refill = argv[4];

        bool Refill = false;
        if(refill!="") Refill=true;

        FileSaver finalHistos;
        finalHistos.setName(OUTPUT.c_str());
        bool checkfile = finalHistos.CheckFile();

	FileSaver finalResults;
        finalResults.setName((OUTPUT+"_Results").c_str());

	TChain * chainDT = InputFileReader(INPUT1.c_str(),"Event");
    	TChain * chainMC = InputFileReader(INPUT2.c_str(),"Event");


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

        Variables * vars = new Variables;

	cout<<"****************************** ANALYIS ******************************************"<<endl;


	//Efficiency * RigBinFullSetEff = new Efficiency(finalHistos,"RigBinFullSetEff","RigBinFullSetEff",PRB,"IsProtonMC","IsProtonMC&IsPreselected&DistanceCut&LikelihoodCut");


	AllRangesEfficiency * Preselections_P = new AllRangesEfficiency(finalHistos,"PresEff_P","PreselectionEfficiency","IsProtonMC","IsProtonMC","IsProtonMC","IsProtonMC&IsPreselected","IsProtonMC&IsPreselected&IsFromNaF","IsProtonMC&IsPreselected&IsFromAgl",Refill);
	AllRangesEfficiency * Quality_P       = new AllRangesEfficiency(finalHistos,"QualEff_P","QualityEfficiency"     ,"IsProtonMC&IsPreselected","IsProtonMC&IsPreselected&DistanceCut&LikelihoodCut&TemplatesMassCut",Refill);
        AllRangesEfficiency * FullSet_P      = new AllRangesEfficiency(finalHistos,"FullsetEff_P","FullsetEfficiency"  ,"","",Refill);


	AllRangesEfficiency * Preselections_D = new AllRangesEfficiency(finalHistos,"PresEff_D","PreselectionEfficiency","IsDeutonMC","IsDeutonMC","IsDeutonMC","IsDeutonMC&IsPreselected","IsDeutonMC&IsPreselected&IsFromNaF","IsDeutonMC&IsPreselected&IsFromAgl",Refill);
	AllRangesEfficiency * Quality_D       = new AllRangesEfficiency(finalHistos,"QualEff_D","QualityEfficiency"     ,"IsDeutonMC&IsPreselected","IsDeutonMC&IsPreselected&DistanceCut&LikelihoodCut&TemplatesMassCut",Refill);
        AllRangesEfficiency * FullSet_D       = new AllRangesEfficiency(finalHistos,"FullsetEff_D","FullsetEfficiency"  ,"","",Refill);


	AllRangesEfficiency * RICH_P = new AllRangesEfficiency(finalHistos,"RICHEff_P","RICHEfficiency","IsProtonMC&IsMinimumBias","IsProtonMC&IsMinimumBias","IsProtonMC&IsMinimumBias","IsProtonMC&IsMinimumBias","IsProtonMC&IsMinimumBias&IsFromNaF","IsProtonMC&IsMinimumBias&IsFromAgl",Refill);
	AllRangesEfficiency * RICH_D = new AllRangesEfficiency(finalHistos,"RICHEff_D","RICHEfficiency","IsDeutonMC&IsMinimumBias","IsDeutonMC&IsMinimumBias","IsDeutonMC&IsMinimumBias","IsDeutonMC&IsMinimumBias","IsDeutonMC&IsMinimumBias&IsFromNaF","IsDeutonMC&IsMinimumBias&IsFromAgl",Refill);


//	RigBinFullSetEff->Fill(treeMC,vars,GetGenMomentum,Refill);

		ParallelFiller<AllRangesEfficiency *> Filler;
		Filler.AddObject2beFilled(Preselections_P,GetBetaGen,GetBetaGen);
		Filler.AddObject2beFilled(Quality_P	 ,GetBetaGen,GetBetaGen);
		Filler.AddObject2beFilled(Preselections_D,GetBetaGen,GetBetaGen);
		Filler.AddObject2beFilled(Quality_D	 ,GetBetaGen,GetBetaGen);
		Filler.AddObject2beFilled(RICH_D	 ,GetBetaGen,GetBetaGen);
		Filler.AddObject2beFilled(RICH_P	 ,GetBetaGen,GetBetaGen);
		Filler.ReinitializeAll(Refill);
		//main loop
		Filler.LoopOnMC(DBarReader(chainMC, true ),vars);

	FullSet_P 	->CloneEfficiency(Preselections_P);
	FullSet_D       ->CloneEfficiency(Preselections_D);


//	RigBinFullSetEff->Save(finalHistos);	
	Preselections_P	->Save(finalHistos);
        Quality_P       ->Save(finalHistos);
	Preselections_D ->Save(finalHistos);
	Quality_D       ->Save(finalHistos);
	RICH_P 	->Save(finalHistos);
	RICH_D       ->Save(finalHistos);
	FullSet_P 	->Save(finalHistos);
	FullSet_D       ->Save(finalHistos);


//	RigBinFullSetEff->Eval_Efficiency();
	Preselections_P	->Eval_Efficiency();
        Quality_P       ->Eval_Efficiency();
        Preselections_D ->Eval_Efficiency();
        Quality_D       ->Eval_Efficiency();
	RICH_P	        ->Eval_Efficiency();
        RICH_D  	->Eval_Efficiency();
	FullSet_P       ->Eval_Efficiency();
        FullSet_D       ->Eval_Efficiency();
	FullSet_P 	->ComposeEfficiency(Quality_P);
	FullSet_D       ->ComposeEfficiency(Quality_D);
	FullSet_P       ->Eval_FittedEfficiency();
        FullSet_D       ->Eval_FittedEfficiency();

//	RigBinFullSetEff->SaveResults(finalResults);
	Preselections_P	->SaveResults(finalResults);
	Quality_P       ->SaveResults(finalResults);
	Preselections_D ->SaveResults(finalResults);
	Quality_D       ->SaveResults(finalResults);
	RICH_P 		->SaveResults(finalResults);
	RICH_D 		->SaveResults(finalResults);
	FullSet_P 	->SaveResults(finalResults);
	FullSet_D 	->SaveResults(finalResults);

	return 0;
}


	

