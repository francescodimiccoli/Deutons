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

	TH1::SetDefaultSumw2();
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


	Efficiency * RigBinFullSetEff = new Efficiency(finalHistos,"RigBinFullSetEff","RigBinFullSetEff",PRB,"IsProtonMC","IsProtonMC&IsPreselected&DistanceCut&LikelihoodCut");

        AllRangesEfficiency * Trigger_P = new AllRangesEfficiency(finalHistos,"Trigger_P","Trigger",
        "IsProtonMC",
        "IsProtonMC",
        "IsProtonMC",
        "IsProtonMC&IsPhysTrig",
        "IsProtonMC&IsPhysTrig",
        "IsProtonMC&IsPhysTrig",
        Refill);

        AllRangesEfficiency * Trigger_D = new AllRangesEfficiency(finalHistos,"Trigger_D","Trigger",
        "IsDeutonMC",
        "IsDeutonMC",
        "IsDeutonMC",
        "IsDeutonMC&IsPhysTrig",
        "IsDeutonMC&IsPhysTrig",
        "IsDeutonMC&IsPhysTrig",
        Refill);

	AllRangesEfficiency * Fragmentation_P = new AllRangesEfficiency(finalHistos,"Fragmentation_P","Fragmentation",
	"IsProtonMC&IsPhysTrig",			      
	"IsProtonMC&IsPhysTrig",			   	      
	"IsProtonMC&IsPhysTrig",
	"IsPurePMC",			      
	"IsPurePMC",			   	      
	"IsPurePMC",
	Refill);

	AllRangesEfficiency * Fragmentation_D = new AllRangesEfficiency(finalHistos,"Fragmentation_D","Fragmentation",
	"IsDeutonMC&IsPhysTrig",			      
	"IsDeutonMC&IsPhysTrig",			   	      
	"IsDeutonMC&IsPhysTrig",
	"IsPureDMC",			      
	"IsPureDMC",			   	      
	"IsPureDMC",
	Refill);

	AllRangesEfficiency * Basic_P = new AllRangesEfficiency(finalHistos,"BasicEff_P","BasicEfficiency",
	"IsPurePMC",			      
	"IsPurePMC",			   	      
	"IsPurePMC",
	"IsPurePMC&IsPositive&IsBasic",
	"IsPurePMC&IsPositive&IsBasic",
	"IsPurePMC&IsPositive&IsBasic",
	Refill);

	AllRangesEfficiency * Preselections_P = new AllRangesEfficiency(finalHistos,"PresEff_P","PreselectionEfficiency",
	"IsPurePMC&IsPositive&IsBasic",			      
	"IsPurePMC&IsPositive&IsBasic",			   	      
	"IsPurePMC&IsPositive&IsBasic",
	"IsPurePMC&IsPositive&IsPreselected",
	"IsPurePMC&IsPositive&IsPreselected&IsFromNaF_nosel",
	"IsPurePMC&IsPositive&IsPreselected&IsFromAgl_nosel",
	Refill);

	AllRangesEfficiency * Quality_P       = new AllRangesEfficiency(finalHistos,"QualEff_P","QualityEfficiency",
	"IsPurePMC&IsPositive&IsPreselected",
	"IsPurePMC&IsPositive&IsPreselected&IsFromNaF_nosel",
	"IsPurePMC&IsPositive&IsPreselected&IsFromAgl_nosel",
	"IsPurePMC&IsPositive&IsPreselected&DistanceCut&LikelihoodCut",
	"IsPurePMC&IsPositive&IsPreselected&IsFromNaF&DistanceCut&LikelihoodCut&RICHBDTCut",
	"IsPurePMC&IsPositive&IsPreselected&IsFromAgl&DistanceCut&LikelihoodCut&RICHBDTCut"
	,Refill);
        
	AllRangesEfficiency * FullSet_P      = new AllRangesEfficiency(finalHistos,"FullsetEff_P","FullsetEfficiency",
	"",
	"",
	Refill);

	AllRangesEfficiency * Basic_D = new AllRangesEfficiency(finalHistos,"BasicEff_D","BasicEfficiency",
	"IsPureDMC",			      
	"IsPureDMC",			   	      
	"IsPureDMC",
	"IsPureDMC&IsPositive&IsBasic",
	"IsPureDMC&IsPositive&IsBasic",
	"IsPureDMC&IsPositive&IsBasic",
	Refill);

	AllRangesEfficiency * Preselections_D = new AllRangesEfficiency(finalHistos,"PresEff_D","PreselectionEfficiency",
	"IsPureDMC&IsPositive&IsBasic",			      
	"IsPureDMC&IsPositive&IsBasic",
	"IsPureDMC&IsPositive&IsBasic",
	"IsPureDMC&IsPositive&IsPreselected",
	"IsPureDMC&IsPositive&IsPreselected&IsFromNaF_nosel",
	"IsPureDMC&IsPositive&IsPreselected&IsFromAgl_nosel",
	Refill);

	AllRangesEfficiency * Quality_D       = new AllRangesEfficiency(finalHistos,"QualEff_D","QualityEfficiency",
	"IsPureDMC&IsPositive&IsPreselected",
	"IsPureDMC&IsPositive&IsPreselected&IsFromNaF_nosel",
	"IsPureDMC&IsPositive&IsPreselected&IsFromAgl_nosel",
	"IsPureDMC&IsPositive&IsPreselected&DistanceCut&LikelihoodCut",
	"IsPureDMC&IsPositive&IsPreselected&IsFromNaF&DistanceCut&LikelihoodCut&RICHBDTCut",
	"IsPureDMC&IsPositive&IsPreselected&IsFromAgl&DistanceCut&LikelihoodCut&RICHBDTCut",
	Refill);
        
	AllRangesEfficiency * FullSet_D       = new AllRangesEfficiency(finalHistos,"FullsetEff_D","FullsetEfficiency"  ,
	"",
	"",
	Refill);

	AllRangesEfficiency * RICH_P = new AllRangesEfficiency(finalHistos,"RICHEff_P","RICHEfficiency",
	"IsPurePMC&IsPositive&IsPreselected",
	"IsPurePMC&IsPositive&IsPreselected&IsFromNaF_nosel",	     
	"IsPurePMC&IsPositive&IsPreselected&IsFromAgl_nosel",
	"IsPurePMC&IsPositive&IsPreselected",
	"IsPurePMC&IsPositive&IsPreselected&IsFromNaF",
	"IsPurePMC&IsPositive&IsPreselected&IsFromAgl",
	Refill);
	
	AllRangesEfficiency * RICH_D = new AllRangesEfficiency(finalHistos,"RICHEff_D","RICHEfficiency",
	"IsPureDM&IsPositive&IsPreselected",
	"IsPureDMC&IsPositive&IsPreselected&IsFromNaF_nosel",
	"IsPureDMC&IsPositive&IsPreselected&IsFromAgl_nosel",
	"IsPureDMC&IsPositive&IsPreselected",
	"IsPureDMC&IsPositive&IsPreselected&IsFromNaF",
	"IsPureDMC&IsPositive&IsPreselected&IsFromAgl",
	Refill);

	AllRangesEfficiency * RICH_PQual = new AllRangesEfficiency(finalHistos,"RICHEff_PQual","RICHQualEfficiency",
	"IsPurePMC&IsPositive&IsPreselected&LikelihoodCut",
	"IsPurePMC&IsPositive&IsPreselected&LikelihoodCut&IsFromNaF",
	"IsPurePMC&IsPositive&IsPreselected&LikelihoodCut&IsFromAgl",
	"IsPurePMC&IsPositive&IsPreselected&LikelihoodCut",
	"IsPurePMC&IsPositive&IsPreselected&LikelihoodCut&IsFromNaF&RICHBDTCut",
	"IsPurePMC&IsPositive&IsPreselected&LikelihoodCut&IsFromAgl&RICHBDTCut",
	Refill);
	
	AllRangesEfficiency * RICH_DQual = new AllRangesEfficiency(finalHistos,"RICHEff_DQual","RICHQualEfficiency",
	"IsPureDMC&IsPositive&IsPreselected&LikelihoodCut",
	"IsPureDMC&IsPositive&IsPreselected&LikelihoodCut&IsFromNaF",
	"IsPureDMC&IsPositive&IsPreselected&LikelihoodCut&IsFromAgl",
	"IsPureDMC&IsPositive&IsPreselected&LikelihoodCut",
	"IsPureDMC&IsPositive&IsPreselected&LikelihoodCut&IsFromNaF&RICHBDTCut",
	"IsPureDMC&IsPositive&IsPreselected&LikelihoodCut&IsFromAgl&RICHBDTCut",
	Refill);

	AllRangesEfficiency * FullSetTOT_P       = new AllRangesEfficiency(finalHistos,"FullsetTOTEff_P","FullsetTOTEfficiency"  ,
	"",
	"",
	Refill);

	AllRangesEfficiency * FullSetTOT_D       = new AllRangesEfficiency(finalHistos,"FullsetTOTEff_D","FullsetTOTEfficiency"  ,
	"",
	"",
	Refill);

	//efficiency cascade
	AllRangesEfficiency * CascadeTrigg = new AllRangesEfficiency(finalHistos,"CascadeTrigg","CascadeTrigg",
        "",
        "",
	Refill);

 	AllRangesEfficiency * CascadeFragm = new AllRangesEfficiency(finalHistos,"CascadeFragm","CascadeFragm",
        "",
        "",
	Refill);

 	AllRangesEfficiency * CascadeBasic = new AllRangesEfficiency(finalHistos,"CascadeBasic","CascadeBasic",
        "",
        "",
        Refill);

 	AllRangesEfficiency * CascadeGolden = new AllRangesEfficiency(finalHistos,"CascadeGolden","CascadeGolden",
        "",
        "",
	Refill);
	
	AllRangesEfficiency * CascadeGoldenRICH = new AllRangesEfficiency(finalHistos,"CascadeGoldenRICH","CascadeGoldenRICH",
        "",
        "",
	Refill);
	
	AllRangesEfficiency * CascadeQuality = new AllRangesEfficiency(finalHistos,"CascadeQuality","CascadeQuality",
        "",
        "",
	Refill);


	//RigBinFullSetEff->Fill(DBarReader(chainMC, true ),vars,GetGenMomentum,Refill);

	ParallelFiller<AllRangesEfficiency *> Filler;
	Filler.AddObject2beFilled(Trigger_P,GetBetaGen,GetBetaGen);
        Filler.AddObject2beFilled(Trigger_D,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(Fragmentation_P,GetBetaGen,GetBetaGen);
        Filler.AddObject2beFilled(Fragmentation_D,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(Basic_P,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(Preselections_P,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(Quality_P	 ,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(Basic_D,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(Preselections_D,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(Quality_D	 ,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(RICH_D	 ,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(RICH_P	 ,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(RICH_DQual	 ,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(RICH_PQual	 ,GetBetaGen,GetBetaGen);
	Filler.ReinitializeAll(Refill);
	//main loop
	Filler.LoopOnMC(DBarReader(chainMC, true ),vars);

	FullSetTOT_P 	->CloneEfficiency(Trigger_P);
	FullSetTOT_D    ->CloneEfficiency(Trigger_D);

	FullSet_P 	->CloneEfficiency(Preselections_P);
	FullSet_D       ->CloneEfficiency(Preselections_D);

	Trigger_D ->Save(finalHistos);	
	Trigger_P ->Save(finalHistos);
	Fragmentation_D ->Save(finalHistos);	
	Fragmentation_P	->Save(finalHistos);
	RigBinFullSetEff->Save(finalHistos);	
	Basic_P		->Save(finalHistos);
	Basic_D         ->Save(finalHistos);
	Preselections_P	->Save(finalHistos);
        Quality_P       ->Save(finalHistos);
	Preselections_D ->Save(finalHistos);
	Quality_D       ->Save(finalHistos);
	RICH_P 		->Save(finalHistos);
	RICH_D       	->Save(finalHistos);
	RICH_PQual 	->Save(finalHistos);
	RICH_DQual      ->Save(finalHistos);
	FullSet_P 	->Save(finalHistos);
	FullSet_D       ->Save(finalHistos);

	Trigger_D ->Eval_Efficiency();	
	Trigger_P	->Eval_Efficiency();
	Fragmentation_D ->Eval_Efficiency();	
	Fragmentation_P	->Eval_Efficiency();
	RigBinFullSetEff->Eval_Efficiency();
	Basic_D         ->Eval_Efficiency();
	Basic_P         ->Eval_Efficiency();
	Preselections_P	->Eval_Efficiency();
        Quality_P       ->Eval_Efficiency();
        Preselections_D ->Eval_Efficiency();
        Quality_D       ->Eval_Efficiency();
	RICH_P	        ->Eval_Efficiency();
        RICH_D  	->Eval_Efficiency();
	RICH_PQual      ->Eval_Efficiency();
        RICH_DQual  	->Eval_Efficiency();
	
	FullSet_P       ->Eval_Efficiency();
        FullSet_D       ->Eval_Efficiency();
	FullSet_P 	->ComposeEfficiency(Quality_P);
	FullSet_D       ->ComposeEfficiency(Quality_D);
	FullSet_P       ->Eval_FittedEfficiency();
        FullSet_D       ->Eval_FittedEfficiency();

	
	FullSetTOT_P      ->Eval_Efficiency();
        FullSetTOT_D      ->Eval_Efficiency();
	FullSetTOT_P      ->ComposeEfficiency(Fragmentation_P);
	FullSetTOT_D      ->ComposeEfficiency(Fragmentation_D);
	FullSetTOT_P      ->ComposeEfficiency(FullSet_P);
	FullSetTOT_D      ->ComposeEfficiency(FullSet_D);


	//cascade
	CascadeTrigg		->CloneEfficiency(Trigger_D);	
	CascadeFragm		->CloneEfficiency(Trigger_D);	
	CascadeBasic		->CloneEfficiency(Trigger_D);	
	CascadeGolden		->CloneEfficiency(Trigger_D);	
	CascadeGoldenRICH	->CloneEfficiency(Trigger_D);	
	CascadeQuality		->CloneEfficiency(Trigger_D);	

	CascadeTrigg		->Eval_Efficiency();	
	CascadeFragm		->Eval_Efficiency();	
	CascadeBasic		->Eval_Efficiency();	
	CascadeGolden		->Eval_Efficiency();	
	CascadeGoldenRICH	->Eval_Efficiency();	
	CascadeQuality		->Eval_Efficiency();	

	CascadeFragm		->ComposeEfficiency(Fragmentation_D);	
	CascadeBasic		->ComposeEfficiency(Fragmentation_D);	
	CascadeGolden		->ComposeEfficiency(Fragmentation_D);	
	CascadeGoldenRICH	->ComposeEfficiency(Fragmentation_D);	
	CascadeQuality		->ComposeEfficiency(Fragmentation_D);	

	CascadeBasic		->ComposeEfficiency(Basic_D);	
	CascadeGolden		->ComposeEfficiency(Basic_D);	
	CascadeGoldenRICH	->ComposeEfficiency(Basic_D);	
	CascadeQuality		->ComposeEfficiency(Basic_D);	

	CascadeGolden		->ComposeEfficiency(Preselections_D);	
	CascadeGoldenRICH	->ComposeEfficiency(Preselections_D);	
	CascadeQuality		->ComposeEfficiency(Preselections_D);	

	CascadeGoldenRICH	->ComposeEfficiency(RICH_D);	

	CascadeQuality		->ComposeEfficiency(Quality_D);	


	Trigger_D	->SaveResults(finalResults);	
	Trigger_P	->SaveResults(finalResults);	
	Fragmentation_D ->SaveResults(finalResults);	
	Fragmentation_P	->SaveResults(finalResults);	
	RigBinFullSetEff->SaveResults(finalResults);
	Basic_P		->SaveResults(finalResults);
	Basic_D		->SaveResults(finalResults);
	Preselections_P	->SaveResults(finalResults);
	Quality_P       ->SaveResults(finalResults);
	Preselections_D ->SaveResults(finalResults);
	Quality_D       ->SaveResults(finalResults);
	RICH_P 		->SaveResults(finalResults);
	RICH_D 		->SaveResults(finalResults);
	RICH_PQual 	->SaveResults(finalResults);
	RICH_DQual      ->SaveResults(finalResults);
	
	FullSet_P 	->SaveResults(finalResults);
	FullSet_D 	->SaveResults(finalResults);
	FullSetTOT_P 	->SaveResults(finalResults);
	FullSetTOT_D 	->SaveResults(finalResults);

	CascadeTrigg		->SaveResults(finalResults);		
	CascadeFragm		->SaveResults(finalResults);
	CascadeBasic		->SaveResults(finalResults);
	CascadeGolden		->SaveResults(finalResults);
	CascadeGoldenRICH	->SaveResults(finalResults);
	CascadeQuality		->SaveResults(finalResults);





	return 0;
}


	

