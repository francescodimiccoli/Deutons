#ifndef COUNTEXTRACTION_CPP
#define COUNTEXTRACTION_CPP


#include "Analyzer.h"
#include "Variables.hpp"
#include "TemplateFITbetasmear.h"
#include "Efficiency.h"

void Analyzer::BookCountsAnalysis(FileSaver finalhistos, FileSaver finalresults, bool refill)
{

	bool checkfile = finalhistos.CheckFile();

	FileSaver LatWeights;
	LatWeights.setName("LatWeights/Weights.root");
	LatReweighter * weighter = new LatReweighter(LatWeights,"LatWeights");	
	
	cout<<"****************************** BINS ***************************************"<<endl;
    	SetUpEffCorrBinning();
 
	cout<<"****************************** Counts ANALYIS ******************************************"<<endl;

	BadEventSimulator * NaFBadEvSimulator= new BadEventSimulator("IsFromNaF",22,0.72,1); 
	BadEventSimulator * AglBadEvSimulator= new BadEventSimulator("IsFromAgl",250,0.95,1); 


	//simple event count
	Efficiency * CountsHE     = new Efficiency(finalhistos    ,"HEPCounts"         ,"HEPCounts"      ,PRB,"IsPositive&IsPrimary&IsMinimumBias&IsLooseCharge1","IsPositive&IsPrimary&IsMinimumBias&IsLooseCharge1");
	Efficiency * CountsQualHE = new Efficiency(finalhistos,"HEPCountsQual"     ,"HEPCountsQual"  ,PRB,"IsPositive&IsPrimary&IsMinimumBias&IsLooseCharge1&IsCleaning","IsPositive&IsPrimary&IsMinimumBias&IsLooseCharge1&IsCleaning");
	Efficiency * CountsTOF    = new Efficiency(finalhistos    ,"TOFPCounts"	   ,"TOFPCounts"	,ToFPB,"IsPositive&IsPrimary&IsMinimumBias&IsLooseCharge1&IsCleaning&IsGoodTime&IsOnlyFromToF","IsPositive&IsPrimary&IsMinimumBias&IsLooseCharge1&IsCleaning&IsGoodTime&IsOnlyFromToF");
	Efficiency * CountsNaF    = new Efficiency(finalhistos    ,"NaFPCounts"	   ,"NaFPCounts"	,NaFPB,"IsPositive&IsPrimary&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut"    ,"IsPositive&IsPrimary&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut"    );
	Efficiency * CountsAgl    = new Efficiency(finalhistos    ,"AglPCounts"	   ,"AglPCounts"	,AglPB,"IsPositive&IsPrimary&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut"    ,"IsPositive&IsPrimary&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut"    );

	// Extraction of counts with Template Fit

	//  TemplateFIT * SmearingCheck = new TemplateFIT("SmearingCheck",PRB,"IsPositive&IsPreselected&LikelihoodCut&DistanceCut&IsOnlyFromToF",60,0.3,1.6);	  
	TemplateFIT * TOFfits= new TemplateFIT("TOFfits",ToFDB,"IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsGoodTime&IsOnlyFromToF"       ,150,0.4,7.5,7);
	TemplateFIT * NaFfits= new TemplateFIT("NaFfits",NaFDB,"IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut"           ,60,0.4,5,true,7,400,200);
	TemplateFIT * Aglfits= new TemplateFIT("Aglfits",AglDB,"IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut"           ,60,0.4,5,true,7,110,80);	

	if(!refill&&checkfile) {	

		CountsHE    ->ReinitializeHistos(refill); 
	        CountsQualHE->ReinitializeHistos(refill); 
                CountsTOF   ->ReinitializeHistos(refill); 
	        CountsNaF   ->ReinitializeHistos(refill); 
                CountsAgl   ->ReinitializeHistos(refill); 


		//TemplateFIT * SmearingCheck = new TemplateFIT(finalhistos,"SmearingCheck",PRB);
		TOFfits= new TemplateFIT(finalhistos,"TOFfits",ToFDB,false,7);
		NaFfits= new TemplateFIT(finalhistos,"NaFfits",NaFDB,true,7,400,200);
		Aglfits= new TemplateFIT(finalhistos,"Aglfits",AglDB,true,7,110,80);
		//	NaFfits->SetFitWithNoiseMode();
		//	Aglfits->SetFitWithNoiseMode();


		//SmearingCheck->SetFitRangeByQuantiles(0.05,0.95);
		//SmearingCheck->SetFitConstraints(0.99,1,0.0001,0.001,0.01,0.001);
		//SmearingCheck->DisableFit();
		//SmearingCheck->ExtractCounts(finalhistos);	
		//SmearingCheck->SaveFitResults(finalresults);

		//TOFfits->DisableFit();
		TOFfits->SetFitRange(0.6,4);
		TOFfits->SetFitConstraints(0.9,1,0.015,0.06,0.005,0.015,true);
		TOFfits->ExtractCounts(finalhistos);	
		TOFfits->SaveFitResults(finalresults);

		NaFfits->SetFitRange(0.6,5);
		//NaFfits->DisableFit();
		NaFfits->SetFitConstraints(0.9,1,0.001,0.1,0.0001,0.0005);
		NaFfits->ExtractCounts(finalhistos);
		NaFfits->SaveFitResults(finalresults);

		Aglfits->SetFitRange(0.6,5);
		//Aglfits->DisableFit();
		Aglfits->SetFitConstraints(0.9,1,0.001,0.1,0.0001,0.0005);
		Aglfits->ExtractCounts(finalhistos);
		Aglfits->SaveFitResults(finalresults);	


		CountsHE	->Eval_Efficiency();	
		CountsQualHE	->Eval_Efficiency();
		CountsTOF	->Eval_Efficiency();
		CountsNaF	->Eval_Efficiency();
		CountsAgl	->Eval_Efficiency();

		CountsHE	->SaveResults(finalresults);   
		CountsQualHE	->SaveResults(finalresults);
		CountsTOF	->SaveResults(finalresults);
		CountsNaF	->SaveResults(finalresults);
		CountsAgl	->SaveResults(finalresults);

	}


	//SmearingCheck -> SetDefaultOutFile(finalhistos);
	TOFfits	    -> SetDefaultOutFile(finalhistos);
	NaFfits	    -> SetDefaultOutFile(finalhistos);
	Aglfits	    -> SetDefaultOutFile(finalhistos);

	CountsHE	-> SetDefaultOutFile(finalhistos);
	CountsQualHE-> SetDefaultOutFile(finalhistos);
	CountsTOF	-> SetDefaultOutFile(finalhistos);
	CountsNaF	-> SetDefaultOutFile(finalhistos);
	CountsAgl	-> SetDefaultOutFile(finalhistos);


	// SmearingCheck->SetLatitudeReweighter(weighter);
	TOFfits	 ->SetLatitudeReweighter(weighter);	
	NaFfits	 ->SetLatitudeReweighter(weighter);	
	Aglfits	 ->SetLatitudeReweighter(weighter);	

	NaFfits->SetUpBadEventSimulator(NaFBadEvSimulator);
	Aglfits->SetUpBadEventSimulator(AglBadEvSimulator);
	NaFfits->SetFitWithNoiseMode();
	Aglfits->SetFitWithNoiseMode();

	//  Filler.AddObject2beFilled(SmearingCheck,GetBetaTOF,GetRigidity);
	Filler.AddObject2beFilled(TOFfits,GetRecMassTOF ,GetBetaTOF);
	Filler.AddObject2beFilled(NaFfits,GetRecMassRICH,GetBetaRICH);
	Filler.AddObject2beFilled(Aglfits,GetRecMassRICH,GetBetaRICH);
	Filler.AddObject2beFilled(CountsHE,GetRigidity,GetRigidity); 
	Filler.AddObject2beFilled(CountsQualHE,GetRigidity,GetRigidity); 
	Filler.AddObject2beFilled(CountsTOF,GetRigidity,GetRigidity); 
	Filler.AddObject2beFilled(CountsNaF,GetRigidity,GetRigidity); 
	Filler.AddObject2beFilled(CountsAgl,GetRigidity,GetRigidity); 
	Filler.ReinitializeAll(refill);


}

#endif

