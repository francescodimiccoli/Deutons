#ifndef COUNTEXTRACTION_CPP
#define COUNTEXTRACTION_CPP


#include "Analyzer.h"
#include "Variables.hpp"
#include "TemplateFITbetasmear.h"
#include "Efficiency.h"

void Analyzer::BookCountsAnalysis(FileSaver finalhistos, FileSaver finalresults, bool refill)
{

	bool checkfile = finalhistos.CheckFile();
	check_file = checkfile;

	FileSaver LatWeights;
	LatWeights.setName("/data1/home/fdimicco/Deutons/DirectAnalysis/LatWeights/Weights.root");
	LatReweighter * weighter = new LatReweighter(LatWeights,"LatWeights");	
	
	cout<<"****************************** BINS ***************************************"<<endl;
    	SetUpUsualBinning();
 
	cout<<"****************************** Counts ANALYIS ******************************************"<<endl;

	BadEventSimulator * NaFBadEvSimulator= new BadEventSimulator("IsFromNaF",22,0.72,1); 
	BadEventSimulator * AglBadEvSimulator= new BadEventSimulator("IsFromAgl",250,0.95,1); 

	//proton flux tests
	Efficiency * CountsTests[9];

	for(int i=0;i<10;i++) {
		CountsTests[i] = new Efficiency(finalhistos,("HEPTests_v"+to_string(i)).c_str(),("HEPTests_v"+to_string(i)).c_str(),PRB,("IsPositive&IsPrimary&IsStandardSel_v"+to_string(i)).c_str(),("IsPositive&IsPrimary&IsStandardSel_v"+to_string(i)).c_str());
	}

	//simple event count
	Efficiency * CountsHE     = new Efficiency(finalhistos    ,"HEPCounts"         ,"HEPCounts"      ,PRB,"IsPositive&IsPrimary&IsBaseline","IsPositive&IsPrimary&IsBaseline");
	Efficiency * CountsL1HE   = new Efficiency(finalhistos    ,"HEPCountsL1"         ,"HEPCountsL1"  ,PRB,"IsPositive&IsPrimary&IsBaseline&L1LooseCharge1","IsPositive&IsPrimary&IsBaseline&L1LooseCharge1");
	Efficiency * CountsQualHE = new Efficiency(finalhistos,"HEPCountsQual"     ,"HEPCountsQual"      ,PRB,"IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning","IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning");

//"IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&IsOnlyFromToF"

	Efficiency * CountsTOF    = new Efficiency(finalhistos    ,"TOFPCounts"	   ,"TOFPCounts"	,ToFRigB,"IsPositive&IsBaseline&L1LooseCharge1&IsCleaning","IsPositive&IsBaseline&L1LooseCharge1&IsCleaning");
	Efficiency * CountsNaF    = new Efficiency(finalhistos    ,"NaFPCounts"	   ,"NaFPCounts"	,NaFRigB,"IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut"    ,"IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut"    );
	Efficiency * CountsAgl    = new Efficiency(finalhistos    ,"AglPCounts"	   ,"AglPCounts"	,AglRigB,"IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut"    ,"IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut"    );

	// Extraction of counts with Template Fit

	//TemplateFIT * SmearingCheck = new TemplateFIT("SmearingCheck",PRB,"IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsOnlyFromToF",60,0.3,1.6);	  
	TemplateFIT * TOFfits= new TemplateFIT("TOFfits",ToFDB,"IsPositive&IsBaseline&L1LooseCharge1&IsCleaning"       ,150,0.4,7.5,11);
	TemplateFIT * NaFfits= new TemplateFIT("NaFfits",NaFDB,"IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut"           ,60,0.4,5,true,11,400,200);
	TemplateFIT * Aglfits= new TemplateFIT("Aglfits",AglDB,"IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut"           ,60,0.4,5,true,11,110,80);	

	if(!refill&&checkfile) {	

		CountsHE    ->ReinitializeHistos(refill); 
	       	CountsL1HE    ->ReinitializeHistos(refill); 
	        CountsQualHE->ReinitializeHistos(refill); 
	        CountsTOF   ->ReinitializeHistos(refill); 
	        CountsNaF   ->ReinitializeHistos(refill); 
                CountsAgl   ->ReinitializeHistos(refill); 
		for(int i=0;i<10;i++)  CountsTests[i] ->ReinitializeHistos(refill);

		//TemplateFIT * SmearingCheck = new TemplateFIT(finalhistos,"SmearingCheck",PRB);
		TOFfits= new TemplateFIT(finalhistos,"TOFfits",ToFDB,false,11);
		NaFfits= new TemplateFIT(finalhistos,"NaFfits",NaFDB,true,11,400,200);
		Aglfits= new TemplateFIT(finalhistos,"Aglfits",AglDB,true,11,110,90);
		//	NaFfits->SetFitWithNoiseMode();
		//	Aglfits->SetFitWithNoiseMode();


		//SmearingCheck->SetFitRangeByQuantiles(0.05,0.95);
		//SmearingCheck->SetFitConstraints(0.99,1,0.0001,0.001,0.01,0.001);
		//SmearingCheck->DisableFit();
		//SmearingCheck->ExtractCounts(finalhistos);	
		//SmearingCheck->SaveFitResults(finalresults);

		TOFfits->SetFitRange(0.6,4);
//		TOFfits->DisableFit();
		TOFfits->SetFitConstraints(0.9,1,0.015,0.16,0.01,0.02,true);
		TOFfits->ExtractCounts(finalhistos);	
		TOFfits->SaveFitResults(finalresults);

		NaFfits->SetFitRange(0.6,5);
//		NaFfits->DisableFit();
		NaFfits->SetFitConstraints(0.9,1,0.001,0.1,0.0005,0.02,true);
		NaFfits->ExtractCounts(finalhistos);
		NaFfits->SaveFitResults(finalresults);

		Aglfits->SetFitRange(0.6,5);
//		Aglfits->DisableFit();
		Aglfits->SetFitConstraints(0.9,1,0.001,0.1,0.0005,0.02,true);
		Aglfits->ExtractCounts(finalhistos);
		Aglfits->SaveFitResults(finalresults);	

		for(int i=0;i<10;i++)  CountsTests[i] ->Eval_Efficiency();
		CountsHE	->Eval_Efficiency();	
		CountsL1HE	->Eval_Efficiency();	
		CountsQualHE	->Eval_Efficiency();
		CountsTOF	->Eval_Efficiency();
		CountsNaF	->Eval_Efficiency();
		CountsAgl	->Eval_Efficiency();

		for(int i=0;i<10;i++)  CountsTests[i] ->SaveResults(finalresults);
		CountsHE	->SaveResults(finalresults);   
		CountsL1HE	->SaveResults(finalresults);   
		CountsQualHE	->SaveResults(finalresults);
		CountsTOF	->SaveResults(finalresults);
		CountsNaF	->SaveResults(finalresults);
		CountsAgl	->SaveResults(finalresults);

	}


	//SmearingCheck -> SetDefaultOutFile(finalhistos);
	TOFfits	    -> SetDefaultOutFile(finalhistos);
	NaFfits	    -> SetDefaultOutFile(finalhistos);
	Aglfits	    -> SetDefaultOutFile(finalhistos);

	for(int i=0;i<10;i++)  CountsTests[i] ->SetDefaultOutFile(finalhistos);
	CountsHE	-> SetDefaultOutFile(finalhistos);
	CountsL1HE	-> SetDefaultOutFile(finalhistos);
	CountsQualHE-> SetDefaultOutFile(finalhistos);
	CountsTOF	-> SetDefaultOutFile(finalhistos);
	CountsNaF	-> SetDefaultOutFile(finalhistos);
	CountsAgl	-> SetDefaultOutFile(finalhistos);


	//SmearingCheck->SetLatitudeReweighter(weighter);
	TOFfits	 ->SetLatitudeReweighter(weighter);	
	NaFfits	 ->SetLatitudeReweighter(weighter);	
	Aglfits	 ->SetLatitudeReweighter(weighter);	

	NaFfits->SetUpBadEventSimulator(NaFBadEvSimulator);
	Aglfits->SetUpBadEventSimulator(AglBadEvSimulator);
	NaFfits->SetFitWithNoiseMode();
	Aglfits->SetFitWithNoiseMode();

	//Filler.AddObject2beFilled(SmearingCheck,GetBetaTOF,GetRigidity);
	Filler.AddObject2beFilled(TOFfits,GetRecMassTOF ,GetBetaTOF);
	Filler.AddObject2beFilled(NaFfits,GetRecMassRICH,GetBetaRICH);
	Filler.AddObject2beFilled(Aglfits,GetRecMassRICH,GetBetaRICH);
	Filler.AddObject2beFilled(CountsHE,GetRigidity,GetRigidity); 
	Filler.AddObject2beFilled(CountsL1HE,GetRigidity,GetRigidity); 
	Filler.AddObject2beFilled(CountsQualHE,GetRigidity,GetRigidity); 
	Filler.AddObject2beFilled(CountsTOF,GetRigidity,GetRigidity); 
	Filler.AddObject2beFilled(CountsNaF,GetRigidity,GetRigidity); 
	Filler.AddObject2beFilled(CountsAgl,GetRigidity,GetRigidity); 
	

	for(int i=0;i<10;i++)  Filler.AddObject2beFilled(CountsTests[i],GetRigidity,GetRigidity);

	Filler.ReinitializeAll(refill);


}

#endif

