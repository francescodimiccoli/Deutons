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

	Efficiency * CountsTOF    = new Efficiency(finalhistos    ,"TOFPCounts"	   ,"TOFPCounts"	,ToFRigB,"IsPositive&IsBaseline&L1LooseCharge1&IsCleaning","IsPositive&IsBaseline&L1LooseCharge1&IsCleaning");
	Efficiency * CountsNaF    = new Efficiency(finalhistos    ,"NaFPCounts"	   ,"NaFPCounts"	,NaFRigB,"IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut"    ,"IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut"    );
	Efficiency * CountsAgl    = new Efficiency(finalhistos    ,"AglPCounts"	   ,"AglPCounts"	,AglRigB,"IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut"    ,"IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut"    );

	// Extraction of counts with Template Fit

	//TemplateFIT * SmearingCheck = new TemplateFIT("SmearingCheck",PRB,"IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsOnlyFromToF",60,0.3,1.6);	  
	TemplateFIT * TOFfits= new TemplateFIT("TOFDfits",ToFDB,"IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsOnlyFromToF"       ,150,0.4,7.5,false,11);
	TemplateFIT * NaFfits= new TemplateFIT("NaFDfits",NaFDB,"IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut"           ,60,0.4,5,true,11,250,150);
	TemplateFIT * Aglfits= new TemplateFIT("AglDfits",AglDB,"IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut"           ,60,0.4,5,true,11,75,55);	

	TemplateFIT * TOFfits_P= new TemplateFIT("TOFPfits",ToFPB,"IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsOnlyFromToF"       ,150,0.4,7.5,false,11);
	TemplateFIT * NaFfits_P= new TemplateFIT("NaFPfits",NaFPB,"IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut"           ,60,0.4,5,true,11,250,150);
	TemplateFIT * Aglfits_P= new TemplateFIT("AglPfits",AglPB,"IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut"           ,60,0.4,5,true,11,75,55);	


	if(!refill&&checkfile) {	

		CountsHE    ->ReinitializeHistos(refill); 
	       	CountsL1HE    ->ReinitializeHistos(refill); 
	        CountsQualHE->ReinitializeHistos(refill); 
	        CountsTOF   ->ReinitializeHistos(refill); 
	        CountsNaF   ->ReinitializeHistos(refill); 
                CountsAgl   ->ReinitializeHistos(refill); 
		for(int i=0;i<10;i++)  CountsTests[i] ->ReinitializeHistos(refill);

		//TemplateFIT * SmearingCheck = new TemplateFIT(finalhistos,"SmearingCheck",PRB);
		TOFfits= new TemplateFIT(finalhistos,"TOFDfits",ToFDB,false,11);
		NaFfits= new TemplateFIT(finalhistos,"NaFDfits",NaFDB,true,11,250,150);
		Aglfits= new TemplateFIT(finalhistos,"AglDfits",AglDB,true,11,75,55);

		TOFfits_P= new TemplateFIT(finalhistos,"TOFPfits",ToFPB,false,11);
		NaFfits_P= new TemplateFIT(finalhistos,"NaFPfits",NaFPB,true,11,250,150);
		Aglfits_P= new TemplateFIT(finalhistos,"AglPfits",AglPB,true,11,75,55);


		//	NaFfits->SetFitWithNoiseMode();
		//	Aglfits->SetFitWithNoiseMode();


		//SmearingCheck->SetFitRangeByQuantiles(0.05,0.95);
		//SmearingCheck->SetFitConstraints(0.99,1,0.0001,0.001,0.01,0.001);
		//SmearingCheck->DisableFit();
		//SmearingCheck->ExtractCounts(finalhistos);	
		//SmearingCheck->SaveFitResults(finalresults);

		TOFfits->SetFitRange(0.6,4);
//		TOFfits->DisableFit();
		TOFfits->SetFitConstraints(0.9,1,0.015,0.16,0.005,0.02,true);
		TOFfits->ExtractCounts(finalhistos);	
		TOFfits->SaveFitResults(finalresults);

		NaFfits->SetFitRange(0.6,6);
//		NaFfits->DisableFit();
		NaFfits->SetFitConstraints(0.9,1,0.001,0.1,0.0005,0.02,true);
		NaFfits->ExtractCounts(finalhistos);
		NaFfits->SaveFitResults(finalresults);

		Aglfits->SetFitRange(0.6,5);
//		Aglfits->DisableFit();
		Aglfits->SetFitConstraints(0.9,1,0.001,0.1,0.0005,0.02,true);
		Aglfits->ExtractCounts(finalhistos);
		Aglfits->SaveFitResults(finalresults);	

		TOFfits_P->SetFitRange(0.6,4);
//		TOFfits_P->DisableFit();
		TOFfits_P->SetFitConstraints(0.9,1,0.015,0.16,0.005,0.02,true);
		TOFfits_P->ExtractCounts(finalhistos);	
		TOFfits_P->SaveFitResults(finalresults);

		NaFfits_P->SetFitRange(0.6,6);
//		NaFfits_P->DisableFit();
		NaFfits_P->SetFitConstraints(0.9,1,0.001,0.1,0.0005,0.02,true);
		NaFfits_P->ExtractCounts(finalhistos);
		NaFfits_P->SaveFitResults(finalresults);

		Aglfits_P->SetFitRange(0.6,5);
//		Aglfits_P->DisableFit();
		Aglfits_P->SetFitConstraints(0.9,1,0.001,0.1,0.0005,0.02,true);
		Aglfits_P->ExtractCounts(finalhistos);
		Aglfits_P->SaveFitResults(finalresults);	


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
	TOFfits_P	    -> SetDefaultOutFile(finalhistos);
	NaFfits_P	    -> SetDefaultOutFile(finalhistos);
	Aglfits_P	    -> SetDefaultOutFile(finalhistos);


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
	TOFfits_P	 ->SetLatitudeReweighter(weighter);	
	NaFfits_P	 ->SetLatitudeReweighter(weighter);	
	Aglfits_P	 ->SetLatitudeReweighter(weighter);	

	NaFfits->SetUpBadEventSimulator(NaFBadEvSimulator);
	Aglfits->SetUpBadEventSimulator(AglBadEvSimulator);
	NaFfits_P->SetUpBadEventSimulator(NaFBadEvSimulator);
	Aglfits_P->SetUpBadEventSimulator(AglBadEvSimulator);


	//Filler.AddObject2beFilled(SmearingCheck,GetBetaTOF,GetRigidity);
	Filler.AddObject2beFilled(TOFfits,GetRecMassTOF ,GetBetaTOF);
	Filler.AddObject2beFilled(NaFfits,GetRecMassRICH,GetBetaRICH);
	Filler.AddObject2beFilled(Aglfits,GetRecMassRICH,GetBetaRICH);
	Filler.AddObject2beFilled(TOFfits_P,GetRecMassTOF ,GetBetaTOF);
	Filler.AddObject2beFilled(NaFfits_P,GetRecMassRICH,GetBetaRICH);
	Filler.AddObject2beFilled(Aglfits_P,GetRecMassRICH,GetBetaRICH);

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

