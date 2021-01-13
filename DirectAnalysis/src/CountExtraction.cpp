#ifndef COUNTEXTRACTION_CPP
#define COUNTEXTRACTION_CPP


#include "Analyzer.h"
#include "Variables.hpp"
#include "TemplateFITbetasmear.h"
#include "Efficiency.h"
#include "GlobalPaths.h"

void Analyzer::BookCountsAnalysis(FileSaver finalhistos, FileSaver finalresults, bool refill)
{

	bool checkfile = finalhistos.CheckFile();
	check_file = checkfile;


	cout<<"****************************** BINS ***************************************"<<endl;
    	SetUpUsualBinning();
 
	cout<<"****************************** Counts ANALYIS ******************************************"<<endl;

	BadEventSimulator * NaFBadEvSimulator= new BadEventSimulator("IsFromNaF",22,0.72,1); 
	BadEventSimulator * AglBadEvSimulator= new BadEventSimulator("IsFromAgl",250,0.95,1); 


	Efficiency * CountsHE     = new Efficiency(finalhistos    ,"HEPCounts"         ,"HEPCounts"      ,PRB,"RigSafetyCut&IsPositive&IsPrimary&IsBaseline","RigSafetyCut&IsPositive&IsPrimary&IsBaseline");
	Efficiency * CountsL1HE   = new Efficiency(finalhistos    ,"HEPCountsL1"       ,"HEPCountsL1"    ,PRB,"RigSafetyCut&IsPositive&IsPrimary&IsBaseline&L1LooseCharge1","RigSafetyCut&IsPositive&IsPrimary&IsBaseline&L1LooseCharge1");
	Efficiency * CountsQualHE = new Efficiency(finalhistos    ,"HEPCountsQual"     ,"HEPCountsQual"  ,PRB,"RigSafetyCut&IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning","RigSafetyCut&IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning");

	Efficiency * CountsTOF    = new Efficiency(finalhistos    ,"TOFPCounts"	   ,"TOFPCounts"	,GlobalRig.GetToFPBins(),"RigSafetyCut&IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning"    ,"RigSafetyCut&IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning");
	Efficiency * CountsNaF    = new Efficiency(finalhistos    ,"NaFPCounts"	   ,"NaFPCounts"	,GlobalRig.GetNaFPBins(),"RigSafetyCut&IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning"    ,"RigSafetyCut&IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning"    );
	Efficiency * CountsAgl    = new Efficiency(finalhistos    ,"AglPCounts"	   ,"AglPCounts"	,GlobalRig.GetAglPBins(),"RigSafetyCut&IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning"    ,"RigSafetyCut&IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning"    );

	Efficiency * CountsBetaTOF    = new Efficiency(finalhistos    ,"TOFPCountsBeta"	   ,"TOFPCountsBeta"	,GlobalRig.GetToFPBins(),"RigSafetyCut&IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&QualityTOF"   ,"RigSafetyCut&IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&QualityTOF");
	Efficiency * CountsBetaNaF    = new Efficiency(finalhistos    ,"NaFPCountsBeta"	   ,"NaFPCountsBeta"	,GlobalRig.GetNaFPBins(),"RigSafetyCut&IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut"    ,"RigSafetyCut&IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut"     );
	Efficiency * CountsBetaAgl    = new Efficiency(finalhistos    ,"AglPCountsBeta"	   ,"AglPCountsBeta"	,GlobalRig.GetAglPBins(),"RigSafetyCut&IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut"    ,"RigSafetyCut&IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut"     );


	Efficiency * CountsFitTOF    = new Efficiency(finalhistos    ,"TOFPCountsFit"	   ,"TOFPCountsFit"	,Global.GetToFPBins(),"RigSafetyCut&IsPositive&IsPrimaryBetaTOFP&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&QualityTOF"   ,"RigSafetyCut&IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&QualityTOF");
	Efficiency * CountsFitNaF    = new Efficiency(finalhistos    ,"NaFPCountsFit"	   ,"NaFPCountsFit"	,Global.GetNaFPBins(),"RigSafetyCut&IsPositive&IsPrimaryBetaRICP&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut"    ,"RigSafetyCut&IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut"     );
	Efficiency * CountsFitAgl    = new Efficiency(finalhistos    ,"AglPCountsFit"	   ,"AglPCountsFit"	,Global.GetAglPBins(),"RigSafetyCut&IsPositive&IsPrimaryBetaRICP&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut"    ,"RigSafetyCut&IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut"     );



	// Extraction of counts with Template Fit

	TemplateFIT * TOFfits= new TemplateFIT("TOFDfits",Global.GetToFDBins(),"IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&QualityTOF","IsPrimaryBetaTOFD"       ,150,0.4,7.5,false,9,50,170,1);
//	TemplateFIT * NaFfits= new TemplateFIT("NaFDfits",Global.GetNaFDBins(),"IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut","IsPrimaryBetaRICD"           ,60,0.4,5,true,11,250,150,0);
//	TemplateFIT * Aglfits= new TemplateFIT("AglDfits",Global.GetAglDBins(),"IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut","IsPrimaryBetaRICD"           ,60,0.4,5,true,11,10,25,1);	
	TemplateFIT * NaFfits= new TemplateFIT("NaFDfits",Global.GetNaFDBins(),"IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut","IsPrimaryBetaRICD"           ,30,0.4,5,true,7,250,200,0);
	TemplateFIT * Aglfits= new TemplateFIT("AglDfits",Global.GetAglDBins(),"IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut","IsPrimaryBetaRICD"           ,60,0.4,5.5,true,11,20,25,1);	


	TemplateFIT * TOFfits_P= new TemplateFIT("TOFPfits",Global.GetToFPBins(),"IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&QualityTOF","IsPrimaryBetaTOFP"       ,150,0.4,7.5,false,9,50,170,1);
//	TemplateFIT * NaFfits_P= new TemplateFIT("NaFPfits",Global.GetNaFPBins(),"IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut","IsPrimaryBetaRICP"           ,60,0.4,5,true,11,250,150,0);
//	TemplateFIT * Aglfits_P= new TemplateFIT("AglPfits",Global.GetAglPBins(),"IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut","IsPrimaryBetaRICP"           ,60,0.4,5,true,11,10,25,1);	
	TemplateFIT * NaFfits_P= new TemplateFIT("NaFPfits",Global.GetNaFPBins(),"IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut","IsPrimaryBetaRICP"           ,30,0.4,5,true,7,250,200,0);
	TemplateFIT * Aglfits_P= new TemplateFIT("AglPfits",Global.GetAglPBins(),"IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut","IsPrimaryBetaRICP"           ,60,0.4,5.5,true,11,20,25,1);	



	if(!refill&&checkfile) {	

		CountsHE    ->ReinitializeHistos(refill); 
	       	CountsL1HE    ->ReinitializeHistos(refill); 
	        CountsQualHE->ReinitializeHistos(refill); 
	        CountsTOF   ->ReinitializeHistos(refill); 
	        CountsNaF   ->ReinitializeHistos(refill); 
                CountsAgl   ->ReinitializeHistos(refill); 
		CountsBetaTOF   ->ReinitializeHistos(refill); 
	        CountsBetaNaF   ->ReinitializeHistos(refill); 
                CountsBetaAgl   ->ReinitializeHistos(refill); 
		CountsFitTOF   ->ReinitializeHistos(refill); 
	        CountsFitNaF   ->ReinitializeHistos(refill); 
                CountsFitAgl   ->ReinitializeHistos(refill); 


	
		TOFfits= new TemplateFIT(finalhistos,"TOFDfits",Global.GetToFDBins(),false,9,50,170,1);
		NaFfits= new TemplateFIT(finalhistos,"NaFDfits",Global.GetNaFDBins(),true,7,250,200,0);
		Aglfits= new TemplateFIT(finalhistos,"AglDfits",Global.GetAglDBins(),true,11,20,25,1);

		TOFfits_P= new TemplateFIT(finalhistos,"TOFPfits",Global.GetToFPBins(),false,9,50,170,1);
		NaFfits_P= new TemplateFIT(finalhistos,"NaFPfits",Global.GetNaFPBins(),true,7,250,200,0);
		Aglfits_P= new TemplateFIT(finalhistos,"AglPfits",Global.GetAglPBins(),true,11,20,25,1);


		//	NaFfits->SetFitWithNoiseMode();
		//	Aglfits->SetFitWithNoiseMode();

		bool disable_fits=false;

		TOFfits->SetFitRange(0.85,5.5);
		if(disable_fits) TOFfits->DisableFit();
		TOFfits->SetFitConstraints(0.9,1,0.015,0.16,0.005,0.02,true);
		TOFfits->ExtractCounts(finalhistos);	
		if(!disable_fits) TOFfits->SaveFitResults(finalresults);

	
		NaFfits->SetFitRange(0.85,6);
		if(disable_fits) NaFfits->DisableFit();
		NaFfits->SetFitConstraints(0.9,1,0.001,0.1,0.000005,0.02,true);
		NaFfits->ExtractCounts(finalhistos);
		if(!disable_fits) NaFfits->SaveFitResults(finalresults);
		

		Aglfits->SetLocalFit();
		Aglfits->SetFitRange(0.5,6);
		if(disable_fits)  Aglfits->DisableFit();
		Aglfits->SetFitConstraints(0.9,1,0.001,0.1,0.000005,0.02,true);
		Aglfits->ExtractCounts(finalhistos);
		if(!disable_fits) 	Aglfits->SaveFitResults(finalresults);	

		disable_fits = true;
/*	
		TOFfits_P->SetFitRange(0.5,5.5);
		if(disable_fits) TOFfits_P->DisableFit();
		TOFfits_P->SetFitConstraints(0.9,1,0.015,0.16,0.005,0.02,true);
		TOFfits_P->ExtractCounts(finalhistos);	
		if(!disable_fits) 	TOFfits_P->SaveFitResults(finalresults);

		NaFfits_P->SetFitRange(0.6,6);
		if(disable_fits) NaFfits_P->DisableFit();
		NaFfits_P->SetFitConstraints(0.9,1,0.001,0.1,0.0005,0.02,true);
		NaFfits_P->ExtractCounts(finalhistos);
		if(!disable_fits) NaFfits_P->SaveFitResults(finalresults);

		Aglfits->SetLocalFit();
		Aglfits_P->SetFitRange(0.6,6);
		if(disable_fits) Aglfits_P->DisableFit();
		Aglfits_P->SetFitConstraints(0.9,1,0.001,0.1,0.0005,0.02,true);
		Aglfits_P->ExtractCounts(finalhistos);
		if(!disable_fits) 	Aglfits_P->SaveFitResults(finalresults);	

*/
		CountsHE	->Eval_Efficiency();	
		CountsL1HE	->Eval_Efficiency();	
		CountsQualHE	->Eval_Efficiency();
		CountsTOF	->Eval_Efficiency();
		CountsNaF	->Eval_Efficiency();
		CountsAgl	->Eval_Efficiency();
		CountsBetaTOF	->Eval_Efficiency();
		CountsBetaNaF	->Eval_Efficiency();
		CountsBetaAgl	->Eval_Efficiency();
		CountsFitTOF	->Eval_Efficiency();
		CountsFitNaF	->Eval_Efficiency();
		CountsFitAgl	->Eval_Efficiency();


		CountsHE	->SaveResults(finalresults);   
		CountsL1HE	->SaveResults(finalresults);   
		CountsQualHE	->SaveResults(finalresults);
		CountsTOF	->SaveResults(finalresults);
		CountsNaF	->SaveResults(finalresults);
		CountsAgl	->SaveResults(finalresults);
		CountsBetaTOF	->SaveResults(finalresults);
		CountsBetaNaF	->SaveResults(finalresults);
		CountsBetaAgl	->SaveResults(finalresults);
		CountsFitTOF	->SaveResults(finalresults);
		CountsFitNaF	->SaveResults(finalresults);
		CountsFitAgl	->SaveResults(finalresults);


	}


	TOFfits	    -> SetDefaultOutFile(finalhistos);
	NaFfits	    -> SetDefaultOutFile(finalhistos);
	Aglfits	    -> SetDefaultOutFile(finalhistos);
	TOFfits_P	    -> SetDefaultOutFile(finalhistos);
	NaFfits_P	    -> SetDefaultOutFile(finalhistos);
	Aglfits_P	    -> SetDefaultOutFile(finalhistos);

	CountsHE	-> SetDefaultOutFile(finalhistos);
	CountsL1HE	-> SetDefaultOutFile(finalhistos);
	CountsQualHE-> SetDefaultOutFile(finalhistos);
	CountsTOF	-> SetDefaultOutFile(finalhistos);
	CountsNaF	-> SetDefaultOutFile(finalhistos);
	CountsAgl	-> SetDefaultOutFile(finalhistos);
	CountsBetaTOF	-> SetDefaultOutFile(finalhistos);
	CountsBetaNaF	-> SetDefaultOutFile(finalhistos);
	CountsBetaAgl	-> SetDefaultOutFile(finalhistos);
	CountsFitTOF	-> SetDefaultOutFile(finalhistos);
	CountsFitNaF	-> SetDefaultOutFile(finalhistos);
	CountsFitAgl	-> SetDefaultOutFile(finalhistos);


	NaFfits->SetUpBadEventSimulator(NaFBadEvSimulator);
	Aglfits->SetUpBadEventSimulator(AglBadEvSimulator);
	NaFfits_P->SetUpBadEventSimulator(NaFBadEvSimulator);
	Aglfits_P->SetUpBadEventSimulator(AglBadEvSimulator);

	Filler.AddObject2beFilled(TOFfits,GetRecMassTOF ,GetBetaTOF);
	Filler.AddObject2beFilled(NaFfits,GetRecMassRICH,GetBetaRICH);
	Filler.AddObject2beFilled(Aglfits,GetRecMassRICH,GetBetaRICH);
	//Filler.AddObject2beFilled(TOFfits_P,GetRecMassTOF ,GetBetaTOF);
	//Filler.AddObject2beFilled(NaFfits_P,GetRecMassRICH,GetBetaRICH);
	//Filler.AddObject2beFilled(Aglfits_P,GetRecMassRICH,GetBetaRICH);
	

	Filler.AddObject2beFilled(CountsHE,GetRigidity,GetRigidity); 
	Filler.AddObject2beFilled(CountsL1HE,GetRigidity,GetRigidity); 
	Filler.AddObject2beFilled(CountsQualHE,GetRigidity,GetRigidity); 
	Filler.AddObject2beFilled(CountsTOF,GetRigidity,GetRigidity); 
	Filler.AddObject2beFilled(CountsNaF,GetRigidity,GetRigidity); 
	Filler.AddObject2beFilled(CountsAgl,GetRigidity,GetRigidity); 
	Filler.AddObject2beFilled(CountsBetaTOF,GetRigidity,GetRigidity); 
	Filler.AddObject2beFilled(CountsBetaNaF,GetRigidity,GetRigidity); 
	Filler.AddObject2beFilled(CountsBetaAgl,GetRigidity,GetRigidity); 

	
	Filler.AddObject2beFilled(CountsFitTOF,GetBetaTOF ,GetBetaTOF ); 
	Filler.AddObject2beFilled(CountsFitNaF,GetBetaRICH,GetBetaRICH); 
	Filler.AddObject2beFilled(CountsFitAgl,GetBetaRICH,GetBetaRICH); 
	


	Filler.ReinitializeAll(refill);


}

#endif

