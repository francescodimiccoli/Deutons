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
 
	cout<<"****************************** Counts ANALYSIS ******************************************"<<endl;

	BadEventSimulator * NaFBadEvSimulator= new BadEventSimulator("IsFromNaF",22,0.72,1); 
	BadEventSimulator * AglBadEvSimulator= new BadEventSimulator("IsFromAgl",250,0.95,1); 

													       
	Efficiency * CountsHE     = new Efficiency(finalhistos    ,"HEPCounts"         ,"HEPCounts"      ,PRB,"RigSafetyCut&IsPositive&IsPrimary&IsBaseline","RigSafetyCut&IsPositive&IsPrimary&IsBaseline");
	Efficiency * CountsL1HE   = new Efficiency(finalhistos    ,"HEPCountsL1"       ,"HEPCountsL1"    ,PRB,"RigSafetyCut&IsPositive&IsPrimary&IsBaseline&L1LooseCharge1","RigSafetyCut&IsPositive&IsPrimary&IsBaseline&L1LooseCharge1");
	Efficiency * CountsQualHE = new Efficiency(finalhistos    ,"HEPCountsQual"     ,"HEPCountsQual"  ,PRB,"RigSafetyCut&IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning","RigSafetyCut&IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning");
	Efficiency * CountsHeQualHE = new Efficiency(finalhistos    ,"HEHeCountsQual"     ,"HEHeCountsQual"  ,HeRB,"RigSafetyCut&IsPositive&IsPrimary&IsBaselineHe&L1LooseCharge2&IsCleaningHe","RigSafetyCut&IsPositive&IsPrimary&IsBaselineHe&L1LooseCharge2&IsCleaningHe");



	Efficiency * CountsTOF    = new Efficiency(finalhistos    ,"TOFPCounts"	   ,"TOFPCounts"	,GlobalRig.GetToFPBins(),"RigSafetyCut&IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning"    ,"RigSafetyCut&IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning");
	Efficiency * CountsNaF    = new Efficiency(finalhistos    ,"NaFPCounts"	   ,"NaFPCounts"	,GlobalRig.GetNaFPBins(),"RigSafetyCut&IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning"    ,"RigSafetyCut&IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning"    );
	Efficiency * CountsAgl    = new Efficiency(finalhistos    ,"AglPCounts"	   ,"AglPCounts"	,GlobalRig.GetAglPBins(),"RigSafetyCut&IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning"    ,"RigSafetyCut&IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning"    );

	Efficiency * CountsBetaTOF    = new Efficiency(finalhistos    ,"TOFPCountsBeta"	   ,"TOFPCountsBeta"	,GlobalRig.GetToFPBins(),"RigSafetyCut&IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&QualityTOF"   ,"RigSafetyCut&IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&QualityTOF");
	Efficiency * CountsBetaNaF    = new Efficiency(finalhistos    ,"NaFPCountsBeta"	   ,"NaFPCountsBeta"	,GlobalRig.GetNaFPBins(),"RigSafetyCut&IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut"    ,"RigSafetyCut&IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut"     );
	Efficiency * CountsBetaAgl    = new Efficiency(finalhistos    ,"AglPCountsBeta"	   ,"AglPCountsBeta"	,GlobalRig.GetAglPBins(),"IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut"    ,"RigSafetyCut&IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut"     );


	Efficiency * CountsFitTOF    = new Efficiency(finalhistos    ,"TOFPCountsFit"	   ,"TOFPCountsFit"	,Global.GetToFPBins(),"RigSafetyCut&IsPositive&IsPrimaryBetaTOFP&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&QualityTOF"   ,"RigSafetyCut&IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&QualityTOF");
	Efficiency * CountsFitNaF    = new Efficiency(finalhistos    ,"NaFPCountsFit"	   ,"NaFPCountsFit"	,Global.GetNaFPBins(),"RigSafetyCut&IsPositive&IsPrimaryBetaRICP&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut"    ,"RigSafetyCut&IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut"     );
	Efficiency * CountsFitAgl    = new Efficiency(finalhistos    ,"AglPCountsFit"	   ,"AglPCountsFit"	,Global.GetAglPBins(),"RigSafetyCut&IsPositive&IsPrimaryBetaRICP&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut"    ,"RigSafetyCut&IsPositive&IsPrimary&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut"     );

	Efficiency * CountsHeTOF    = new Efficiency(finalhistos    ,"TOFHeCounts"	   ,"TOFHeCounts"	,Global_HeRig.GetToFDBins(),"RigSafetyCut&IsPositive&IsPrimary&IsBaselineHe&L1LooseCharge2&IsCleaningHe"    ,"RigSafetyCut&IsPositive&IsPrimary&IsBaselineHe&L1LooseCharge2&IsCleaningHe");
	Efficiency * CountsHeNaF    = new Efficiency(finalhistos    ,"NaFHeCounts"	   ,"NaFHeCounts"	,Global_HeRig.GetNaFDBins(),"RigSafetyCut&IsPositive&IsPrimary&IsBaselineHe&L1LooseCharge2&IsCleaningHe"    ,"RigSafetyCut&IsPositive&IsPrimary&IsBaselineHe&L1LooseCharge2&IsCleaningHe"    );
	Efficiency * CountsHeAgl    = new Efficiency(finalhistos    ,"AglHeCounts"	   ,"AglHeCounts"	,Global_HeRig.GetAglDBins(),"RigSafetyCut&IsPositive&IsPrimary&IsBaselineHe&L1LooseCharge2&IsCleaningHe"    ,"RigSafetyCut&IsPositive&IsPrimary&IsBaselineHe&L1LooseCharge2&IsCleaningHe"    );

	Efficiency * CountsHeBetaTOF    = new Efficiency(finalhistos    ,"TOFHeBetaCounts"	   ,"TOFHeBetaCounts"	,Global_HeRig.GetToFDBins(),"RigSafetyCut&IsPositive&IsPrimary&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsGoodTimeHe"     	        ,"RigSafetyCut&IsPositive&IsPrimary&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsGoodTimeHe");
	Efficiency * CountsHeBetaNaF    = new Efficiency(finalhistos    ,"NaFHeBetaCounts"	   ,"NaFHeBetaCounts"	,Global_HeRig.GetNaFDBins(),"RigSafetyCut&IsPositive&IsPrimary&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromNaF&RICHHeCutNaF"    ,"RigSafetyCut&IsPositive&IsPrimary&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromNaF&RICHHeCutNaF");
	Efficiency * CountsHeBetaAgl    = new Efficiency(finalhistos    ,"AglHeBetaCounts"	   ,"AglHeBetaCounts"	,Global_HeRig.GetAglDBins(),"RigSafetyCut&IsPositive&IsPrimary&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromAgl&RICHHeCutAgl"    ,"RigSafetyCut&IsPositive&IsPrimary&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromAgl&RICHHeCutAgl");

	Efficiency * CountsHeFitTOF    = new Efficiency(finalhistos    ,"TOFHeFitCounts"	   ,"TOFHeFitCounts"	,Global_He.GetToFDBins(),"RigSafetyCut&IsPositive&IsPrimaryBetaTOF4He&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsGoodTimeHe"     	   ,"RigSafetyCut&IsPositive&IsPrimary&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsGoodTimeHe");
	Efficiency * CountsHeFitNaF    = new Efficiency(finalhistos    ,"NaFHeFitCounts"	   ,"NaFHeFitCounts"	,Global_He.GetNaFDBins(),"RigSafetyCut&IsPositive&IsPrimaryBetaRIC4He&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromNaF&RICHHeCutNaF"    ,"RigSafetyCut&IsPositive&IsPrimary&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromNaF&RICHHeCutNaF");
	Efficiency * CountsHeFitAgl    = new Efficiency(finalhistos    ,"AglHeFitCounts"	   ,"AglHeFitCounts"	,Global_He.GetAglDBins(),"RigSafetyCut&IsPositive&IsPrimaryBetaRIC4He&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromAgl&RICHHeCutAgl"    ,"RigSafetyCut&IsPositive&IsPrimary&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromAgl&RICHHeCutAgl");



	// Extraction of counts with Template Fit

	
	//deuterons
	TemplateFIT * TOFfits= new TemplateFIT("TOFDfits",Global.GetToFDBins()  ,"IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&QualityTOF","&IsProtonMC","&IsDeutonMC","&IsProtonMC","IsPrimaryBetaTOFD"       ,150,0.4,7.5,false,9,50,60,1);
//	TemplateFIT * NaFfits= new TemplateFIT("NaFDfits",Global.GetNaFDBins(),"IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut","&IsProtonMC","&IsDeutonMC","&IsProtonMC","IsPrimaryBetaRICD"           ,60,0.4,5,true,9,450,300,1);
	TemplateFIT * NaFfits= new TemplateFIT("NaFDfits",Global.GetNaFDBins(),"IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut","&IsProtonMC","&IsProtonMC","&IsProtonMC","IsPrimaryBetaRICD"           ,90,0.4,5,true,9,40,45,1);
	TemplateFIT * Aglfits= new TemplateFIT("AglDfits",Global.GetAglDBins(),"IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut","&IsProtonMC","&IsDeutonMC","&IsProtonMC","IsPrimaryBetaRICD"           ,90,0.4,7,true,11,40,45,1);	

	//protons
	TemplateFIT * TOFfits_P= new TemplateFIT("TOFPfits",Global.GetToFPBins(),"IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&QualityTOF","&IsProtonMC","&IsDeutonMC","&IsProtonMC","IsPrimaryBetaTOFP"       ,150,0.4,7.5,false,9,50,60,1);
	//TemplateFIT * NaFfits_P= new TemplateFIT("NaFPfits",Global.GetNaFPBins(),"IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut","&IsProtonMC","&IsDeutonMC","&IsProtonMC","IsPrimaryBetaRICP"           ,60,0.4,5,true,9,450,300,0);
	TemplateFIT * NaFfits_P= new TemplateFIT("NaFPfits",Global.GetNaFPBins(),"IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut","&IsProtonMC","&IsProtonMC","&IsProtonMC","IsPrimaryBetaRICP"           ,90,0.4,5,true,9,20,25,1);
	TemplateFIT * Aglfits_P= new TemplateFIT("AglPfits",Global.GetAglPBins(),"IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut","&IsProtonMC","&IsDeutonMC","&IsProtonMC","IsPrimaryBetaRICP"           ,90,0.4,7,true,11,20,25,1);	

	//Helium4
	TemplateFIT * TOFfits_He= new TemplateFIT("TOFHefits",Global_He.GetToFDBins(),"IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsGoodTimeHe" ,"&IsHeliumMC","&IsHeliumMC","&IsHeliumMC","IsPrimaryBetaTOF4He"       ,150,0.4,7.5,false,9,90,90,1);
	TemplateFIT * NaFfits_He= new TemplateFIT("NaFHefits",Global_He.GetNaFDBins(),"IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromNaF&RICHHeCutNaF","&IsHeliumMC","&IsHeliumMC","&IsHeliumMC","IsPrimaryBetaRIC4He"           ,90,0.4,5,true,7,250,200,1);
	TemplateFIT * Aglfits_He= new TemplateFIT("AglHEfits",Global_He.GetAglDBins(),"IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromAgl&RICHHeCutAgl","&IsHeliumMC","&IsHeliumMC","&IsHeliumMC","IsPrimaryBetaRIC4He"           ,90,0.4,5.5,true,11,75,55,1);	

	//Helium3
	TemplateFIT * TOFfits_He3= new TemplateFIT("TOFHe3fits",Global_He.GetToFPBins(),"IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsGoodTimeHe" ,"&IsHeliumMC","&IsHeliumMC","&IsHeliumMC","IsPrimaryBetaTOF3He"       ,150,0.4,7.5,false,9,90,90,1);
	TemplateFIT * NaFfits_He3= new TemplateFIT("NaFHe3fits",Global_He.GetNaFPBins(),"IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromNaF&RICHHeCutNaF","&IsHeliumMC","&IsHeliumMC","&IsHeliumMC","IsPrimaryBetaRIC3He"           ,90,0.4,5,true,7,250,200,1);
	TemplateFIT * Aglfits_He3= new TemplateFIT("AglHE3fits",Global_He.GetAglPBins(),"IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromAgl&RICHHeCutAgl","&IsHeliumMC","&IsHeliumMC","&IsHeliumMC","IsPrimaryBetaRIC3He"           ,90,0.4,5.5,true,11,75,55,1);	

/*
	TOFfits	 ->SetAsExtern();
	NaFfits	 ->SetAsExtern();
	Aglfits	 ->SetAsExtern();

	TOFfits_P->SetAsExtern();
	NaFfits_P->SetAsExtern();
	Aglfits_P->SetAsExtern();

	TOFfits_He->SetAsExtern();
	NaFfits_He->SetAsExtern();
	Aglfits_He->SetAsExtern();
*/	
/*
//	TOFfits->SetNotWeightedMC();
	//TOFfits->SetNotTunedMC();
	NaFfits->SetNotWeightedMC();
	//NaFfits->SetNotTunedMC();
	Aglfits->SetNotWeightedMC();
	//Aglfits->SetNotTunedMC();
//	TOFfits_P->SetNotWeightedMC();
	//TOFfits_P->SetNotTunedMC();
	NaFfits_P->SetNotWeightedMC();
	//NaFfits_P->SetNotTunedMC();
	Aglfits_P->SetNotWeightedMC();
	//Aglfits_P->SetNotTunedMC();

	TOFfits_He->SetNotWeightedMC();
//	TOFfits_He->SetNotTunedMC();
	NaFfits_He->SetNotWeightedMC();
//	NaFfits_He->SetNotTunedMC();
	Aglfits_He->SetNotWeightedMC();
//	Aglfits_He->SetNotTunedMC();
	TOFfits_He3->SetNotWeightedMC();
//	TOFfits_He3->SetNotTunedMC();
	NaFfits_He3->SetNotWeightedMC();
//	NaFfits_He3->SetNotTunedMC();
	Aglfits_He3->SetNotWeightedMC();
//	Aglfits_He3->SetNotTunedMC();
*/	
	if(!refill&&checkfile) {	

		CountsHE    ->ReinitializeHistos(refill); 
	       	CountsL1HE    ->ReinitializeHistos(refill); 
	        CountsQualHE->ReinitializeHistos(refill); 
	        CountsHeQualHE->ReinitializeHistos(refill); 
	        CountsTOF   ->ReinitializeHistos(refill); 
	        CountsNaF   ->ReinitializeHistos(refill); 
                CountsAgl   ->ReinitializeHistos(refill); 
		CountsBetaTOF   ->ReinitializeHistos(refill); 
	        CountsBetaNaF   ->ReinitializeHistos(refill); 
                CountsBetaAgl   ->ReinitializeHistos(refill); 
		CountsFitTOF   ->ReinitializeHistos(refill); 
	        CountsFitNaF   ->ReinitializeHistos(refill); 
                CountsFitAgl   ->ReinitializeHistos(refill); 
	        CountsHeTOF   ->ReinitializeHistos(refill); 
	        CountsHeNaF   ->ReinitializeHistos(refill); 
                CountsHeAgl   ->ReinitializeHistos(refill); 
		CountsHeBetaTOF   ->ReinitializeHistos(refill); 
	        CountsHeBetaNaF   ->ReinitializeHistos(refill); 
                CountsHeBetaAgl   ->ReinitializeHistos(refill); 
		CountsHeFitTOF   ->ReinitializeHistos(refill); 
	        CountsHeFitNaF   ->ReinitializeHistos(refill); 
                CountsHeFitAgl   ->ReinitializeHistos(refill); 
		
	
		TOFfits	 = new TemplateFIT(finalhistos,"TOFDfits",Global.GetToFDBins(),false,9,50,60,1);
		NaFfits	 = new TemplateFIT(finalhistos,"NaFDfits",Global.GetNaFDBins(),true,9,40,45,1,1);
		Aglfits	 = new TemplateFIT(finalhistos,"AglDfits",Global.GetAglDBins(),true,11,40,45,1,1);

		TOFfits_P= new TemplateFIT(finalhistos,"TOFPfits",Global.GetToFPBins(),false,9,50,60,1,1);
		NaFfits_P= new TemplateFIT(finalhistos,"NaFPfits",Global.GetNaFPBins(),true,9,20,25,1,1);
		Aglfits_P= new TemplateFIT(finalhistos,"AglPfits",Global.GetAglPBins(),true,11,20,25,1,1);

		TOFfits_He= new TemplateFIT(finalhistos,"TOFHefits",Global_He.GetToFDBins(),false,9,90,90,1);
		NaFfits_He= new TemplateFIT(finalhistos,"NaFHefits",Global_He.GetNaFDBins(),true,7,250,200,1,1);
		Aglfits_He= new TemplateFIT(finalhistos,"AglHEfits",Global_He.GetAglDBins(),true,11,75,55,1,1);

		TOFfits_He3= new TemplateFIT(finalhistos,"TOFHe3fits",Global_He.GetToFPBins(),false,9,90,90,1);
		NaFfits_He3= new TemplateFIT(finalhistos,"NaFHe3fits",Global_He.GetNaFPBins(),true,7,250,200,1,1);
		Aglfits_He3= new TemplateFIT(finalhistos,"AglHE3fits",Global_He.GetAglPBins(),true,11,75,55,1,1);

		//	NaFfits->SetFitWithNoiseMode();
		//	Aglfits->SetFitWithNoiseMode();

		TFile * RegularizationFile = TFile::Open("/afs/cern.ch/work/f/fdimicco/private/Deutons/DirectAnalysis/EffSyst/Time.root");
		
		bool disable_fits=false;
/*
		TOFfits->SetLocalFit();
	//	TOFfits->SetRegularizationFile(RegularizationFile,filelistDT.c_str());
		TOFfits->SetFitRange(0.6,5);
		if(disable_fits) TOFfits->DisableFit();
		TOFfits->SetFitConstraints(0.94,1,0.0001,0.16,0.0005,0.0025,false);
		TOFfits->SetAdjoustTailMode(0.6,2.3,11,25);
	//	TOFfits->UseCentroid();
		TOFfits->ExtractCounts(finalhistos);	
		if(!disable_fits) TOFfits->SaveFitResults(finalresults);
*/	
//		disable_fits=true;


		NaFfits->SetLocalFit();
 		//NaFfits->SetFitWithNoiseMode();
	//	NaFfits->SetRegularizationFile(RegularizationFile,filelistDT.c_str());
		NaFfits->SetFitRange(0.85,6);
		if(disable_fits) NaFfits->DisableFit();
		NaFfits->SetFitConstraints(0.9,1,0.001,0.1,0.0005,0.01,false);
		NaFfits->SetFixModTailMode(-1.5,1.7,6,8);
		NaFfits->SetAdjoustPeakMode(1,-0.02);
//		NaFfits->SetLowStatDMode();
		NaFfits->ExtractCounts(finalhistos);
		if(!disable_fits) NaFfits->SaveFitResults(finalresults);
/*
		disable_fits=false;

		Aglfits->SetLocalFit();
		Aglfits->SetFitRange(0.5,7);
		if(disable_fits)  Aglfits->DisableFit();
		Aglfits->SetFitConstraints(0.9,1,0.00072,0.1,0.0001,0.02,false);
                Aglfits->SetFixModTailMode(-0.6,1.8,6,5);
		
		Aglfits->ExtractCounts(finalhistos);
		if(!disable_fits) Aglfits->SaveFitResults(finalresults);	

		disable_fits=false;

		TOFfits_P->SetLocalFit();
		TOFfits_P->SetFitRange(0.6,5);
		if(disable_fits) TOFfits_P->DisableFit();
		TOFfits_P->SetFitConstraints(0.94,1,0.0001,0.16,0.0005,0.025,false);
		TOFfits_P->SetAdjoustTailMode(0.6,2.3,3,25);
		TOFfits_P->ExtractCounts(finalhistos);	
		if(!disable_fits) 	TOFfits_P->SaveFitResults(finalresults);
*/
	//	disable_fits=false;

		NaFfits_P->SetLocalFit();
		NaFfits_P->SetFitRange(0.85,6);
		if(disable_fits) NaFfits_P->DisableFit();
		NaFfits_P->SetFitConstraints(0.9,1,0.001,0.1,0.0005,0.02,false);
		NaFfits_P->SetFixModTailMode(-1.5,1.7,6,8);
		NaFfits_P->SetAdjoustPeakMode(1,-0.02);
		//NaFfits_P->SetLowStatDMode();
		NaFfits_P->ExtractCounts(finalhistos);
		if(!disable_fits) NaFfits_P->SaveFitResults(finalresults);
/*
		disable_fits=false;


		Aglfits_P->SetLocalFit();
		Aglfits_P->SetFitRange(0.6,7);
		if(disable_fits) Aglfits_P->DisableFit();
		Aglfits_P->SetFitConstraints(0.9,1,0.00072,0.1,0.0001,0.02,false);
                Aglfits_P->SetFixModTailMode(-0.6,1.8,5,5);
	//	Aglfits_P->SetUseBestOnly();
		Aglfits_P->ExtractCounts(finalhistos);
		if(!disable_fits) 	Aglfits_P->SaveFitResults(finalresults);	


		//disable_fits=true;

		TOFfits_He->SetLocalFit();
//		TOFfits_He->SetRegularizationFile(RegularizationFile,filelistDT.c_str());
		TOFfits_He->SetFitRange(0.65,5.5);
		TOFfits_He->SetNotWeightedMC();
		TOFfits_He->SetNotTunedMC();
		TOFfits_He->SetSharpenMode(0.045);
		if(disable_fits) TOFfits_He->DisableFit();
		TOFfits_He->SetFitConstraints(0.6,1,0.015,0.6,0.000001,0.02,false);
		TOFfits_He->SetUseBestOnly();
		TOFfits_He->ExtractCounts(finalhistos);	
		if(!disable_fits) 	TOFfits_He->SaveFitResults(finalresults);
		
//		disable_fits=false;


		NaFfits_He->SetLocalFit();
	//	NaFfits_He->SetRegularizationFile(RegularizationFile,filelistDT.c_str());
		NaFfits_He->SetFitRange(0.6,6);
		NaFfits_He->SetNotWeightedMC();
		NaFfits_He->SetNotTunedMC();
		NaFfits_He->SetAdjoustPeakMode(1,0.02);
		if(disable_fits) NaFfits_He->DisableFit();
		NaFfits_He->SetFitConstraints(0.6,1,0.001,0.6,0.0005,0.02,false);
		NaFfits_He->SetSharpenMode(0.045);
		NaFfits_He->UseCentroid();
		NaFfits_He->ExtractCounts(finalhistos);
		if(!disable_fits) NaFfits_He->SaveFitResults(finalresults);

//		disable_fits=false;

		Aglfits_He->SetLocalFit();
	//	Aglfits_He->SetRegularizationFile(RegularizationFile,filelistDT.c_str());
		Aglfits_He->SetFitRange(0.5,6);
		Aglfits_He->SetNotWeightedMC();
		Aglfits_He->SetNotTunedMC();
		if(disable_fits) Aglfits_He->DisableFit();
		Aglfits_He->SetFitConstraints(0.6,1,0.001,0.6,0.0005,0.02,false);
		Aglfits_He->SetSharpenMode(0.015);
		Aglfits_He->SetAdjoustPeakMode(0,0.01);
		Aglfits_He->SetAdjoustPeakMode(1,0.01);
		Aglfits_He->ExtractCounts(finalhistos);
		if(!disable_fits) 	Aglfits_He->SaveFitResults(finalresults);	


		//disable_fits=false;


		TOFfits_He3->SetLocalFit();
//		TOFfits_He3->SetRegularizationFile(RegularizationFile,filelistDT.c_str());
		TOFfits_He3->SetFitRange(0.65,5.5);
		TOFfits_He3->SetNotWeightedMC();
		TOFfits_He3->SetNotTunedMC();
		TOFfits_He3->SetSharpenMode(0.045);
//		TOFfits_He3->SetAdjoustPeakMode(0,-0.02);
		if(disable_fits) TOFfits_He3->DisableFit();
		TOFfits_He3->SetFitConstraints(0.6,1,0.015,0.6,0.000001,0.02,false);
		TOFfits_He3->SetUseBestOnly();
		TOFfits_He3->ExtractCounts(finalhistos);	
		if(!disable_fits) 	TOFfits_He3->SaveFitResults(finalresults);
		
//		disable_fits=false;


		NaFfits_He3->SetLocalFit();
	//	NaFfits_He3->SetRegularizationFile(RegularizationFile,filelistDT.c_str());
		NaFfits_He3->SetFitRange(0.6,6);
		NaFfits_He3->SetNotWeightedMC();
		NaFfits_He3->SetNotTunedMC();
		//NaFfits_He3->SetAdjoustPeakMode(0,0.01);
		NaFfits_He3->SetAdjoustPeakMode(1,0.02);
		if(disable_fits) NaFfits_He3->DisableFit();
		NaFfits_He3->SetFitConstraints(0.6,1,0.001,0.6,0.0005,0.02,false);
		NaFfits_He3->SetSharpenMode(0.045);
		NaFfits_He3->UseCentroid();
		NaFfits_He3->ExtractCounts(finalhistos);
		if(!disable_fits) NaFfits_He3->SaveFitResults(finalresults);

		//disable_fits=false;

		Aglfits_He3->SetLocalFit();
	//	Aglfits_He3->SetRegularizationFile(RegularizationFile,filelistDT.c_str());
		Aglfits_He3->SetFitRange(0.5,6);
		Aglfits_He3->SetNotWeightedMC();
		Aglfits_He3->SetNotTunedMC();
		if(disable_fits) Aglfits_He3->DisableFit();
		Aglfits_He3->SetFitConstraints(0.6,1,0.001,0.6,0.0005,0.02,false);
		Aglfits_He3->SetSharpenMode(0.015);
		Aglfits_He3->SetAdjoustPeakMode(0,0.01);
		Aglfits_He3->SetAdjoustPeakMode(1,0.01);
		Aglfits_He3->ExtractCounts(finalhistos);
		if(!disable_fits) 	Aglfits_He3->SaveFitResults(finalresults);	
	
	*/
		CountsHE	->Eval_Efficiency();	
		CountsL1HE	->Eval_Efficiency();	
		CountsQualHE	->Eval_Efficiency();
		CountsHeQualHE	->Eval_Efficiency();
		CountsTOF	->Eval_Efficiency();
		CountsNaF	->Eval_Efficiency();
		CountsAgl	->Eval_Efficiency();
		CountsBetaTOF	->Eval_Efficiency();
		CountsBetaNaF	->Eval_Efficiency();
		CountsBetaAgl	->Eval_Efficiency();
		CountsFitTOF	->Eval_Efficiency();
		CountsFitNaF	->Eval_Efficiency();
		CountsFitAgl	->Eval_Efficiency();
		CountsHeTOF	->Eval_Efficiency();
		CountsHeNaF	->Eval_Efficiency();
		CountsHeAgl	->Eval_Efficiency();
		CountsHeBetaTOF	->Eval_Efficiency();
		CountsHeBetaNaF	->Eval_Efficiency();
		CountsHeBetaAgl	->Eval_Efficiency();
		CountsHeFitTOF	->Eval_Efficiency();
		CountsHeFitNaF	->Eval_Efficiency();
		CountsHeFitAgl	->Eval_Efficiency();
		

		CountsHE	->SaveResults(finalresults);   
		CountsL1HE	->SaveResults(finalresults);   
		CountsQualHE	->SaveResults(finalresults);
		CountsHeQualHE	->SaveResults(finalresults);
		CountsTOF	->SaveResults(finalresults);
		CountsNaF	->SaveResults(finalresults);
		CountsAgl	->SaveResults(finalresults);
		CountsBetaTOF	->SaveResults(finalresults);
		CountsBetaNaF	->SaveResults(finalresults);
		CountsBetaAgl	->SaveResults(finalresults);
		CountsFitTOF	->SaveResults(finalresults);
		CountsFitNaF	->SaveResults(finalresults);
		CountsFitAgl	->SaveResults(finalresults);
		CountsHeTOF	->SaveResults(finalresults);
		CountsHeNaF	->SaveResults(finalresults);
		CountsHeAgl	->SaveResults(finalresults);
		CountsHeBetaTOF	->SaveResults(finalresults);
		CountsHeBetaNaF	->SaveResults(finalresults);
		CountsHeBetaAgl	->SaveResults(finalresults);
		CountsHeFitTOF	->SaveResults(finalresults);
		CountsHeFitNaF	->SaveResults(finalresults);
		CountsHeFitAgl	->SaveResults(finalresults);
		

	}


	TOFfits	    	    -> SetDefaultOutFile(finalhistos);
	NaFfits		    -> SetDefaultOutFile(finalhistos);
	Aglfits	   	    -> SetDefaultOutFile(finalhistos);
	TOFfits_P	    -> SetDefaultOutFile(finalhistos);
	NaFfits_P	    -> SetDefaultOutFile(finalhistos);
	Aglfits_P	    -> SetDefaultOutFile(finalhistos);
	TOFfits_He	    -> SetDefaultOutFile(finalhistos);
	NaFfits_He	    -> SetDefaultOutFile(finalhistos);
	Aglfits_He	    -> SetDefaultOutFile(finalhistos);
	TOFfits_He3	    -> SetDefaultOutFile(finalhistos);
	NaFfits_He3	    -> SetDefaultOutFile(finalhistos);
	Aglfits_He3	    -> SetDefaultOutFile(finalhistos);



	TOFfits	    	    ->SetTemplateScaleFactor(1,1,3.); 
	NaFfits		    ->SetTemplateScaleFactor(1,2.,3.); 
	Aglfits	   	    ->SetTemplateScaleFactor(1,1,3.); 
	TOFfits_P	    ->SetTemplateScaleFactor(1,1,3.); 
	NaFfits_P	    ->SetTemplateScaleFactor(1,2,3.); 
	Aglfits_P	    ->SetTemplateScaleFactor(1,1,3.); 
	TOFfits_He	    ->SetTemplateScaleFactor(1,3/4.,0); 
	NaFfits_He	    ->SetTemplateScaleFactor(1,3/4.,0); 
	Aglfits_He	    ->SetTemplateScaleFactor(1,3/4.,0); 
	TOFfits_He3	    ->SetTemplateScaleFactor(1,3/4.,0); 
	NaFfits_He3	    ->SetTemplateScaleFactor(1,3/4.,0); 
	Aglfits_He3	    ->SetTemplateScaleFactor(1,3/4.,0); 



	CountsHE	-> SetDefaultOutFile(finalhistos);
	CountsL1HE	-> SetDefaultOutFile(finalhistos);
	CountsQualHE-> SetDefaultOutFile(finalhistos);
	CountsHeQualHE-> SetDefaultOutFile(finalhistos);
	CountsTOF	-> SetDefaultOutFile(finalhistos);
	CountsNaF	-> SetDefaultOutFile(finalhistos);
	CountsAgl	-> SetDefaultOutFile(finalhistos);
	CountsBetaTOF	-> SetDefaultOutFile(finalhistos);
	CountsBetaNaF	-> SetDefaultOutFile(finalhistos);
	CountsBetaAgl	-> SetDefaultOutFile(finalhistos);
	CountsFitTOF	-> SetDefaultOutFile(finalhistos);
	CountsFitNaF	-> SetDefaultOutFile(finalhistos);
	CountsFitAgl	-> SetDefaultOutFile(finalhistos);
	CountsHeTOF	-> SetDefaultOutFile(finalhistos);
	CountsHeNaF	-> SetDefaultOutFile(finalhistos);
	CountsHeAgl	-> SetDefaultOutFile(finalhistos);
	CountsHeBetaTOF	-> SetDefaultOutFile(finalhistos);
	CountsHeBetaNaF	-> SetDefaultOutFile(finalhistos);
	CountsHeBetaAgl	-> SetDefaultOutFile(finalhistos);
	CountsHeFitTOF	-> SetDefaultOutFile(finalhistos);
	CountsHeFitNaF	-> SetDefaultOutFile(finalhistos);
	CountsHeFitAgl	-> SetDefaultOutFile(finalhistos);
	

	NaFfits->SetUpBadEventSimulator(NaFBadEvSimulator);
	Aglfits->SetUpBadEventSimulator(AglBadEvSimulator);
	NaFfits_P->SetUpBadEventSimulator(NaFBadEvSimulator);
	Aglfits_P->SetUpBadEventSimulator(AglBadEvSimulator);

	Filler.AddObject2beFilled(TOFfits,GetRecMassTOF ,GetBetaTOF);
	Filler.AddObject2beFilled(NaFfits,GetRecMassRICH,GetBetaRICH);
	Filler.AddObject2beFilled(Aglfits,GetRecMassRICH,GetBetaRICH);
	Filler.AddObject2beFilled(TOFfits_P,GetRecMassTOF ,GetBetaTOF);
	Filler.AddObject2beFilled(NaFfits_P,GetRecMassRICH,GetBetaRICH);
	Filler.AddObject2beFilled(Aglfits_P,GetRecMassRICH,GetBetaRICH);
	Filler.AddObject2beFilled(TOFfits_He,GetRecMassTOF ,GetBetaTOF);
	Filler.AddObject2beFilled(NaFfits_He,GetRecMassRICH,GetBetaRICH);
	Filler.AddObject2beFilled(Aglfits_He,GetRecMassRICH,GetBetaRICH);
	Filler.AddObject2beFilled(TOFfits_He3,GetRecMassTOF ,GetBetaTOF);
	Filler.AddObject2beFilled(NaFfits_He3,GetRecMassRICH,GetBetaRICH);
	Filler.AddObject2beFilled(Aglfits_He3,GetRecMassRICH,GetBetaRICH);
	

	Filler.AddObject2beFilled(CountsHE,GetRigidity,GetRigidity); 
	Filler.AddObject2beFilled(CountsL1HE,GetRigidity,GetRigidity); 
	Filler.AddObject2beFilled(CountsQualHE,GetRigidity,GetRigidity); 
	Filler.AddObject2beFilled(CountsHeQualHE,GetRigidity,GetRigidity); 
	Filler.AddObject2beFilled(CountsTOF,GetRigidity,GetRigidity); 
	Filler.AddObject2beFilled(CountsNaF,GetRigidity,GetRigidity); 
	Filler.AddObject2beFilled(CountsAgl,GetRigidity,GetRigidity); 
	Filler.AddObject2beFilled(CountsBetaTOF,GetRigidity,GetRigidity); 
	Filler.AddObject2beFilled(CountsBetaNaF,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(CountsBetaAgl,GetRigidity,GetRigidity); 
	Filler.AddObject2beFilled(CountsFitTOF,GetBetaTOF ,GetBetaTOF ); 
	Filler.AddObject2beFilled(CountsFitNaF,GetBetaRICH,GetBetaRICH); 
	Filler.AddObject2beFilled(CountsFitAgl,GetBetaRICH,GetBetaRICH); 
	Filler.AddObject2beFilled(CountsHeTOF,GetRigidity,GetRigidity); 
	Filler.AddObject2beFilled(CountsHeNaF,GetRigidity,GetRigidity); 
	Filler.AddObject2beFilled(CountsHeAgl,GetRigidity,GetRigidity); 
	Filler.AddObject2beFilled(CountsHeBetaTOF,GetRigidity,GetRigidity); 
	Filler.AddObject2beFilled(CountsHeBetaNaF,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(CountsHeBetaAgl,GetRigidity,GetRigidity); 
	Filler.AddObject2beFilled(CountsHeFitTOF,GetBetaTOF ,GetBetaTOF ); 
	Filler.AddObject2beFilled(CountsHeFitNaF,GetBetaRICH,GetBetaRICH); 
	Filler.AddObject2beFilled(CountsHeFitAgl,GetBetaRICH,GetBetaRICH); 
	

	Filler.ReinitializeAll(refill);


}

#endif

