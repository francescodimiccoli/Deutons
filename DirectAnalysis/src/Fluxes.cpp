#ifndef FLUXES_H
#define FLUXES_H

#include "Analyzer.h"
#include "Efficiency.h"
#include "AllRangesEfficiency.h"
#include "Flux.h"
#include "ResultMerger.h"


//MC parameters taken from /cvmfs/ams.cern.ch/Offline/AMSDataDir/DataManagement/DataSetsDesc
TH1F * DoRatio(TH1F* numerator, TH1F * denominator,std::string name);
TH1F * DoRatio_Extended_subtraction(TH1F* numerator, TH1F * denominator, TH1F* extended,std::string name);



void Analyzer::BookFluxAnalysis(FileSaver finalhistos, FileSaver finalresults, bool refill){

	bool checkfile = finalhistos.CheckFile();
	check_file = checkfile;
	cout<<"****************************** BINS ***************************************"<<endl;
        SetUpUsualBinning();
 
	cout<<"****************************** FLUXES EVALUATION ******************************************"<<endl;

	Flux * DFluxTOF  = new Flux(finalhistos,finalresults, "DFluxTOF", "Acceptance_DTOF","Acceptance","TOFDfits/Fit Results/Deuteron Counts","TOFDfits/Fit Results/StatErrorD","ExposureTOF",Global.GetToFDBins());
	Flux * DFluxNaF  = new Flux(finalhistos,finalresults, "DFluxNaF", "Acceptance_DNaF","Acceptance","NaFDfits/Fit Results/Deuteron Counts","NaFDfits/Fit Results/StatErrorD","ExposureNaF",Global.GetNaFDBins());
	Flux * DFluxAgl  = new Flux(finalhistos,finalresults, "DFluxAgl", "Acceptance_DAgl","Acceptance","AglDfits/Fit Results/Deuteron Counts","AglDfits/Fit Results/StatErrorD","ExposureAgl",Global.GetAglDBins());


	Flux * HEPFlux   = new Flux(finalhistos,finalresults,"PFluxHE"  ,"Acceptance_HE"    ,"Acceptance","HEPCounts/HEPCounts/HEPCounts_before","","HEExposure"		,PRB);
	Flux * HEPFluxL1 = new Flux(finalhistos,finalresults,"PFluxL1HE","Acceptance_L1HE"  ,"Acceptance","HEPCountsL1/HEPCountsL1/HEPCountsL1_before","","HEExposure"	,PRB);
	Flux * HEPFluxQ  = new Flux(finalhistos,finalresults,"PFluxQHE" ,"Acceptance_QualHE","Acceptance","HEPCountsQual/HEPCountsQual/HEPCountsQual_before","","HEExposure",PRB);
	Flux * HEHeFluxQ  = new Flux(finalhistos,finalresults,"HeFluxQHE" ,"Acceptance_HeQualHE","Acceptance","HEHeCountsQual/HEHeCountsQual/HEHeCountsQual_before","","HEExposure",HeRB);


	Flux * PFluxTOF  = new Flux(finalhistos,finalresults, "PFluxTOF", "Acceptance_PTOF","Acceptance","TOFPfits/Fit Results/Proton Counts","TOFPfits/Fit Results/StatErrorP","ExposureTOF",Global.GetToFPBins());
	Flux * PFluxNaF  = new Flux(finalhistos,finalresults, "PFluxNaF", "Acceptance_PNaF","Acceptance","NaFPfits/Fit Results/Proton Counts","NaFPfits/Fit Results/StatErrorP","ExposureNaF",Global.GetNaFPBins());
	Flux * PFluxAgl  = new Flux(finalhistos,finalresults, "PFluxAgl", "Acceptance_PAgl","Acceptance","AglPfits/Fit Results/Proton Counts","AglPfits/Fit Results/StatErrorP","ExposureAgl",Global.GetAglPBins());
	
	Flux * RigPTOF = new Flux(finalhistos,finalresults, "RigPTOF", "Acceptance_RigPTOF","Acceptance","TOFPCounts/TOFPCounts/TOFPCounts_before","","ExposureTOF",GlobalRig.GetToFPBins());
	Flux * RigPNaF = new Flux(finalhistos,finalresults, "RigPNaF", "Acceptance_RigPNaF","Acceptance","NaFPCounts/NaFPCounts/NaFPCounts_before","","ExposureNaF",GlobalRig.GetNaFPBins());
	Flux * RigPAgl = new Flux(finalhistos,finalresults, "RigPAgl", "Acceptance_RigPAgl","Acceptance","AglPCounts/AglPCounts/AglPCounts_before","","ExposureAgl",GlobalRig.GetAglPBins());

	Flux * RigBetaPTOF = new Flux(finalhistos,finalresults, "RigBetaPTOF", "Acceptance_RigBetaPTOF","Acceptance","TOFPCountsBeta/TOFPCountsBeta/TOFPCountsBeta_before","","ExposureTOF",GlobalRig.GetToFPBins());
	Flux * RigBetaPNaF = new Flux(finalhistos,finalresults, "RigBetaPNaF", "Acceptance_RigBetaPNaF","Acceptance","NaFPCountsBeta/NaFPCountsBeta/NaFPCountsBeta_before","","ExposureNaF",GlobalRig.GetNaFPBins());
	Flux * RigBetaPAgl = new Flux(finalhistos,finalresults, "RigBetaPAgl", "Acceptance_RigBetaPAgl","Acceptance","AglPCountsBeta/AglPCountsBeta/AglPCountsBeta_before","","ExposureAgl",GlobalRig.GetAglPBins());

	Flux * FitBetaPTOF = new Flux(finalhistos,finalresults, "FitBetaPTOF", "Acceptance_PTOF","Acceptance","TOFPCountsFit/TOFPCountsFit/TOFPCountsFit_before","","ExposureTOF",Global.GetToFPBins());
	Flux * FitBetaPNaF = new Flux(finalhistos,finalresults, "FitBetaPNaF", "Acceptance_PNaF","Acceptance","NaFPCountsFit/NaFPCountsFit/NaFPCountsFit_before","","ExposureNaF",Global.GetNaFPBins());
	Flux * FitBetaPAgl = new Flux(finalhistos,finalresults, "FitBetaPAgl", "Acceptance_PAgl","Acceptance","AglPCountsFit/AglPCountsFit/AglPCountsFit_before","","ExposureAgl",Global.GetAglPBins());


	Flux * RigHeTOF = new Flux(finalhistos,finalresults, "RigHeTOF", "Acceptance_RigHeTOF","Acceptance","TOFHeCounts/TOFHeCounts/TOFHeCounts_before","","ExposureTOF",Global_HeRig.GetToFDBins());
	Flux * RigHeNaF = new Flux(finalhistos,finalresults, "RigHeNaF", "Acceptance_RigHeNaF","Acceptance","NaFHeCounts/NaFHeCounts/NaFHeCounts_before","","ExposureNaF",Global_HeRig.GetNaFDBins());
	Flux * RigHeAgl = new Flux(finalhistos,finalresults, "RigHeAgl", "Acceptance_RigHeAgl","Acceptance","AglHeCounts/AglHeCounts/AglHeCounts_before","","ExposureAgl",Global_HeRig.GetAglDBins());

	Flux * RigBetaHeTOF = new Flux(finalhistos,finalresults, "RigBetaHeTOF", "Acceptance_HeTOF","Acceptance","TOFHeBetaCountsBeta/TOFHeBetaCounts/TOFHeBetaCounts_before","","ExposureTOF",Global_HeRig.GetToFDBins());
	Flux * RigBetaHeNaF = new Flux(finalhistos,finalresults, "RigBetaHeNaF", "Acceptance_HeNaF","Acceptance","NaFHeBetaCountsBeta/NaFHeBetaCounts/NaFHeBetaCounts_before","","ExposureNaF",Global_HeRig.GetNaFDBins());
	Flux * RigBetaHeAgl = new Flux(finalhistos,finalresults, "RigBetaHeAgl", "Acceptance_HeAgl","Acceptance","AglHeBetaCountsBeta/AglHeBetaCounts/AglHeBetaCounts_before","","ExposureAgl",Global_HeRig.GetAglDBins());

	Flux * FitBetaHeTOF = new Flux(finalhistos,finalresults, "FitBetaHeTOF", "Acceptance_HeTOF","Acceptance","TOFHeFitCountsBeta/TOFHeFitCounts/TOFHeFitCounts_before","","ExposureTOF",Global_He.GetToFDBins());
	Flux * FitBetaHeNaF = new Flux(finalhistos,finalresults, "FitBetaHeNaF", "Acceptance_HeNaF","Acceptance","NaFHeFitCountsBeta/NaFHeFitCounts/NaFHeFitCounts_before","","ExposureNaF",Global_He.GetNaFDBins());
	Flux * FitBetaHeAgl = new Flux(finalhistos,finalresults, "FitBetaHeAgl", "Acceptance_HeAgl","Acceptance","AglHeFitCountsBeta/AglHeFitCounts/AglHeFitCounts_before","","ExposureAgl",Global_He.GetAglDBins());


	Flux * He3FluxTOF  = new Flux(finalhistos,finalresults, "He3FluxTOF", "Acceptance_HeTOF","Acceptance","TOFHefits/Fit Results/Deuteron Counts","TOFHefits/Fit Results/StatErrorD","ExposureTOF",Global_He.GetToFDBins());
	Flux * He3FluxNaF  = new Flux(finalhistos,finalresults, "He3FluxNaF", "Acceptance_HeNaF","Acceptance","NaFHefits/Fit Results/Deuteron Counts","NaFHefits/Fit Results/StatErrorD","ExposureNaF",Global_He.GetNaFDBins());
	Flux * He3FluxAgl  = new Flux(finalhistos,finalresults, "He3FluxAgl", "Acceptance_HeAgl","Acceptance","AglHEfits/Fit Results/Deuteron Counts","AglHEfits/Fit Results/StatErrorD","ExposureAgl",Global_He.GetAglDBins());

	Flux * He4FluxTOF  = new Flux(finalhistos,finalresults, "He4FluxTOF", "Acceptance_HeTOF","Acceptance","TOFHefits/Fit Results/Proton Counts","TOFHefits/Fit Results/StatErrorP","ExposureTOF",Global_He.GetToFDBins());
	Flux * He4FluxNaF  = new Flux(finalhistos,finalresults, "He4FluxNaF", "Acceptance_HeNaF","Acceptance","NaFHefits/Fit Results/Proton Counts","NaFHefits/Fit Results/StatErrorP","ExposureNaF",Global_He.GetNaFDBins());
	Flux * He4FluxAgl  = new Flux(finalhistos,finalresults, "He4FluxAgl", "Acceptance_HeAgl","Acceptance","AglHEfits/Fit Results/Proton Counts","AglHEfits/Fit Results/StatErrorP","ExposureAgl",Global_He.GetAglDBins());


	DFluxTOF ->SetDefaultOutFile(finalhistos); 
	DFluxNaF ->SetDefaultOutFile(finalhistos);
	DFluxAgl ->SetDefaultOutFile(finalhistos);

	HEPFlux  ->SetDefaultOutFile(finalhistos);
	HEPFluxL1 ->SetDefaultOutFile(finalhistos);
	HEPFluxQ ->SetDefaultOutFile(finalhistos);
	HEHeFluxQ ->SetDefaultOutFile(finalhistos);


	PFluxTOF ->SetDefaultOutFile(finalhistos);
	PFluxNaF ->SetDefaultOutFile(finalhistos);
	PFluxAgl ->SetDefaultOutFile(finalhistos);

	RigPTOF->SetDefaultOutFile(finalhistos);
	RigPNaF->SetDefaultOutFile(finalhistos);
	RigPAgl->SetDefaultOutFile(finalhistos);

	RigBetaPTOF->SetDefaultOutFile(finalhistos);
	RigBetaPNaF->SetDefaultOutFile(finalhistos);
	RigBetaPAgl->SetDefaultOutFile(finalhistos);

	FitBetaPTOF->SetDefaultOutFile(finalhistos);
	FitBetaPNaF->SetDefaultOutFile(finalhistos);
	FitBetaPAgl->SetDefaultOutFile(finalhistos);

	RigHeTOF->SetDefaultOutFile(finalhistos);
	RigHeNaF->SetDefaultOutFile(finalhistos);
	RigHeAgl->SetDefaultOutFile(finalhistos);

	RigBetaHeTOF->SetDefaultOutFile(finalhistos);
	RigBetaHeNaF->SetDefaultOutFile(finalhistos);
	RigBetaHeAgl->SetDefaultOutFile(finalhistos);

	FitBetaHeTOF->SetDefaultOutFile(finalhistos);
	FitBetaHeNaF->SetDefaultOutFile(finalhistos);
	FitBetaHeAgl->SetDefaultOutFile(finalhistos);


	He3FluxTOF ->SetDefaultOutFile(finalhistos); 
	He3FluxNaF ->SetDefaultOutFile(finalhistos);
	He3FluxAgl ->SetDefaultOutFile(finalhistos);

	He4FluxTOF ->SetDefaultOutFile(finalhistos); 
	He4FluxNaF ->SetDefaultOutFile(finalhistos);
	He4FluxAgl ->SetDefaultOutFile(finalhistos);


	cout<<"********** EXPOSURE TIME & GEOM. ACCEPT. ********"<<endl;

	Filler_RTI.AddObject2beFilled(DFluxTOF,GetBetaGen,GetBetaGen,"IsDeutonMC");
	Filler_RTI.AddObject2beFilled(DFluxNaF,GetBetaGen,GetBetaGen,"IsDeutonMC");
	Filler_RTI.AddObject2beFilled(DFluxAgl,GetBetaGen,GetBetaGen,"IsDeutonMC");
	Filler_RTI.AddObject2beFilled(PFluxTOF,GetBetaGen,GetBetaGen,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(PFluxNaF,GetBetaGen,GetBetaGen,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(PFluxAgl,GetBetaGen,GetBetaGen,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(HEPFlux,GetGenMomentum,GetGenMomentum,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(HEPFluxL1,GetGenMomentum,GetGenMomentum,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(HEPFluxQ,GetGenMomentum,GetGenMomentum,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(HEHeFluxQ,GetGenMomentum,GetGenMomentum,"IsHeliumMC");
	Filler_RTI.AddObject2beFilled(RigPTOF,GetBetaGen,GetBetaGen,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(RigPNaF,GetBetaGen,GetBetaGen,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(RigPAgl,GetBetaGen,GetBetaGen,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(RigBetaPTOF,GetBetaGen,GetBetaGen,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(RigBetaPNaF,GetBetaGen,GetBetaGen,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(RigBetaPAgl,GetBetaGen,GetBetaGen,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(FitBetaPTOF,GetBetaGen,GetBetaGen,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(FitBetaPNaF,GetBetaGen,GetBetaGen,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(FitBetaPAgl,GetBetaGen,GetBetaGen,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(RigHeTOF,GetBetaGen,GetBetaGen,"IsHeliumMC");
	Filler_RTI.AddObject2beFilled(RigHeNaF,GetBetaGen,GetBetaGen,"IsHeliumMC");
	Filler_RTI.AddObject2beFilled(RigHeAgl,GetBetaGen,GetBetaGen,"IsHeliumMC");
	Filler_RTI.AddObject2beFilled(RigBetaHeTOF,GetBetaGen,GetBetaGen,"IsHeliumMC");
	Filler_RTI.AddObject2beFilled(RigBetaHeNaF,GetBetaGen,GetBetaGen,"IsHeliumMC");
	Filler_RTI.AddObject2beFilled(RigBetaHeAgl,GetBetaGen,GetBetaGen,"IsHeliumMC");
	Filler_RTI.AddObject2beFilled(FitBetaHeTOF,GetBetaGen,GetBetaGen,"IsHeliumMC");
	Filler_RTI.AddObject2beFilled(FitBetaHeNaF,GetBetaGen,GetBetaGen,"IsHeliumMC");
	Filler_RTI.AddObject2beFilled(FitBetaHeAgl,GetBetaGen,GetBetaGen,"IsHeliumMC");
	Filler_RTI.AddObject2beFilled(He3FluxTOF,GetBetaGen,GetBetaGen,"IsHeliumMC");
	Filler_RTI.AddObject2beFilled(He3FluxNaF,GetBetaGen,GetBetaGen,"IsHeliumMC");
	Filler_RTI.AddObject2beFilled(He3FluxAgl,GetBetaGen,GetBetaGen,"IsHeliumMC");
	Filler_RTI.AddObject2beFilled(He4FluxTOF,GetBetaGen,GetBetaGen,"IsHeliumMC");
	Filler_RTI.AddObject2beFilled(He4FluxNaF,GetBetaGen,GetBetaGen,"IsHeliumMC");
	Filler_RTI.AddObject2beFilled(He4FluxAgl,GetBetaGen,GetBetaGen,"IsHeliumMC");
	Filler_RTI.ReinitializeAll(refill);


	if(!refill&&checkfile) {
	float corr_D  = 1;//0.8618;
	float corr_p  = 1.16;
	float corr_he = 1;


	TFile * UnfoldingFile = TFile::Open("/afs/cern.ch/work/f/fdimicco/private/Deutons/DirectAnalysis/EffSyst/Time.root");
	TH2F* avg; 	

	cout<<"************* D FLUX ************"<<endl;
		
	if(UnfoldingFile){	
		avg  = (TH2F*)UnfoldingFile->Get("/Unfolding/Deuteron/sTOF/unfoldingTOFD_timeavg");
		DFluxTOF->Set_UnfoldingTime(avg,filelistDT.c_str());	
		avg  = (TH2F*)UnfoldingFile->Get("/Unfolding/Deuteron/sNaF/unfoldingNaFD_timeavg");
		DFluxNaF->Set_UnfoldingTime(avg,filelistDT.c_str());
		avg  = (TH2F*)UnfoldingFile->Get("/Unfolding/Deuteron/sAgl/unfoldingAglD_timeavg");
		DFluxAgl->Set_UnfoldingTime(avg,filelistDT.c_str());
	}	
	
	//DFluxNaF->SetForceFolded(); 
        //DFluxAgl->SetForceFolded();  //sistemare unfolding agl 
      
	DFluxTOF-> Eval_Flux(corr_D,1.25,3.2,6,0.01);
	DFluxNaF-> Eval_Flux(corr_D,3.4,8.2,6,0.01);
	DFluxAgl-> Eval_Flux(corr_D,8,17.8,4,0.1,true);

	DFluxTOF->SaveResults(finalresults);
	DFluxNaF->SaveResults(finalresults);
	DFluxAgl->SaveResults(finalresults);

		cout<<"************* P FLUX ************"<<endl;

		avg  = (TH2F*)UnfoldingFile->Get("/Unfolding/Proton/sTOF/unfoldingTOFP_timeavg");
		PFluxTOF->Set_UnfoldingTime(avg,filelistDT.c_str());	
		avg  = (TH2F*)UnfoldingFile->Get("/Unfolding/Proton/sNaF/unfoldingNaFP_timeavg");
		PFluxNaF->Set_UnfoldingTime(avg,filelistDT.c_str());
        	avg  = (TH2F*)UnfoldingFile->Get("/Unfolding/Proton/sAgl/unfoldingAglP_timeavg");
		PFluxAgl->Set_UnfoldingTime(avg,filelistDT.c_str());

		PFluxTOF-> Eval_Flux(1);//,1.2,1.5,3,0.01);
		PFluxNaF-> Eval_Flux(1);//,1.5,4.0,7,0.1);
		PFluxAgl-> Eval_Flux(1);//,4.0,8.9,7,0.1);

		PFluxTOF->SaveResults(finalresults);
		PFluxNaF->SaveResults(finalresults);
		PFluxAgl->SaveResults(finalresults);

		cout<<"************* P+D FLUX ************"<<endl;

		if(UnfoldingFile){	
			avg  = (TH2F*)UnfoldingFile->Get("/Unfolding/Proton/sTOF/unfoldingTOFP_timeavg");
			FitBetaPTOF->Set_UnfoldingTime(avg,filelistDT.c_str());	
			avg  = (TH2F*)UnfoldingFile->Get("/Unfolding/Proton/sNaF/unfoldingNaFP_timeavg");
			FitBetaPNaF->Set_UnfoldingTime(avg,filelistDT.c_str());
			avg  = (TH2F*)UnfoldingFile->Get("/Unfolding/Proton/sAgl/unfoldingAglP_timeavg");
			FitBetaPAgl->Set_UnfoldingTime(avg,filelistDT.c_str());
		}

		FitBetaPTOF-> Eval_Flux(corr_p,1.2,1.57,3,0.01);
		FitBetaPNaF-> Eval_Flux(corr_p,1.5,4.0,7,0.1);
		FitBetaPAgl-> Eval_Flux(corr_p,4.0,8.9,7,0.1);

		FitBetaPTOF->ChangeName("FitBetaPTOF");
		FitBetaPNaF->ChangeName("FitBetaPNaF");
		FitBetaPAgl->ChangeName("FitBetaPAgl");

		RigBetaPTOF-> Eval_Flux(corr_p,1.2,1.57,3,0.01);
		RigBetaPNaF-> Eval_Flux(corr_p,1.5,4.0,7,0.1);
		RigBetaPAgl-> Eval_Flux(corr_p,4.0,8.9,7,0.1);

		RigBetaPTOF->ChangeName("RigBetaPTOF");
		RigBetaPNaF->ChangeName("RigBetaPNaF");
		RigBetaPAgl->ChangeName("RigBetaPAgl");

		RigPTOF-> Eval_Flux(corr_p,1.2,1.57,3,0.01);
		RigPNaF-> Eval_Flux(corr_p,1.5,4.0,7,0.1); ;
		RigPAgl-> Eval_Flux(corr_p,4.0,8.9,7,0.1); ;

		RigPTOF->ChangeName("RigPTOF");
		RigPNaF->ChangeName("RigPNaF");
		RigPAgl->ChangeName("RigPAgl");

		RigPTOF->SaveResults(finalresults);
		RigPNaF->SaveResults(finalresults);
		RigPAgl->SaveResults(finalresults);

		RigBetaPTOF->SaveResults(finalresults);
		RigBetaPNaF->SaveResults(finalresults);
		RigBetaPAgl->SaveResults(finalresults);
	
		FitBetaPTOF->SaveResults(finalresults);
		FitBetaPNaF->SaveResults(finalresults);
		FitBetaPAgl->SaveResults(finalresults);
	

		cout<<"************* HE P+D FLUX ************"<<endl;

		HEPFlux -> Eval_Flux(corr_p);
		HEPFluxL1 -> Eval_Flux(corr_p);
		HEPFluxQ -> Eval_Flux(corr_p);
		HEHeFluxQ -> Eval_Flux();


		HEPFlux ->SaveResults(finalresults);
		HEPFluxL1 ->SaveResults(finalresults);
		HEPFluxQ->SaveResults(finalresults);
		HEHeFluxQ->SaveResults(finalresults);


		cout<<"************* P FLUX ************"<<endl;

		PFluxTOF-> Eval_Flux(corr_p,1.2,1.57,3,0.01);
		PFluxNaF-> Eval_Flux(corr_p,1.5,4.0,7,0.1);
		PFluxAgl-> Eval_Flux(corr_p,4.0,8.9,7,0.1);


		PFluxTOF->SaveResults(finalresults);
		PFluxNaF->SaveResults(finalresults);
		PFluxAgl->SaveResults(finalresults);


		cout<<"************* AllHe FLUX **********"<<endl;
		RigHeTOF->Eval_Flux(corr_he);//,1.9,3.2,3,0.01);
		RigHeNaF->Eval_Flux(corr_he);//,3.4,8.2,6,0.01);
		RigHeAgl->Eval_Flux(corr_he);//,8,17.5,5,0.01);

		RigHeTOF->SaveResults(finalresults);
		RigHeNaF->SaveResults(finalresults);
		RigHeAgl->SaveResults(finalresults);

		cout<<"************* 4He FLUX ************"<<endl;

		He4FluxTOF-> Eval_Flux(corr_he);//,1.9,3.2,3,0.01);
		He4FluxNaF-> Eval_Flux(corr_he);//,3.4,8.2,6,0.01);
		He4FluxAgl-> Eval_Flux(corr_he);//,8,17.5,5,0.01);

		He4FluxTOF->SaveResults(finalresults);
		He4FluxNaF->SaveResults(finalresults);
		He4FluxAgl->SaveResults(finalresults);

		cout<<"************* 3He FLUX ************"<<endl;
		
		He3FluxTOF-> Eval_Flux(corr_he);//,1.9,3.2,3,0.01);
		He3FluxNaF-> Eval_Flux(corr_he);//,3.4,8.2,6,0.01);
		He3FluxAgl-> Eval_Flux(corr_he);//,8,17.5,5,0.01);

		He3FluxTOF->SaveResults(finalresults);
		He3FluxNaF->SaveResults(finalresults);
		He3FluxAgl->SaveResults(finalresults);


		cout<<"************** BINNING *****************"<<endl;

		cout<<"************* MERGING RANGES ************"<<endl;

		Particle proton(0.9382720813, 1, 1);  // proton mass 938 MeV
		Particle deuton(1.8756129   , 1, 2);  // deuterium mass 1876 MeV, Z=1, A=2

		ResultMerger*	DeutonRes   = new ResultMerger(finalresults,"DFlux",Global,DFluxTOF,DFluxNaF,DFluxAgl,deuton);
		ResultMerger*	RigPRes     = new ResultMerger(finalresults,"RigPFlux",Global,RigPTOF,RigPNaF,RigPAgl,proton);
		ResultMerger*	RigBetaPRes = new ResultMerger(finalresults,"RigBetaPFlux",Global,RigBetaPTOF,RigBetaPNaF,RigBetaPAgl,proton);
		ResultMerger*	FitBetaPRes = new ResultMerger(finalresults,"FitBetaPFlux",Global,FitBetaPTOF,FitBetaPNaF,FitBetaPAgl,proton);
		ResultMerger*	ProtonRes   = new ResultMerger(finalresults,"PFlux",Global,PFluxTOF,PFluxNaF,PFluxAgl,proton);
	
		ResultMerger*	He4Res 	    = new ResultMerger(finalresults,"He4Flux",Global_He,He4FluxTOF,He4FluxNaF,He4FluxAgl,deuton);
		ResultMerger*	He3Res 	    = new ResultMerger(finalresults,"He3Flux",Global_He,He3FluxTOF,He3FluxNaF,He3FluxAgl,deuton);
		ResultMerger*	RigHeRes    = new ResultMerger(finalresults,"RigHeFlux",Global_HeRig,RigHeTOF,RigHeNaF,RigHeAgl,deuton);
		

		DeutonRes->SaveResults(finalresults);
		ProtonRes->SaveResults(finalresults);
		RigPRes->SaveResults(finalresults);
		RigBetaPRes->SaveResults(finalresults);
		FitBetaPRes->SaveResults(finalresults);
		He4Res->SaveResults(finalresults);
		He3Res->SaveResults(finalresults);
		RigHeRes->SaveResults(finalresults);


	
		//////////////// RIG RATIO ////////////////////////

	
		//Counts
		TH1F * ratiocountsDP_ = DoRatio(DeutonRes->counts,ProtonRes->counts,"DPratiocounts");
		TH1F * ratiocountsDHe4_ = DoRatio(DeutonRes->counts,He4Res->counts,"DHe4ratiocounts");


		//Acceptance
		TH1F * ratioaccDP_ = DoRatio(DeutonRes->effAcc,ProtonRes->effAcc,"DPratioacc");
		TH1F * ratioaccDHe4_ = DoRatio(DeutonRes->effAcc,He4Res->effAcc,"DHe4ratioacc");


		//Unfolding
		TH1F * ratiounfDP_ = DoRatio(DeutonRes->unfolding,ProtonRes->unfolding,"DPratiounf");


		//Flux
		TH1F * ratioDP_ = DoRatio_Extended_subtraction(DeutonRes->flux,ProtonRes->flux,HEPFluxL1->GetFlux_rig(),"DPratio");


		//Flux_Unfolded
		TH1F * ratioDP_unf_ = DoRatio_Extended_subtraction(DeutonRes->flux_unf,ProtonRes->flux_unf,HEPFluxL1->GetFlux_unf(),"DPratio_unf");


		//Flux_Folded_stat

		TH1F * ratioDP_rig_stat_ = DoRatio_Extended_subtraction(DeutonRes->flux_stat,ProtonRes->flux_stat,HEPFluxL1->GetFlux_rig_stat(),"DPratio_rig_stat");


		//Flux_Unfolded_stat
		TH1F * ratioDP_unf_stat_ = DoRatio_Extended_subtraction(DeutonRes->flux_unf_stat,ProtonRes->flux_unf_stat,HEPFluxL1->GetFlux_unf_stat(),"DPratio_unf_stat");



		TH1F * ratioDHe4 	  = DoRatio(DeutonRes->flux,He4Res->flux,"DHe4ratio");
		TH1F * ratioDHe4_stat	  = DoRatio(DeutonRes->flux_stat,He4Res->flux_stat,"DHe4ratio_stat");
		TH1F * ratioDHe4_unf 	  = DoRatio(DeutonRes->flux_unf,He4Res->flux_unf,"DHe4ratio_unf");
		TH1F * ratioDHe4_unf_stat = DoRatio(DeutonRes->flux_unf_stat,He4Res->flux_unf_stat,"DHe4ratio_unf_stat");

		TH1F * ratioHe3He4 	  = DoRatio(He3Res->flux_ekin,He4Res->flux_ekin,"He3He4ratio");
		TH1F * ratioHe3He4_ekin_unf 	  = DoRatio(He3Res->flux_ekin_unf,He4Res->flux_ekin_unf,"He3He4ratio_ekin_unf");
		TH1F * ratioHe3He4_stat	  = DoRatio(He3Res->flux_stat,He4Res->flux_stat,"He3He4ratio_stat");
		TH1F * ratioHe3He4_unf 	  = DoRatio(He3Res->flux_unf,He4Res->flux_unf,"He3He4ratio_unf");
		TH1F * ratioHe3He4_unf_stat = DoRatio(He3Res->flux_unf_stat,He4Res->flux_unf_stat,"He3He4ratio_unf_stat");

		TH1F * ratioDAllHe 	  = DoRatio(DeutonRes->flux,RigHeRes->flux,"DAllHeratio");
		TH1F * ratioDAllHe_stat	  = DoRatio(DeutonRes->flux_stat,RigHeRes->flux_stat,"DAllHeratio_stat");
		TH1F * ratioDAllHe_unf 	  = DoRatio(DeutonRes->flux_unf,RigHeRes->flux_unf,"DAllHeratio_unf");
		TH1F * ratioDAllHe_unf_stat = DoRatio(DeutonRes->flux_unf_stat,RigHeRes->flux_unf_stat,"DAllHeratio_unf_stat");

  		//Flux Unfolding Test	
		TH1F * ratioPP = DoRatio(FitBetaPRes->flux,RigBetaPRes->flux,"unfoldingtest_ratiofolded");
		TH1F * ratioPP_stat = DoRatio(FitBetaPRes->flux_stat,RigBetaPRes->flux_stat,"unfoldingtest_ratiofolded_stat");
		TH1F * ratioPP_unf = DoRatio(FitBetaPRes->flux_unf,RigBetaPRes->flux_unf,"unfoldingtest_ratiounfolded");
		TH1F * ratioPP_unf_stat = DoRatio(FitBetaPRes->flux_unf_stat,RigBetaPRes->flux_unf_stat,"unfoldingtest_ratiounfolded_stat");





		finalresults.Add(ratioDP_);	
		finalresults.Add(ratioDP_unf_);	
		finalresults.Add(ratioDP_unf_stat_);	
		finalresults.Add(ratioDP_rig_stat_);	

		finalresults.Add(ratiocountsDP_);	
		finalresults.Add(ratioaccDP_);	
		finalresults.Add(ratiounfDP_);	
	
		finalresults.Add(ratiocountsDHe4_);	
		finalresults.Add(ratioaccDHe4_);	
	
		finalresults.Add(ratioDHe4);	
		finalresults.Add(ratioDHe4_unf);	
		finalresults.Add(ratioDHe4_stat);	
		finalresults.Add(ratioDHe4_unf_stat);	

		finalresults.Add(ratioHe3He4);	
		finalresults.Add(ratioHe3He4_ekin_unf);	
		finalresults.Add(ratioHe3He4_unf);	
		finalresults.Add(ratioHe3He4_stat);	
		finalresults.Add(ratioHe3He4_unf_stat);	
	
		finalresults.Add(ratioDAllHe);	
		finalresults.Add(ratioDAllHe_unf);	
		finalresults.Add(ratioDAllHe_stat);	
		finalresults.Add(ratioDAllHe_unf_stat);	

		finalresults.Add(ratioPP);
		finalresults.Add(ratioPP_stat);
		finalresults.Add(ratioPP_unf);
		finalresults.Add(ratioPP_unf_stat);



		finalresults.writeObjsInFolder("MergedRatios/");	



	}
	return;
}





















TH1F * DoRatio(TH1F* numerator, TH1F * denominator,std::string name){

	TH1F * ratio = (TH1F*) numerator->Clone(name.c_str());
	TH1F * rebinned_denominator = (TH1F*) numerator->Clone("rebinned_denom");
	rebinned_denominator->Reset();
	cout<<"RATIO CALC "<<name<<endl;
	for(int i=0;i<rebinned_denominator->GetNbinsX();i++) {
		rebinned_denominator->SetBinContent(i+1,denominator->GetBinContent(denominator->FindBin(numerator->GetBinCenter(i+1))));
		rebinned_denominator->SetBinError(i+1,denominator->GetBinError(denominator->FindBin(numerator->GetBinCenter(i+1))));
		cout<<i<<" "<<numerator->GetBinLowEdge(i+1)<<" "<<numerator->GetBinLowEdge(i+2)<<" "<<numerator->GetBinContent(i+1)<<" "<<denominator->GetBinContent(denominator->FindBin(numerator->GetBinCenter(i+1)))<<endl;	

	}
	ratio->Divide(rebinned_denominator);

	for(int i=0;i<rebinned_denominator->GetNbinsX();i++) {

		float error = pow(rebinned_denominator->GetBinError(i+1)/rebinned_denominator->GetBinContent(i+1),2);
		error+= pow(numerator->GetBinError(i+1)/numerator->GetBinContent(i+1),2);

	}


	return ratio;	


}


TH1F * DoRatio_Extended_subtraction(TH1F* numerator, TH1F * denominator, TH1F* extended,std::string name){

	TH1F * ratio = (TH1F*) numerator->Clone(name.c_str());
	TH1F * rebinned_denominator = (TH1F*) numerator->Clone("rebinned_denom");
	rebinned_denominator->Reset();
	for(int i=0;i<rebinned_denominator->GetNbinsX();i++) {

		 if(rebinned_denominator->GetBinLowEdge(i+1)<=denominator->GetBinLowEdge(denominator->GetNbinsX())){
			rebinned_denominator->SetBinContent(i+1,denominator->GetBinContent(denominator->FindBin(numerator->GetBinCenter(i+1)))-numerator->GetBinContent(i+1));
			rebinned_denominator->SetBinError(i+1,denominator->GetBinError(denominator->FindBin(numerator->GetBinCenter(i+1))));
		}
		else{
			rebinned_denominator->SetBinContent(i+1,extended->GetBinContent(extended->FindBin(rebinned_denominator->GetBinCenter(i+1)))-numerator->GetBinContent(i+1));
			rebinned_denominator->SetBinError(i+1,extended->GetBinError(extended->FindBin(rebinned_denominator->GetBinCenter(i+1))));
		}
	}
	ratio->Divide(rebinned_denominator);

	for(int i=0;i<rebinned_denominator->GetNbinsX();i++) {

		float error = pow(rebinned_denominator->GetBinError(i+1)/rebinned_denominator->GetBinContent(i+1),2);
		error+= pow(numerator->GetBinError(i+1)/numerator->GetBinContent(i+1),2);

	}


	return ratio;	
}


#endif
