#ifndef FLUXES_H
#define FLUXES_H

#include "Analyzer.h"
#include "Efficiency.h"
#include "AllRangesEfficiency.h"
#include "Flux.h"

//MC parameters taken from /cvmfs/ams.cern.ch/Offline/AMSDataDir/DataManagement/DataSetsDesc

void Analyzer::BookFluxAnalysis(FileSaver finalhistos, FileSaver finalresults, bool refill){

	bool checkfile = finalhistos.CheckFile();
	check_file = checkfile;
	cout<<"****************************** BINS ***************************************"<<endl;
        SetUpUsualBinning();
 
	cout<<"****************************** FLUXES EVALUATION ******************************************"<<endl;

	Flux * DFluxTOF  = new Flux(finalhistos,finalresults, "DFluxTOF", "Acceptance_DTOF","Acceptance","TOFDfits/Fit Results/Deuteron Counts","TOFDfits/Fit Results/StatErrorD","ExposureTOF",Global.GetToFDBins());
	Flux * DFluxNaF  = new Flux(finalhistos,finalresults, "DFluxNaF", "Acceptance_DNaF","Acceptance","NaFDfits/Fit Results/Deuteron Counts","NaFDfits/Fit Results/StatErrorD","ExposureNaF",Global.GetNaFDBins());
	Flux * DFluxAgl  = new Flux(finalhistos,finalresults, "DFluxAgl", "Acceptance_DAgl","Acceptance","AglDfits/Fit Results/Deuteron Counts","AglDfits/Fit Results/StatErrorD","ExposureAgl",Global.GetAglDBins());

	Flux * HEPFlux   = new Flux(finalhistos,finalresults,"PFluxHE"  ,"Acceptance_HE"    ,"Acceptance","HEPCounts/HEPCounts/HEPCounts_after","","HEExposure"		,PRB);
	Flux * HEPFluxL1 = new Flux(finalhistos,finalresults,"PFluxL1HE","Acceptance_L1HE"  ,"Acceptance","HEPCountsL1/HEPCountsL1/HEPCountsL1_before","","HEExposure"	,PRB);
	Flux * HEPFluxQ  = new Flux(finalhistos,finalresults,"PFluxQHE" ,"Acceptance_QualHE","Acceptance","HEPCountsQual/HEPCountsQual/HEPCountsQual_after","","HEExposure",PRB);

	Flux * PFluxTOF  = new Flux(finalhistos,finalresults, "PFluxTOF", "Acceptance_PTOF","Acceptance","TOFPfits/Fit Results/Proton Counts","TOFPfits/Fit Results/StatErrorP","ExposureTOF",Global.GetToFPBins());
	Flux * PFluxNaF  = new Flux(finalhistos,finalresults, "PFluxNaF", "Acceptance_PNaF","Acceptance","NaFPfits/Fit Results/Proton Counts","NaFPfits/Fit Results/StatErrorP","ExposureNaF",Global.GetNaFPBins());
	Flux * PFluxAgl  = new Flux(finalhistos,finalresults, "PFluxAgl", "Acceptance_PAgl","Acceptance","AglPfits/Fit Results/Proton Counts","AglPfits/Fit Results/StatErrorP","ExposureAgl",Global.GetAglPBins());
	
	Flux * RigPTOF = new Flux(finalhistos,finalresults, "RigPTOF", "Acceptance_RigPTOF","Acceptance","TOFPCounts/TOFPCounts/TOFPCounts_before","","ExposureTOF",GlobalRig.GetToFPBins());
	Flux * RigPNaF = new Flux(finalhistos,finalresults, "RigPNaF", "Acceptance_RigPNaF","Acceptance","NaFPCounts/NaFPCounts/NaFPCounts_before","","ExposureNaF",GlobalRig.GetNaFPBins());
	Flux * RigPAgl = new Flux(finalhistos,finalresults, "RigPAgl", "Acceptance_RigPAgl","Acceptance","AglPCounts/AglPCounts/AglPCounts_before","","ExposureAgl",GlobalRig.GetAglPBins());

	Flux * RigBetaPTOF = new Flux(finalhistos,finalresults, "RigBetaPTOF", "Acceptance_PTOF","Acceptance","TOFPCountsBeta/TOFPCountsBeta/TOFPCountsBeta_before","","ExposureTOF",GlobalRig.GetToFPBins());
	Flux * RigBetaPNaF = new Flux(finalhistos,finalresults, "RigBetaPNaF", "Acceptance_PNaF","Acceptance","NaFPCountsBeta/NaFPCountsBeta/NaFPCountsBeta_before","","ExposureNaF",GlobalRig.GetNaFPBins());
	Flux * RigBetaPAgl = new Flux(finalhistos,finalresults, "RigBetaPAgl", "Acceptance_PAgl","Acceptance","AglPCountsBeta/AglPCountsBeta/AglPCountsBeta_before","","ExposureAgl",GlobalRig.GetAglPBins());

	Flux * FitBetaPTOF = new Flux(finalhistos,finalresults, "FitBetaPTOF", "Acceptance_PTOF","Acceptance","TOFPCountsFit/TOFPCountsFit/TOFPCountsFit_before","","ExposureTOF",Global.GetToFPBins());
	Flux * FitBetaPNaF = new Flux(finalhistos,finalresults, "FitBetaPNaF", "Acceptance_PNaF","Acceptance","NaFPCountsFit/NaFPCountsFit/NaFPCountsFit_before","","ExposureNaF",Global.GetNaFPBins());
	Flux * FitBetaPAgl = new Flux(finalhistos,finalresults, "FitBetaPAgl", "Acceptance_PAgl","Acceptance","AglPCountsFit/AglPCountsFit/AglPCountsFit_before","","ExposureAgl",Global.GetAglPBins());



	DFluxTOF ->SetDefaultOutFile(finalhistos); 
	DFluxNaF ->SetDefaultOutFile(finalhistos);
	DFluxAgl ->SetDefaultOutFile(finalhistos);

	HEPFlux  ->SetDefaultOutFile(finalhistos);
	HEPFluxL1 ->SetDefaultOutFile(finalhistos);
	HEPFluxQ ->SetDefaultOutFile(finalhistos);

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
	Filler_RTI.AddObject2beFilled(RigPTOF,GetBetaGen,GetBetaGen,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(RigPNaF,GetBetaGen,GetBetaGen,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(RigPAgl,GetBetaGen,GetBetaGen,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(RigBetaPTOF,GetBetaGen,GetBetaGen,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(RigBetaPNaF,GetBetaGen,GetBetaGen,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(RigBetaPAgl,GetBetaGen,GetBetaGen,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(FitBetaPTOF,GetBetaGen,GetBetaGen,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(FitBetaPNaF,GetBetaGen,GetBetaGen,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(FitBetaPAgl,GetBetaGen,GetBetaGen,"IsProtonMC");
	Filler_RTI.ReinitializeAll(refill);


	if(!refill&&checkfile) {
	float corr_D = 0.8618;
	float corr_p = 1/0.89510935;
	cout<<"************* D FLUX ************"<<endl;
		
		DFluxTOF-> Eval_Flux(corr_D,1.25,3.2,6,0.01);
		DFluxNaF-> Eval_Flux(corr_D,3.4,8.2,6,0.01);
		DFluxAgl-> Eval_Flux(corr_D);//,8,17.5,6,0.01);


		DFluxTOF->SaveResults(finalresults);
		DFluxNaF->SaveResults(finalresults);
		DFluxAgl->SaveResults(finalresults);

		cout<<"************* P FLUX ************"<<endl;

	/*	PFluxTOF-> Eval_Flux(1,1.2,1.5,3,0.01);
		PFluxNaF-> Eval_Flux(1,1.5,4.0,7,0.1);
		PFluxAgl-> Eval_Flux(1,4.0,8.9,7,0.1);

		PFluxTOF->SaveResults(finalresults);
		PFluxNaF->SaveResults(finalresults);
		PFluxAgl->SaveResults(finalresults);
*/
		cout<<"************* P+D FLUX ************"<<endl;

		FitBetaPTOF-> Eval_Flux(corr_p,1.2,1.5,3,0.01);
		FitBetaPNaF-> Eval_Flux(corr_p,1.5,4.0,7,0.1);
		FitBetaPAgl-> Eval_Flux(corr_p,4.0,8.9,7,0.1,5);

		FitBetaPTOF->ChangeName("FitBetaPTOF");
		FitBetaPNaF->ChangeName("FitBetaPNaF");
		FitBetaPAgl->ChangeName("FitBetaPAgl");

		RigBetaPTOF-> Eval_Flux(corr_p,1.2,1.5,3,0.01);
		RigBetaPNaF-> Eval_Flux(corr_p,1.5,4.0,7,0.1);
		RigBetaPAgl-> Eval_Flux(corr_p,4.0,8.9,7,0.1,5);

		RigBetaPTOF->ChangeName("RigBetaPTOF");
		RigBetaPNaF->ChangeName("RigBetaPNaF");
		RigBetaPAgl->ChangeName("RigBetaPAgl");

		RigPTOF-> Eval_Flux(corr_p);
		RigPNaF-> Eval_Flux(corr_p);
		RigPAgl-> Eval_Flux(corr_p);

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

		HEPFlux ->SaveResults(finalresults);
		HEPFluxL1 ->SaveResults(finalresults);
		HEPFluxQ->SaveResults(finalresults);


		cout<<"************** BINNING *****************"<<endl;

		cout<<"************* MERGING RANGES ************"<<endl;

		TH1F * RigP         = Global.MergeSubDResult_P(RigPTOF->GetFlux_rig(),RigPNaF->GetFlux_rig(),RigPAgl->GetFlux_rig());
		TH1F * RigBetaP     = Global.MergeSubDResult_P(RigBetaPTOF->GetFlux_rig(),RigBetaPNaF->GetFlux_rig(),RigBetaPAgl->GetFlux_rig());
		TH1F * FitBetaP     = Global.MergeSubDResult_P(FitBetaPTOF->GetFlux_rig(),FitBetaPNaF->GetFlux_rig(),FitBetaPAgl->GetFlux_rig());
		//TH1F * PFlux        = Global.MergeSubDResult_P(PFluxTOF->GetFlux_rig(),PFluxNaF->GetFlux_rig(),PFluxAgl->GetFlux_rig());
		TH1F * DFlux        = Global.MergeSubDResult_D(DFluxTOF->GetFlux_rig(),DFluxNaF->GetFlux_rig(),DFluxAgl->GetFlux_rig());
	
		TH1F * RigP_stat         = Global.MergeSubDResult_P(RigPTOF->GetStatError(),RigPNaF->GetStatError(),RigPAgl->GetStatError());
		TH1F * RigBetaP_stat     = Global.MergeSubDResult_P(RigBetaPTOF->GetStatError(),RigBetaPNaF->GetStatError(),RigBetaPAgl->GetStatError());
		TH1F * FitBetaP_stat     = Global.MergeSubDResult_P(FitBetaPTOF->GetStatError(),FitBetaPNaF->GetStatError(),FitBetaPAgl->GetStatError());
	//	TH1F * PFlux_stat        = Global.MergeSubDResult_P(PFluxTOF->GetStatError(),PFluxNaF->GetStatError(),PFluxAgl->GetStatError());
		TH1F * DFlux_stat        = Global.MergeSubDResult_D(DFluxTOF->GetStatError(),DFluxNaF->GetStatError(),DFluxAgl->GetStatError());
		
		TH1F * RigBetaP_unf = Global.MergeSubDResult_P(RigBetaPTOF->GetFlux_unf(),RigBetaPNaF->GetFlux_unf(),RigBetaPAgl->GetFlux_unf());
		TH1F * FitBetaP_unf = Global.MergeSubDResult_P(FitBetaPTOF->GetFlux_unf(),FitBetaPNaF->GetFlux_unf(),FitBetaPAgl->GetFlux_unf());
		//TH1F * PFlux_unf    = Global.MergeSubDResult_P(PFluxTOF->GetFlux_unf(),PFluxNaF->GetFlux_unf(),PFluxAgl->GetFlux_unf());
		TH1F * DFlux_unf    = Global.MergeSubDResult_D(DFluxTOF->GetFlux_unf(),DFluxNaF->GetFlux_unf(),DFluxAgl->GetFlux_unf());
		

		RigP_stat         = ConvertBinnedHisto( RigP_stat     , "RigP_stat"     ,Global.GetGlobalPBins(),false);   
               	RigBetaP_stat     = ConvertBinnedHisto( RigBetaP_stat     , "RigBetaP_stat"     ,Global.GetGlobalPBins(),false);   
                FitBetaP_stat     = ConvertBinnedHisto( FitBetaP_stat     , "FitBetaP_stat"     ,Global.GetGlobalPBins(),false);   
              //  PFlux_stat        = ConvertBinnedHisto( PFlux_stat        , "PFlux_stat"        ,Global.GetGlobalPBins(),false);   
                DFlux_stat        = ConvertBinnedHisto( DFlux_stat        , "DFlux_stat"        ,Global.GetGlobalDBins(),false);   
        
		RigP          = ConvertBinnedHisto( RigP         , "RigP"     ,Global.GetGlobalPBins(),false);   
        	RigBetaP      = ConvertBinnedHisto( RigBetaP     , "RigBetaP"     ,Global.GetGlobalPBins(),false);   
                FitBetaP      = ConvertBinnedHisto( FitBetaP     , "FitBetaP"     ,Global.GetGlobalPBins(),false);   
               // PFlux         = ConvertBinnedHisto( PFlux        , "PFlux"        ,Global.GetGlobalPBins(),false);   
                DFlux         = ConvertBinnedHisto( DFlux        , "DFlux"        ,Global.GetGlobalDBins(),false);   
                                                                               
                RigBetaP_unf  = ConvertBinnedHisto( RigBetaP_unf , "RigBetaP_unf"  ,Global.GetGlobalPBins(),false);   
                FitBetaP_unf  = ConvertBinnedHisto( FitBetaP_unf , "FitBetaP_unf"  ,Global.GetGlobalPBins(),false);   
               // PFlux_unf     = ConvertBinnedHisto( PFlux_unf    , "PFlux_unf"     ,Global.GetGlobalPBins(),false);   
 		DFlux_unf     = ConvertBinnedHisto( DFlux_unf    , "DFlux_unf"     ,Global.GetGlobalDBins(),false);   
               

		finalresults.Add(RigP);
		finalresults.Add(RigP_stat);
		finalresults.Add(RigBetaP);
		finalresults.Add(RigBetaP_unf);
		finalresults.Add(RigBetaP_stat);
		finalresults.Add(FitBetaP);
		finalresults.Add(FitBetaP_unf);
		finalresults.Add(FitBetaP_stat);
	//	finalresults.Add(PFlux);
	//	finalresults.Add(PFlux_unf);
	//	finalresults.Add(PFlux_stat);
		finalresults.Add(DFlux);
		finalresults.Add(DFlux_unf);
		finalresults.Add(DFlux_stat);
	




		finalresults.writeObjsInFolder("Fluxes/");	


/*		cout<<"************* MERGING RANGES ************"<<endl;
		
		TH1F * MergedRange_D = Global.MergeSubDResult_D(DFluxTOF->GetFlux(),DFluxNaF->GetFlux(),DFluxAgl->GetFlux());
		TH1F * MergedRange_P = Global.MergeSubDResult_P(PFluxTOF->GetFlux(),PFluxNaF->GetFlux(),PFluxAgl->GetFlux());

		SetUpRigTOIBinning();		
		TH1F * MergedRange_rig_D = Global.MergeSubDResult_D(DFluxTOF->GetFlux_rig(),DFluxNaF->GetFlux_rig(),DFluxAgl->GetFlux_rig());
		TH1F * MergedRange_rig_P = Global.MergeSubDResult_P(PFluxTOF->GetFlux_rig(),PFluxNaF->GetFlux_rig(),PFluxAgl->GetFlux_rig());

		MergedRange_D-> SetTitle("MergedRange_D_Ekin");
		MergedRange_P-> SetTitle("MergedRange_P_Ekin");
		MergedRange_D-> SetName("MergedRange_D_Ekin");
		MergedRange_P-> SetName("MergedRange_P_Ekin");

		TH1F * Ratio_R    = Global.MergedRatio(MergedRange_rig_D,MergedRange_rig_P); 
		TH1F * Ratio_Ekin = Global.MergedRatio_Ekin(MergedRange_D,MergedRange_P); 
	
		finalresults.Add(MergedRange_D);
		finalresults.Add(MergedRange_rig_D);
		finalresults.Add(MergedRange_P);
		finalresults.Add(MergedRange_rig_P);
		finalresults.Add(Ratio_R);
		finalresults.Add(Ratio_Ekin);
*/
		finalresults.writeObjsInFolder("Fluxes/");	
	}
	return;
}

#endif
