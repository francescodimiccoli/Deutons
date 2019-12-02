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

//	Flux * FluxTests[10];

//	for(int i=0;i<10;i++) FluxTests[i] = new Flux(finalhistos,finalresults,("PFluxTests_v"+to_string(i)).c_str(),("CountsTestEff_rig_v" + to_string(i)).c_str(),("CountsTestEff_rig_v" + to_string(i)).c_str(),("HEPTests_v"+to_string(i) + "/HEPTests_v"+to_string(i)+"/HEPTests_v"+to_string(i)+"_before").c_str(),"HEExposure","For_Acceptance_P",PRB);

	Flux * DFluxTOF  = new Flux(finalhistos,finalresults, "DFluxTOF", "Acceptance_DTOF","Acceptance","TOFDfits/Fit Results/Primary Deuteron Counts","ExposureTOF",Global.GetToFDBins());
	Flux * DFluxNaF  = new Flux(finalhistos,finalresults, "DFluxNaF", "Acceptance_DNaF","Acceptance","NaFDfits/Fit Results/Primary Deuteron Counts","ExposureNaF",Global.GetNaFDBins());
	Flux * DFluxAgl  = new Flux(finalhistos,finalresults, "DFluxAgl", "Acceptance_DAgl","Acceptance","AglDfits/Fit Results/Primary Deuteron Counts","ExposureAgl",Global.GetAglDBins());

	Flux * HEPFlux   = new Flux(finalhistos,finalresults,"PFluxHE"  ,"Acceptance_HE"    ,"Acceptance","HEPCounts/HEPCounts/HEPCounts_before","HEExposure"		,PRB);
	Flux * HEPFluxL1 = new Flux(finalhistos,finalresults,"PFluxL1HE","Acceptance_L1HE"  ,"Acceptance","HEPCountsL1/HEPCountsL1/HEPCountsL1_before","HEExposure"	,PRB);
	Flux * HEPFluxQ  = new Flux(finalhistos,finalresults,"PFluxQHE" ,"Acceptance_QualHE","Acceptance","HEPCountsQual/HEPCountsQual/HEPCountsQual_before","HEExposure",PRB);

	Flux * PFluxTOF  = new Flux(finalhistos,finalresults, "PFluxTOF", "Acceptance_PTOF","Acceptance","TOFPfits/Fit Results/Primary Proton Counts","ExposureTOF",Global.GetToFPBins());
	Flux * PFluxNaF  = new Flux(finalhistos,finalresults, "PFluxNaF", "Acceptance_PNaF","Acceptance","NaFPfits/Fit Results/Primary Proton Counts","ExposureNaF",Global.GetNaFPBins());
	Flux * PFluxAgl  = new Flux(finalhistos,finalresults, "PFluxAgl", "Acceptance_PAgl","Acceptance","AglPfits/Fit Results/Primary Proton Counts","ExposureAgl",Global.GetAglPBins());
	
	Flux * RigPTOF = new Flux(finalhistos,finalresults, "RigPTOF", "Acceptance_PTOF","Acceptance","TOFPCounts/TOFPCounts/TOFPCounts_before","ExposureTOF",GlobalRig.GetToFPBins());
	Flux * RigPNaF = new Flux(finalhistos,finalresults, "RigPNaF", "Acceptance_PNaF","Acceptance","NaFPCounts/NaFPCounts/NaFPCounts_before","ExposureNaF",GlobalRig.GetNaFPBins());
	Flux * RigPAgl = new Flux(finalhistos,finalresults, "RigPAgl", "Acceptance_PAgl","Acceptance","AglPCounts/AglPCounts/AglPCounts_before","ExposureAgl",GlobalRig.GetAglPBins());


	DFluxTOF ->SetDefaultOutFile(finalhistos); 
	DFluxNaF ->SetDefaultOutFile(finalhistos);
	DFluxAgl ->SetDefaultOutFile(finalhistos);

	HEPFlux  ->SetDefaultOutFile(finalhistos);
	HEPFluxL1 ->SetDefaultOutFile(finalhistos);
	HEPFluxQ ->SetDefaultOutFile(finalhistos);
//	for(int i=0;i<10;i++) FluxTests[i] ->SetDefaultOutFile(finalhistos);

	PFluxTOF ->SetDefaultOutFile(finalhistos);
	PFluxNaF ->SetDefaultOutFile(finalhistos);
	PFluxAgl ->SetDefaultOutFile(finalhistos);


	RigPTOF->SetDefaultOutFile(finalhistos);
	RigPNaF->SetDefaultOutFile(finalhistos);
	RigPAgl->SetDefaultOutFile(finalhistos);



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
//	for(int i=0;i<10;i++) Filler_RTI.AddObject2beFilled(FluxTests[i],GetGenMomentum,GetGenMomentum,"IsProtonMC");
	Filler_RTI.ReinitializeAll(refill);


	if(!refill&&checkfile) {
		cout<<"************* D FLUX ************"<<endl;

		DFluxTOF-> Eval_Flux();
		DFluxNaF-> Eval_Flux();
		DFluxAgl-> Eval_Flux();


		DFluxTOF->SaveResults(finalresults);
		DFluxNaF->SaveResults(finalresults);
		DFluxAgl->SaveResults(finalresults);

		cout<<"************* P FLUX ************"<<endl;

		HEPFlux -> Eval_Flux();
		HEPFluxL1 -> Eval_Flux();
		HEPFluxQ -> Eval_Flux();
//		for(int i=0;i<10;i++) FluxTests[i] ->Eval_Flux();

		PFluxTOF-> Eval_Flux();
		PFluxNaF-> Eval_Flux();
		PFluxAgl-> Eval_Flux();
		RigPTOF-> Eval_Flux();
		RigPNaF-> Eval_Flux();
		RigPAgl-> Eval_Flux();

		RigPTOF->ChangeName("RigPTOF");
		RigPNaF->ChangeName("RigPNaF");
		RigPAgl->ChangeName("RigPAgl");

		HEPFlux ->SaveResults(finalresults);
		HEPFluxL1 ->SaveResults(finalresults);
		HEPFluxQ->SaveResults(finalresults);
	//	for(int i=0;i<10;i++) FluxTests[i] ->SaveResults(finalresults);
		PFluxTOF->SaveResults(finalresults);
		PFluxNaF->SaveResults(finalresults);
		PFluxAgl->SaveResults(finalresults);
		RigPTOF->SaveResults(finalresults);
		RigPNaF->SaveResults(finalresults);
		RigPAgl->SaveResults(finalresults);

		cout<<"************* MERGING RANGES ************"<<endl;
		
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

		finalresults.writeObjsInFolder("Fluxes/");	
	}
	return;
}

#endif
