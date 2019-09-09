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

//	for(int i=0;i<10;i++) FluxTests[i] = new Flux(finalhistos,finalresults,("PFluxTests_v"+to_string(i)).c_str(),("CountsTestEff_rig_v" + to_string(i)).c_str(),("CountsTestEff_rig_v" + to_string(i)).c_str(),("HEPTests_v"+to_string(i) + "/HEPTests_v"+to_string(i)+"/HEPTests_v"+to_string(i)+"_after").c_str(),"HEExposure","For_Acceptance_P",PRB);

	Flux * DFluxTOF  = new Flux(finalhistos,finalresults, "DFluxTOF", "FullSetTOT_D_TOF","FullSetTOT","TOFDfits/Fit Results/Primary Deuteron Counts","ExposureTOF",Global.GetToFDBins());
	Flux * DFluxNaF  = new Flux(finalhistos,finalresults, "DFluxNaF", "FullSetTOT_D_NaF","FullSetTOT","NaFDfits/Fit Results/Primary Deuteron Counts","ExposureNaF",Global.GetNaFDBins());
	Flux * DFluxAgl  = new Flux(finalhistos,finalresults, "DFluxAgl", "FullSetTOT_D_Agl","FullSetTOT","AglDfits/Fit Results/Primary Deuteron Counts","ExposureAgl",Global.GetAglDBins());

	Flux * HEPFlux   = new Flux(finalhistos,finalresults,"PFluxHE","RigBinBaselineEff_Trig","RigBinBaselineEff_Trig","HEPCounts/HEPCounts/HEPCounts_after","HEExposure",PRB);
	Flux * HEPFluxL1 = new Flux(finalhistos,finalresults,"PFluxL1HE","RigBinBaselineL1Eff_Trig","RigBinBaselineL1Eff_Trig","HEPCountsL1/HEPCountsL1/HEPCountsL1_after","HEExposure",PRB);
	Flux * HEPFluxQ  = new Flux(finalhistos,finalresults,"PFluxQHE","RigBinQualEff","RigBinQualEff","HEPCountsQual/HEPCountsQual/HEPCountsQual_after","HEExposure",PRB);

	Flux * PFluxTOF  = new Flux(finalhistos,finalresults, "PFluxTOF", "FullSetTOT_P_TOF","FullSetTOT","TOFPfits/Fit Results/Primary Proton Counts","ExposureTOF",Global.GetToFPBins());
	Flux * PFluxNaF  = new Flux(finalhistos,finalresults, "PFluxNaF", "FullSetTOT_P_NaF","FullSetTOT","NaFPfits/Fit Results/Primary Proton Counts","ExposureNaF",Global.GetNaFPBins());
	Flux * PFluxAgl  = new Flux(finalhistos,finalresults, "PFluxAgl", "FullSetTOT_P_Agl","FullSetTOT","AglPfits/Fit Results/Primary Proton Counts","ExposureAgl",Global.GetAglPBins());
	
	Flux * PFluxTOF_Ekin  = new Flux(finalhistos,finalresults, "PFluxTOF_Ekin", "FullSetTOT_P_TOF","FullSetTOT","TOFDfits/Fit Results/Primary Proton Counts","ExposureTOF",Global.GetToFDBins());
	Flux * PFluxNaF_Ekin  = new Flux(finalhistos,finalresults, "PFluxNaF_Ekin", "FullSetTOT_P_NaF","FullSetTOT","NaFDfits/Fit Results/Primary Proton Counts","ExposureNaF",Global.GetNaFDBins());
	Flux * PFluxAgl_Ekin  = new Flux(finalhistos,finalresults, "PFluxAgl_Ekin", "FullSetTOT_P_Agl","FullSetTOT","AglDfits/Fit Results/Primary Proton Counts","ExposureAgl",Global.GetAglDBins());

	Flux * DummyDTOF = new Flux(finalhistos,finalresults, "DummyDTOF", "FullSetTOT_D_TOF","FullSetTOT","TOFDfits/Fit Results/Primary Deuteron Counts","ExposureTOF",Global.GetToFDBins());
	Flux * DummyDNaF = new Flux(finalhistos,finalresults, "DummyDNaF", "FullSetTOT_D_NaF","FullSetTOT","NaFDfits/Fit Results/Primary Deuteron Counts","ExposureNaF",Global.GetNaFDBins());
	Flux * DummyDAgl = new Flux(finalhistos,finalresults, "DummyDAgl", "FullSetTOT_D_Agl","FullSetTOT","AglDfits/Fit Results/Primary Deuteron Counts","ExposureAgl",Global.GetAglDBins());

	Flux * DummyPTOF = new Flux(finalhistos,finalresults, "DummyPTOF", "FullSetTOT_P_TOF","FullSetTOT","TOFPCounts/TOFPCounts/TOFPCounts_after","ExposureTOF",Global.GetToFPBins());
	Flux * DummyPNaF = new Flux(finalhistos,finalresults, "DummyPNaF", "FullSetTOT_P_NaF","FullSetTOT","NaFPCounts/NaFPCounts/NaFPCounts_after","ExposureNaF",Global.GetNaFPBins());
	Flux * DummyPAgl = new Flux(finalhistos,finalresults, "DummyPAgl", "FullSetTOT_P_Agl","FullSetTOT","AglPCounts/AglPCounts/AglPCounts_after","ExposureAgl",Global.GetAglPBins());


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

	PFluxTOF_Ekin ->SetDefaultOutFile(finalhistos);
	PFluxNaF_Ekin ->SetDefaultOutFile(finalhistos);
	PFluxAgl_Ekin ->SetDefaultOutFile(finalhistos);

	DummyDTOF->SetDefaultOutFile(finalhistos);
	DummyDNaF->SetDefaultOutFile(finalhistos);
	DummyDAgl->SetDefaultOutFile(finalhistos);

	DummyPTOF->SetDefaultOutFile(finalhistos);
	DummyPNaF->SetDefaultOutFile(finalhistos);
	DummyPAgl->SetDefaultOutFile(finalhistos);



	cout<<"********** EXPOSURE TIME & GEOM. ACCEPT. ********"<<endl;

	Filler_RTI.AddObject2beFilled(DFluxTOF,GetBetaGen,GetBetaGen,"IsDeutonMC");
	Filler_RTI.AddObject2beFilled(DFluxNaF,GetBetaGen,GetBetaGen,"IsDeutonMC");
	Filler_RTI.AddObject2beFilled(DFluxAgl,GetBetaGen,GetBetaGen,"IsDeutonMC");
	Filler_RTI.AddObject2beFilled(PFluxTOF,GetBetaGen,GetBetaGen,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(PFluxNaF,GetBetaGen,GetBetaGen,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(PFluxAgl,GetBetaGen,GetBetaGen,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(PFluxTOF_Ekin,GetBetaGen,GetBetaGen,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(PFluxNaF_Ekin,GetBetaGen,GetBetaGen,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(PFluxAgl_Ekin,GetBetaGen,GetBetaGen,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(HEPFlux,GetGenMomentum,GetGenMomentum,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(HEPFluxL1,GetGenMomentum,GetGenMomentum,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(HEPFluxQ,GetGenMomentum,GetGenMomentum,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(DummyPTOF,GetBetaGen,GetBetaGen,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(DummyPNaF,GetBetaGen,GetBetaGen,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(DummyPAgl,GetBetaGen,GetBetaGen,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(DummyDTOF,GetBetaGen,GetBetaGen,"IsDeutonMC");
	Filler_RTI.AddObject2beFilled(DummyDNaF,GetBetaGen,GetBetaGen,"IsDeutonMC");
	Filler_RTI.AddObject2beFilled(DummyDAgl,GetBetaGen,GetBetaGen,"IsDeutonMC");
//	for(int i=0;i<10;i++) Filler_RTI.AddObject2beFilled(FluxTests[i],GetGenMomentum,GetGenMomentum,"IsProtonMC");
	Filler_RTI.ReinitializeAll(refill);


	if(!refill&&checkfile) {
		cout<<"************* D FLUX ************"<<endl;

		DFluxTOF-> Eval_Flux();
		DFluxNaF-> Eval_Flux();
		DFluxAgl-> Eval_Flux();
		DummyDTOF-> Eval_Flux();
		DummyDNaF-> Eval_Flux();
		DummyDAgl-> Eval_Flux();

		DummyDTOF->ChangeName("DummyDTOF");
		DummyDNaF->ChangeName("DummyDNaF");
		DummyDAgl->ChangeName("DummyDAgl");

		DFluxTOF->SaveResults(finalresults);
		DFluxNaF->SaveResults(finalresults);
		DFluxAgl->SaveResults(finalresults);
		DummyDTOF->SaveResults(finalresults);
		DummyDNaF->SaveResults(finalresults);
		DummyDAgl->SaveResults(finalresults);

		cout<<"************* P FLUX ************"<<endl;

		HEPFlux -> Eval_Flux();
		HEPFluxL1 -> Eval_Flux();
		HEPFluxQ -> Eval_Flux();
//		for(int i=0;i<10;i++) FluxTests[i] ->Eval_Flux();

		PFluxTOF-> Eval_Flux();
		PFluxNaF-> Eval_Flux();
		PFluxAgl-> Eval_Flux();
		PFluxTOF_Ekin-> Eval_Flux();
		PFluxNaF_Ekin-> Eval_Flux();
		PFluxAgl_Ekin-> Eval_Flux();
		DummyPTOF-> Eval_Flux();
		DummyPNaF-> Eval_Flux();
		DummyPAgl-> Eval_Flux();

		DummyPTOF->ChangeName("DummyPTOF");
		DummyPNaF->ChangeName("DummyPNaF");
		DummyPAgl->ChangeName("DummyPAgl");

		HEPFlux ->SaveResults(finalresults);
		HEPFluxL1 ->SaveResults(finalresults);
		HEPFluxQ->SaveResults(finalresults);
	//	for(int i=0;i<10;i++) FluxTests[i] ->SaveResults(finalresults);
		PFluxTOF->SaveResults(finalresults);
		PFluxNaF->SaveResults(finalresults);
		PFluxAgl->SaveResults(finalresults);
		PFluxTOF_Ekin->SaveResults(finalresults);
		PFluxNaF_Ekin->SaveResults(finalresults);
		PFluxAgl_Ekin->SaveResults(finalresults);
		DummyPTOF->SaveResults(finalresults);
		DummyPNaF->SaveResults(finalresults);
		DummyPAgl->SaveResults(finalresults);

		cout<<"************* D/P RATIO ************"<<endl;

		TH1F * DPRatioTOF = DFluxTOF->Eval_FluxRatio(PFluxTOF_Ekin,"DP ratio TOF");
		TH1F * DPRatioNaF = DFluxNaF->Eval_FluxRatio(PFluxNaF_Ekin,"DP ratio NaF");
		TH1F * DPRatioAgl = DFluxAgl->Eval_FluxRatio(PFluxAgl_Ekin,"DP ratio Agl");

		finalresults.Add(DPRatioTOF);
		finalresults.Add(DPRatioNaF);
		finalresults.Add(DPRatioAgl);

		finalresults.writeObjsInFolder("Fluxes/");	
	}
	return;
}

#endif
