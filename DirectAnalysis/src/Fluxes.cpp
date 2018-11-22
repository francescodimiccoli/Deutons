#ifndef FLUXES_H
#define FLUXES_H

#include "Analyzer.h"
#include "../include/Efficiency.h"
#include "../include/AllRangesEfficiency.h"
#include "../include/Flux.h"
#include "../include/EffCorr.h"
//#include "../include/EffCorrTemplate.h"


//MC parameters taken from /cvmfs/ams.cern.ch/Offline/AMSDataDir/DataManagement/DataSetsDesc

void Analyzer::BookFluxAnalysis(FileSaver finalhistos, FileSaver finalresults, bool refill){

	bool checkfile = finalhistos.CheckFile();
	check_file = checkfile;
	cout<<"****************************** BINS ***************************************"<<endl;
        SetUpUsualBinning();
 
	cout<<"****************************** FLUXES EVALUATION ******************************************"<<endl;


	Flux * DFluxTOF  = new Flux(finalhistos,finalresults, "DFluxTOF", "FullSetTOT_D_TOF","FullSetTOT","TOFfits/Fit Results/Primary Deuteron Counts","ExposureTOF",ToFDB);
	Flux * DFluxNaF  = new Flux(finalhistos,finalresults, "DFluxNaF", "FullSetTOT_D_NaF","FullSetTOT","NaFfits/Fit Results/Primary Deuteron Counts","ExposureNaF",NaFDB);
	Flux * DFluxAgl  = new Flux(finalhistos,finalresults, "DFluxAgl", "FullSetTOT_D_Agl","FullSetTOT","Aglfits/Fit Results/Primary Deuteron Counts","ExposureAgl",AglDB);

	Flux * HEPFlux   = new Flux(finalhistos,finalresults,"PFluxHE","RigBinBaselineEff","RigBinBaselineEff","HEPCounts/HEPCounts/HEPCounts_before","HEExposure",PRB);
	Flux * HEPFluxQ  = new Flux(finalhistos,finalresults,"PFluxQHE","RigBinQualEff","RigBinQualEff","HEPCountsQual/HEPCountsQual/HEPCountsQual_before","HEExposure",PRB);
	Flux * PFluxTOF  = new Flux(finalhistos,finalresults, "PFluxTOF", "FullSetTOT_P_TOF","FullSetTOT","TOFfits/Fit Results/Primary Proton Counts","ExposureTOF",ToFPB);
	Flux * PFluxNaF  = new Flux(finalhistos,finalresults, "PFluxNaF", "FullSetTOT_P_NaF","FullSetTOT","NaFfits/Fit Results/Primary Proton Counts","ExposureNaF",NaFPB);
	Flux * PFluxAgl  = new Flux(finalhistos,finalresults, "PFluxAgl", "FullSetTOT_P_Agl","FullSetTOT","Aglfits/Fit Results/Primary Proton Counts","ExposureAgl",AglPB);

	Flux * DummyDTOF = new Flux(finalhistos,finalresults, "DummyDTOF", "FullSetTOT_D_TOF","FullSetTOT","TOFfits/Fit Results/Primary Deuteron Counts","ExposureTOF",ToFDB);
	Flux * DummyDNaF = new Flux(finalhistos,finalresults, "DummyDNaF", "FullSetTOT_D_NaF","FullSetTOT","NaFfits/Fit Results/Primary Deuteron Counts","ExposureNaF",NaFDB);
	Flux * DummyDAgl = new Flux(finalhistos,finalresults, "DummyDAgl", "FullSetTOT_D_Agl","FullSetTOT","Aglfits/Fit Results/Primary Deuteron Counts","ExposureAgl",AglDB);

	Flux * DummyPTOF = new Flux(finalhistos,finalresults, "DummyPTOF", "FullSetTOT_P_TOF","FullSetTOT","TOFPCounts/TOFPCounts/TOFPCounts_before","ExposureTOF",ToFRigB);
	Flux * DummyPNaF = new Flux(finalhistos,finalresults, "DummyPNaF", "FullSetTOT_P_NaF","FullSetTOT","NaFPCounts/NaFPCounts/NaFPCounts_before","ExposureNaF",NaFRigB);
	Flux * DummyPAgl = new Flux(finalhistos,finalresults, "DummyPAgl", "FullSetTOT_P_Agl","FullSetTOT","AglPCounts/AglPCounts/AglPCounts_before","ExposureAgl",AglRigB);


	DFluxTOF ->SetDefaultOutFile(finalhistos); 
	DFluxNaF ->SetDefaultOutFile(finalhistos);
	DFluxAgl ->SetDefaultOutFile(finalhistos);

	HEPFlux  ->SetDefaultOutFile(finalhistos);
	HEPFluxQ ->SetDefaultOutFile(finalhistos);
	PFluxTOF ->SetDefaultOutFile(finalhistos);
	PFluxNaF ->SetDefaultOutFile(finalhistos);
	PFluxAgl ->SetDefaultOutFile(finalhistos);

	DummyDTOF->SetDefaultOutFile(finalhistos);
	DummyDNaF->SetDefaultOutFile(finalhistos);
	DummyDAgl->SetDefaultOutFile(finalhistos);

	DummyPTOF->SetDefaultOutFile(finalhistos);
	DummyPNaF->SetDefaultOutFile(finalhistos);
	DummyPAgl->SetDefaultOutFile(finalhistos);




	cout<<"****************************** EFFICIENCY CORR. READING  ******************************************"<<endl;
	std::string before;
        std::string after;
        before = "";
        after  = ""; 

	//baseline eff. corr.
	EffCorr * TriggerEffCorr_HE = new EffCorr(finalresults,"TriggerEffCorr_HE","Trigger Eff. Corr",PRB,before,after,"",    "IsPurePMC","IsPureDMC","IsDeutonMC"); 
	EffCorr * TriggerEffCorr_TOF = new EffCorr(finalresults,"TriggerEffCorr_TOF","Trigger Eff. Corr",ToFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * TriggerEffCorr_NaF = new EffCorr(finalresults,"TriggerEffCorr_NaF","Trigger Eff. Corr",NaFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * TriggerEffCorr_Agl = new EffCorr(finalresults,"TriggerEffCorr_Agl","Trigger Eff. Corr",AglPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");

	EffCorr * L1PickUpEffCorr_HE = new EffCorr(finalresults,"L1PickUpEffCorr_HE","L1PickUp Eff. Corr",PRB,before,after,"",    "IsPurePMC","IsPureDMC","IsDeutonMC"); 
	EffCorr * L1PickUpEffCorr_TOF = new EffCorr(finalresults,"L1PickUpEffCorr_TOF","L1PickUp Eff. Corr",ToFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * L1PickUpEffCorr_NaF = new EffCorr(finalresults,"L1PickUpEffCorr_NaF","L1PickUp Eff. Corr",NaFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * L1PickUpEffCorr_Agl = new EffCorr(finalresults,"L1PickUpEffCorr_Agl","L1PickUp Eff. Corr",AglPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");

	EffCorr * MinTOFEffCorr_HE  = new EffCorr(finalresults,"MinTOFEffCorr_HE" ,"Min TOF Eff. Corr",PRB,before,after,  "","IsPurePMC","IsPureDMC","IsDeutonMC"); 
	EffCorr * MinTOFEffCorr_TOF = new EffCorr(finalresults,"MinTOFEffCorr_TOF","Min TOF Eff. Corr",ToFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * MinTOFEffCorr_NaF = new EffCorr(finalresults,"MinTOFEffCorr_NaF","Min TOF Eff. Corr",NaFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * MinTOFEffCorr_Agl = new EffCorr(finalresults,"MinTOFEffCorr_Agl","Min TOF Eff. Corr",AglPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");

	EffCorr * TrackerEffCorr_HE  = new EffCorr(finalresults,"TrackerEffCorr_HE","Tracker Eff. Corr",PRB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * TrackerEffCorr_TOF = new EffCorr(finalresults,"TrackerEffCorr_TOF","Tracker Eff. Corr",ToFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * TrackerEffCorr_NaF = new EffCorr(finalresults,"TrackerEffCorr_NaF","Tracker Eff. Corr",NaFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * TrackerEffCorr_Agl = new EffCorr(finalresults,"TrackerEffCorr_Agl","Tracker Eff. Corr",AglPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");

	//Other selections eff. corr.
	EffCorr * GoodChi_HE 		= new EffCorr(finalresults,"GoodChiEffCorr_HE","GoodChi Eff. Corr",PRB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodChi_TOF 		= new EffCorr(finalresults,"GoodChiEffCorr_TOF","GoodChi Eff. Corr",ToFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodChi_NaF 		= new EffCorr(finalresults,"GoodChiEffCorr_NaF","GoodChi Eff. Corr",NaFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodChi_Agl 		= new EffCorr(finalresults,"GoodChiEffCorr_Agl","GoodChi Eff. Corr",AglPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");

	EffCorr * GoodUtof_HE 		= new EffCorr(finalresults,"GoodUtofEffCorr_HE","GoodUtof Eff. Corr",PRB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodUtof_TOF 		= new EffCorr(finalresults,"GoodUtofEffCorr_TOF","GoodUtof Eff. Corr",ToFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodUtof_NaF 		= new EffCorr(finalresults,"GoodUtofEffCorr_NaF","GoodUtof Eff. Corr",NaFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodUtof_Agl 		= new EffCorr(finalresults,"GoodUtofEffCorr_Agl","GoodUtof Eff. Corr",AglPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");

        EffCorr * GoodLtof_HE 		= new EffCorr(finalresults,"GoodLTOFEffCorr_HE","GoodQTrack Eff. Corr",PRB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
        EffCorr * GoodLtof_TOF 		= new EffCorr(finalresults,"GoodLTOFEffCorr_TOF","GoodQTrack Eff. Corr",ToFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
        EffCorr * GoodLtof_NaF 		= new EffCorr(finalresults,"GoodLTOFEffCorr_NaF","GoodQTrack Eff. Corr",NaFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
        EffCorr * GoodLtof_Agl 		= new EffCorr(finalresults,"GoodLTOFEffCorr_Agl","GoodQTrack Eff. Corr",AglPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");

	EffCorr * Good1Track_HE 	= new EffCorr(finalresults,"Good1TrackEffCorr_HE","Good1Track Eff. Corr",PRB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * Good1Track_TOF 	= new EffCorr(finalresults,"Good1TrackEffCorr_TOF","Good1Track Eff. Corr",ToFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * Good1Track_NaF 	= new EffCorr(finalresults,"Good1TackEffCorr_NaF", "Good1Track Eff. Corr",NaFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * Good1Track_Agl 	= new EffCorr(finalresults,"Good1TrackEffCorr_Agl","Good1Track Eff. Corr",AglPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");

	EffCorr * GoodQTrack_HE 	= new EffCorr(finalresults,"GoodQTrackEffCorr_HE","GoodQTrack Eff. Corr",PRB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodQTrack_TOF 	= new EffCorr(finalresults,"GoodQTrackEffCorr_TOF","GoodQTrack Eff. Corr",ToFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodQTrack_NaF 	= new EffCorr(finalresults,"GoodQTrackEffCorr_NaF","GoodQTrack Eff. Corr",NaFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodQTrack_Agl 	= new EffCorr(finalresults,"GoodQTrackEffCorr_Agl","GoodQTrack Eff. Corr",AglPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");

	EffCorr * GoodTime_TOF 		= new EffCorr(finalresults,"GoodTimeEffCorr_TOF","GoodTime Eff. Corr",ToFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * RICHEffCorr_NaF 	= new EffCorr(finalresults,"RICHCorrection_NaF","RICH Eff. Corr",NaFPB,before,(after+"&IsFromNaF").c_str(),"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * RICHEffCorr_Agl 	= new EffCorr(finalresults,"RICHCorrection_Agl","RICH Eff. Corr",AglPB,before,(after+"&IsFromAgl").c_str(),"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * RICHQualEffCorr_NaF 	= new EffCorr(finalresults,"RICHQualCorrection_NaF","RICH Qual Eff. Corr",NaFPB,(before+"&IsFromNaF").c_str(),(after+"&IsFromNaF&RICHBDTCut").c_str(),"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * RICHQualEffCorr_Agl 	= new EffCorr(finalresults,"RICHqualCorrection_Agl","RICH Qual. Eff. Corr",AglPB,(before+"&IsFromNaF").c_str(),(after+"&IsFromAgl&RICHBDTCut").c_str(),"","IsPurePMC","IsPureDMC","IsDeutonMC");

	TrackerEffCorr_TOF -> Eval_Efficiencies();
	TrackerEffCorr_NaF -> Eval_Efficiencies();
	TrackerEffCorr_Agl -> Eval_Efficiencies();


	TrackerEffCorr_TOF -> Eval_Corrections();
	TrackerEffCorr_NaF -> Eval_Corrections();
	TrackerEffCorr_Agl -> Eval_Corrections();

	TrackerEffCorr_HE -> SetToConstantValue(1.02);
	TrackerEffCorr_TOF -> SetToConstantValue(1.02);
	TrackerEffCorr_NaF -> SetToConstantValue(1.02);
	TrackerEffCorr_Agl -> SetToConstantValue(1.02);


	cout<<"********** EXPOSURE TIME & GEOM. ACCEPT. ********"<<endl;

	Filler_RTI.AddObject2beFilled(DFluxTOF,GetBetaGen,GetBetaGen,"IsDeutonMC");
	Filler_RTI.AddObject2beFilled(DFluxNaF,GetBetaGen,GetBetaGen,"IsDeutonMC");
	Filler_RTI.AddObject2beFilled(DFluxAgl,GetBetaGen,GetBetaGen,"IsDeutonMC");
	Filler_RTI.AddObject2beFilled(PFluxTOF,GetBetaGen,GetBetaGen,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(PFluxNaF,GetBetaGen,GetBetaGen,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(PFluxAgl,GetBetaGen,GetBetaGen,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(HEPFlux,GetGenMomentum,GetGenMomentum,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(HEPFluxQ,GetGenMomentum,GetGenMomentum,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(DummyPTOF,GetBetaGen,GetBetaGen,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(DummyPNaF,GetBetaGen,GetBetaGen,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(DummyPAgl,GetBetaGen,GetBetaGen,"IsProtonMC");
	Filler_RTI.AddObject2beFilled(DummyDTOF,GetBetaGen,GetBetaGen,"IsDeutonMC");
	Filler_RTI.AddObject2beFilled(DummyDNaF,GetBetaGen,GetBetaGen,"IsDeutonMC");
	Filler_RTI.AddObject2beFilled(DummyDAgl,GetBetaGen,GetBetaGen,"IsDeutonMC");
	Filler_RTI.ReinitializeAll(refill);

	if(!refill&&checkfile) {
		cout<<"************* D FLUX ************"<<endl;

		DFluxTOF->Set_MCPar(0.5,200,1,"D.B1081/d.pl1.l1.0_5200.2_01.info");	
		DFluxNaF->Set_MCPar(0.5,200,1,"D.B1081/d.pl1.l1.0_5200.2_01.info");	
		DFluxAgl->Set_MCPar(0.5,200,1,"D.B1081/d.pl1.l1.0_5200.2_01.info");	
		DummyDTOF->Set_MCPar(0.5,200,1,"D.B1081/d.pl1.l1.0_5200.2_01.info");	
		DummyDNaF->Set_MCPar(0.5,200,1,"D.B1081/d.pl1.l1.0_5200.2_01.info");	
		DummyDAgl->Set_MCPar(0.5,200,1,"D.B1081/d.pl1.l1.0_5200.2_01.info");	


		DFluxTOF->ApplyEfficCorr(GoodChi_TOF 	     ->GetGlobCorrection());
		DFluxNaF->ApplyEfficCorr(GoodChi_NaF 	     ->GetGlobCorrection());
		DFluxAgl->ApplyEfficCorr(GoodChi_Agl 	     ->GetGlobCorrection());
		DFluxTOF->ApplyEfficCorr(GoodUtof_TOF 	     ->GetGlobCorrection());
		DFluxNaF->ApplyEfficCorr(GoodUtof_NaF 	     ->GetGlobCorrection());
		DFluxAgl->ApplyEfficCorr(GoodUtof_Agl 	     ->GetGlobCorrection());
		DFluxTOF->ApplyEfficCorr(GoodLtof_TOF 	     ->GetGlobCorrection());
		DFluxNaF->ApplyEfficCorr(GoodLtof_NaF 	     ->GetGlobCorrection());
		DFluxAgl->ApplyEfficCorr(GoodLtof_Agl 	     ->GetGlobCorrection());
		DFluxTOF->ApplyEfficCorr(Good1Track_TOF     ->GetGlobCorrection());
		DFluxNaF->ApplyEfficCorr(Good1Track_NaF     ->GetGlobCorrection());
		DFluxAgl->ApplyEfficCorr(Good1Track_Agl     ->GetGlobCorrection());
		DFluxTOF->ApplyEfficCorr(GoodQTrack_TOF     ->GetGlobCorrection());
		DFluxNaF->ApplyEfficCorr(GoodQTrack_NaF     ->GetGlobCorrection());
		DFluxAgl->ApplyEfficCorr(GoodQTrack_Agl     ->GetGlobCorrection());
		DFluxTOF->ApplyEfficCorr(GoodTime_TOF 	     ->GetGlobCorrection());
		DFluxNaF->ApplyEfficCorr(RICHEffCorr_NaF    ->GetGlobCorrection());
		DFluxAgl->ApplyEfficCorr(RICHEffCorr_Agl    ->GetGlobCorrection());
		DFluxNaF->ApplyEfficCorr(RICHQualEffCorr_NaF->GetGlobCorrection());
		DFluxAgl->ApplyEfficCorr(RICHQualEffCorr_Agl->GetGlobCorrection());


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

		HEPFlux ->Set_MCPar(0.5,100,1,"Pr.B1119/pr.pl1.05100.3_01.info");	
		HEPFluxQ->Set_MCPar(0.5,100,1,"Pr.B1119/pr.pl1.05100.3_01.info");	
		PFluxTOF->Set_MCPar(0.5,100,1,"Pr.B1119/pr.pl1.05100.3_01.info");	
		PFluxNaF->Set_MCPar(0.5,100,1,"Pr.B1119/pr.pl1.05100.3_01.info");	
		PFluxAgl->Set_MCPar(0.5,100,1,"Pr.B1119/pr.pl1.05100.3_01.info");	
		DummyPTOF->Set_MCPar(0.5,100,1,"Pr.B1119/pr.pl1.05100.3_01.info");	
		DummyPNaF->Set_MCPar(0.5,100,1,"Pr.B1119/pr.pl1.05100.3_01.info");	
		DummyPAgl->Set_MCPar(0.5,100,1,"Pr.B1119/pr.pl1.05100.3_01.info");	


		//baseline flux corrections
		HEPFlux->ApplyEfficCorr(TriggerEffCorr_HE   ->GetGlobCorrection());
		HEPFlux->ApplyEfficCorr(L1PickUpEffCorr_HE  ->GetGlobCorrection());
		HEPFlux->ApplyEfficCorr(TrackerEffCorr_HE   ->GetGlobCorrection());
		HEPFlux->ApplyEfficCorr(MinTOFEffCorr_HE    ->GetGlobCorrection());

		//quality flux corrections
		HEPFluxQ->ApplyEfficCorr(TriggerEffCorr_HE   ->GetGlobCorrection());
		HEPFluxQ->ApplyEfficCorr(L1PickUpEffCorr_HE  ->GetGlobCorrection());
		HEPFluxQ->ApplyEfficCorr(TrackerEffCorr_HE   ->GetGlobCorrection());
		HEPFluxQ->ApplyEfficCorr(MinTOFEffCorr_HE    ->GetGlobCorrection());
		HEPFluxQ->ApplyEfficCorr(GoodChi_HE 	     ->GetGlobCorrection());
		HEPFluxQ->ApplyEfficCorr(GoodUtof_HE 	     ->GetGlobCorrection());
		HEPFluxQ->ApplyEfficCorr(GoodLtof_HE 	     ->GetGlobCorrection());
		HEPFluxQ->ApplyEfficCorr(Good1Track_HE     ->GetGlobCorrection());
		HEPFluxQ->ApplyEfficCorr(GoodQTrack_HE     ->GetGlobCorrection());

		//corrections for flux with R edges
		DummyPTOF->ApplyEfficCorr(TriggerEffCorr_TOF ->GetGlobCorrection());
		DummyPNaF->ApplyEfficCorr(TriggerEffCorr_NaF ->GetGlobCorrection());
		DummyPAgl->ApplyEfficCorr(TriggerEffCorr_Agl ->GetGlobCorrection());
		DummyPTOF->ApplyEfficCorr(L1PickUpEffCorr_TOF ->GetGlobCorrection());
		DummyPNaF->ApplyEfficCorr(L1PickUpEffCorr_NaF ->GetGlobCorrection());
		DummyPAgl->ApplyEfficCorr(L1PickUpEffCorr_Agl ->GetGlobCorrection());
		DummyPTOF->ApplyEfficCorr(TrackerEffCorr_TOF ->GetGlobCorrection());
		DummyPNaF->ApplyEfficCorr(TrackerEffCorr_NaF ->GetGlobCorrection());
		DummyPAgl->ApplyEfficCorr(TrackerEffCorr_Agl ->GetGlobCorrection());
		DummyPTOF->ApplyEfficCorr(MinTOFEffCorr_TOF ->GetGlobCorrection());
		DummyPNaF->ApplyEfficCorr(MinTOFEffCorr_NaF ->GetGlobCorrection());
		DummyPAgl->ApplyEfficCorr(MinTOFEffCorr_Agl ->GetGlobCorrection());

		DummyPTOF->ApplyEfficCorr(GoodChi_TOF 	     ->GetGlobCorrection());
		DummyPNaF->ApplyEfficCorr(GoodChi_NaF 	     ->GetGlobCorrection());
		DummyPAgl->ApplyEfficCorr(GoodChi_Agl 	     ->GetGlobCorrection());
		DummyPTOF->ApplyEfficCorr(GoodUtof_TOF 	     ->GetGlobCorrection());
		DummyPNaF->ApplyEfficCorr(GoodUtof_NaF 	     ->GetGlobCorrection());
		DummyPAgl->ApplyEfficCorr(GoodUtof_Agl 	     ->GetGlobCorrection());
		DummyPTOF->ApplyEfficCorr(GoodLtof_TOF 	     ->GetGlobCorrection());
		DummyPNaF->ApplyEfficCorr(GoodLtof_NaF 	     ->GetGlobCorrection());
		DummyPAgl->ApplyEfficCorr(GoodLtof_Agl 	     ->GetGlobCorrection());
		DummyPTOF->ApplyEfficCorr(Good1Track_TOF     ->GetGlobCorrection());
		DummyPNaF->ApplyEfficCorr(Good1Track_NaF     ->GetGlobCorrection());
		DummyPAgl->ApplyEfficCorr(Good1Track_Agl     ->GetGlobCorrection());
		DummyPTOF->ApplyEfficCorr(GoodQTrack_TOF     ->GetGlobCorrection());
		DummyPNaF->ApplyEfficCorr(GoodQTrack_NaF     ->GetGlobCorrection());
		DummyPAgl->ApplyEfficCorr(GoodQTrack_Agl     ->GetGlobCorrection());
		DummyPTOF->ApplyEfficCorr(GoodTime_TOF 	     ->GetGlobCorrection());
		DummyPNaF->ApplyEfficCorr(RICHEffCorr_NaF    ->GetGlobCorrection());
		DummyPAgl->ApplyEfficCorr(RICHEffCorr_Agl    ->GetGlobCorrection());
		DummyPNaF->ApplyEfficCorr(RICHQualEffCorr_NaF->GetGlobCorrection());
		DummyPAgl->ApplyEfficCorr(RICHQualEffCorr_Agl->GetGlobCorrection());


		//corrections for Analysis flux	
		PFluxTOF->ApplyEfficCorr(TriggerEffCorr_TOF ->GetGlobCorrection());
		PFluxNaF->ApplyEfficCorr(TriggerEffCorr_NaF ->GetGlobCorrection());
		PFluxAgl->ApplyEfficCorr(TriggerEffCorr_Agl ->GetGlobCorrection());
		PFluxTOF->ApplyEfficCorr(L1PickUpEffCorr_TOF ->GetGlobCorrection());
		PFluxNaF->ApplyEfficCorr(L1PickUpEffCorr_NaF ->GetGlobCorrection());
		PFluxAgl->ApplyEfficCorr(L1PickUpEffCorr_Agl ->GetGlobCorrection());
		PFluxTOF->ApplyEfficCorr(TrackerEffCorr_TOF ->GetGlobCorrection());
		PFluxNaF->ApplyEfficCorr(TrackerEffCorr_NaF ->GetGlobCorrection());
		PFluxAgl->ApplyEfficCorr(TrackerEffCorr_Agl ->GetGlobCorrection());
		PFluxTOF->ApplyEfficCorr(MinTOFEffCorr_TOF ->GetGlobCorrection());
		PFluxNaF->ApplyEfficCorr(MinTOFEffCorr_NaF ->GetGlobCorrection());
		PFluxAgl->ApplyEfficCorr(MinTOFEffCorr_Agl ->GetGlobCorrection());

		PFluxTOF->ApplyEfficCorr(GoodChi_TOF 	     ->GetGlobCorrection());
		PFluxNaF->ApplyEfficCorr(GoodChi_NaF 	     ->GetGlobCorrection());
		PFluxAgl->ApplyEfficCorr(GoodChi_Agl 	     ->GetGlobCorrection());
		PFluxTOF->ApplyEfficCorr(GoodUtof_TOF 	     ->GetGlobCorrection());
		PFluxNaF->ApplyEfficCorr(GoodUtof_NaF 	     ->GetGlobCorrection());
		PFluxAgl->ApplyEfficCorr(GoodUtof_Agl 	     ->GetGlobCorrection());
		PFluxTOF->ApplyEfficCorr(GoodLtof_TOF 	     ->GetGlobCorrection());
		PFluxNaF->ApplyEfficCorr(GoodLtof_NaF 	     ->GetGlobCorrection());
		PFluxAgl->ApplyEfficCorr(GoodLtof_Agl 	     ->GetGlobCorrection());
		PFluxTOF->ApplyEfficCorr(Good1Track_TOF     ->GetGlobCorrection());
		PFluxNaF->ApplyEfficCorr(Good1Track_NaF     ->GetGlobCorrection());
		PFluxAgl->ApplyEfficCorr(Good1Track_Agl     ->GetGlobCorrection());
		PFluxTOF->ApplyEfficCorr(GoodQTrack_TOF     ->GetGlobCorrection());
		PFluxNaF->ApplyEfficCorr(GoodQTrack_NaF     ->GetGlobCorrection());
		PFluxAgl->ApplyEfficCorr(GoodQTrack_Agl     ->GetGlobCorrection());
		PFluxTOF->ApplyEfficCorr(GoodTime_TOF 	     ->GetGlobCorrection());
		PFluxNaF->ApplyEfficCorr(RICHEffCorr_NaF    ->GetGlobCorrection());
		PFluxAgl->ApplyEfficCorr(RICHEffCorr_Agl    ->GetGlobCorrection());
		PFluxNaF->ApplyEfficCorr(RICHQualEffCorr_NaF->GetGlobCorrection());
		PFluxAgl->ApplyEfficCorr(RICHQualEffCorr_Agl->GetGlobCorrection());

		HEPFlux -> Eval_Flux();
		HEPFluxQ -> Eval_Flux();
		PFluxTOF-> Eval_Flux();
		PFluxNaF-> Eval_Flux();
		PFluxAgl-> Eval_Flux();
		DummyPTOF-> Eval_Flux();
		DummyPNaF-> Eval_Flux();
		DummyPAgl-> Eval_Flux();

		DummyPTOF->ChangeName("DummyPTOF");
		DummyPNaF->ChangeName("DummyPNaF");
		DummyPAgl->ChangeName("DummyPAgl");

		HEPFlux ->SaveResults(finalresults);
		HEPFluxQ->SaveResults(finalresults);
		PFluxTOF->SaveResults(finalresults);
		PFluxNaF->SaveResults(finalresults);
		PFluxAgl->SaveResults(finalresults);
		DummyPTOF->SaveResults(finalresults);
		DummyPNaF->SaveResults(finalresults);
		DummyPAgl->SaveResults(finalresults);


		cout<<"************* D/P RATIO ************"<<endl;


		TH1F * DPRatioTOF = DFluxTOF->Eval_FluxRatio(PFluxTOF,"DP ratio TOF");
		TH1F * DPRatioNaF = DFluxNaF->Eval_FluxRatio(PFluxNaF,"DP ratio NaF");
		TH1F * DPRatioAgl = DFluxAgl->Eval_FluxRatio(PFluxAgl,"DP ratio Agl");

		finalresults.Add(DPRatioTOF);
		finalresults.Add(DPRatioNaF);
		finalresults.Add(DPRatioAgl);

		finalresults.writeObjsInFolder("Fluxes/");	
	}
	return;
}

#endif
