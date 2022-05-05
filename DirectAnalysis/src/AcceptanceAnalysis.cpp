#include "Analyzer.h"
#include "../include/filesaver.h"
#include "../include/Efficiency.h"
#include "../include/EffCorr.h"
#include "../include/EffCorrTemplate.h"
#include "Acceptance.h"
#include "TStyle.h"

void AnalyzeEffCorr(EffCorr * Correction, FileSaver  finalhistos, FileSaver  finalresults, TFile * systfile, std::string systidirname,std::string systerrname, std::string avgname,std::string tima,float shift =0, bool IsTrig=false);

void DrawBlockCascade(std::vector<EffCorr *> block, std::string Blockname,FileSaver Plots,float rangemin, float rangemax);
void DrawCorrectionBlock(std::vector<EffCorr *> block, std::string Blockname,FileSaver Plots);
void DrawTotalUncertainty(std::vector<std::vector<EffCorr *>> Total,FileSaver Plots, std::string name,std::vector<std::string> Names);
	

void Analyzer::BookAcceptanceAnalysis(FileSaver finalhistos, FileSaver finalresults, bool refill)
{

	bool checkfile = finalhistos.CheckFile();
        check_file = checkfile;
	cout<<"****************************** BINS ***************************************"<<endl;
        SetUpTOIBinning();
 
	cout<<"****************************** Efficiency corr. ANALYIS ******************************************"<<endl;

	BadEventSimulator * NaFBadEvSimulator= new BadEventSimulator("IsFromNaF",30,0.82,1);
        BadEventSimulator * AglBadEvSimulator= new BadEventSimulator("IsFromAgl",120,0.96,1);


	//Baseline efficiency corrections
	std::string before;
        std::string after;
	before = "IsPhysTrig&IsDownGoing&IsGoodTOFStandaloneQ1&IsGoodTRDStandaloneQ1&HasECAL";
        after  = "IsPhysTrig&IsDownGoing&IsGoodTOFStandaloneQ1&IsGoodTRDStandaloneQ1&HasECAL&IsGoodTrack"; 
	EffCorr * TrackerEffCorr_HE = new EffCorr(finalhistos,"TrackerEffCorr_HE","Tracker Eff. Corr",false,1,before,after,"",true);
	before = "IsPhysTrig&IsDownGoing&IsGoodTOFStandaloneQ2&IsGoodTRDStandaloneQ2&HasECAL";
        after  = "IsPhysTrig&IsDownGoing&IsGoodTOFStandaloneQ2&IsGoodTRDStandaloneQ2&HasECAL&IsGoodTrack"; 
	EffCorr * TrackerEffCorr_Z2 = new EffCorr(finalhistos,"TrackerEffCorr_Z2","Tracker Eff. Corr",false,2,before,after,"",true,"IsHeliumMC");

	before = "IsPhysTrig&IsDownGoing&IsGoodTOFStandaloneQ1&IsGoodTRDStandaloneQ1&HasECAL&IsGoodTrack";
        after  = "IsPhysTrig&IsDownGoing&IsGoodTOFStandaloneQ1&IsGoodTRDStandaloneQ1&HasECAL&IsGoodTrack&IsCharge1Track"; 
	EffCorr * GoodQTrack_HE  = new EffCorr(finalhistos,"GoodQTrackEffCorr_HE","GoodQTrack Eff. Corr",false,1,before,after,"",true);
	before = "IsPhysTrig&IsDownGoing&IsGoodTOFStandaloneQ2&IsGoodTRDStandaloneQ2&HasECAL&IsGoodTrack";
        after  = "IsPhysTrig&IsDownGoing&IsGoodTOFStandaloneQ2&IsGoodTRDStandaloneQ2&HasECAL&IsGoodTrack&IsCharge2Track"; 
	EffCorr * GoodQTrack_Z2  = new EffCorr(finalhistos,"GoodQTrackEffCorr_Z2","GoodQTrack Eff. Corr",false,2,before,after,"",true,"IsHeliumMC");

	before = "IsPhysTrig&IsDownGoing&IsGoodTOFStandaloneQ1&IsGoodTRDStandaloneQ1&HasECAL&IsGoodTrack&IsCharge1Track";
        after  = "IsPhysTrig&IsDownGoing&IsGoodTOFStandaloneQ1&IsGoodTRDStandaloneQ1&HasECAL&IsGoodTrack&IsCharge1Track&IsGoodChi2"; 
	EffCorr * GoodChi_HE  = new EffCorr(finalhistos,"GoodChiEffCorr_HE","GoodChi Eff. Corr",false,1,before,after,"",true);
	before = "IsPhysTrig&IsDownGoing&IsGoodTOFStandaloneQ2&IsGoodTRDStandaloneQ2&HasECAL&IsGoodTrack&IsCharge2Track";
        after  = "IsPhysTrig&IsDownGoing&IsGoodTOFStandaloneQ2&IsGoodTRDStandaloneQ2&HasECAL&IsGoodTrack&IsCharge2Track&IsGoodChi2"; 
	EffCorr * GoodChi_Z2  = new EffCorr(finalhistos,"GoodChiEffCorr_Z2","GoodChi Eff. Corr",false,2,before,after,"",true,"IsHeliumMC");

	before = "IsPhysTrig&IsDownGoing&IsGoodTOFStandaloneQ1&IsGoodTRDStandaloneQ1&HasECAL&IsGoodTrack&IsCharge1Track&IsGoodChi2";
        after  = "IsPhysTrig&IsDownGoing&IsGoodTOFStandaloneQ1&IsGoodTRDStandaloneQ1&HasECAL&IsGoodTrack&IsCharge1Track&IsGoodChi2&IsGoodKalman"; 
	EffCorr * KalmanEffCorr_HE = new EffCorr(finalhistos,"KalmanEffCorr_HE","Kalman Eff. Corr",false,1,before,after,"",true);
	before = "IsPhysTrig&IsDownGoing&IsGoodTOFStandaloneQ2&IsGoodTRDStandaloneQ2&HasECAL&IsGoodTrack&IsCharge2Track&IsGoodChi2";
        after  = "IsPhysTrig&IsDownGoing&IsGoodTOFStandaloneQ2&IsGoodTRDStandaloneQ2&HasECAL&IsGoodTrack&IsCharge2Track&IsGoodChi2&IsGoodKalman"; 
	EffCorr * KalmanEffCorr_Z2 = new EffCorr(finalhistos,"KalmanEffCorr_Z2","Kalman Eff. Corr",false,2,before,after,"",true,"IsHeliumMC");

	// L1 Pickup efficiency corrections
	before = "IsDownGoing&IsPhysTrig&IsGoodTOFStandaloneQ1&IsGoodTRDStandaloneQ1&IsGoodTrack&IsGoodChi2&IsCharge1Track&Is1TrTrack&IsMinTOF&IsL1Fiducial&IsGoodTime";
        after  = "IsDownGoing&IsPhysTrig&IsGoodTOFStandaloneQ1&IsGoodTRDStandaloneQ1&IsGoodTrack&IsGoodChi2&IsCharge1Track&Is1TrTrack&IsMinTOF&IsL1Fiducial&IsGoodTime&IsL1HitNearExtrapol";
	EffCorr * L1PickUpGeom_HE = new EffCorr(finalhistos,"L1PickUpGeom_HE","L1PickUp Geom.",false,1,before,after,"IsPrimaryInner",true); 
	before = "IsDownGoing&IsPhysTrig&IsGoodTOFStandaloneQ2&IsGoodTRDStandaloneQ2&IsGoodTrack&IsGoodChi2&IsCharge2Track&Is1TrTrack&IsMinTOF&IsL1Fiducial&IsGoodTimeHe";
        after  = "IsDownGoing&IsPhysTrig&IsGoodTOFStandaloneQ2&IsGoodTRDStandaloneQ2&IsGoodTrack&IsGoodChi2&IsCharge2Track&Is1TrTrack&IsMinTOF&IsL1Fiducial&IsGoodTimeHe&IsL1HitNearExtrapol";
	EffCorr * L1PickUpGeom_Z2 = new EffCorr(finalhistos,"L1PickUpGeom_Z2","L1PickUp Geom.",false,2,before,after,"IsPrimaryInner",true,"IsHeliumMC"); 

	before = "IsDownGoing&IsPhysTrig&IsGoodTOFStandaloneQ1&IsGoodTrack&IsGoodChi2&IsCharge1Track&IsL1HitNearExtrapol&IsMinTOF&IsL1Fiducial&IsGoodTime";
        after  = "IsDownGoing&IsPhysTrig&IsGoodTOFStandaloneQ1&IsGoodTrack&IsGoodChi2&IsCharge1Track&IsL1HitNearExtrapol&IsMinTOF&IsL1Fiducial&IsGoodTime&L1LooseCharge1";
	EffCorr * L1PickUpEffCorr_HE = new EffCorr(finalhistos,"L1PickUpEffCorr_HE","L1PickUp Eff. Corr",false,1,before,after,"IsPrimaryInner",true); 
	before = "IsDownGoing&IsPhysTrig&IsGoodTOFStandaloneQ2&IsGoodTrack&IsGoodChi2&IsCharge2Track&IsL1HitNearExtrapol&IsMinTOF&IsL1Fiducial&IsGoodTimeHe";
        after  = "IsDownGoing&IsPhysTrig&IsGoodTOFStandaloneQ2&IsGoodTrack&IsGoodChi2&IsCharge2Track&IsL1HitNearExtrapol&IsMinTOF&IsL1Fiducial&IsGoodTimeHe&L1LooseCharge2";
	EffCorr * L1PickUpEffCorr_Z2 = new EffCorr(finalhistos,"L1PickUpEffCorr_Z2","L1PickUp Eff. Corr",false,2,before,after,"IsPrimaryInner",true,"IsHeliumMC"); 

	// Trigger Eff
	before = "IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge1Track&IsGoodKalman&L1LooseCharge1&IsLUT2";
        after  = "IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge1Track&IsGoodKalman&L1LooseCharge1&IsPhysTrig";
	EffCorr * TriggerEffCorr_HE  = new EffCorr(finalhistos,"TriggerEffCorr_HE" ,"Trigger Eff. Corr",false,1  ,before,after,"IsPrimary"); 
	before = "IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge2Track&IsGoodKalman&L1LooseCharge2&IsLUT2";
        after  = "IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge2Track&IsGoodKalman&L1LooseCharge2&IsPhysTrig";
	EffCorr * TriggerEffCorr_Z2  = new EffCorr(finalhistos,"TriggerEffCorr_Z2" ,"Trigger Eff. Corr",false,2  ,before,after,"IsPrimary",false,"IsHeliumMC"); 


        before = "IsDownGoin&IsGoodTrack&IsGoodChi2&IsCharge1Track&IsGoodKalman&L1LooseCharge1&HasL9&IsLUT2";
        after  = "IsDownGoin&IsGoodTrack&IsGoodChi2&IsCharge1Track&IsGoodKalman&L1LooseCharge1&HasL9&IsPhysTrig";
	EffCorr * TriggerFullSpan_HE  = new EffCorr(finalhistos,"TriggerFullSpan_HE" ,"Trigger Full Span",false,2  ,before,after,"IsPrimary"); 


	// Good Z=1 and Golden Efficiency corrections
	before = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1";
        after  = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1&Is1TrTrack"; 
	EffCorr * Good1Track_HE  = new EffCorr(finalhistos,"Good1TrackEffCorr_HE","Good1Track Eff. Corr",true,1,before,after,"IsPrimary");
	before = "IsPositive&IsPhysTrig&IsBaselineHe&L1LooseCharge2";
        after  = "IsPositive&IsPhysTrig&IsBaselineHe&L1LooseCharge2&Is1TrTrack"; 
	EffCorr * Good1Track_Z2  = new EffCorr(finalhistos,"Good1TrackEffCorr_Z2","Good1Track Eff. Corr",true,2,before,after,"IsPrimary",false,"isHeliumMC");

	before = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1&Is1TrTrack";
	after  = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1&Is1TrTrack&IsCharge1LTOF";
        EffCorr * GoodLtof_HE  = new EffCorr(finalhistos,"GoodLTOFEffCorr_HE" ,"GoodLtof Eff. Corr",true,1,before,after,"IsPrimary");

	before = "IsPositive&IsPhysTrig&IsBaselineHe&L1LooseCharge2&Is1TrTrack";
	after  = "IsPositive&IsPhysTrig&IsBaselineHe&L1LooseCharge2&Is1TrTrack&IsCharge2LTOF";
        EffCorr * GoodLtof_Z2  = new EffCorr(finalhistos,"GoodLTOFEffCorr_Z2" ,"GoodLtof Eff. Corr",true,2,before,after,"IsPrimary",false,"IsHeliumMC");

	before = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1&Is1TrTrack&IsCharge1LTOF";
        after  = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1&Is1TrTrack&IsCharge1LTOF&IsCharge1UTOF"; 
	EffCorr * GoodUtof_HE  = new EffCorr(finalhistos,"GoodUtofEffCorr_HE" ,"GoodUtof Eff. Corr",true,1,before,after,"IsPrimary");

	before = "IsPositive&IsPhysTrig&IsBaselineHe&L1LooseCharge2&Is1TrTrack&IsCharge2LTOF";
        after  = "IsPositive&IsPhysTrig&IsBaselineHe&L1LooseCharge2&Is1TrTrack&IsCharge2LTOF&IsCharge2UTOF"; 
	EffCorr * GoodUtof_Z2  = new EffCorr(finalhistos,"GoodUtofEffCorr_Z2" ,"GoodUtof Eff. Corr",true,2,before,after,"IsPrimary",false,"IsHeliumMC");

	before = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1&IsCleaning";
        after  = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime"; 
	EffCorr * GoodTime_TOF = new EffCorr(finalhistos,"GoodTimeEffCorr_TOF","GoodTime Eff. Corr",true,1,before,after,"IsPrimary");
	before = "IsPositive&IsPhysTrig&IsBaselineHe&L1LooseCharge2&IsCleaningHe";
        after  = "IsPositive&IsPhysTrig&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsGoodTimeHe"; 
	EffCorr * GoodTimeZ2_TOF = new EffCorr(finalhistos,"GoodTimeEffCorrZ2_TOF","GoodTime Eff. Corr",true,2,before,after,"IsPrimary",false,"IsHeliumMC");

	before = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime";
        after  = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&QualityTOF"; 
	EffCorr * Quality_TOF = new EffCorr(finalhistos,"QualityEffCorr_TOF","Quality TOF Eff. Corr",true,1,before,after,"IsPrimary");
	
	before = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1&IsCleaning";
        after  = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1&IsCleaning";
	EffCorr * RICHEffCorr_NaF = new EffCorr(finalhistos,"RICHCorrection_NaF","RICH Eff. Corr",true,1,before,(after+"&IsFromNaF").c_str(),"IsPrimary");
	EffCorr * RICHEffCorr_Agl = new EffCorr(finalhistos,"RICHCorrection_Agl","RICH Eff. Corr",true,1,before,(after+"&IsFromAgl").c_str(),"IsPrimary");
	before = "IsPositive&IsPhysTrig&IsBaselineHe&L1LooseCharge2&IsCleaningHe";
        after  = "IsPositive&IsPhysTrig&IsBaselineHe&L1LooseCharge2&IsCleaningHe";
	EffCorr * RICHEffCorrZ2_NaF = new EffCorr(finalhistos,"RICHCorrectionZ2_NaF","RICH Eff. Corr",true,2,before,(after+"&IsFromNaF").c_str(),"IsPrimary",false,"IsHeliumMC");
	EffCorr * RICHEffCorrZ2_Agl = new EffCorr(finalhistos,"RICHCorrectionZ2_Agl","RICH Eff. Corr",true,2,before,(after+"&IsFromAgl").c_str(),"IsPrimary",false,"IsHeliumMC");

	before = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1&IsCleaning";
        after  = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1&IsCleaning";
	EffCorr * RICHQualEffCorr_NaF = new EffCorr(finalhistos,"RICHQualCorrection_NaF","RICH Qual. Eff. Corr",true,1,(before+"&IsFromNaF").c_str(),(after+"&IsFromNaF&RICHBDTCut").c_str(),"IsPrimary");
	EffCorr * RICHQualEffCorr_Agl = new EffCorr(finalhistos,"RICHqualCorrection_Agl","RICH Qual. Eff. Corr",true,1,(before+"&IsFromAgl").c_str(),(after+"&IsFromAgl&RICHBDTCut").c_str(),"IsPrimary");
	before = "IsPositive&IsPhysTrig&IsBaselineHe&L1LooseCharge2&IsCleaningHe";
        after  = "IsPositive&IsPhysTrig&IsBaselineHe&L1LooseCharge2&IsCleaningHe";
	EffCorr * RICHQualEffCorrZ2_NaF = new EffCorr(finalhistos,"RICHQualCorrectionZ2_NaF","RICH Qual. Eff. Corr",true,2,(before+"&IsFromNaF").c_str(),(after+"&IsFromNaF&RICHHeCutNaF").c_str(),"IsPrimary",false,"IsHeliumMC");
	EffCorr * RICHQualEffCorrZ2_Agl = new EffCorr(finalhistos,"RICHqualCorrectionZ2_Agl","RICH Qual. Eff. Corr",true,2,(before+"&IsFromAgl").c_str(),(after+"&IsFromAgl&RICHHeCutAgl").c_str(),"IsPrimary",false,"IsHeliumMC");


	//Fragmentation
	Efficiency * Fragmentation_P = new Efficiency(finalhistos,"Fragmentation_P","Fragmentation_P",Quality_TOF->GetBins(),"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&QualityTOF","IsPurePMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&QualityTOF");
	Efficiency * Fragmentation_D = new Efficiency(finalhistos,"Fragmentation_D","Fragmentation_D",Quality_TOF->GetBins(),"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&QualityTOF","IsPureDMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&QualityTOF");


	//Acceptance


	Acceptance * Acceptance_HE     = new Acceptance(finalhistos,"Acceptance_HE"	,"Acceptance","IsProtonMC","IsProtonMC&IsPositive&IsBaseline"	    ,PRB,UnfoldingToF);
	Acceptance * Acceptance_L1HE   = new Acceptance(finalhistos,"Acceptance_L1HE"	,"Acceptance","IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1"	    ,PRB,UnfoldingToF);
	Acceptance * Acceptance_QualHE = new Acceptance(finalhistos,"Acceptance_QualHE" ,"Acceptance","IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",PRB,UnfoldingToF);

	Acceptance * Acceptance_HeQualHE = new Acceptance(finalhistos,"Acceptance_HeQualHE" ,"Acceptance","IsHeliumMC","IsHeliumMC&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe",HeRB,UnfoldingToF_He);

	Acceptance * Acceptance_RigPTOF = new Acceptance(finalhistos,"Acceptance_RigPTOF","Acceptance","IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",GlobalRig.GetToFPBins(),UnfoldingToF);
	Acceptance * Acceptance_RigPNaF = new Acceptance(finalhistos,"Acceptance_RigPNaF","Acceptance","IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",GlobalRig.GetNaFPBins(),UnfoldingNaF);
	Acceptance * Acceptance_RigPAgl = new Acceptance(finalhistos,"Acceptance_RigPAgl","Acceptance","IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",GlobalRig.GetAglPBins(),UnfoldingAgl);

	Acceptance * Acceptance_RigBetaPTOF = new Acceptance(finalhistos,"Acceptance_RigBetaPTOF","Acceptance","IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&QualityTOF",GlobalRig.GetToFPBins(),UnfoldingToF);
	Acceptance * Acceptance_RigBetaPNaF = new Acceptance(finalhistos,"Acceptance_RigBetaPNaF","Acceptance","IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut",GlobalRig.GetNaFPBins(),UnfoldingNaF);
	Acceptance * Acceptance_RigBetaPAgl = new Acceptance(finalhistos,"Acceptance_RigBetaPAgl","Acceptance","IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut",GlobalRig.GetAglPBins(),UnfoldingAgl);

	Acceptance * Acceptance_RigHeTOF = new Acceptance(finalhistos,"Acceptance_RigHeTOF","Acceptance","IsHeliumMC","IsHeliumMC&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe",Global_HeRig.GetToFDBins(),UnfoldingToF_He);
	Acceptance * Acceptance_RigHeNaF = new Acceptance(finalhistos,"Acceptance_RigHeNaF","Acceptance","IsHeliumMC","IsHeliumMC&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe",Global_HeRig.GetNaFDBins(),UnfoldingNaF_He);
	Acceptance * Acceptance_RigHeAgl = new Acceptance(finalhistos,"Acceptance_RigHeAgl","Acceptance","IsHeliumMC","IsHeliumMC&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe",Global_HeRig.GetAglDBins(),UnfoldingAgl_He);

	Acceptance * Acceptance_PTOF = new Acceptance(finalhistos,"Acceptance_PTOF","Acceptance","IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&QualityTOF",Global.GetToFPBins(),UnfoldingToF);
	Acceptance * Acceptance_PNaF = new Acceptance(finalhistos,"Acceptance_PNaF","Acceptance","IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut",Global.GetNaFPBins(),UnfoldingNaF);
	Acceptance * Acceptance_PAgl = new Acceptance(finalhistos,"Acceptance_PAgl","Acceptance","IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut",Global.GetAglPBins(),UnfoldingAgl);

	Acceptance * Acceptance_DTOF = new Acceptance(finalhistos,"Acceptance_DTOF","Acceptance","IsDeutonMC","IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&QualityTOF",Global.GetToFDBins(),UnfoldingToF_D);
	Acceptance * Acceptance_DNaF = new Acceptance(finalhistos,"Acceptance_DNaF","Acceptance","IsDeutonMC","IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut",Global.GetNaFDBins(),UnfoldingNaF_D);
	Acceptance * Acceptance_DAgl = new Acceptance(finalhistos,"Acceptance_DAgl","Acceptance","IsDeutonMC","IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut",Global.GetAglDBins(),UnfoldingAgl_D);

	Acceptance * Acceptance_HeTOF = new Acceptance(finalhistos,"Acceptance_HeTOF","Acceptance","IsHeliumMC","IsHeliumMC&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsGoodTimeHe",Global_He.GetToFDBins(),UnfoldingToF_He);
	Acceptance * Acceptance_HeNaF = new Acceptance(finalhistos,"Acceptance_HeNaF","Acceptance","IsHeliumMC","IsHeliumMC&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromNaF&RICHHeCutNaF",Global_He.GetNaFDBins(),UnfoldingNaF_He);
	Acceptance * Acceptance_HeAgl = new Acceptance(finalhistos,"Acceptance_HeAgl","Acceptance","IsHeliumMC","IsHeliumMC&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromAgl&RICHHeCutAgl",Global_He.GetAglDBins(),UnfoldingAgl_He);

	Acceptance * Acceptance_He3TOF = new Acceptance(finalhistos,"Acceptance_He3TOF","Acceptance","IsHeliumMC","IsHeliumMC&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsGoodTimeHe",Global_He.GetToFPBins(),UnfoldingToF_He);
	Acceptance * Acceptance_He3NaF = new Acceptance(finalhistos,"Acceptance_He3NaF","Acceptance","IsHeliumMC","IsHeliumMC&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromNaF&RICHHeCutNaF",Global_He.GetNaFPBins(),UnfoldingNaF_He);
	Acceptance * Acceptance_He3Agl = new Acceptance(finalhistos,"Acceptance_He3Agl","Acceptance","IsHeliumMC","IsHeliumMC&IsPositive&IsBaselineHe&L1LooseCharge2&IsCleaningHe&IsFromAgl&RICHHeCutAgl",Global_He.GetAglPBins(),UnfoldingAgl_He);


	TriggerEffCorr_HE	->SetDefaultOutFile(finalhistos); 
	TriggerEffCorr_Z2	->SetDefaultOutFile(finalhistos); 
	TriggerFullSpan_HE	->SetDefaultOutFile(finalhistos); 
	L1PickUpEffCorr_HE ->SetDefaultOutFile(finalhistos); 
	L1PickUpGeom_HE ->SetDefaultOutFile(finalhistos); 
	TrackerEffCorr_HE ->SetDefaultOutFile(finalhistos); 
	KalmanEffCorr_HE ->SetDefaultOutFile(finalhistos); 
	GoodChi_HE   ->SetDefaultOutFile(finalhistos); 
	GoodUtof_HE ->SetDefaultOutFile(finalhistos); 
	GoodLtof_HE  ->SetDefaultOutFile(finalhistos); 
	GoodQTrack_HE  ->SetDefaultOutFile(finalhistos); 
	Good1Track_HE  ->SetDefaultOutFile(finalhistos); 
	GoodTime_TOF 	->SetDefaultOutFile(finalhistos); 
	RICHEffCorr_NaF 	->SetDefaultOutFile(finalhistos); 
	RICHEffCorr_Agl 	->SetDefaultOutFile(finalhistos); 
	RICHQualEffCorr_NaF 	->SetDefaultOutFile(finalhistos); 
	RICHQualEffCorr_Agl 	->SetDefaultOutFile(finalhistos); 
	Quality_TOF      	->SetDefaultOutFile(finalhistos); 
	Fragmentation_P	->SetDefaultOutFile(finalhistos); 
	Fragmentation_D	->SetDefaultOutFile(finalhistos); 
	L1PickUpEffCorr_Z2 ->SetDefaultOutFile(finalhistos); 
	L1PickUpGeom_Z2 ->SetDefaultOutFile(finalhistos); 
	TrackerEffCorr_Z2 ->SetDefaultOutFile(finalhistos); 
	KalmanEffCorr_Z2 ->SetDefaultOutFile(finalhistos); 
	GoodChi_Z2   ->SetDefaultOutFile(finalhistos); 
	GoodUtof_Z2 ->SetDefaultOutFile(finalhistos); 
	GoodLtof_Z2  ->SetDefaultOutFile(finalhistos); 
	GoodQTrack_Z2  ->SetDefaultOutFile(finalhistos); 
	Good1Track_Z2  ->SetDefaultOutFile(finalhistos); 
	GoodTimeZ2_TOF 	->SetDefaultOutFile(finalhistos); 
	RICHEffCorrZ2_NaF 	->SetDefaultOutFile(finalhistos); 
	RICHEffCorrZ2_Agl 	->SetDefaultOutFile(finalhistos); 
	RICHQualEffCorrZ2_NaF 	->SetDefaultOutFile(finalhistos); 
	RICHQualEffCorrZ2_Agl 	->SetDefaultOutFile(finalhistos); 
	Quality_TOF      	->SetDefaultOutFile(finalhistos); 
	Acceptance_PTOF 	->SetDefaultOutFile(finalhistos); 
 	Acceptance_PNaF 	->SetDefaultOutFile(finalhistos); 
 	Acceptance_PAgl 	->SetDefaultOutFile(finalhistos); 
	Acceptance_DTOF 	->SetDefaultOutFile(finalhistos); 
 	Acceptance_DNaF 	->SetDefaultOutFile(finalhistos); 
 	Acceptance_DAgl 	->SetDefaultOutFile(finalhistos); 
	Acceptance_HeTOF 	->SetDefaultOutFile(finalhistos); 
 	Acceptance_HeNaF 	->SetDefaultOutFile(finalhistos); 
 	Acceptance_HeAgl 	->SetDefaultOutFile(finalhistos); 
	Acceptance_He3TOF 	->SetDefaultOutFile(finalhistos); 
 	Acceptance_He3NaF 	->SetDefaultOutFile(finalhistos); 
 	Acceptance_He3Agl 	->SetDefaultOutFile(finalhistos); 
	Acceptance_RigPTOF 	->SetDefaultOutFile(finalhistos); 
 	Acceptance_RigPNaF 	->SetDefaultOutFile(finalhistos); 
 	Acceptance_RigPAgl 	->SetDefaultOutFile(finalhistos); 
	Acceptance_RigBetaPTOF 	->SetDefaultOutFile(finalhistos); 
 	Acceptance_RigBetaPNaF 	->SetDefaultOutFile(finalhistos); 
 	Acceptance_RigBetaPAgl 	->SetDefaultOutFile(finalhistos); 
	Acceptance_RigHeTOF 	->SetDefaultOutFile(finalhistos); 
 	Acceptance_RigHeNaF 	->SetDefaultOutFile(finalhistos); 
 	Acceptance_RigHeAgl 	->SetDefaultOutFile(finalhistos); 
	Acceptance_HE    	->SetDefaultOutFile(finalhistos); 
        Acceptance_L1HE  	->SetDefaultOutFile(finalhistos); 
	Acceptance_QualHE	->SetDefaultOutFile(finalhistos); 
	Acceptance_HeQualHE	->SetDefaultOutFile(finalhistos); 



	Filler.AddObject2beFilled(TriggerEffCorr_HE,GetRigidityInner,GetRigidityInner);
	Filler.AddObject2beFilled(TriggerEffCorr_Z2,GetRigidityInner,GetRigidityInner);
	Filler.AddObject2beFilled(TriggerFullSpan_HE,GetRigidityInner,GetRigidityInner);
	Filler.AddObject2beFilled(TrackerEffCorr_HE,GetMomentumProxy,GetMomentumProxy);
	Filler.AddObject2beFilled(GoodQTrack_HE,GetMomentumProxy,GetMomentumProxy);
	Filler.AddObject2beFilled(GoodChi_HE   ,GetMomentumProxy,GetMomentumProxy);	
	Filler.AddObject2beFilled(KalmanEffCorr_HE,GetMomentumProxy,GetMomentumProxy);
	Filler.AddObject2beFilled(L1PickUpGeom_HE,GetRigidityInner,GetRigidityInner);
	Filler.AddObject2beFilled(L1PickUpEffCorr_HE,GetRigidityInner,GetRigidityInner);
	Filler.AddObject2beFilled(GoodUtof_HE  ,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(GoodLtof_HE  ,GetRigidity,GetRigidity);	
	Filler.AddObject2beFilled(Good1Track_HE,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(GoodTime_TOF,GetRigidity,GetRigidity);	
	Filler.AddObject2beFilled(Quality_TOF,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(RICHEffCorr_NaF,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(RICHEffCorr_Agl,GetRigidity,GetRigidity);	
	Filler.AddObject2beFilled(RICHQualEffCorr_NaF,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(RICHQualEffCorr_Agl,GetRigidity,GetRigidity);	

	Filler.AddObject2beFilled(TrackerEffCorr_Z2,GetMomentumProxyHe,GetMomentumProxyHe);
	Filler.AddObject2beFilled(GoodQTrack_Z2,GetMomentumProxyHe,GetMomentumProxyHe);
	Filler.AddObject2beFilled(GoodChi_Z2   ,GetMomentumProxyHe,GetMomentumProxyHe);	
	Filler.AddObject2beFilled(KalmanEffCorr_Z2,GetMomentumProxyHe,GetMomentumProxyHe);
	Filler.AddObject2beFilled(L1PickUpGeom_Z2,GetRigidityInner,GetRigidityInner);
	Filler.AddObject2beFilled(L1PickUpEffCorr_Z2,GetRigidityInner,GetRigidityInner);
	Filler.AddObject2beFilled(GoodUtof_Z2  ,GetRigidity,GetRigidity);
	
	Filler.AddObject2beFilled(GoodLtof_Z2  ,GetRigidity,GetRigidity);	
	Filler.AddObject2beFilled(Good1Track_Z2,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(GoodTimeZ2_TOF,GetRigidity,GetRigidity);	
	Filler.AddObject2beFilled(RICHEffCorrZ2_NaF,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(RICHEffCorrZ2_Agl,GetRigidity,GetRigidity);	
	Filler.AddObject2beFilled(RICHQualEffCorrZ2_NaF,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(RICHQualEffCorrZ2_Agl,GetRigidity,GetRigidity);
	
	Filler.AddObject2beFilled(Acceptance_PTOF,GetSmearedBetaTOF ,GetSmearedBetaTOF );
        Filler.AddObject2beFilled(Acceptance_PNaF,GetSmearedBetaRICH,GetSmearedBetaRICH);
        Filler.AddObject2beFilled(Acceptance_PAgl,GetSmearedBetaRICH,GetSmearedBetaRICH);
	Filler.AddObject2beFilled(Acceptance_DTOF,GetSmearedBetaTOF ,GetSmearedBetaTOF );
        Filler.AddObject2beFilled(Acceptance_DNaF,GetSmearedBetaRICH,GetSmearedBetaRICH);
        Filler.AddObject2beFilled(Acceptance_DAgl,GetSmearedBetaRICH,GetSmearedBetaRICH);
	Filler.AddObject2beFilled(Acceptance_HeTOF,GetSmearedBetaTOF ,GetSmearedBetaTOF );
        Filler.AddObject2beFilled(Acceptance_HeNaF,GetSmearedBetaRICH,GetSmearedBetaRICH);
        Filler.AddObject2beFilled(Acceptance_HeAgl,GetSmearedBetaRICH,GetSmearedBetaRICH);
	Filler.AddObject2beFilled(Acceptance_He3TOF,GetSmearedBetaTOF ,GetSmearedBetaTOF );
        Filler.AddObject2beFilled(Acceptance_He3NaF,GetSmearedBetaRICH,GetSmearedBetaRICH);
        Filler.AddObject2beFilled(Acceptance_He3Agl,GetSmearedBetaRICH,GetSmearedBetaRICH);


	Filler.AddObject2beFilled(Acceptance_RigPTOF, 	GetRigidity,GetRigidity);
        Filler.AddObject2beFilled(Acceptance_RigPNaF,	GetRigidity,GetRigidity);
        Filler.AddObject2beFilled(Acceptance_RigPAgl, 	GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(Acceptance_RigHeTOF, 	GetRigidity,GetRigidity);
        Filler.AddObject2beFilled(Acceptance_RigHeNaF,	GetRigidity,GetRigidity);
        Filler.AddObject2beFilled(Acceptance_RigHeAgl, 	GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(Acceptance_RigBetaPTOF, 	GetRigidity,GetRigidity);
        Filler.AddObject2beFilled(Acceptance_RigBetaPNaF,	GetRigidity,GetRigidity);
        Filler.AddObject2beFilled(Acceptance_RigBetaPAgl, 	GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(Acceptance_HE,    	GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(Acceptance_L1HE , 	GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(Acceptance_QualHE,	GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(Acceptance_HeQualHE,	GetRigidity,GetRigidity);

	Filler.AddObject2beFilled(Fragmentation_P,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(Fragmentation_D,GetGenMomentum,GetGenMomentum);

	Filler.ReinitializeAll(refill);

/*
	FillerAcc.AddObject2beFilled(Acceptance_PTOF,GetSmearedBetaTOF ,GetSmearedBetaTOF );
        FillerAcc.AddObject2beFilled(Acceptance_PNaF,GetSmearedBetaRICH,GetSmearedBetaRICH);
        FillerAcc.AddObject2beFilled(Acceptance_PAgl,GetSmearedBetaRICH,GetSmearedBetaRICH);
	FillerAcc.AddObject2beFilled(Acceptance_DTOF,GetSmearedBetaTOF ,GetSmearedBetaTOF );
        FillerAcc.AddObject2beFilled(Acceptance_DNaF,GetSmearedBetaRICH,GetSmearedBetaRICH);
        FillerAcc.AddObject2beFilled(Acceptance_DAgl,GetSmearedBetaRICH,GetSmearedBetaRICH);
	FillerAcc.AddObject2beFilled(Acceptance_RigPTOF, 	GetSmearedBetaTOF ,GetSmearedBetaTOF );
        FillerAcc.AddObject2beFilled(Acceptance_RigPNaF,	GetSmearedBetaRICH,GetSmearedBetaRICH);
        FillerAcc.AddObject2beFilled(Acceptance_RigPAgl, 	GetSmearedBetaRICH,GetSmearedBetaRICH);
	FillerAcc.AddObject2beFilled(Acceptance_HE,    	GetSmearedBetaTOF ,GetSmearedBetaTOF );
	FillerAcc.AddObject2beFilled(Acceptance_L1HE , 	GetSmearedBetaRICH,GetSmearedBetaRICH);
	FillerAcc.AddObject2beFilled(Acceptance_QualHE,	GetSmearedBetaRICH,GetSmearedBetaRICH);
	
*/

	if(!refill&&checkfile) {	


	
		Acceptance_PTOF->Set_MCPar(0.5,100,1,"Pr.B1200/pr.pl1.05100.4_00.info",filelistMC1.c_str(),0.14699260);	
		Acceptance_PNaF->Set_MCPar(0.5,100,1,"Pr.B1200/pr.pl1.05100.4_00.info",filelistMC1.c_str(),0.14699260);	
		Acceptance_PAgl->Set_MCPar(0.5,100,1,"Pr.B1200/pr.pl1.05100.4_00.info",filelistMC1.c_str(),0.14699260);	

		Acceptance_DTOF->Set_MCPar(0.5,100,1,"D.B1220/d.pl1.05100.info",filelistMC2.c_str(),0.15262600,0.05);	
		Acceptance_DNaF->Set_MCPar(0.5,100,1,"D.B1220/d.pl1.05100.info",filelistMC2.c_str(),0.15262600,0.05);	
		Acceptance_DAgl->Set_MCPar(0.5,100,1,"D.B1220/d.pl1.05100.info",filelistMC2.c_str(),0.15262600,0.05);	

		Acceptance_HeTOF->Set_MCPar(1,1500,1,"He.B1200/he4.pl1.21000.4_00.info",filelistMC3.c_str(),0.1381);	
		Acceptance_HeNaF->Set_MCPar(1,1500,1,"He.B1200/he4.pl1.21000.4_00.info",filelistMC3.c_str(),0.1381);	
		Acceptance_HeAgl->Set_MCPar(1,1500,1,"He.B1200/he4.pl1.21000.4_00.info",filelistMC3.c_str(),0.1381);	

		Acceptance_He3TOF->Set_MCPar(1,1500,1,"He.B1200/he4.pl1.21000.4_00.info",filelistMC3.c_str(),0.1381);	
		Acceptance_He3NaF->Set_MCPar(1,1500,1,"He.B1200/he4.pl1.21000.4_00.info",filelistMC3.c_str(),0.1381);	
		Acceptance_He3Agl->Set_MCPar(1,1500,1,"He.B1200/he4.pl1.21000.4_00.info",filelistMC3.c_str(),0.1381);	


		Acceptance_RigPTOF->Set_MCPar(0.5,100,1,"Pr.B1200/pr.pl1.05100.4_00.info",filelistMC1.c_str(),0.14699260);	
		Acceptance_RigPNaF->Set_MCPar(0.5,100,1,"Pr.B1200/pr.pl1.05100.4_00.info",filelistMC1.c_str(),0.14699260);	
		Acceptance_RigPAgl->Set_MCPar(0.5,100,1,"Pr.B1200/pr.pl1.05100.4_00.info",filelistMC1.c_str(),0.14699260);	
	
		Acceptance_RigBetaPTOF->Set_MCPar(0.5,100,1,"Pr.B1200/pr.pl1.05100.4_00.info",filelistMC1.c_str(),0.14699260);	
		Acceptance_RigBetaPNaF->Set_MCPar(0.5,100,1,"Pr.B1200/pr.pl1.05100.4_00.info",filelistMC1.c_str(),0.14699260);	
		Acceptance_RigBetaPAgl->Set_MCPar(0.5,100,1,"Pr.B1200/pr.pl1.05100.4_00.info",filelistMC1.c_str(),0.14699260);	
	
		Acceptance_RigHeTOF->Set_MCPar(1,500,1,"He.B1200/he4.pl1.21000.4_00.info",filelistMC3.c_str(),0.1381);	
		Acceptance_RigHeNaF->Set_MCPar(1,500,1,"He.B1200/he4.pl1.21000.4_00.info",filelistMC3.c_str(),0.1381);	
		Acceptance_RigHeAgl->Set_MCPar(1,500,1,"He.B1200/he4.pl1.21000.4_00.info",filelistMC3.c_str(),0.1381);	

		Acceptance_HE    ->Set_MCPar(0.5,100,1,"Pr.B1200/pr.pl1.05100.4_00.info",filelistMC1.c_str(),0.14699260);	
		Acceptance_L1HE  ->Set_MCPar(0.5,100,1,"Pr.B1200/pr.pl1.05100.4_00.info",filelistMC1.c_str(),0.14699260);	
		Acceptance_QualHE->Set_MCPar(0.5,100,1,"Pr.B1200/pr.pl1.05100.4_00.info",filelistMC1.c_str(),0.14699260);	
		Acceptance_HeQualHE->Set_MCPar(1,500,1,"He.B1200/he4.pl1.21000.4_00.info",filelistMC3.c_str(),0.1381);	


		Acceptance_PTOF->SetModeled();
		Acceptance_PNaF->SetModeled();
		Acceptance_PAgl->SetModeled();

		Acceptance_DTOF->SetModeled();
		Acceptance_DNaF->SetModeled();
		Acceptance_DAgl->SetModeled();

		Acceptance_HeTOF->SetModeled();
		Acceptance_HeNaF->SetModeled();
		Acceptance_HeAgl->SetModeled();

		Acceptance_He3TOF->SetModeled();
		Acceptance_He3NaF->SetModeled();
		Acceptance_He3Agl->SetModeled();

		Acceptance_RigPTOF->SetModeled();
		Acceptance_RigPNaF->SetModeled();
		Acceptance_RigPAgl->SetModeled();

		Acceptance_RigBetaPTOF->SetModeled();
		Acceptance_RigBetaPNaF->SetModeled();
		Acceptance_RigBetaPAgl->SetModeled();

		Fragmentation_P->Eval_Efficiency();
		Fragmentation_D->Eval_Efficiency();

		//Fragmentation_P->SaveResults(finalresults);
	        //Fragmentation_D->SaveResults(finalresults);

		TFile * EffSystFile = TFile::Open("/afs/cern.ch/work/f/fdimicco/private/Deutons/DirectAnalysis/EffSyst/Time.root");
			
	
		AnalyzeEffCorr(	TriggerEffCorr_HE  , finalhistos, finalresults,EffSystFile,"","","",filelistDT.c_str(),1,true);
		AnalyzeEffCorr(	TriggerEffCorr_Z2  , finalhistos, finalresults,EffSystFile,"","","",filelistDT.c_str(),1.45,true);
		AnalyzeEffCorr(	TriggerFullSpan_HE  , finalhistos, finalresults,EffSystFile,"","","",filelistDT.c_str(),1,true);
		AnalyzeEffCorr(	L1PickUpEffCorr_HE , finalhistos, finalresults,EffSystFile,"","","",filelistDT.c_str(),1);
		AnalyzeEffCorr(	L1PickUpGeom_HE , finalhistos, finalresults,EffSystFile,"","","",filelistDT.c_str(),1);
		AnalyzeEffCorr(	TrackerEffCorr_HE , finalhistos, finalresults,EffSystFile,"","","",filelistDT.c_str(),3.);
		AnalyzeEffCorr(	KalmanEffCorr_HE , finalhistos, finalresults,EffSystFile,"","","",filelistDT.c_str(),1);
		AnalyzeEffCorr(	GoodChi_HE  , finalhistos, finalresults,EffSystFile,"Eff. Corrections Sys/Track Quality/","GoodChi_Err"		,"GoodChi_timeavg"	,filelistDT.c_str(),0.8);  
		AnalyzeEffCorr(	GoodUtof_HE , finalhistos, finalresults,EffSystFile,"Eff. Corrections Sys/NoInteraction/","GoodUtof_Err"	,"GoodUtof_timeavg"	,filelistDT.c_str(),0.65);
		AnalyzeEffCorr(	GoodLtof_HE  , finalhistos, finalresults ,EffSystFile,"Eff. Corrections Sys/NoInteraction/","GoodLtof_Err"	,"GoodLtof_timeavg"	,filelistDT.c_str(),0.65); 
		AnalyzeEffCorr(	GoodQTrack_HE , finalhistos, finalresults,EffSystFile,"Eff. Corrections Sys/Track Quality/","GoodQTrack_Err"	,"GoodQTrack_timeavg"	,filelistDT.c_str(),0.8);
		AnalyzeEffCorr(	Good1Track_HE , finalhistos, finalresults,EffSystFile,"Eff. Corrections Sys/NoInteraction/","Good1Track_Err"	,"Good1Track_timeavg"	,filelistDT.c_str(),0.65); 
		AnalyzeEffCorr(	GoodTime_TOF , finalhistos, finalresults,EffSystFile, "Eff. Corrections Sys/TOF/","GoodTime_Err"		,"GoodTime_timeavg"	,filelistDT.c_str(),0.65);
		AnalyzeEffCorr(	Quality_TOF , finalhistos, finalresults,EffSystFile,    "Eff. Corrections Sys/TOF/","QualityTime_Err"		,"Quality_timeavg"	,filelistDT.c_str(),0.65);
		AnalyzeEffCorr(	RICHEffCorr_NaF , finalhistos, finalresults,EffSystFile,"Eff. Corrections Sys/RICH CIEMAT/","RICHNTime_Err"	,"RICHTF1_Ntimeavg"	,filelistDT.c_str(),1);
		AnalyzeEffCorr(	RICHEffCorr_Agl , finalhistos, finalresults,EffSystFile,"Eff. Corrections Sys/RICH CIEMAT/","RICHATime_Err"	,"RICHTF1_Atimeavg"	,filelistDT.c_str(),3.2);
		AnalyzeEffCorr(	RICHQualEffCorr_NaF , finalhistos, finalresults,EffSystFile,"Eff. Corrections Sys/RICH BDT/","RICHQNTime_Err"	,"RICHQualTF1_Ntimeavg",filelistDT.c_str(),1);
		AnalyzeEffCorr(	RICHQualEffCorr_Agl , finalhistos, finalresults,EffSystFile,"Eff. Corrections Sys/RICH BDT/","RICHQATime_Err"	,"RICHQualTF1_Atimeavg",filelistDT.c_str(),3.2);

		AnalyzeEffCorr(	L1PickUpEffCorr_Z2 , finalhistos, finalresults,EffSystFile,"","","",filelistDT.c_str(),1.05);
		AnalyzeEffCorr(	L1PickUpGeom_Z2 , finalhistos, finalresults,EffSystFile,"","","",filelistDT.c_str(),1.05);
		AnalyzeEffCorr(	TrackerEffCorr_Z2 , finalhistos, finalresults,EffSystFile,"","","",filelistDT.c_str(),1.25);
		AnalyzeEffCorr(	KalmanEffCorr_Z2 , finalhistos, finalresults,EffSystFile,"","","",filelistDT.c_str(),1.25);
		AnalyzeEffCorr(	GoodChi_Z2  , finalhistos, finalresults,EffSystFile,"Eff. Corrections Sys. He/Track Quality/","GoodChi_Err"		,"GoodChi_timeavg"	,filelistDT.c_str(),1.25);  
		AnalyzeEffCorr(	GoodUtof_Z2 , finalhistos, finalresults,EffSystFile,"Eff. Corrections Sys. He/NoInteraction/","GoodUtof_Err"		,"GoodUtof_timeavg"	,filelistDT.c_str(),0.35);
		AnalyzeEffCorr(	GoodLtof_Z2  , finalhistos, finalresults ,EffSystFile,"Eff. Corrections Sys. He/NoInteraction/","GoodLtof_Err"		,"GoodLtof_timeavg"	,filelistDT.c_str(),0.35); 
		AnalyzeEffCorr(	GoodQTrack_Z2 , finalhistos, finalresults,EffSystFile,"Eff. Corrections Sys. He/Track Quality/","GoodQTrack_Err"	,"GoodQTrack_timeavg"	,filelistDT.c_str(),1.25);
		AnalyzeEffCorr(	Good1Track_Z2 , finalhistos, finalresults,EffSystFile,"Eff. Corrections Sys. He/NoInteraction/","Good1Track_Err"	,"Good1Track_timeavg"	,filelistDT.c_str(),0.23); 
		AnalyzeEffCorr(	GoodTimeZ2_TOF , finalhistos, finalresults,EffSystFile, "Eff. Corrections Sys. He/TOF/","GoodTime_Err"			,"GoodTime_timeavg"	,filelistDT.c_str(),0.3);
		AnalyzeEffCorr(	RICHEffCorrZ2_NaF , finalhistos, finalresults,EffSystFile,"Eff. Corrections Sys. He/RICH CIEMAT/","RICHNTime_Err"	,"RICHTF1_Ntimeavg"    ,filelistDT.c_str(),1.05);
		AnalyzeEffCorr(	RICHEffCorrZ2_Agl , finalhistos, finalresults,EffSystFile,"Eff. Corrections Sys. He/RICH CIEMAT/","RICHATime_Err"	,"RICHTF1_Atimeavg"    ,filelistDT.c_str(),3.05);
		AnalyzeEffCorr(	RICHQualEffCorrZ2_NaF , finalhistos, finalresults,EffSystFile,"Eff. Corrections Sys. He/RICH BDT/","RICHQNTime_Err"	,"RICHQualTF1_Ntimeavg",filelistDT.c_str(),0.9);
		AnalyzeEffCorr(	RICHQualEffCorrZ2_Agl , finalhistos, finalresults,EffSystFile,"Eff. Corrections Sys. He/RICH BDT/","RICHQATime_Err"	,"RICHQualTF1_Atimeavg",filelistDT.c_str(),3.05);


		//Baseline effcorr
		Acceptance_HE->ApplyEfficCorr(TrackerEffCorr_HE);
		Acceptance_HE->ApplyEfficCorr(GoodQTrack_HE);
		Acceptance_HE->ApplyEfficCorr(GoodChi_HE);
		Acceptance_HE->ApplyEfficCorr(KalmanEffCorr_HE);
	
		Acceptance_L1HE->ApplyEfficCorr(TrackerEffCorr_HE);
		Acceptance_L1HE->ApplyEfficCorr(GoodQTrack_HE);
		Acceptance_L1HE->ApplyEfficCorr(GoodChi_HE);
		Acceptance_L1HE->ApplyEfficCorr(KalmanEffCorr_HE);
	
		Acceptance_QualHE->ApplyEfficCorr(TrackerEffCorr_HE);
		Acceptance_QualHE->ApplyEfficCorr(GoodQTrack_HE);
		Acceptance_QualHE->ApplyEfficCorr(GoodChi_HE);
		Acceptance_QualHE->ApplyEfficCorr(KalmanEffCorr_HE);

		Acceptance_RigPTOF->ApplyEfficCorr(TrackerEffCorr_HE);
		Acceptance_RigPTOF->ApplyEfficCorr(GoodQTrack_HE);
		Acceptance_RigPTOF->ApplyEfficCorr(GoodChi_HE);
		Acceptance_RigPTOF->ApplyEfficCorr(KalmanEffCorr_HE);

		Acceptance_RigPNaF->ApplyEfficCorr(TrackerEffCorr_HE);
		Acceptance_RigPNaF->ApplyEfficCorr(GoodQTrack_HE);
		Acceptance_RigPNaF->ApplyEfficCorr(GoodChi_HE);
		Acceptance_RigPNaF->ApplyEfficCorr(KalmanEffCorr_HE);

		Acceptance_RigPAgl->ApplyEfficCorr(TrackerEffCorr_HE);
		Acceptance_RigPAgl->ApplyEfficCorr(GoodQTrack_HE);
		Acceptance_RigPAgl->ApplyEfficCorr(GoodChi_HE);
		Acceptance_RigPAgl->ApplyEfficCorr(KalmanEffCorr_HE);

		Acceptance_RigBetaPTOF->ApplyEfficCorr(TrackerEffCorr_HE);
		Acceptance_RigBetaPTOF->ApplyEfficCorr(GoodQTrack_HE);
		Acceptance_RigBetaPTOF->ApplyEfficCorr(GoodChi_HE);
		Acceptance_RigBetaPTOF->ApplyEfficCorr(KalmanEffCorr_HE);

		Acceptance_RigBetaPNaF->ApplyEfficCorr(TrackerEffCorr_HE);
		Acceptance_RigBetaPNaF->ApplyEfficCorr(GoodQTrack_HE);
		Acceptance_RigBetaPNaF->ApplyEfficCorr(GoodChi_HE);
		Acceptance_RigBetaPNaF->ApplyEfficCorr(KalmanEffCorr_HE);

		Acceptance_RigBetaPAgl->ApplyEfficCorr(TrackerEffCorr_HE);
		Acceptance_RigBetaPAgl->ApplyEfficCorr(GoodQTrack_HE);
		Acceptance_RigBetaPAgl->ApplyEfficCorr(GoodChi_HE);
		Acceptance_RigBetaPAgl->ApplyEfficCorr(KalmanEffCorr_HE);

		Acceptance_PTOF->ApplyEfficCorr(TrackerEffCorr_HE);
		Acceptance_PTOF->ApplyEfficCorr(GoodQTrack_HE);
		Acceptance_PTOF->ApplyEfficCorr(GoodChi_HE);
		Acceptance_PTOF->ApplyEfficCorr(KalmanEffCorr_HE);
	
		Acceptance_PNaF->ApplyEfficCorr(TrackerEffCorr_HE);
		Acceptance_PNaF->ApplyEfficCorr(GoodQTrack_HE);
		Acceptance_PNaF->ApplyEfficCorr(GoodChi_HE);
		Acceptance_PNaF->ApplyEfficCorr(KalmanEffCorr_HE);
	
		Acceptance_PAgl->ApplyEfficCorr(TrackerEffCorr_HE);
		Acceptance_PAgl->ApplyEfficCorr(GoodQTrack_HE);
		Acceptance_PAgl->ApplyEfficCorr(GoodChi_HE);
		Acceptance_PAgl->ApplyEfficCorr(KalmanEffCorr_HE);

		Acceptance_DTOF->ApplyEfficCorr(TrackerEffCorr_HE);
		Acceptance_DTOF->ApplyEfficCorr(GoodQTrack_HE);
		Acceptance_DTOF->ApplyEfficCorr(GoodChi_HE);
		Acceptance_DTOF->ApplyEfficCorr(KalmanEffCorr_HE);
	
		Acceptance_DNaF->ApplyEfficCorr(TrackerEffCorr_HE);
		Acceptance_DNaF->ApplyEfficCorr(GoodQTrack_HE);
		Acceptance_DNaF->ApplyEfficCorr(GoodChi_HE);
		Acceptance_DNaF->ApplyEfficCorr(KalmanEffCorr_HE);
	
		Acceptance_DAgl->ApplyEfficCorr(TrackerEffCorr_HE);
		Acceptance_DAgl->ApplyEfficCorr(GoodQTrack_HE);
		Acceptance_DAgl->ApplyEfficCorr(GoodChi_HE);
		Acceptance_DAgl->ApplyEfficCorr(KalmanEffCorr_HE);

		Acceptance_HeTOF->ApplyEfficCorr(TrackerEffCorr_Z2);
		Acceptance_HeTOF->ApplyEfficCorr(GoodQTrack_Z2);
		Acceptance_HeTOF->ApplyEfficCorr(GoodChi_Z2);
		Acceptance_HeTOF->ApplyEfficCorr(KalmanEffCorr_Z2);
		Acceptance_HeNaF->ApplyEfficCorr(TrackerEffCorr_Z2);
		Acceptance_HeNaF->ApplyEfficCorr(GoodQTrack_Z2);
		Acceptance_HeNaF->ApplyEfficCorr(GoodChi_Z2);
		Acceptance_HeNaF->ApplyEfficCorr(KalmanEffCorr_Z2);
		Acceptance_HeAgl->ApplyEfficCorr(TrackerEffCorr_Z2);
		Acceptance_HeAgl->ApplyEfficCorr(GoodQTrack_Z2);
		Acceptance_HeAgl->ApplyEfficCorr(GoodChi_Z2);
		Acceptance_HeAgl->ApplyEfficCorr(KalmanEffCorr_Z2);

		Acceptance_He3TOF->ApplyEfficCorr(TrackerEffCorr_Z2);
		Acceptance_He3TOF->ApplyEfficCorr(GoodQTrack_Z2);
		Acceptance_He3TOF->ApplyEfficCorr(GoodChi_Z2);
		Acceptance_He3TOF->ApplyEfficCorr(KalmanEffCorr_Z2);
		Acceptance_He3NaF->ApplyEfficCorr(TrackerEffCorr_Z2);
		Acceptance_He3NaF->ApplyEfficCorr(GoodQTrack_Z2);
		Acceptance_He3NaF->ApplyEfficCorr(GoodChi_Z2);
		Acceptance_He3NaF->ApplyEfficCorr(KalmanEffCorr_Z2);
		Acceptance_He3Agl->ApplyEfficCorr(TrackerEffCorr_Z2);
		Acceptance_He3Agl->ApplyEfficCorr(GoodQTrack_Z2);
		Acceptance_He3Agl->ApplyEfficCorr(GoodChi_Z2);
		Acceptance_He3Agl->ApplyEfficCorr(KalmanEffCorr_Z2);


		Acceptance_RigHeTOF->ApplyEfficCorr(TrackerEffCorr_Z2);
		Acceptance_RigHeTOF->ApplyEfficCorr(GoodQTrack_Z2);
		Acceptance_RigHeTOF->ApplyEfficCorr(GoodChi_Z2);
		Acceptance_RigHeTOF->ApplyEfficCorr(KalmanEffCorr_Z2);
		Acceptance_RigHeNaF->ApplyEfficCorr(TrackerEffCorr_Z2);
		Acceptance_RigHeNaF->ApplyEfficCorr(GoodQTrack_Z2);
		Acceptance_RigHeNaF->ApplyEfficCorr(GoodChi_Z2);
		Acceptance_RigHeNaF->ApplyEfficCorr(KalmanEffCorr_Z2);
		Acceptance_RigHeAgl->ApplyEfficCorr(TrackerEffCorr_Z2);
		Acceptance_RigHeAgl->ApplyEfficCorr(GoodQTrack_Z2);
		Acceptance_RigHeAgl->ApplyEfficCorr(GoodChi_Z2);
		Acceptance_RigHeAgl->ApplyEfficCorr(KalmanEffCorr_Z2);


		//L1 effcorr
		Acceptance_L1HE->ApplyEfficCorr(L1PickUpGeom_HE);
		Acceptance_L1HE->ApplyEfficCorr(L1PickUpEffCorr_HE);
		Acceptance_QualHE->ApplyEfficCorr(L1PickUpGeom_HE);
		Acceptance_QualHE->ApplyEfficCorr(L1PickUpEffCorr_HE);

		Acceptance_RigPTOF->ApplyEfficCorr(L1PickUpGeom_HE);
		Acceptance_RigPTOF->ApplyEfficCorr(L1PickUpEffCorr_HE);
		Acceptance_RigPNaF->ApplyEfficCorr(L1PickUpGeom_HE);
		Acceptance_RigPNaF->ApplyEfficCorr(L1PickUpEffCorr_HE);
		Acceptance_RigPAgl->ApplyEfficCorr(L1PickUpGeom_HE);
		Acceptance_RigPAgl->ApplyEfficCorr(L1PickUpEffCorr_HE);

		Acceptance_RigBetaPTOF->ApplyEfficCorr(L1PickUpGeom_HE);
		Acceptance_RigBetaPTOF->ApplyEfficCorr(L1PickUpEffCorr_HE);
		Acceptance_RigBetaPNaF->ApplyEfficCorr(L1PickUpGeom_HE);
		Acceptance_RigBetaPNaF->ApplyEfficCorr(L1PickUpEffCorr_HE);
		Acceptance_RigBetaPAgl->ApplyEfficCorr(L1PickUpGeom_HE);
		Acceptance_RigBetaPAgl->ApplyEfficCorr(L1PickUpEffCorr_HE);

		Acceptance_PTOF->ApplyEfficCorr(L1PickUpGeom_HE);
		Acceptance_PTOF->ApplyEfficCorr(L1PickUpEffCorr_HE);
		Acceptance_PNaF->ApplyEfficCorr(L1PickUpGeom_HE);
		Acceptance_PNaF->ApplyEfficCorr(L1PickUpEffCorr_HE);
		Acceptance_PAgl->ApplyEfficCorr(L1PickUpGeom_HE);
		Acceptance_PAgl->ApplyEfficCorr(L1PickUpEffCorr_HE);

		Acceptance_DTOF->ApplyEfficCorr(L1PickUpGeom_HE);
		Acceptance_DTOF->ApplyEfficCorr(L1PickUpEffCorr_HE);
		Acceptance_DNaF->ApplyEfficCorr(L1PickUpGeom_HE);
		Acceptance_DNaF->ApplyEfficCorr(L1PickUpEffCorr_HE);
		Acceptance_DAgl->ApplyEfficCorr(L1PickUpGeom_HE);
		Acceptance_DAgl->ApplyEfficCorr(L1PickUpEffCorr_HE);

		Acceptance_HeTOF->ApplyEfficCorr(L1PickUpGeom_Z2);
		Acceptance_HeTOF->ApplyEfficCorr(L1PickUpEffCorr_Z2);
		Acceptance_HeNaF->ApplyEfficCorr(L1PickUpGeom_Z2);
		Acceptance_HeNaF->ApplyEfficCorr(L1PickUpEffCorr_Z2);
		Acceptance_HeAgl->ApplyEfficCorr(L1PickUpGeom_Z2);
		Acceptance_HeAgl->ApplyEfficCorr(L1PickUpEffCorr_Z2);

		Acceptance_He3TOF->ApplyEfficCorr(L1PickUpGeom_Z2);
		Acceptance_He3TOF->ApplyEfficCorr(L1PickUpEffCorr_Z2);
		Acceptance_He3NaF->ApplyEfficCorr(L1PickUpGeom_Z2);
		Acceptance_He3NaF->ApplyEfficCorr(L1PickUpEffCorr_Z2);
		Acceptance_He3Agl->ApplyEfficCorr(L1PickUpGeom_Z2);
		Acceptance_He3Agl->ApplyEfficCorr(L1PickUpEffCorr_Z2);


		Acceptance_RigHeTOF->ApplyEfficCorr(L1PickUpGeom_Z2);
		Acceptance_RigHeTOF->ApplyEfficCorr(L1PickUpEffCorr_Z2);
		Acceptance_RigHeNaF->ApplyEfficCorr(L1PickUpGeom_Z2);
		Acceptance_RigHeNaF->ApplyEfficCorr(L1PickUpEffCorr_Z2);
		Acceptance_RigHeAgl->ApplyEfficCorr(L1PickUpGeom_Z2);
		Acceptance_RigHeAgl->ApplyEfficCorr(L1PickUpEffCorr_Z2);


		//Quality effcorr
		Acceptance_QualHE->ApplyEfficCorr(Good1Track_HE);	
		Acceptance_QualHE->ApplyEfficCorr(GoodUtof_HE);	
		Acceptance_QualHE->ApplyEfficCorr(GoodLtof_HE);	
		
		Acceptance_RigPTOF->ApplyEfficCorr(GoodUtof_HE);
		Acceptance_RigPNaF->ApplyEfficCorr(GoodUtof_HE);
		Acceptance_RigPAgl->ApplyEfficCorr(GoodUtof_HE);
		Acceptance_RigPTOF->ApplyEfficCorr(GoodLtof_HE);
		Acceptance_RigPNaF->ApplyEfficCorr(GoodLtof_HE);
		Acceptance_RigPAgl->ApplyEfficCorr(GoodLtof_HE);
		Acceptance_RigPTOF->ApplyEfficCorr(Good1Track_HE);
		Acceptance_RigPNaF->ApplyEfficCorr(Good1Track_HE);
		Acceptance_RigPAgl->ApplyEfficCorr(Good1Track_HE);
		
		Acceptance_RigBetaPTOF->ApplyEfficCorr(GoodUtof_HE);
		Acceptance_RigBetaPNaF->ApplyEfficCorr(GoodUtof_HE);
		Acceptance_RigBetaPAgl->ApplyEfficCorr(GoodUtof_HE);
		Acceptance_RigBetaPTOF->ApplyEfficCorr(GoodLtof_HE);
		Acceptance_RigBetaPNaF->ApplyEfficCorr(GoodLtof_HE);
		Acceptance_RigBetaPAgl->ApplyEfficCorr(GoodLtof_HE);
		Acceptance_RigBetaPTOF->ApplyEfficCorr(Good1Track_HE);
		Acceptance_RigBetaPNaF->ApplyEfficCorr(Good1Track_HE);
		Acceptance_RigBetaPAgl->ApplyEfficCorr(Good1Track_HE);

		Acceptance_PTOF->ApplyEfficCorr(GoodUtof_HE);
		Acceptance_PNaF->ApplyEfficCorr(GoodUtof_HE);
		Acceptance_PAgl->ApplyEfficCorr(GoodUtof_HE);
		Acceptance_PTOF->ApplyEfficCorr(GoodLtof_HE);
		Acceptance_PNaF->ApplyEfficCorr(GoodLtof_HE);
		Acceptance_PAgl->ApplyEfficCorr(GoodLtof_HE);
		Acceptance_PTOF->ApplyEfficCorr(Good1Track_HE);
		Acceptance_PNaF->ApplyEfficCorr(Good1Track_HE);
		Acceptance_PAgl->ApplyEfficCorr(Good1Track_HE);
	
		Acceptance_DTOF->ApplyEfficCorr(GoodUtof_HE);
		Acceptance_DNaF->ApplyEfficCorr(GoodUtof_HE);
		Acceptance_DAgl->ApplyEfficCorr(GoodUtof_HE);
		Acceptance_DTOF->ApplyEfficCorr(GoodLtof_HE);
		Acceptance_DNaF->ApplyEfficCorr(GoodLtof_HE);
		Acceptance_DAgl->ApplyEfficCorr(GoodLtof_HE);
		Acceptance_DTOF->ApplyEfficCorr(Good1Track_HE);
		Acceptance_DNaF->ApplyEfficCorr(Good1Track_HE);
		Acceptance_DAgl->ApplyEfficCorr(Good1Track_HE);
	
		Acceptance_HeTOF->ApplyEfficCorr(GoodUtof_Z2);
		Acceptance_HeNaF->ApplyEfficCorr(GoodUtof_Z2);
		Acceptance_HeAgl->ApplyEfficCorr(GoodUtof_Z2);
		Acceptance_HeTOF->ApplyEfficCorr(GoodLtof_Z2);
		Acceptance_HeNaF->ApplyEfficCorr(GoodLtof_Z2);
		Acceptance_HeAgl->ApplyEfficCorr(GoodLtof_Z2);
		Acceptance_HeTOF->ApplyEfficCorr(Good1Track_Z2);
		Acceptance_HeNaF->ApplyEfficCorr(Good1Track_Z2);
		Acceptance_HeAgl->ApplyEfficCorr(Good1Track_Z2);
	
		Acceptance_He3TOF->ApplyEfficCorr(GoodUtof_Z2);
		Acceptance_He3NaF->ApplyEfficCorr(GoodUtof_Z2);
		Acceptance_He3Agl->ApplyEfficCorr(GoodUtof_Z2);
		Acceptance_He3TOF->ApplyEfficCorr(GoodLtof_Z2);
		Acceptance_He3NaF->ApplyEfficCorr(GoodLtof_Z2);
		Acceptance_He3Agl->ApplyEfficCorr(GoodLtof_Z2);
		Acceptance_He3TOF->ApplyEfficCorr(Good1Track_Z2);
		Acceptance_He3NaF->ApplyEfficCorr(Good1Track_Z2);
		Acceptance_He3Agl->ApplyEfficCorr(Good1Track_Z2);

		Acceptance_RigHeTOF->ApplyEfficCorr(GoodUtof_Z2);
		Acceptance_RigHeNaF->ApplyEfficCorr(GoodUtof_Z2);
		Acceptance_RigHeAgl->ApplyEfficCorr(GoodUtof_Z2);
		Acceptance_RigHeTOF->ApplyEfficCorr(GoodLtof_Z2);
		Acceptance_RigHeNaF->ApplyEfficCorr(GoodLtof_Z2);
		Acceptance_RigHeAgl->ApplyEfficCorr(GoodLtof_Z2);
		Acceptance_RigHeTOF->ApplyEfficCorr(Good1Track_Z2);
		Acceptance_RigHeNaF->ApplyEfficCorr(Good1Track_Z2);
		Acceptance_RigHeAgl->ApplyEfficCorr(Good1Track_Z2);
		
			
		//Trigger effcorr
		Acceptance_L1HE	 	->ApplyEfficCorr(TriggerEffCorr_HE);
		Acceptance_QualHE	->ApplyEfficCorr(TriggerEffCorr_HE);

		Acceptance_RigPTOF	->ApplyEfficCorr(TriggerEffCorr_HE);
		Acceptance_RigPNaF	->ApplyEfficCorr(TriggerEffCorr_HE);
		Acceptance_RigPAgl	->ApplyEfficCorr(TriggerEffCorr_HE);
		Acceptance_RigBetaPTOF	->ApplyEfficCorr(TriggerEffCorr_HE);
		Acceptance_RigBetaPNaF	->ApplyEfficCorr(TriggerEffCorr_HE);
		Acceptance_RigBetaPAgl	->ApplyEfficCorr(TriggerEffCorr_HE);
		Acceptance_PTOF	->ApplyEfficCorr(TriggerEffCorr_HE);
		Acceptance_PNaF	->ApplyEfficCorr(TriggerEffCorr_HE);
		Acceptance_PAgl	->ApplyEfficCorr(TriggerEffCorr_HE);
		Acceptance_DTOF	->ApplyEfficCorr(TriggerEffCorr_HE);
		Acceptance_DNaF	->ApplyEfficCorr(TriggerEffCorr_HE);
		Acceptance_DAgl	->ApplyEfficCorr(TriggerEffCorr_HE);
		Acceptance_HeTOF	->ApplyEfficCorr(TriggerEffCorr_Z2);
		Acceptance_HeNaF	->ApplyEfficCorr(TriggerEffCorr_Z2);
		Acceptance_HeAgl	->ApplyEfficCorr(TriggerEffCorr_Z2);
		Acceptance_He3TOF	->ApplyEfficCorr(TriggerEffCorr_Z2);
		Acceptance_He3NaF	->ApplyEfficCorr(TriggerEffCorr_Z2);
		Acceptance_He3Agl	->ApplyEfficCorr(TriggerEffCorr_Z2);
		Acceptance_RigHeTOF	->ApplyEfficCorr(TriggerEffCorr_Z2);
		Acceptance_RigHeNaF	->ApplyEfficCorr(TriggerEffCorr_Z2);
		Acceptance_RigHeAgl	->ApplyEfficCorr(TriggerEffCorr_Z2);
				
		//velocity effcorr
	
		Acceptance_PTOF->ApplyEfficCorr(GoodTime_TOF);
		Acceptance_PTOF->ApplyEfficCorr(Quality_TOF);
	        Acceptance_PNaF->ApplyEfficCorr(RICHEffCorr_NaF);
                Acceptance_PAgl->ApplyEfficCorr(RICHEffCorr_Agl);
	        Acceptance_PNaF->ApplyEfficCorr(RICHQualEffCorr_NaF);
                Acceptance_PAgl->ApplyEfficCorr(RICHQualEffCorr_Agl);
	
		Acceptance_RigBetaPTOF->ApplyEfficCorr(GoodTime_TOF);
		Acceptance_RigBetaPTOF->ApplyEfficCorr(Quality_TOF);
	        Acceptance_RigBetaPNaF->ApplyEfficCorr(RICHEffCorr_NaF);
                Acceptance_RigBetaPAgl->ApplyEfficCorr(RICHEffCorr_Agl);
	        Acceptance_RigBetaPNaF->ApplyEfficCorr(RICHQualEffCorr_NaF);
                Acceptance_RigBetaPAgl->ApplyEfficCorr(RICHQualEffCorr_Agl);

		Acceptance_DTOF->ApplyEfficCorr(GoodTime_TOF);
		Acceptance_DTOF->ApplyEfficCorr(Quality_TOF);

	        Acceptance_DNaF->ApplyEfficCorr(RICHEffCorr_NaF);
                Acceptance_DAgl->ApplyEfficCorr(RICHEffCorr_Agl);//,0.2);

	        Acceptance_DNaF->ApplyEfficCorr(RICHQualEffCorr_NaF);
          	Acceptance_DAgl->ApplyEfficCorr(RICHQualEffCorr_Agl);//,0.2);
	

		Acceptance_HeTOF->ApplyEfficCorr(GoodTimeZ2_TOF);
	        Acceptance_HeNaF->ApplyEfficCorr(RICHEffCorrZ2_NaF);
                Acceptance_HeAgl->ApplyEfficCorr(RICHEffCorrZ2_Agl);
	        Acceptance_HeNaF->ApplyEfficCorr(RICHQualEffCorrZ2_NaF);
                Acceptance_HeAgl->ApplyEfficCorr(RICHQualEffCorrZ2_Agl);
	
		Acceptance_He3TOF->ApplyEfficCorr(GoodTimeZ2_TOF);
	        Acceptance_He3NaF->ApplyEfficCorr(RICHEffCorrZ2_NaF);
                Acceptance_He3Agl->ApplyEfficCorr(RICHEffCorrZ2_Agl);
	        Acceptance_He3NaF->ApplyEfficCorr(RICHQualEffCorrZ2_NaF);
                Acceptance_He3Agl->ApplyEfficCorr(RICHQualEffCorrZ2_Agl);
	


		Acceptance_PTOF 	->  EvalEffAcc(timeindex,1.2);
		Acceptance_PNaF 	->  EvalEffAcc(timeindex,1.2);
		Acceptance_PAgl 	->  EvalEffAcc(timeindex,1.2);
		
	Acceptance_PTOF 	->SaveResults(finalresults);
        Acceptance_PNaF 	->SaveResults(finalresults);
        Acceptance_PAgl 	->SaveResults(finalresults);

	Acceptance_DTOF 	->  EvalEffAcc(timeindex,1.2);
        Acceptance_DNaF 	->  EvalEffAcc(timeindex,1.2);
        Acceptance_DAgl 	->  EvalEffAcc(timeindex,1.2);

	Acceptance_DTOF 	->SaveResults(finalresults);
        Acceptance_DNaF 	->SaveResults(finalresults);
        Acceptance_DAgl 	->SaveResults(finalresults);

	Acceptance_HeTOF 	->  EvalEffAcc(timeindex,1.2,true);
        Acceptance_HeNaF 	->  EvalEffAcc(timeindex,1.2,true);
        Acceptance_HeAgl 	->  EvalEffAcc(timeindex,1.2,true);
		
	Acceptance_HeTOF 	->SaveResults(finalresults);
        Acceptance_HeNaF 	->SaveResults(finalresults);
        Acceptance_HeAgl 	->SaveResults(finalresults);

	Acceptance_He3TOF 	->  EvalEffAcc(timeindex,1.2,true);
        Acceptance_He3NaF 	->  EvalEffAcc(timeindex,1.2,true);
        Acceptance_He3Agl 	->  EvalEffAcc(timeindex,1.2,true);
		
	Acceptance_He3TOF 	->SaveResults(finalresults);
        Acceptance_He3NaF 	->SaveResults(finalresults);
        Acceptance_He3Agl 	->SaveResults(finalresults);

	Acceptance_RigPTOF 	->  EvalEffAcc(timeindex,1.2);
        Acceptance_RigPNaF 	->  EvalEffAcc(timeindex,1.2);
        Acceptance_RigPAgl 	->  EvalEffAcc(timeindex,1.2);
	
	Acceptance_RigPTOF 	->SaveResults(finalresults);
        Acceptance_RigPNaF 	->SaveResults(finalresults);
        Acceptance_RigPAgl 	->SaveResults(finalresults);

	Acceptance_RigBetaPTOF 	->  EvalEffAcc(timeindex,1.2);
        Acceptance_RigBetaPNaF 	->  EvalEffAcc(timeindex,1.2);
        Acceptance_RigBetaPAgl 	->  EvalEffAcc(timeindex,1.2);
	
	Acceptance_RigBetaPTOF 	->SaveResults(finalresults);
        Acceptance_RigBetaPNaF 	->SaveResults(finalresults);
        Acceptance_RigBetaPAgl 	->SaveResults(finalresults);


	Acceptance_RigHeTOF 	->  EvalEffAcc(timeindex,1.2);
        Acceptance_RigHeNaF 	->  EvalEffAcc(timeindex,1.2);
        Acceptance_RigHeAgl 	->  EvalEffAcc(timeindex,1.2);
	
	Acceptance_RigHeTOF 	->SaveResults(finalresults);
        Acceptance_RigHeNaF 	->SaveResults(finalresults);
        Acceptance_RigHeAgl 	->SaveResults(finalresults);


	Acceptance_HE     	->  EvalEffAcc(timeindex,1.2);
        Acceptance_L1HE   	->  EvalEffAcc(timeindex,1.2);
        Acceptance_QualHE 	->  EvalEffAcc(timeindex,1.2);
	
	Acceptance_HE     	->SaveResults(finalresults);
        Acceptance_L1HE   	->SaveResults(finalresults);
        Acceptance_QualHE 	->SaveResults(finalresults);




	//Plotting EffCorr
		
	//Z1
	std::vector<EffCorr *> Trigg;
	Trigg.push_back(TriggerEffCorr_HE);
	DrawBlockCascade(Trigg,"Z1/Trigger",finalresults,0.1,40);
	DrawCorrectionBlock(Trigg,"Z1/Trigger",finalresults);

	std::vector<EffCorr *> TriggFS;
	TriggFS.push_back(TriggerFullSpan_HE);
	DrawBlockCascade(TriggFS,"Z1/TriggerFS",finalresults,0.1,40);
	DrawCorrectionBlock(TriggFS,"Z1/TriggerFS",finalresults);

	std::vector<EffCorr *> Track;
	Track.push_back(TrackerEffCorr_HE);
	Track.push_back(GoodQTrack_HE);
	Track.push_back(GoodChi_HE);
	Track.push_back(KalmanEffCorr_HE);
	DrawBlockCascade(Track	   ,"Z1/Tracking",finalresults,0.1,40);
	DrawCorrectionBlock(Track  ,"Z1/Tracking",finalresults);

	std::vector<EffCorr *> L1PickUp;
	L1PickUp.push_back(L1PickUpGeom_HE);
	L1PickUp.push_back(L1PickUpEffCorr_HE);
	DrawBlockCascade(L1PickUp,"Z1/L1PickUp",finalresults,0.1,40);       
	DrawCorrectionBlock(L1PickUp,"Z1/L1PickUp",finalresults); 


	std::vector<EffCorr *> NoInteractions;
	NoInteractions.push_back(Good1Track_HE);
	NoInteractions.push_back(GoodUtof_HE);
	NoInteractions.push_back(GoodLtof_HE);
	DrawBlockCascade(NoInteractions,"Z1/NoInteractions",finalresults,0.1,40);	
	DrawCorrectionBlock(NoInteractions,"Z1/NoInteractions",finalresults);	

	std::vector<EffCorr *> GoodBetaTOF;
	GoodBetaTOF.push_back(GoodTime_TOF);
	GoodBetaTOF.push_back(Quality_TOF);
	DrawBlockCascade(GoodBetaTOF,"Z1/GoodBetaTOF",finalresults,0.1,40);	
	DrawCorrectionBlock(GoodBetaTOF,"Z1/GoodBetaTOF",finalresults);	


	std::vector<EffCorr *> GoodBetaNaF;
	GoodBetaNaF.push_back(RICHEffCorr_NaF);
	GoodBetaNaF.push_back(RICHQualEffCorr_NaF);
	DrawBlockCascade(GoodBetaNaF,"Z1/GoodBetaNaF",finalresults,0.1,40);	
	DrawCorrectionBlock(GoodBetaNaF,"Z1/GoodBetaNaF",finalresults);	


	std::vector<EffCorr *> GoodBetaAgl;
	GoodBetaAgl.push_back(RICHEffCorr_Agl);
	GoodBetaAgl.push_back(RICHQualEffCorr_Agl);
	DrawBlockCascade(GoodBetaAgl,"Z1/GoodBetaAgl",finalresults,0.1,40);	
	DrawCorrectionBlock(GoodBetaAgl,"Z1/GoodBetaAgl",finalresults);	


	std::vector<std::string > Names{"Trigger","L1 PickUp (Geom)","L1 PickUp (Association)","Tracking","NoInteractions","GoodVelocity"};

	std::vector<std::vector<EffCorr *>> TotalTOF;
	TotalTOF.push_back(Trigg);
	TotalTOF.push_back(L1PickUp);
	TotalTOF.push_back(Track);
	TotalTOF.push_back(NoInteractions);
	TotalTOF.push_back(GoodBetaTOF);
	DrawTotalUncertainty(TotalTOF,finalresults,"Z1_TOF",Names);

	std::vector<std::vector<EffCorr *>> TotalNaF;
	TotalNaF.push_back(Trigg);
	TotalNaF.push_back(L1PickUp);
	TotalNaF.push_back(Track);
	TotalNaF.push_back(NoInteractions);
	TotalNaF.push_back(GoodBetaNaF);
	DrawTotalUncertainty(TotalNaF,finalresults,"Z1_NaF",Names);

	std::vector<std::vector<EffCorr *>> TotalAgl;
	TotalAgl.push_back(Trigg);
	TotalAgl.push_back(L1PickUp);
	TotalAgl.push_back(Track);
	TotalAgl.push_back(NoInteractions);
	TotalAgl.push_back(GoodBetaAgl);
	DrawTotalUncertainty(TotalAgl,finalresults,"Z1_Agl",Names);

	//Z2
	std::vector<EffCorr *> TriggHe;
	TriggHe.push_back(TriggerEffCorr_Z2);
	DrawBlockCascade(TriggHe,"Z2/Trigger",finalresults,0.1,40);
	DrawCorrectionBlock(TriggHe,"Z2/Trigger",finalresults);

	std::vector<EffCorr *> TrackHe;
	TrackHe.push_back(TrackerEffCorr_Z2);
	TrackHe.push_back(GoodQTrack_Z2);
	TrackHe.push_back(GoodChi_Z2);
	TrackHe.push_back(KalmanEffCorr_Z2);
	DrawBlockCascade(TrackHe	   ,"Z2/Tracking",finalresults,0.1,40);
	DrawCorrectionBlock(TrackHe  ,"Z2/Tracking",finalresults);

	std::vector<EffCorr *> L1PickUpHe;
	L1PickUpHe.push_back(L1PickUpGeom_Z2);
	L1PickUpHe.push_back(L1PickUpEffCorr_Z2);
	DrawBlockCascade(L1PickUpHe,"Z2/L1PickUp",finalresults,0.1,40);       
	DrawCorrectionBlock(L1PickUpHe,"Z2/L1PickUp",finalresults); 


	std::vector<EffCorr *> NoInteractionsHe;
	NoInteractionsHe.push_back(Good1Track_Z2);
	NoInteractionsHe.push_back(GoodUtof_Z2);
	NoInteractionsHe.push_back(GoodLtof_Z2);
	DrawBlockCascade(NoInteractionsHe,"Z2/NoInteractions",finalresults,0.1,40);	
	DrawCorrectionBlock(NoInteractionsHe,"Z2/NoInteractions",finalresults);	

	std::vector<EffCorr *> GoodBetaTOFHe;
	GoodBetaTOFHe.push_back(GoodTimeZ2_TOF);
	DrawBlockCascade(GoodBetaTOFHe,"Z2/GoodBetaTOF",finalresults,0.1,40);	
	DrawCorrectionBlock(GoodBetaTOFHe,"Z2/GoodBetaTOF",finalresults);	


	std::vector<EffCorr *> GoodBetaNaFHe;
	GoodBetaNaFHe.push_back(RICHEffCorrZ2_NaF);
	GoodBetaNaFHe.push_back(RICHQualEffCorrZ2_NaF);
	DrawBlockCascade(GoodBetaNaFHe,"Z2/GoodBetaNaF",finalresults,0.1,40);	
	DrawCorrectionBlock(GoodBetaNaFHe,"Z2/GoodBetaNaF",finalresults);	


	std::vector<EffCorr *> GoodBetaAglHe;
	GoodBetaAglHe.push_back(RICHEffCorrZ2_Agl);
	GoodBetaAglHe.push_back(RICHQualEffCorrZ2_Agl);
	DrawBlockCascade(GoodBetaAglHe,"Z2/GoodBetaAgl",finalresults,0.1,40);	
	DrawCorrectionBlock(GoodBetaAglHe,"Z2/GoodBetaAgl",finalresults);	


	std::vector<std::string > NamesHe{"Trigger","L1 PickUp (Geom)","L1 PickUp (Association)","Tracking","NoInteractions","GoodVelocity"};

	std::vector<std::vector<EffCorr *>> TotalTOFHe;
	TotalTOFHe.push_back(TriggHe);
	TotalTOFHe.push_back(L1PickUpHe);
	TotalTOFHe.push_back(TrackHe);
	TotalTOFHe.push_back(NoInteractionsHe);
	TotalTOFHe.push_back(GoodBetaTOFHe);
	DrawTotalUncertainty(TotalTOFHe,finalresults,"Z2_TOF",NamesHe);

	std::vector<std::vector<EffCorr *>> TotalNaFHe;
	TotalNaFHe.push_back(TriggHe);
	TotalNaFHe.push_back(L1PickUpHe);
		TotalNaFHe.push_back(TrackHe);
		TotalNaFHe.push_back(NoInteractionsHe);
		TotalNaFHe.push_back(GoodBetaNaFHe);
		DrawTotalUncertainty(TotalNaFHe,finalresults,"Z2_NaF",NamesHe);

		std::vector<std::vector<EffCorr *>> TotalAglHe;
		TotalAglHe.push_back(TriggHe);
		TotalAglHe.push_back(L1PickUpHe);
		TotalAglHe.push_back(TrackHe);
		TotalAglHe.push_back(NoInteractionsHe);
		TotalAglHe.push_back(GoodBetaAglHe);
		DrawTotalUncertainty(TotalAglHe,finalresults,"Z2_Agl",NamesHe);






	}

	return ;
}


void AnalyzeEffCorr(EffCorr * Correction, FileSaver  finalhistos, FileSaver  finalresults, TFile * systfile, std::string systdirname, std::string systerrname, std::string avgname, std::string tima,float shift, bool IsTrig){

	TH1F * syst=0x0;
        TH2F * avg=0x0;
 
/*	if(systfile) {
		syst = (TH1F*)systfile->Get((systdirname+systerrname).c_str());
	//	avg  = (TH2F*)systfile->Get((systdirname+avgname).c_str());
		cout<<"AVERAGE EFFCORR: "<<systdirname<<" "<<syst<<" "<<avg<<endl;
	}
*/
	if(IsTrig) Correction   -> SetAsTrigEffCorr();
	Correction   -> Set_SystStat(syst,avg,tima.c_str());
	Correction   -> Eval_Efficiencies();
	Correction   -> Eval_Corrections(shift);
	Correction   -> SaveResults(finalresults);
}	

	
void DrawBlockCascade(std::vector<EffCorr *> block, std::string Blockname,FileSaver Plots,float rangemin, float rangemax){

	TH1F * Efficiencies[block.size()];
	TH1F * EfficienciesMC[block.size()];



	for(int i=0;i<block.size();i++){
		Efficiencies[i] = block[i]->GetGlobEfficiency();
		EfficienciesMC[i] = block[i]->GetMCEfficiency();
		if(i!=0) Efficiencies[i]->Multiply(Efficiencies[i-1]);
		if(i!=0) EfficienciesMC[i]->Multiply(EfficienciesMC[i-1]);
	}	
	
	TCanvas * c1 = new TCanvas("Cascade");
	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetLogx();
	


	TH1F * Efficiencies_plot[block.size()];
	TH1F * EfficienciesMC_plot[block.size()];
	TGraphErrors * Efficiencies_Plot[block.size()];
	TGraphErrors * EfficienciesMC_Plot[block.size()];
	for(int i=0;i<block.size();i++){
		Efficiencies_plot[i] =ConvertBinnedHisto(Efficiencies[i],Efficiencies[i]->GetName(),block[0]->GetBins(),block[0]->IsEkin());
		EfficienciesMC_plot[i] =ConvertBinnedHisto(EfficienciesMC[i],Efficiencies[i]->GetName(),block[0]->GetBins(),block[0]->IsEkin());

		Efficiencies_Plot[i] = new TGraphErrors(Efficiencies_plot[i]);
		EfficienciesMC_Plot[i] = new TGraphErrors(EfficienciesMC_plot[i]);

		Efficiencies_Plot[i]->SetMarkerColor(i+2);
		Efficiencies_Plot[i]->SetMarkerSize(1.2);
		Efficiencies_Plot[i]->SetMarkerStyle(8);
		EfficienciesMC_Plot[i]->SetMarkerStyle(0);
		Efficiencies_Plot[i]->SetLineWidth(0);
	
		Efficiencies_Plot[i]->GetYaxis()->SetRangeUser(0,1.1);
		if(i==0) Efficiencies_Plot[i]->Draw("AP");
		else	 Efficiencies_Plot[i]->Draw("Psame");
		EfficienciesMC_Plot[i]->SetLineWidth(0);
		EfficienciesMC_Plot[i]->SetMarkerColor(i+2);
		EfficienciesMC_Plot[i]->SetMarkerSize(1.4);
		EfficienciesMC_Plot[i]->SetMarkerStyle(4);
		EfficienciesMC_Plot[i]->Draw("Psame");
	}


	Plots.Add(c1);
	Plots.writeObjsInFolder(("AcceptancePlots/Blocks/"+Blockname).c_str());	
}


void DrawCorrectionBlock(std::vector<EffCorr *> block, std::string Blockname,FileSaver Plots){
	TStyle* m_gStyle= new TStyle();
        m_gStyle->SetPalette(55);
	int nColors = m_gStyle->GetNumberOfColors();


	TGraphErrors * TotalCorrection = new TGraphErrors();
	float y=1;
	for(int bin=0;bin<block[0]->GetBins().size();bin++){
		y=1;
		for(int i=0;i<block.size();i++){
			if(block[i]->IsEkin()) y*=block[i]->GetCorrectionModel()->Eval(block[i]->GetBins().EkPerMassBinCent(bin));
			else y*=block[i]->GetCorrectionModel()->Eval(block[i]->GetBins().RigBinCent(bin));
		}
		if(block[0]->IsEkin()) TotalCorrection->SetPoint(bin,block[0]->GetBins().EkPerMassBinCent(bin),y);
		else TotalCorrection->SetPoint(bin,block[0]->GetBins().RigBinCent(bin),y); 	
		float toterr=0;
		for(int i=0;i<block.size();i++){
			toterr+=pow(block[i]->GetGlobCorrection()->GetBinError(bin+1),2);
			toterr+=pow(block[i]->GetSyst_Stat()->GetBinContent(bin+1),2);
		}
		TotalCorrection->SetPointError(bin,0,pow(toterr,0.5));
	}
	TotalCorrection->SetLineColor(1);
	TotalCorrection->SetLineWidth(2);	
	TotalCorrection->SetLineStyle(2);	
	TotalCorrection->SetFillColor(1);
	TotalCorrection->SetFillStyle(3002);


	TotalCorrection->SetTitle((Blockname+ " Corrections").c_str());
	TotalCorrection->GetYaxis()->SetRangeUser(0.8*y,1.1*block[0]->GetCorrectionModel()->Eval(5));
	TotalCorrection->GetYaxis()->SetTitle("Correction");
	if(block[0]->IsEkin())	TotalCorrection->GetXaxis()->SetTitle("Ekin [GeV/n]");
	else TotalCorrection->GetXaxis()->SetTitle("R [GV]");
	
	TCanvas * c1 = new TCanvas("Correction Block");
	if(block.size()>1) TotalCorrection->Draw("Ae4C");
	for(int i=0;i<block.size();i++) block[i]->GetCorrectionModel()->SetLineColor(i+2);	
	for(int i=0;i<block.size();i++) block[i]->GetCorrectionModel()->Draw("same");

	TGraphErrors * raw_corrections[block.size()];
	for(int i=0;i<block.size();i++){
		raw_corrections[i]=new TGraphErrors(block[i]->GetGlobCorrection());
		for(int j=0;j<block[i]->GetGlobCorrection()->GetNbinsX();j++) {
			if(block[i]->GetGlobCorrection()->GetBinError(j+1)>3*block[i]->GetGlobCorrection()->GetBinError(j+2)
				||block[i]->GetGlobCorrection()->GetBinError(j+1)>0.1)
				raw_corrections[i]->RemovePoint(j);
			raw_corrections[i]->SetPointError(j,0,raw_corrections[i]->GetErrorY(j));
		}
		raw_corrections[i]->SetMarkerColor(i+2);
		raw_corrections[i]->SetMarkerStyle(8);
		raw_corrections[i]->SetLineColor(i+2);
		raw_corrections[i]->SetLineWidth(2);
	if(block.size()==1)	raw_corrections[i]->Draw("APsame");
	else	raw_corrections[i]->Draw("Psame");
	}

	TGraphErrors * correction_err[block.size()];
	for(int i=0;i<block.size();i++){
		correction_err[i]=new TGraphErrors();
	
		for(int j=0;j<block[i]->GetGlobCorrection()->GetNbinsX();j++) {
			if(block[i]->IsEkin()) correction_err[i]->SetPoint(j,block[i]->GetBins().EkPerMassBinCent(j),block[i]->GetCorrectionModel()->Eval(block[i]->GetBins().EkPerMassBinCent(j)));
			else  correction_err[i]->SetPoint(j,block[i]->GetBins().RigBinCent(j),block[i]->GetCorrectionModel()->Eval(block[i]->GetBins().RigBinCent(j)));
			TH1F* state=block[i]->GetStat_Err();
			TH1F* syste=block[i]->GetSyst_Stat();
			correction_err[i]->SetPointError(j,0,pow(pow(state->GetBinContent(j+1),2)+pow(syste->GetBinContent(j+1),2),0.5));
		}	
		correction_err[i]->SetFillColor(i+2);
		correction_err[i]->SetFillStyle(3002);
		correction_err[i]->SetLineColor(0);
		correction_err[i]->SetLineWidth(0);
		correction_err[i]->Draw("e4same");
	}


	Plots.Add(c1);
	Plots.writeObjsInFolder(("AcceptancePlots/Blocks/"+Blockname).c_str());	


//	TH1F * Total_StatError = (TH1F *) block[0]->GetStat_Err()->Clone("Total_StatError");
//	TH1F * Total_SystError = (TH1F *) block[0]->GetSyst_Err()->Clone("Total_SystError");

	
	TH1F * Total_StatError = CreateHisto("Total_StatError",block[0]->GetBins(),block[0]->IsEkin()); 
	TH1F * Total_SystError = CreateHisto("Total_SystError",block[0]->GetBins(),block[0]->IsEkin());
	TH1F * Total_SystStat  = CreateHisto("Total_SystStat",block[0]->GetBins(),block[0]->IsEkin());



	for(int bin=0;bin<block[0]->GetBins().size();bin++){
		float errstat=0;
		float errsyst=0;
		float errsyststat=0;
		
		for(int i=0;i<block.size();i++){
		errstat+=pow(block[i]->GetStat_Err()->GetBinContent(bin+1),2);			
			errsyst+=pow(block[i]->GetSyst_Err()->GetBinContent(bin+1),2);			
			errsyststat+=pow(block[i]->GetSyst_Stat()->GetBinContent(bin+1),2);			
		}	
		cout<<"systatat: "<<errsyststat<<endl;
		Total_StatError->SetBinContent(bin+1,pow(errstat,0.5));
		Total_SystError->SetBinContent(bin+1,pow(errsyst,0.5));
		Total_SystStat->SetBinContent(bin+1,pow(errsyststat,0.5));
	}

	TCanvas * c2 = new TCanvas("Block Error");
	
	Total_SystError->SetLineColor(2);
        Total_StatError->SetLineColor(1);
	Total_SystError->SetLineWidth(2);
        Total_StatError->SetLineWidth(1);
	Total_SystStat->SetLineColor(4);
        Total_SystStat->SetLineWidth(2);
	Total_SystError->Draw("hist");
	Total_StatError->Draw("hist,same");	
	Total_SystStat->Draw("hist,same");	


	Plots.Add(c2);
	Plots.writeObjsInFolder(("AcceptancePlots/Blocks/"+Blockname).c_str());	


}

void DrawTotalUncertainty(std::vector<std::vector<EffCorr *>> Total,FileSaver Plots, std::string name, std::vector<std::string> Names){
	TH1F * Total_Uncertainty[Total.size()];

	for(int i =0;i<Total.size();i++) 
		Total_Uncertainty[i]=CreateHisto(Names[i].c_str(),Total[0][0]->GetBins(),Total[0][0]->IsEkin());

			
	for(int bin=0;bin<Total[0][0]->GetBins().size();bin++){
		for(int i =0;i<Total.size();i++){
				float errcorrection=0;
				for(int j=0;j<Total[i].size();j++){
					errcorrection+=pow(Total[i][j]->GetStat_Err()->GetBinContent(bin+1),2);
					errcorrection+=pow(Total[i][j]->GetSyst_Err()->GetBinContent(bin+1),2);
					errcorrection+=pow(Total[i][j]->GetSyst_Stat()->GetBinContent(bin+1),2);
				}
				Total_Uncertainty[i]->SetBinContent(bin+1,pow(errcorrection,0.5));
			}			
	  }

	TH1F * TOTAL = CreateHisto("Total Uncertainty",Total[0][0]->GetBins(),Total[0][0]->IsEkin());
	for(int bin=0;bin<Total_Uncertainty[0]->GetNbinsX();bin++){
		float content=0;
		for(int i =0;i<Total.size();i++){
			content += pow(Total_Uncertainty[i]->GetBinContent(bin+1),2);
			}
		TOTAL->SetBinContent(bin+1,pow(content,0.5) );
	}


	TCanvas * c2 = new TCanvas(("Global Error "+name).c_str());
	TOTAL->SetLineColor(1);
	TOTAL->SetLineWidth(3);
	TOTAL->Draw("hist");
	
	for(int i =0;i<Total.size();i++){

		Total_Uncertainty[i]->SetLineColor(i+2);
		Total_Uncertainty[i]->SetLineWidth(2);
		Total_Uncertainty[i]->Draw("hist,same");
	}		
	Plots.Add(c2);
        Plots.writeObjsInFolder("Global Acceptance Errors");
}




