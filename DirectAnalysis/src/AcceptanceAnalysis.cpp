#include "Analyzer.h"
#include "../include/filesaver.h"
#include "../include/Efficiency.h"
#include "../include/EffCorr.h"
#include "../include/EffCorrTemplate.h"
#include "Acceptance.h"


void AnalyzeEffCorr(EffCorr * Correction, FileSaver  finalhistos, FileSaver  finalresults, TFile * systfile, std::string systname, float shift =0, bool IsTrig=false);
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
	before = "IsDownGoing&IsGoodTOFStandaloneQ1&IsGoodTRDStandaloneQ1&HasECAL";
        after  = "IsDownGoing&IsGoodTOFStandaloneQ1&IsGoodTRDStandaloneQ1&HasECAL&IsGoodTrack"; 
	EffCorr * TrackerEffCorr_HE = new EffCorr(finalhistos,"TrackerEffCorr_HE","Tracker Eff. Corr",false,before,after,"",true);

	before = "IsDownGoing&IsGoodTOFStandaloneQ1&IsGoodTRDStandaloneQ1&HasECAL&IsGoodTrack";
        after  = "IsDownGoing&IsGoodTOFStandaloneQ1&IsGoodTRDStandaloneQ1&HasECAL&IsGoodTrack&IsCharge1Track"; 
	EffCorr * GoodQTrack_HE  = new EffCorr(finalhistos,"GoodQTrackEffCorr_HE","GoodQTrack Eff. Corr",false,before,after,"",true);

	before = "IsDownGoing&IsGoodTOFStandaloneQ1&IsGoodTRDStandaloneQ1&HasECAL&IsGoodTrack&IsCharge1Track";
        after  = "IsDownGoing&IsGoodTOFStandaloneQ1&IsGoodTRDStandaloneQ1&HasECAL&IsGoodTrack&IsCharge1Track&IsGoodChi2"; 
	EffCorr * GoodChi_HE  = new EffCorr(finalhistos,"GoodChiEffCorr_HE","GoodChi Eff. Corr",false,before,after,"",true);

	before = "IsDownGoing&IsGoodTOFStandaloneQ1&IsGoodTRDStandaloneQ1&HasECAL&IsGoodTrack&IsCharge1Track&IsGoodChi2";
        after  = "IsDownGoing&IsGoodTOFStandaloneQ1&IsGoodTRDStandaloneQ1&HasECAL&IsGoodTrack&IsCharge1Track&IsGoodChi2&IsGoodKalman"; 
	EffCorr * KalmanEffCorr_HE = new EffCorr(finalhistos,"KalmanEffCorr_HE","Kalman Eff. Corr",false,before,after,"",true);

	// L1 Pickup efficiency corrections
	before = "IsDownGoing&IsPhysTrig&IsGoodTOFStandaloneQ1&IsGoodTRDStandaloneQ1&IsGoodTrack&IsGoodChi2&IsCharge1Track&Is1TrTrack&IsMinTOF&IsL1Fiducial&IsGoodTime";
        after  = "IsDownGoing&IsPhysTrig&IsGoodTOFStandaloneQ1&IsGoodTRDStandaloneQ1&IsGoodTrack&IsGoodChi2&IsCharge1Track&Is1TrTrack&IsMinTOF&IsL1Fiducial&IsGoodTime&IsL1HitNearExtrapol";
	EffCorr * L1PickUpGeom_HE = new EffCorr(finalhistos,"L1PickUpGeom_HE","L1PickUp Geom.",false,before,after,"IsPrimaryInner",true); 

	before = "IsDownGoing&IsPhysTrig&IsGoodTOFStandaloneQ1&IsGoodTrack&IsGoodChi2&IsCharge1Track&IsL1HitNearExtrapol&IsMinTOF&IsL1Fiducial&IsGoodTime";
        after  = "IsDownGoing&IsPhysTrig&IsGoodTOFStandaloneQ1&IsGoodTrack&IsGoodChi2&IsCharge1Track&IsL1HitNearExtrapol&IsMinTOF&IsL1Fiducial&IsGoodTime&L1LooseCharge1";
	EffCorr * L1PickUpEffCorr_HE = new EffCorr(finalhistos,"L1PickUpEffCorr_HE","L1PickUp Eff. Corr",false,before,after,"IsPrimaryInner",true); 

	// Trigger Eff
	before = "IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge1Track&IsGoodKalman&L1LooseCharge1&IsLUT2";
        after  = "IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge1Track&IsGoodKalman&L1LooseCharge1&IsPhysTrig";
	EffCorr * TriggerEffCorr_HE  = new EffCorr(finalhistos,"TriggerEffCorr_HE" ,"Trigger Eff. Corr",false  ,before,after,"IsPrimary"); 

        before = "IsDownGoin&IsGoodTrack&IsGoodChi2&IsCharge1Track&IsGoodKalman&L1LooseCharge1&HasL9&IsLUT2";
        after  = "IsDownGoin&IsGoodTrack&IsGoodChi2&IsCharge1Track&IsGoodKalman&L1LooseCharge1&HasL9&IsPhysTrig";
	EffCorr * TriggerFullSpan_HE  = new EffCorr(finalhistos,"TriggerFullSpan_HE" ,"Trigger Full Span",false  ,before,after,"IsPrimary"); 


	// Good Z=1 and Golden Efficiency corrections
	before = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1";
        after  = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1&Is1TrTrack"; 
	EffCorr * Good1Track_HE  = new EffCorr(finalhistos,"Good1TrackEffCorr_HE","Good1Track Eff. Corr",true,before,after,"");

	before = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1&Is1TrTrack";
	after  = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1&Is1TrTrack&IsCharge1LTOF";
        EffCorr * GoodLtof_HE  = new EffCorr(finalhistos,"GoodLTOFEffCorr_HE" ,"GoodLtof Eff. Corr",true,before,after,"");

	before = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1&Is1TrTrack&IsCharge1LTOF";
        after  = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1&Is1TrTrack&IsCharge1LTOF&IsCharge1UTOF"; 
	EffCorr * GoodUtof_HE  = new EffCorr(finalhistos,"GoodUtofEffCorr_HE" ,"GoodUtof Eff. Corr",true,before,after,"");

	before = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1&IsCleaning";
        after  = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime"; 
	EffCorr * GoodTime_TOF = new EffCorr(finalhistos,"GoodTimeEffCorr_TOF","GoodTime Eff. Corr",true,before,after,"");

	before = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime";
        after  = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&QualityTOF"; 
	EffCorr * Quality_TOF = new EffCorr(finalhistos,"QualityEffCorr_TOF","Quality TOF Eff. Corr",true,before,after,"");
	
	before = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1&IsCleaning";
        after  = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1&IsCleaning";
	EffCorr * RICHEffCorr_NaF = new EffCorr(finalhistos,"RICHCorrection_NaF","RICH Eff. Corr",true,before,(after+"&IsFromNaF").c_str(),"IsPrimary");
	EffCorr * RICHEffCorr_Agl = new EffCorr(finalhistos,"RICHCorrection_Agl","RICH Eff. Corr",true,before,(after+"&IsFromAgl").c_str(),"IsPrimary");

	before = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1&IsCleaning";
        after  = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1&IsCleaning";
	EffCorr * RICHQualEffCorr_NaF = new EffCorr(finalhistos,"RICHQualCorrection_NaF","RICH Qual. Eff. Corr",true,(before+"&IsFromNaF").c_str(),(after+"&IsFromNaF&RICHBDTCut").c_str(),"IsPrimary");
	EffCorr * RICHQualEffCorr_Agl = new EffCorr(finalhistos,"RICHqualCorrection_Agl","RICH Qual. Eff. Corr",true,(before+"&IsFromAgl").c_str(),(after+"&IsFromAgl&RICHBDTCut").c_str(),"IsPrimary");

	//Fragmentation
	Efficiency * Fragmentation_P = new Efficiency(finalhistos,"Fragmentation_P","Fragmentation_P",Quality_TOF->GetBins(),"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&QualityTOF","IsPurePMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&QualityTOF");
	Efficiency * Fragmentation_D = new Efficiency(finalhistos,"Fragmentation_D","Fragmentation_D",Quality_TOF->GetBins(),"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&QualityTOF","IsPureDMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&QualityTOF");


	//Acceptance

	Acceptance * Cascade0 = new Acceptance(finalhistos,"Cascade0"   ,"Acceptance","IsProtonMC","IsProtonMC&IsPositive&IsPhysTrig"                                                                          ,PRB,UnfoldingToF);
	Acceptance * Cascade1 = new Acceptance(finalhistos,"Cascade1"   ,"Acceptance","IsProtonMC","IsProtonMC&IsPositive&IsPhysTrig&IsDownGoing"                                                                ,PRB,UnfoldingToF);
	Acceptance * Cascade2 = new Acceptance(finalhistos,"Cascade2"   ,"Acceptance","IsProtonMC","IsProtonMC&IsPositive&IsPhysTrig&IsDownGoing&IsGoodTrack"                                                  ,PRB,UnfoldingToF);
	Acceptance * Cascade3 = new Acceptance(finalhistos,"Cascade3"   ,"Acceptance","IsProtonMC","IsProtonMC&IsPositive&IsPhysTrig&IsDownGoing&IsGoodTrack&IsGoodChi2"                                               ,PRB,UnfoldingToF);
	Acceptance * Cascade4 = new Acceptance(finalhistos,"Cascade4"   ,"Acceptance","IsProtonMC","IsProtonMC&IsPositive&IsPhysTrig&IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge1Track"                        ,PRB,UnfoldingToF);
	Acceptance * Cascade5 = new Acceptance(finalhistos,"Cascade5"   ,"Acceptance","IsProtonMC","IsProtonMC&IsPositive&IsPhysTrig&IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge1Track&IsGoodKalman"           ,PRB,UnfoldingToF);
	Acceptance * Cascade6 = new Acceptance(finalhistos,"Cascade6"   ,"Acceptance","IsProtonMC","IsProtonMC&IsPositive&IsPhysTrig&IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge1Track&IsGoodKalman&MassSafetyCut",PRB,UnfoldingToF);

	Acceptance * Acceptance_HE     = new Acceptance(finalhistos,"Acceptance_HE"	,"Acceptance","IsProtonMC","IsProtonMC&IsPositive&IsBaseline"	    ,PRB,UnfoldingToF);
	Acceptance * Acceptance_L1HE   = new Acceptance(finalhistos,"Acceptance_L1HE"	,"Acceptance","IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1"	    ,PRB,UnfoldingToF);
	Acceptance * Acceptance_QualHE = new Acceptance(finalhistos,"Acceptance_QualHE" ,"Acceptance","IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",PRB,UnfoldingToF);

	Acceptance * Acceptance_RigPTOF = new Acceptance(finalhistos,"Acceptance_RigPTOF","Acceptance","IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",GlobalRig.GetToFPBins(),UnfoldingToF);
	Acceptance * Acceptance_RigPNaF = new Acceptance(finalhistos,"Acceptance_RigPNaF","Acceptance","IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",GlobalRig.GetNaFPBins(),UnfoldingNaF);
	Acceptance * Acceptance_RigPAgl = new Acceptance(finalhistos,"Acceptance_RigPAgl","Acceptance","IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",GlobalRig.GetAglPBins(),UnfoldingAgl);

	Acceptance * Acceptance_PTOF = new Acceptance(finalhistos,"Acceptance_PTOF","Acceptance","IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&QualityTOF",Global.GetToFPBins(),UnfoldingToF);
	Acceptance * Acceptance_PNaF = new Acceptance(finalhistos,"Acceptance_PNaF","Acceptance","IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut",Global.GetNaFPBins(),UnfoldingNaF);
	Acceptance * Acceptance_PAgl = new Acceptance(finalhistos,"Acceptance_PAgl","Acceptance","IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut",Global.GetAglPBins(),UnfoldingAgl);

	Acceptance * Acceptance_DTOF = new Acceptance(finalhistos,"Acceptance_DTOF","Acceptance","IsDeutonMC","IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&QualityTOF",Global.GetToFDBins(),UnfoldingToF_D);
	Acceptance * Acceptance_DNaF = new Acceptance(finalhistos,"Acceptance_DNaF","Acceptance","IsDeutonMC","IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut",Global.GetNaFDBins(),UnfoldingNaF_D);
	Acceptance * Acceptance_DAgl = new Acceptance(finalhistos,"Acceptance_DAgl","Acceptance","IsDeutonMC","IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut",Global.GetAglDBins(),UnfoldingAgl_D);

	Cascade0         ->SetDefaultOutFile(finalhistos);
	Cascade1  	 ->SetDefaultOutFile(finalhistos);
	Cascade2  	 ->SetDefaultOutFile(finalhistos);
	Cascade3 	 ->SetDefaultOutFile(finalhistos);
	Cascade4 	 ->SetDefaultOutFile(finalhistos);
	Cascade5 	 ->SetDefaultOutFile(finalhistos);
	Cascade6 	 ->SetDefaultOutFile(finalhistos);
	TriggerEffCorr_HE	->SetDefaultOutFile(finalhistos); 
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
	Acceptance_PTOF 	->SetDefaultOutFile(finalhistos); 
 	Acceptance_PNaF 	->SetDefaultOutFile(finalhistos); 
 	Acceptance_PAgl 	->SetDefaultOutFile(finalhistos); 
	Acceptance_DTOF 	->SetDefaultOutFile(finalhistos); 
 	Acceptance_DNaF 	->SetDefaultOutFile(finalhistos); 
 	Acceptance_DAgl 	->SetDefaultOutFile(finalhistos); 
	Acceptance_RigPTOF 	->SetDefaultOutFile(finalhistos); 
 	Acceptance_RigPNaF 	->SetDefaultOutFile(finalhistos); 
 	Acceptance_RigPAgl 	->SetDefaultOutFile(finalhistos); 
	Acceptance_HE    	->SetDefaultOutFile(finalhistos); 
        Acceptance_L1HE  	->SetDefaultOutFile(finalhistos); 
	Acceptance_QualHE	->SetDefaultOutFile(finalhistos); 


	Filler.AddObject2beFilled(Cascade0,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(Cascade1,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(Cascade2,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(Cascade3,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(Cascade4,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(Cascade5,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(Cascade6,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(TriggerEffCorr_HE,GetRigidityInner,GetRigidityInner);
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

	Filler.AddObject2beFilled(Acceptance_PTOF,GetSmearedBetaTOF ,GetSmearedBetaTOF );
        Filler.AddObject2beFilled(Acceptance_PNaF,GetSmearedBetaRICH,GetSmearedBetaRICH);
        Filler.AddObject2beFilled(Acceptance_PAgl,GetSmearedBetaRICH,GetSmearedBetaRICH);
	Filler.AddObject2beFilled(Acceptance_DTOF,GetSmearedBetaTOF ,GetSmearedBetaTOF );
        Filler.AddObject2beFilled(Acceptance_DNaF,GetSmearedBetaRICH,GetSmearedBetaRICH);
        Filler.AddObject2beFilled(Acceptance_DAgl,GetSmearedBetaRICH,GetSmearedBetaRICH);
	Filler.AddObject2beFilled(Acceptance_RigPTOF, 	GetSmearedBetaTOF ,GetSmearedBetaTOF );
        Filler.AddObject2beFilled(Acceptance_RigPNaF,	GetSmearedBetaRICH,GetSmearedBetaRICH);
        Filler.AddObject2beFilled(Acceptance_RigPAgl, 	GetSmearedBetaRICH,GetSmearedBetaRICH);
	Filler.AddObject2beFilled(Acceptance_HE,    	GetSmearedBetaTOF ,GetSmearedBetaTOF );
	Filler.AddObject2beFilled(Acceptance_L1HE , 	GetSmearedBetaRICH,GetSmearedBetaRICH);
	Filler.AddObject2beFilled(Acceptance_QualHE,	GetSmearedBetaRICH,GetSmearedBetaRICH);
	Filler.AddObject2beFilled(Fragmentation_P,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(Fragmentation_D,GetGenMomentum,GetGenMomentum);

	Filler.ReinitializeAll(refill);


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
	


	if(!refill&&checkfile) {	


		Cascade0        ->Set_MCPar(0.5,100,1,"Pr.B1200/pr.pl1.05100.4_00.info",filelistMC.c_str());	
                Cascade1  	->Set_MCPar(0.5,100,1,"Pr.B1200/pr.pl1.05100.4_00.info",filelistMC.c_str());	
                Cascade2  	->Set_MCPar(0.5,100,1,"Pr.B1200/pr.pl1.05100.4_00.info",filelistMC.c_str());	
                Cascade3 	->Set_MCPar(0.5,100,1,"Pr.B1200/pr.pl1.05100.4_00.info",filelistMC.c_str());	
                Cascade4 	->Set_MCPar(0.5,100,1,"Pr.B1200/pr.pl1.05100.4_00.info",filelistMC.c_str());	
                Cascade5 	->Set_MCPar(0.5,100,1,"Pr.B1200/pr.pl1.05100.4_00.info",filelistMC.c_str());	
                Cascade6 	->Set_MCPar(0.5,100,1,"Pr.B1200/pr.pl1.05100.4_00.info",filelistMC.c_str());	

	
		Acceptance_PTOF->Set_MCPar(0.5,100,1,"Pr.B1200/pr.pl1.05100.4_00.info",filelistMC.c_str());	
		Acceptance_PNaF->Set_MCPar(0.5,100,1,"Pr.B1200/pr.pl1.05100.4_00.info",filelistMC.c_str());	
		Acceptance_PAgl->Set_MCPar(0.5,100,1,"Pr.B1200/pr.pl1.05100.4_00.info",filelistMC.c_str());	

		Acceptance_DTOF->Set_MCPar(0.5,100,1,"D.B1220/d.pl1.05100.info",filelistMC.c_str(),0.05);	
		Acceptance_DNaF->Set_MCPar(0.5,100,1,"D.B1220/d.pl1.05100.info",filelistMC.c_str(),0.05);	
		Acceptance_DAgl->Set_MCPar(0.5,100,1,"D.B1220/d.pl1.05100.info",filelistMC.c_str(),0.05);	

		Acceptance_RigPTOF->Set_MCPar(0.5,100,1,"Pr.B1200/pr.pl1.05100.4_00.info",filelistMC.c_str());	
		Acceptance_RigPNaF->Set_MCPar(0.5,100,1,"Pr.B1200/pr.pl1.05100.4_00.info",filelistMC.c_str());	
		Acceptance_RigPAgl->Set_MCPar(0.5,100,1,"Pr.B1200/pr.pl1.05100.4_00.info",filelistMC.c_str());	
	
		Acceptance_HE    ->Set_MCPar(0.5,100,1,"Pr.B1200/pr.pl1.05100.4_00.info",filelistMC.c_str());	
		Acceptance_L1HE  ->Set_MCPar(0.5,100,1,"Pr.B1200/pr.pl1.05100.4_00.info",filelistMC.c_str());	
		Acceptance_QualHE->Set_MCPar(0.5,100,1,"Pr.B1200/pr.pl1.05100.4_00.info",filelistMC.c_str());	

	
		Fragmentation_P->Eval_Efficiency();
	        Fragmentation_D->Eval_Efficiency();

		Fragmentation_P->SaveResults(finalresults);
	        Fragmentation_D->SaveResults(finalresults);

		TFile * EffSystFile = TFile::Open("/afs/cern.ch/work/f/fdimicco/private/Deutons/DirectAnalysis/EffSyst/Time.root");
		
	
		AnalyzeEffCorr(	TriggerEffCorr_HE  , finalhistos, finalresults,EffSystFile,"",1,true);
		AnalyzeEffCorr(	TriggerFullSpan_HE  , finalhistos, finalresults,EffSystFile,"",1,true);
		AnalyzeEffCorr(	L1PickUpEffCorr_HE , finalhistos, finalresults,EffSystFile,"",1);
		AnalyzeEffCorr(	L1PickUpGeom_HE , finalhistos, finalresults,EffSystFile,"",1);
		AnalyzeEffCorr(	TrackerEffCorr_HE , finalhistos, finalresults,EffSystFile,"",1);
		AnalyzeEffCorr(	KalmanEffCorr_HE , finalhistos, finalresults,EffSystFile,"",1);
		AnalyzeEffCorr(	GoodChi_HE  , finalhistos, finalresults,EffSystFile,"Eff. Corrections Sys/Track Quality/GoodChi_Err",0.8);  
		AnalyzeEffCorr(	GoodUtof_HE , finalhistos, finalresults,EffSystFile,"Eff. Corrections Sys/NoInteraction/GoodUtof_Err",0.2);
		AnalyzeEffCorr(	GoodLtof_HE  , finalhistos, finalresults ,EffSystFile,"Eff. Corrections Sys/NoInteraction/GoodLtof_Err",0.2); 
		AnalyzeEffCorr(	GoodQTrack_HE , finalhistos, finalresults,EffSystFile,"Eff. Corrections Sys/Track Quality/GoodQTrack_Err",0.8);
		AnalyzeEffCorr(	Good1Track_HE , finalhistos, finalresults,EffSystFile,"Eff. Corrections Sys/NoInteraction/Good1Track_Err",0.2); 
		AnalyzeEffCorr(	GoodTime_TOF , finalhistos, finalresults,EffSystFile, "Eff. Corrections Sys/TOF/GoodTime_Err",0.3);
		AnalyzeEffCorr(	Quality_TOF , finalhistos, finalresults,EffSystFile,    "Eff. Corrections Sys/TOF/QualityTime_Err",0.2);
		AnalyzeEffCorr(	RICHEffCorr_NaF , finalhistos, finalresults,EffSystFile,"Eff. Corrections Sys/RICH CIEMAT/RICHNTime_Err",1);
		AnalyzeEffCorr(	RICHEffCorr_Agl , finalhistos, finalresults,EffSystFile,"Eff. Corrections Sys/RICH CIEMAT/RICHATime_Err",3);
		AnalyzeEffCorr(	RICHQualEffCorr_NaF , finalhistos, finalresults,EffSystFile,"Eff. Corrections Sys/RICH BDT/RICHQNTime_Err",0.8);
		AnalyzeEffCorr(	RICHQualEffCorr_Agl , finalhistos, finalresults,EffSystFile,"Eff. Corrections Sys/RICH BDT/RICHQATime_Err",1);


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
		
		//Trigger effcorr
		//Acceptance_L1HE->ApplyEfficCorr(TriggerEffCorr_HE);
		//Acceptance_QualHE->ApplyEfficCorr(TriggerEffCorr_HE);
		
		//velocity effcorr
		Acceptance_PTOF->ApplyEfficCorr(GoodTime_TOF);
		Acceptance_PTOF->ApplyEfficCorr(Quality_TOF);
	        Acceptance_PNaF->ApplyEfficCorr(RICHEffCorr_NaF);
                Acceptance_PAgl->ApplyEfficCorr(RICHEffCorr_Agl);
	        Acceptance_PNaF->ApplyEfficCorr(RICHQualEffCorr_NaF);
                Acceptance_PAgl->ApplyEfficCorr(RICHQualEffCorr_Agl);
	
		Acceptance_DTOF->ApplyEfficCorr(GoodTime_TOF);
		Acceptance_DTOF->ApplyEfficCorr(Quality_TOF);
	        Acceptance_DNaF->ApplyEfficCorr(RICHEffCorr_NaF);
                Acceptance_DAgl->ApplyEfficCorr(RICHEffCorr_Agl);
	        Acceptance_DNaF->ApplyEfficCorr(RICHQualEffCorr_NaF);
                Acceptance_DAgl->ApplyEfficCorr(RICHQualEffCorr_Agl);
	
		Cascade0        ->EvalEffAcc(timeindex,1.2);
		Cascade1  	->EvalEffAcc(timeindex,1.2);
		Cascade2  	->EvalEffAcc(timeindex,1.2);
		Cascade3 	->EvalEffAcc(timeindex,1.2);
		Cascade4 	->EvalEffAcc(timeindex,1.2);
		Cascade5 	->EvalEffAcc(timeindex,1.2);
		Cascade6 	->EvalEffAcc(timeindex,1.2);

		Cascade0        ->SaveResults(finalresults);
		Cascade1  	->SaveResults(finalresults);
		Cascade2  	->SaveResults(finalresults);
		Cascade3 	->SaveResults(finalresults);
		Cascade4 	->SaveResults(finalresults);
		Cascade5 	->SaveResults(finalresults);
		Cascade6 	->SaveResults(finalresults);


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

		Acceptance_RigPTOF 	->  EvalEffAcc(timeindex,1.2);
	        Acceptance_RigPNaF 	->  EvalEffAcc(timeindex,1.2);
                Acceptance_RigPAgl 	->  EvalEffAcc(timeindex,1.2);
		
		Acceptance_RigPTOF 	->SaveResults(finalresults);
                Acceptance_RigPNaF 	->SaveResults(finalresults);
                Acceptance_RigPAgl 	->SaveResults(finalresults);

		Acceptance_HE     	->  EvalEffAcc(timeindex,1.2);
	        Acceptance_L1HE   	->  EvalEffAcc(timeindex,1.2);
                Acceptance_QualHE 	->  EvalEffAcc(timeindex,1.2);
		
		Acceptance_HE     	->SaveResults(finalresults);
                Acceptance_L1HE   	->SaveResults(finalresults);
                Acceptance_QualHE 	->SaveResults(finalresults);

	}

	return ;
}


void AnalyzeEffCorr(EffCorr * Correction, FileSaver  finalhistos, FileSaver  finalresults, TFile * systfile, std::string systname, float shift, bool IsTrig){

	TH1F * syst = (TH1F*)systfile->Get(systname.c_str());
	cout<<"SISTEMATIC ERROR: "<<systfile<<" "<<syst<<" "<<systname.c_str()<<endl;
	if(IsTrig) Correction   -> SetAsTrigEffCorr();
	Correction   -> Set_SystStat(syst);
	Correction   -> Eval_Efficiencies();
	Correction   -> Eval_Corrections(shift);
	Correction   -> SaveResults(finalresults);
}	

	


