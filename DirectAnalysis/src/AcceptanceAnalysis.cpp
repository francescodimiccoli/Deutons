#include "Analyzer.h"
#include "../include/filesaver.h"
#include "../include/Efficiency.h"
#include "../include/EffCorr.h"
#include "../include/EffCorrTemplate.h"
#include "Acceptance.h"


void AnalyzeEffCorr(EffCorr * Correction, FileSaver  finalhistos, FileSaver  finalresults,bool IsTrig=false);

void Analyzer::BookAcceptanceAnalysis(FileSaver finalhistos, FileSaver finalresults, bool refill)
{

	bool checkfile = finalhistos.CheckFile();
        check_file = checkfile;
	cout<<"****************************** BINS ***************************************"<<endl;
        SetUpUsualBinning();
 
	cout<<"****************************** Efficiency corr. ANALYIS ******************************************"<<endl;

	BadEventSimulator * NaFBadEvSimulator= new BadEventSimulator("IsFromNaF",30,0.82,1);
        BadEventSimulator * AglBadEvSimulator= new BadEventSimulator("IsFromAgl",120,0.96,1);


	//MC Efficiencies
	Efficiency * Cascade0 = new Efficiency(finalhistos,"Cascade1","Cascade1",PRB,"IsProtonMC","IsProtonMC&IsPositive&IsPhysTrig",true);
	Efficiency * Cascade1 = new Efficiency(finalhistos,"Cascade1","Cascade1",PRB,"IsProtonMC","IsProtonMC&IsPositive&IsBaseline",true);
	Efficiency * Cascade2 = new Efficiency(finalhistos,"Cascade2","Cascade2",PRB,"IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1",true);
	Efficiency * Cascade3 = new Efficiency(finalhistos,"Cascade3","Cascade3",PRB,"IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",true);
	Efficiency * Cascade4 = new Efficiency(finalhistos,"Cascade4","Cascade4",PRB,"IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime",true);
	Efficiency * Cascade5 = new Efficiency(finalhistos,"Cascade5","Cascade5",PRB,"IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF",true);
	Efficiency * Cascade6 = new Efficiency(finalhistos,"Cascade6","Cascade6",PRB,"IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl",true);
	Efficiency * Cascade7 = new Efficiency(finalhistos,"Cascade7","Cascade7",PRB,"IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut",true);
	Efficiency * Cascade8 = new Efficiency(finalhistos,"Cascade8","Cascade8",PRB,"IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut",true);



	//Baseline efficiency corrections
	std::string before;
        std::string after;
        before = "IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge1Track&L1LooseCharge1&IsLUT2";
        after  = "IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge1Track&L1LooseCharge1&IsPhysTrig";
	EffCorr * TriggerEffCorr_HE  = new EffCorr(finalhistos,"TriggerEffCorr_HE" ,"Trigger Eff. Corr",true  ,before,after,"IsPrimary","IsProtonMC","IsPureDMC","IsPurePMC"); 

        before = "IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge1Track&L1LooseCharge1&IsLUT2&HasL9";
        after  = "IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge1Track&L1LooseCharge1&IsPhysTrigi&HasL9";
	//EffCorr * TriggerFullSpan_HE  = new EffCorr(finalhistos,"TriggerFullSpan_HE" ,"Trigger Full Span",true  ,before,after,"IsPrimary","IsProtonMC","IsPureDMC","IsPurePMC"); 


	before = "IsDownGoing&IsPhysTrig&IsL1HitNearExtrapol&IsCleanL1Hit&IsGoodTOFStandaloneQ1&IsGoodTrackPattern";
        after  = "IsDownGoing&IsPhysTrig&IsL1HitNearExtrapol&IsCleanL1Hit&IsGoodTOFStandaloneQ1&IsGoodTrackPattern&HasL1";
	EffCorr * L1PickUpEffCorr_HE = new EffCorr(finalhistos,"L1PickUpEffCorr_HE","L1PickUp Eff. Corr",false,before,after,"IsPrimary",    "IsProtonMC","IsPureDMC","IsPurePMC"); 

	before = "IsDownGoing&IsPhysTrig&IsGoodTrack&L1LooseCharge1&IsCharge1UTOF";
        after  = "IsDownGoing&IsPhysTrig&IsGoodTrack&L1LooseCharge1&IsCharge1UTOF&IsCharge1Track"; 
	EffCorr * GoodQTrack_HE  = new EffCorr(finalhistos,"GoodQTrackEffCorr_HE","GoodQTrack Eff. Corr",true,before,after,"IsPrimary","IsProtonMC","IsPureDMC","IsPurePMC");

	before = "IsDownGoing&IsPhysTrig&IsGoodTrack&IsCharge1Track&L1LooseCharge1";
        after  = "IsDownGoing&IsPhysTrig&IsGoodTrack&IsCharge1Track&L1LooseCharge1&IsGoodChi2"; 
	EffCorr * GoodChi_HE  = new EffCorr(finalhistos,"GoodChiEffCorr_HE","GoodChi Eff. Corr",false,before,after,"IsPrimary","IsProtonMC","IsPureDMC","IsPurePMC");

//        before = "IsDownGoing&IsPhysTrig&IsCompact&IsGoodTOFStandaloneQ1";
//        after  = "IsDownGoing&IsPhysTrig&IsCompact&IsGoodTOFStandaloneQ1&IsGoodTrack"; 
	before = "IsPhysTrig&IsCompact";
        after  = "IsPhysTrig&IsCompact&IsCompact_An"; 
	EffCorr * TrackerEffCorr_HE = new EffCorr(finalhistos,"TrackerEffCorr_HE","Tracker Eff. Corr",false,before,after,"IsPrimary","IsProtonMC","IsPureDMC","IsPurePMC");

	before = "IsDownGoing&IsPhysTrig&IsGoodTrack&IsCharge1Track&L1LooseCharge1";
        after  = "IsDownGoing&IsPhysTrig&IsGoodTrack&IsCharge1Track&L1LooseCharge1&IsGoodL1Status"; 
	EffCorr * StatusL1Check_HE = new EffCorr(finalhistos,"StatusL1Check_HE","StatusL1Check Eff. Corr",true,before,after,"IsPrimary","IsProtonMC","IsPureDMC","IsPurePMC");


	// Good Z=1 and Golden Efficiency corrections
	before = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1";
        after  = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1&Is1TrTrack"; 
	EffCorr * Good1Track_HE  = new EffCorr(finalhistos,"Good1TrackEffCorr_HE","Good1Track Eff. Corr",true,before,after,"IsPrimary","IsProtonMC","IsPureDMC","IsPurePMC");

	before = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1&Is1TrTrack";
	after  = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1&Is1TrTrack&IsCharge1LTOF";
        EffCorr * GoodLtof_HE  = new EffCorr(finalhistos,"GoodLTOFEffCorr_HE" ,"GoodLtof Eff. Corr",true,before,after,"IsPrimary","IsProtonMC","IsPureDMC","IsPurePMC");

	before = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1&Is1TrTrack&IsCharge1LTOF";
        after  = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1&Is1TrTrack&IsCharge1LTOF&IsCharge1UTOF"; 
	EffCorr * GoodUtof_HE  = new EffCorr(finalhistos,"GoodUtofEffCorr_HE" ,"GoodUtof Eff. Corr",true,before,after,"IsPrimary","IsProtonMC","IsPureDMC","IsPurePMC");

	before = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1&IsCleaning";
        after  = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime"; 
	EffCorr * GoodTime_TOF = new EffCorr(finalhistos,"GoodTimeEffCorr_TOF","GoodTime Eff. Corr",true,before,after,"IsPrimary","IsProtonMC","IsPureDMC","IsPurePMC");
	
	before = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1&IsCleaning";
        after  = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1&IsCleaning";
	EffCorr * RICHEffCorr_NaF = new EffCorr(finalhistos,"RICHCorrection_NaF","RICH Eff. Corr",true,before,(after+"&IsFromNaF").c_str(),"IsPrimary","IsProtonMC","IsPureDMC","IsPurePMC");
	EffCorr * RICHEffCorr_Agl = new EffCorr(finalhistos,"RICHCorrection_Agl","RICH Eff. Corr",true,before,(after+"&IsFromAgl").c_str(),"IsPrimary","IsProtonMC","IsPureDMC","IsPurePMC");

	before = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1&IsCleaning";
        after  = "IsPositive&IsPhysTrig&IsBaseline&L1LooseCharge1&IsCleaning";
	EffCorr * RICHQualEffCorr_NaF = new EffCorr(finalhistos,"RICHQualCorrection_NaF","RICH Qual Eff. Corr",true,(before+"&IsFromNaF").c_str(),(after+"&IsFromNaF&RICHBDTCut").c_str(),"IsPrimary","IsProtonMC","IsPureDMC","IsPurePMC");
	EffCorr * RICHQualEffCorr_Agl = new EffCorr(finalhistos,"RICHqualCorrection_Agl","RICH Qual. Eff. Corr",true,(before+"&IsFromAgl").c_str(),(after+"&IsFromAgl&RICHBDTCut").c_str(),"IsPrimary","IsProtonMC","IsPureDMC","IsPurePMC");

	//Acceptance
	Acceptance * Acceptance_HE     = new Acceptance(finalhistos,"Acceptance_HE"	,"Acceptance","IsProtonMC","IsProtonMC&IsPositive&IsBaseline&RigSafetyCut"	    ,PRB);
	Acceptance * Acceptance_L1HE   = new Acceptance(finalhistos,"Acceptance_L1HE"	,"Acceptance","IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&RigSafetyCut"	    ,PRB);
	Acceptance * Acceptance_QualHE = new Acceptance(finalhistos,"Acceptance_QualHE" ,"Acceptance","IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&RigSafetyCut",PRB);

	Acceptance * Acceptance_RigPTOF = new Acceptance(finalhistos,"Acceptance_RigPTOF","Acceptance","IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&RigSafetyCut",GlobalRig.GetToFPBins());
	Acceptance * Acceptance_RigPNaF = new Acceptance(finalhistos,"Acceptance_RigPNaF","Acceptance","IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut",GlobalRig.GetNaFPBins());
	Acceptance * Acceptance_RigPAgl = new Acceptance(finalhistos,"Acceptance_RigPAgl","Acceptance","IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut",GlobalRig.GetAglPBins());

	Acceptance * Acceptance_PTOF = new Acceptance(finalhistos,"Acceptance_PTOF","Acceptance","IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime&RigSafetyCut",Global.GetToFPBins());
	Acceptance * Acceptance_PNaF = new Acceptance(finalhistos,"Acceptance_PNaF","Acceptance","IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut",Global.GetNaFPBins());
	Acceptance * Acceptance_PAgl = new Acceptance(finalhistos,"Acceptance_PAgl","Acceptance","IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut",Global.GetAglPBins());

	Acceptance * Acceptance_DTOF = new Acceptance(finalhistos,"Acceptance_DTOF","Acceptance","IsDeutonMC","IsPureDMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime",Global.GetToFDBins());
	Acceptance * Acceptance_DNaF = new Acceptance(finalhistos,"Acceptance_DNaF","Acceptance","IsDeutonMC","IsPureDMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut",Global.GetNaFDBins());
	Acceptance * Acceptance_DAgl = new Acceptance(finalhistos,"Acceptance_DAgl","Acceptance","IsDeutonMC","IsPureDMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut",Global.GetAglDBins());


	
	Cascade0         ->SetDefaultOutFile(finalhistos);
	Cascade1  	 ->SetDefaultOutFile(finalhistos);
	Cascade2  	 ->SetDefaultOutFile(finalhistos);
	Cascade3 	 ->SetDefaultOutFile(finalhistos);
	Cascade4 	 ->SetDefaultOutFile(finalhistos);
	Cascade5 	 ->SetDefaultOutFile(finalhistos);
	Cascade6 	 ->SetDefaultOutFile(finalhistos);
	Cascade7 	 ->SetDefaultOutFile(finalhistos);
	Cascade8 	 ->SetDefaultOutFile(finalhistos);
	TriggerEffCorr_HE	->SetDefaultOutFile(finalhistos); 
	//TriggerFullSpan_HE	->SetDefaultOutFile(finalhistos); 
	L1PickUpEffCorr_HE ->SetDefaultOutFile(finalhistos); 
	TrackerEffCorr_HE ->SetDefaultOutFile(finalhistos); 
	StatusL1Check_HE ->SetDefaultOutFile(finalhistos); 
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


	Filler.AddObject2beFilled(Cascade0,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(Cascade1,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(Cascade2,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(Cascade3,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(Cascade4,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(Cascade5,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(Cascade6,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(Cascade7,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(Cascade8,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(TriggerEffCorr_HE,GetRigidity,GetRigidity);
//	Filler.AddObject2beFilled(TriggerFullSpan_HE,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(L1PickUpEffCorr_HE,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(TrackerEffCorr_HE,GetMomentumProxy,GetMomentumProxy);
	Filler.AddObject2beFilled(StatusL1Check_HE,GetMomentumProxy,GetMomentumProxy);
	Filler.AddObject2beFilled(GoodChi_HE   ,GetRigidity,GetRigidity);	
	Filler.AddObject2beFilled(GoodUtof_HE  ,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(GoodLtof_HE  ,GetRigidity,GetRigidity);	
	Filler.AddObject2beFilled(GoodQTrack_HE,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(Good1Track_HE,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(GoodTime_TOF,GetRigidity,GetRigidity);	
	Filler.AddObject2beFilled(RICHEffCorr_NaF,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(RICHEffCorr_Agl,GetRigidity,GetRigidity);	
	Filler.AddObject2beFilled(RICHQualEffCorr_NaF,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(RICHQualEffCorr_Agl,GetRigidity,GetRigidity);	

	Filler.AddObject2beFilled(Acceptance_PTOF,GetBetaTOF,GetBetaTOF);
        Filler.AddObject2beFilled(Acceptance_PNaF,GetBetaRICH,GetBetaRICH);
        Filler.AddObject2beFilled(Acceptance_PAgl,GetBetaRICH,GetBetaRICH);
	Filler.AddObject2beFilled(Acceptance_DTOF,GetBetaTOF,GetBetaTOF);
        Filler.AddObject2beFilled(Acceptance_DNaF,GetBetaRICH,GetBetaRICH);
        Filler.AddObject2beFilled(Acceptance_DAgl,GetBetaRICH,GetBetaRICH);
	Filler.AddObject2beFilled(Acceptance_RigPTOF 	,GetRigidity,GetRigidity);
        Filler.AddObject2beFilled(Acceptance_RigPNaF 	,GetRigidity,GetRigidity);
        Filler.AddObject2beFilled(Acceptance_RigPAgl 	,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(Acceptance_HE    	,GetRigidity,GetRigidity);
        Filler.AddObject2beFilled(Acceptance_L1HE  	,GetRigidity,GetRigidity);
        Filler.AddObject2beFilled(Acceptance_QualHE	,GetRigidity,GetRigidity);










	Filler.ReinitializeAll(refill);

	if(!refill&&checkfile) {	
	
		Acceptance_PTOF->Set_MCPar(0.5,100,2,"Pr.B1200/pr.pl1.05100.4_00.info");	
		Acceptance_PNaF->Set_MCPar(0.5,100,2,"Pr.B1200/pr.pl1.05100.4_00.info");	
		Acceptance_PAgl->Set_MCPar(0.5,100,2,"Pr.B1200/pr.pl1.05100.4_00.info");	

		Acceptance_DTOF->Set_MCPar(0.2,1000,2,"D.B1128/d.pl1ph.021000.info",0.02);	
		Acceptance_DNaF->Set_MCPar(0.2,1000,2,"D.B1128/d.pl1ph.021000.info",0.02);	
		Acceptance_DAgl->Set_MCPar(0.2,1000,2,"D.B1128/d.pl1ph.021000.info",0.021);	

		Acceptance_RigPTOF->Set_MCPar(0.5,100,2,"Pr.B1200/pr.pl1.05100.4_00.info");	
		Acceptance_RigPNaF->Set_MCPar(0.5,100,2,"Pr.B1200/pr.pl1.05100.4_00.info");	
		Acceptance_RigPAgl->Set_MCPar(0.5,100,2,"Pr.B1200/pr.pl1.05100.4_00.info");	
	
		Acceptance_HE    ->Set_MCPar(0.5,100,2,"Pr.B1200/pr.pl1.05100.4_00.info");	
		Acceptance_L1HE  ->Set_MCPar(0.5,100,2,"Pr.B1200/pr.pl1.05100.4_00.info");	
		Acceptance_QualHE->Set_MCPar(0.5,100,2,"Pr.B1200/pr.pl1.05100.4_00.info");	


		Cascade0  	->Eval_Efficiency();
		Cascade1  	->Eval_Efficiency();
		Cascade2  	->Eval_Efficiency();
		Cascade3 	->Eval_Efficiency();
		Cascade4 	->Eval_Efficiency();
		Cascade5 	->Eval_Efficiency();
		Cascade6 	->Eval_Efficiency();
		Cascade7 	->Eval_Efficiency();
		Cascade8 	->Eval_Efficiency();

		Cascade1        ->SaveResults(finalresults);
		Cascade1  	->SaveResults(finalresults);
		Cascade2  	->SaveResults(finalresults);
		Cascade3 	->SaveResults(finalresults);
		Cascade4 	->SaveResults(finalresults);
		Cascade5 	->SaveResults(finalresults);
		Cascade6 	->SaveResults(finalresults);
		Cascade7 	->SaveResults(finalresults);
		Cascade8 	->SaveResults(finalresults);
		
		AnalyzeEffCorr(	TriggerEffCorr_HE  , finalhistos, finalresults,true);
		//AnalyzeEffCorr(	TriggerFullSpan_HE  , finalhistos, finalresults,true);
		AnalyzeEffCorr(	L1PickUpEffCorr_HE , finalhistos, finalresults);
		AnalyzeEffCorr(	TrackerEffCorr_HE , finalhistos, finalresults);
		AnalyzeEffCorr(	StatusL1Check_HE , finalhistos, finalresults);
		AnalyzeEffCorr(	GoodChi_HE  , finalhistos, finalresults);  
		AnalyzeEffCorr(	GoodUtof_HE , finalhistos, finalresults);
		AnalyzeEffCorr(	GoodLtof_HE  , finalhistos, finalresults);
		AnalyzeEffCorr(	GoodQTrack_HE , finalhistos, finalresults);
		AnalyzeEffCorr(	Good1Track_HE , finalhistos, finalresults);
		AnalyzeEffCorr(	GoodTime_TOF , finalhistos, finalresults);
		AnalyzeEffCorr(	RICHEffCorr_NaF , finalhistos, finalresults);
		AnalyzeEffCorr(	RICHEffCorr_Agl , finalhistos, finalresults);
		AnalyzeEffCorr(	RICHQualEffCorr_NaF , finalhistos, finalresults);
		AnalyzeEffCorr(	RICHQualEffCorr_Agl , finalhistos, finalresults);

		//quality effcorr
		
		Acceptance_PTOF->ApplyEfficCorr(GoodUtof_HE);
		Acceptance_PNaF->ApplyEfficCorr(GoodUtof_HE);
		Acceptance_PAgl->ApplyEfficCorr(GoodUtof_HE);

		Acceptance_PTOF->ApplyEfficCorr(GoodLtof_HE);
		Acceptance_PNaF->ApplyEfficCorr(GoodLtof_HE);
		Acceptance_PAgl->ApplyEfficCorr(GoodLtof_HE);

		Acceptance_PTOF->ApplyEfficCorr(Good1Track_HE);
		Acceptance_PNaF->ApplyEfficCorr(Good1Track_HE);
		Acceptance_PAgl->ApplyEfficCorr(Good1Track_HE);

		Acceptance_QualHE->ApplyEfficCorr(GoodUtof_HE);
		Acceptance_QualHE->ApplyEfficCorr(GoodLtof_HE);
		Acceptance_QualHE->ApplyEfficCorr(Good1Track_HE);
		
		//velocity effcorr
		Acceptance_PTOF->ApplyEfficCorr(GoodTime_TOF);
		Acceptance_DTOF->ApplyEfficCorr(GoodTime_TOF);

	        Acceptance_PNaF->ApplyEfficCorr(RICHEffCorr_NaF);
                Acceptance_PAgl->ApplyEfficCorr(RICHEffCorr_Agl);
                Acceptance_DNaF->ApplyEfficCorr(RICHEffCorr_NaF);
	        Acceptance_DAgl->ApplyEfficCorr(RICHEffCorr_Agl);

	        Acceptance_PNaF->ApplyEfficCorr(RICHQualEffCorr_NaF);
                Acceptance_PAgl->ApplyEfficCorr(RICHQualEffCorr_Agl);
                Acceptance_DNaF->ApplyEfficCorr(RICHQualEffCorr_NaF);
	        Acceptance_DAgl->ApplyEfficCorr(RICHQualEffCorr_Agl);
		
		Acceptance_PTOF 	->  EvalEffAcc();
	        Acceptance_PNaF 	->  EvalEffAcc();
                Acceptance_PAgl 	->  EvalEffAcc();
		
		Acceptance_PTOF 	->SaveResults(finalresults);
                Acceptance_PNaF 	->SaveResults(finalresults);
                Acceptance_PAgl 	->SaveResults(finalresults);

		Acceptance_DTOF 	->  EvalEffAcc();
	        Acceptance_DNaF 	->  EvalEffAcc();
                Acceptance_DAgl 	->  EvalEffAcc();
		
		Acceptance_DTOF 	->SaveResults(finalresults);
                Acceptance_DNaF 	->SaveResults(finalresults);
                Acceptance_DAgl 	->SaveResults(finalresults);

		Acceptance_RigPTOF 	->  EvalEffAcc();
	        Acceptance_RigPNaF 	->  EvalEffAcc();
                Acceptance_RigPAgl 	->  EvalEffAcc();
		
		Acceptance_RigPTOF 	->SaveResults(finalresults);
                Acceptance_RigPNaF 	->SaveResults(finalresults);
                Acceptance_RigPAgl 	->SaveResults(finalresults);

		Acceptance_HE     	->  EvalEffAcc();
	        Acceptance_L1HE   	->  EvalEffAcc();
                Acceptance_QualHE 	->  EvalEffAcc();
		
		Acceptance_HE     	->SaveResults(finalresults);
                Acceptance_L1HE   	->SaveResults(finalresults);
                Acceptance_QualHE 	->SaveResults(finalresults);

	}

	return ;
}


void AnalyzeEffCorr(EffCorr * Correction, FileSaver  finalhistos, FileSaver  finalresults, bool IsTrig){
	if(IsTrig) Correction   -> SetAsTrigEffCorr();
	Correction   -> Eval_Efficiencies();
	Correction   -> Eval_Corrections();
	Correction   -> SaveResults(finalresults);
}	


