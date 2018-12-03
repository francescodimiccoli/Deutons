#include "Analyzer.h"
#include "../include/filesaver.h"
#include "../include/Efficiency.h"
#include "../include/EffCorr.h"
#include "../include/EffCorrTemplate.h"

void AnalyzeEffCorr(EffCorr * Correction, FileSaver  finalhistos, FileSaver  finalresults);

void Analyzer::BookEffCorrAnalysis(FileSaver finalhistos, FileSaver finalresults, bool refill)
{

	bool checkfile = finalhistos.CheckFile();
        check_file = checkfile;
	cout<<"****************************** BINS ***************************************"<<endl;
        SetUpEffCorrBinning();
 
	cout<<"****************************** Efficiency corr. ANALYIS ******************************************"<<endl;

	BadEventSimulator * NaFBadEvSimulator= new BadEventSimulator("IsFromNaF",30,0.82,1);
        BadEventSimulator * AglBadEvSimulator= new BadEventSimulator("IsFromAgl",120,0.96,1);



	//Baseline efficiency corrections
	std::string before;
        std::string after;
        before = "IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge1Track&IsLooseCharge1";
        after  = "IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge1Track&IsLooseCharge1&IsPhysTrig";
	EffCorr * TriggerEffCorr_HE  = new EffCorr(finalhistos,"TriggerEffCorr_HE" ,"Trigger Eff. Corr",PRB  ,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC"); 
	EffCorr * TriggerEffCorr_TOF = new EffCorr(finalhistos,"TriggerEffCorr_TOF","Trigger Eff. Corr",ToFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * TriggerEffCorr_NaF = new EffCorr(finalhistos,"TriggerEffCorr_NaF","Trigger Eff. Corr",NaFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * TriggerEffCorr_Agl = new EffCorr(finalhistos,"TriggerEffCorr_Agl","Trigger Eff. Corr",AglPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");

	before = "IsDownGoing&IsPhysTrig&IsGoodTOFStandaloneQ1&IsL1HitNearExtrapol&IsCleanL1Hit";
        after  = "IsDownGoing&IsPhysTrig&IsGoodTOFStandaloneQ1&IsL1HitNearExtrapol&IsCleanL1Hit&L1LooseCharge1";
	EffCorr * L1PickUpEffCorr_HE = new EffCorr(finalhistos,"L1PickUpEffCorr_HE","L1PickUp Eff. Corr",PRB,before,after,"IsPrimary",    "IsPurePMC","IsPureDMC","IsDeutonMC"); 
	EffCorr * L1PickUpEffCorr_TOF = new EffCorr(finalhistos,"L1PickUpEffCorr_TOF","L1PickUp Eff. Corr",ToFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * L1PickUpEffCorr_NaF = new EffCorr(finalhistos,"L1PickUpEffCorr_NaF","L1PickUp Eff. Corr",NaFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * L1PickUpEffCorr_Agl = new EffCorr(finalhistos,"L1PickUpEffCorr_Agl","L1PickUp Eff. Corr",AglPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");

	before = "IsDownGoing&IsGoodTrack&IsLooseCharge1";
        after  = "IsDownGoing&IsGoodTrack&IsLooseCharge1&IsCharge1Track"; 
	EffCorr * GoodQTrack_HE  = new EffCorr(finalhistos,"GoodQTrackEffCorr_HE","GoodQTrack Eff. Corr",PRB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodQTrack_TOF = new EffCorr(finalhistos,"GoodQTrackEffCorr_TOF","GoodQTrack Eff. Corr",ToFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodQTrack_NaF = new EffCorr(finalhistos,"GoodQTrackEffCorr_NaF","GoodQTrack Eff. Corr",NaFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodQTrack_Agl = new EffCorr(finalhistos,"GoodQTrackEffCorr_Agl","GoodQTrack Eff. Corr",AglPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");


	before = "IsDownGoing&IsGoodTrack&IsCharge1Track&IsLooseCharge1";
        after  = "IsDownGoing&IsGoodTrack&IsCharge1Track&IsLooseCharge1&IsGoodChi2"; 
	EffCorr * GoodChi_HE  = new EffCorr(finalhistos,"GoodChiEffCorr_HE","GoodChi Eff. Corr",PRB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodChi_TOF = new EffCorr(finalhistos,"GoodChiEffCorr_TOF","GoodChi Eff. Corr",ToFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodChi_NaF = new EffCorr(finalhistos,"GoodChiEffCorr_NaF","GoodChi Eff. Corr",NaFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodChi_Agl = new EffCorr(finalhistos,"GoodChiEffCorr_Agl","GoodChi Eff. Corr",AglPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");

        before = "IsDownGoing&IsPhysTrig&IsGoodTOFStandaloneQ1&IsExtrapolInsideL8";
        after  = "IsDownGoing&IsPhysTrig&IsGoodTrack&IsGoodTOFStandaloneQ1&IsExtrapolInsideL8"; 
	EffCorr * TrackerEffCorr_HE = new EffCorr(finalhistos,"TrackerEffCorr_HE","Tracker Eff. Corr",PRB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");

	// Good Z=1 and Golden Efficiency corrections
		  
	before = "IsPositive&IsMinimumBias&IsLooseCharge1";
        after  = "IsPositive&IsMinimumBias&IsLooseCharge1&IsMinTOF";
	EffCorr * MinTOFEffCorr_HE  = new EffCorr(finalhistos,"MinTOFEffCorr_HE" ,"Min TOF Eff. Corr",PRB,before,after,  "IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC"); 
	EffCorr * MinTOFEffCorr_TOF = new EffCorr(finalhistos,"MinTOFEffCorr_TOF","Min TOF Eff. Corr",ToFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * MinTOFEffCorr_NaF = new EffCorr(finalhistos,"MinTOFEffCorr_NaF","Min TOF Eff. Corr",NaFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * MinTOFEffCorr_Agl = new EffCorr(finalhistos,"MinTOFEffCorr_Agl","Min TOF Eff. Corr",AglPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");

	before = "IsPositive&IsMinimumBias&IsLooseCharge1&IsMinTOF";
        after  = "IsPositive&IsMinimumBias&IsLooseCharge1&IsMinTOF&Is1TrTrack"; 
	EffCorr * Good1Track_HE  = new EffCorr(finalhistos,"Good1TrackEffCorr_HE","Good1Track Eff. Corr",PRB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * Good1Track_TOF = new EffCorr(finalhistos,"Good1TrackEffCorr_TOF","Good1Track Eff. Corr",ToFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * Good1Track_NaF = new EffCorr(finalhistos,"Good1TackEffCorr_NaF", "Good1Track Eff. Corr",NaFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * Good1Track_Agl = new EffCorr(finalhistos,"Good1TrackEffCorr_Agl","Good1Track Eff. Corr",AglPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");

	before = "IsPositive&IsMinimumBias&IsLooseCharge1&Is1TrTrack&IsMinTOF";
        after  = "IsPositive&IsMinimumBias&IsLooseCharge1&Is1TrTrack&IsMinTOF&IsCharge1UTOF"; 
	EffCorr * GoodUtof_HE  = new EffCorr(finalhistos,"GoodUtofEffCorr_HE" ,"GoodUtof Eff. Corr",PRB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodUtof_TOF = new EffCorr(finalhistos,"GoodUtofEffCorr_TOF","GoodUtof Eff. Corr",ToFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodUtof_NaF = new EffCorr(finalhistos,"GoodUtofEffCorr_NaF","GoodUtof Eff. Corr",NaFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodUtof_Agl = new EffCorr(finalhistos,"GoodUtofEffCorr_Agl","GoodUtof Eff. Corr",AglPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");

	before = "IsPositive&IsMinimumBias&IsLooseCharge1&Is1TrTrack&IsMinTOF&IsCharge1UTOF";
        after  = "IsPositive&IsMinimumBias&IsLooseCharge1&Is1TrTrack&IsMinTOF&IsCharge1UTOF&IsCharge1LTOF";
        EffCorr * GoodLtof_HE  = new EffCorr(finalhistos,"GoodLTOFEffCorr_HE" ,"GoodLtof Eff. Corr",PRB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
        EffCorr * GoodLtof_TOF = new EffCorr(finalhistos,"GoodLTOFEffCorr_TOF","GoodLtof Eff. Corr",ToFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
        EffCorr * GoodLtof_NaF = new EffCorr(finalhistos,"GoodLTOFEffCorr_NaF","GoodLtof Eff. Corr",NaFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
        EffCorr * GoodLtof_Agl = new EffCorr(finalhistos,"GoodLTOFEffCorr_Agl","GoodLtof Eff. Corr",AglPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");

	before = "IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning";
        after  = "IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsGoodTime"; 
	EffCorr * GoodTime_TOF = new EffCorr(finalhistos,"GoodTimeEffCorr_TOF","GoodTime Eff. Corr",ToFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	
	before = "IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning";
        after  = "IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning";
	EffCorr * RICHEffCorr_NaF = new EffCorr(finalhistos,"RICHCorrection_NaF","RICH Eff. Corr",NaFPB,before,(after+"&IsFromNaF").c_str(),"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * RICHEffCorr_Agl = new EffCorr(finalhistos,"RICHCorrection_Agl","RICH Eff. Corr",AglPB,before,(after+"&IsFromAgl").c_str(),"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");

	before = "IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning";
        after  = "IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning";
	EffCorr * RICHQualEffCorr_NaF = new EffCorr(finalhistos,"RICHQualCorrection_NaF","RICH Qual Eff. Corr",NaFPB,(before+"&IsFromNaF").c_str(),(after+"&IsFromNaF&RICHBDTCut").c_str(),"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * RICHQualEffCorr_Agl = new EffCorr(finalhistos,"RICHqualCorrection_Agl","RICH Qual. Eff. Corr",AglPB,(before+"&IsFromAgl").c_str(),(after+"&IsFromAgl&RICHBDTCut").c_str(),"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");

	TriggerEffCorr_HE 	->SetDefaultOutFile(finalhistos); 
	TriggerEffCorr_TOF ->SetDefaultOutFile(finalhistos); 
	TriggerEffCorr_NaF ->SetDefaultOutFile(finalhistos); 
	TriggerEffCorr_Agl ->SetDefaultOutFile(finalhistos); 
	
	L1PickUpEffCorr_HE ->SetDefaultOutFile(finalhistos); 
	L1PickUpEffCorr_TOF->SetDefaultOutFile(finalhistos); 
	L1PickUpEffCorr_NaF->SetDefaultOutFile(finalhistos); 
	L1PickUpEffCorr_Agl->SetDefaultOutFile(finalhistos); 
	
	MinTOFEffCorr_HE  ->SetDefaultOutFile(finalhistos); 
	MinTOFEffCorr_TOF ->SetDefaultOutFile(finalhistos); 
	MinTOFEffCorr_NaF ->SetDefaultOutFile(finalhistos); 
	MinTOFEffCorr_Agl ->SetDefaultOutFile(finalhistos); 
	
	TrackerEffCorr_HE ->SetDefaultOutFile(finalhistos); 
	
	GoodChi_HE   ->SetDefaultOutFile(finalhistos); 
	GoodChi_TOF  ->SetDefaultOutFile(finalhistos); 
	GoodChi_NaF    ->SetDefaultOutFile(finalhistos); 
	GoodChi_Agl ->SetDefaultOutFile(finalhistos); 
	
	GoodUtof_HE ->SetDefaultOutFile(finalhistos); 
	GoodUtof_TOF ->SetDefaultOutFile(finalhistos); 
	GoodUtof_NaF ->SetDefaultOutFile(finalhistos); 
	GoodUtof_Agl ->SetDefaultOutFile(finalhistos); 
	
	GoodLtof_HE  ->SetDefaultOutFile(finalhistos); 
	GoodLtof_TOF ->SetDefaultOutFile(finalhistos); 
	GoodLtof_NaF ->SetDefaultOutFile(finalhistos); 
	GoodLtof_Agl ->SetDefaultOutFile(finalhistos); 
	
	GoodQTrack_HE  ->SetDefaultOutFile(finalhistos); 
	GoodQTrack_TOF ->SetDefaultOutFile(finalhistos); 
	GoodQTrack_NaF ->SetDefaultOutFile(finalhistos); 
	GoodQTrack_Agl ->SetDefaultOutFile(finalhistos); 
	
	Good1Track_HE  ->SetDefaultOutFile(finalhistos); 
	Good1Track_TOF ->SetDefaultOutFile(finalhistos); 
	Good1Track_NaF ->SetDefaultOutFile(finalhistos); 
	Good1Track_Agl ->SetDefaultOutFile(finalhistos); 
	
	GoodTime_TOF 	->SetDefaultOutFile(finalhistos); 
	RICHEffCorr_NaF 	->SetDefaultOutFile(finalhistos); 
	RICHEffCorr_Agl 	->SetDefaultOutFile(finalhistos); 
	RICHQualEffCorr_NaF 	->SetDefaultOutFile(finalhistos); 
	RICHQualEffCorr_Agl 	->SetDefaultOutFile(finalhistos); 




	Filler.AddObject2beFilled(TriggerEffCorr_HE,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(TriggerEffCorr_TOF,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(TriggerEffCorr_NaF,GetRigidity,GetRigidity);	
	Filler.AddObject2beFilled(TriggerEffCorr_Agl,GetRigidity,GetRigidity);

	Filler.AddObject2beFilled(L1PickUpEffCorr_HE,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(L1PickUpEffCorr_TOF,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(L1PickUpEffCorr_NaF,GetRigidity,GetRigidity);	
	Filler.AddObject2beFilled(L1PickUpEffCorr_Agl,GetRigidity,GetRigidity);

	Filler.AddObject2beFilled(MinTOFEffCorr_HE , GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(MinTOFEffCorr_TOF,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(MinTOFEffCorr_NaF,GetRigidity,GetRigidity);	
	Filler.AddObject2beFilled(MinTOFEffCorr_Agl,GetRigidity,GetRigidity);

	Filler.AddObject2beFilled(TrackerEffCorr_HE,GetMomentumProxy,GetMomentumProxy);

	Filler.AddObject2beFilled(GoodChi_HE,GetRigidity,GetRigidity);	
	Filler.AddObject2beFilled(GoodChi_TOF,GetRigidity,GetRigidity);	
	Filler.AddObject2beFilled(GoodChi_NaF,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(GoodChi_Agl,GetRigidity,GetRigidity);

	Filler.AddObject2beFilled(GoodUtof_HE,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(GoodUtof_TOF,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(GoodUtof_NaF,GetRigidity,GetRigidity);	
	Filler.AddObject2beFilled(GoodUtof_Agl,GetRigidity,GetRigidity);

	Filler.AddObject2beFilled(GoodLtof_HE,GetRigidity,GetRigidity);	
	Filler.AddObject2beFilled(GoodLtof_TOF,GetRigidity,GetRigidity);	
	Filler.AddObject2beFilled(GoodLtof_NaF,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(GoodLtof_Agl,GetRigidity,GetRigidity);	

	Filler.AddObject2beFilled(GoodQTrack_HE,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(GoodQTrack_TOF,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(GoodQTrack_NaF,GetRigidity,GetRigidity);	
	Filler.AddObject2beFilled(GoodQTrack_Agl,GetRigidity,GetRigidity);

	Filler.AddObject2beFilled(Good1Track_HE,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(Good1Track_TOF,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(Good1Track_NaF,GetRigidity,GetRigidity);	
	Filler.AddObject2beFilled(Good1Track_Agl,GetRigidity,GetRigidity);

	Filler.AddObject2beFilled(GoodTime_TOF,GetRigidity,GetRigidity);	
	Filler.AddObject2beFilled(RICHEffCorr_NaF,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(RICHEffCorr_Agl,GetRigidity,GetRigidity);	
	Filler.AddObject2beFilled(RICHQualEffCorr_NaF,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(RICHQualEffCorr_Agl,GetRigidity,GetRigidity);	

	Filler.ReinitializeAll(refill);

	if(!refill&&checkfile) {	

		AnalyzeEffCorr(	TriggerEffCorr_HE  , finalhistos, finalresults);
		AnalyzeEffCorr(	TriggerEffCorr_TOF , finalhistos, finalresults);
		AnalyzeEffCorr(	TriggerEffCorr_NaF , finalhistos, finalresults);
		AnalyzeEffCorr(	TriggerEffCorr_Agl , finalhistos, finalresults);

		AnalyzeEffCorr(	L1PickUpEffCorr_HE , finalhistos, finalresults);
		AnalyzeEffCorr(	L1PickUpEffCorr_TOF, finalhistos, finalresults);
		AnalyzeEffCorr(	L1PickUpEffCorr_NaF, finalhistos, finalresults);
		AnalyzeEffCorr(	L1PickUpEffCorr_Agl, finalhistos, finalresults);

		AnalyzeEffCorr(	MinTOFEffCorr_HE  , finalhistos, finalresults);
		AnalyzeEffCorr(	MinTOFEffCorr_TOF , finalhistos, finalresults);
		AnalyzeEffCorr(	MinTOFEffCorr_NaF , finalhistos, finalresults);
		AnalyzeEffCorr(	MinTOFEffCorr_Agl , finalhistos, finalresults);

		AnalyzeEffCorr(	TrackerEffCorr_HE , finalhistos, finalresults);

		AnalyzeEffCorr(	GoodChi_HE  , finalhistos, finalresults);  
		AnalyzeEffCorr(	GoodChi_TOF , finalhistos, finalresults);  
		AnalyzeEffCorr(	GoodChi_NaF , finalhistos, finalresults);    
		AnalyzeEffCorr(	GoodChi_Agl , finalhistos, finalresults);

		AnalyzeEffCorr(	GoodUtof_HE , finalhistos, finalresults);
		AnalyzeEffCorr(	GoodUtof_TOF , finalhistos, finalresults);
		AnalyzeEffCorr(	GoodUtof_NaF , finalhistos, finalresults);
		AnalyzeEffCorr(	GoodUtof_Agl , finalhistos, finalresults);

		AnalyzeEffCorr(	GoodLtof_HE  , finalhistos, finalresults);
		AnalyzeEffCorr(	GoodLtof_TOF , finalhistos, finalresults);
		AnalyzeEffCorr(	GoodLtof_NaF , finalhistos, finalresults);
		AnalyzeEffCorr(	GoodLtof_Agl , finalhistos, finalresults);

		AnalyzeEffCorr(	GoodQTrack_HE , finalhistos, finalresults);
		AnalyzeEffCorr(	GoodQTrack_TOF , finalhistos, finalresults);
		AnalyzeEffCorr(	GoodQTrack_NaF , finalhistos, finalresults);
		AnalyzeEffCorr(	GoodQTrack_Agl , finalhistos, finalresults);

		AnalyzeEffCorr(	Good1Track_HE , finalhistos, finalresults);
		AnalyzeEffCorr(	Good1Track_TOF , finalhistos, finalresults);
		AnalyzeEffCorr(	Good1Track_NaF , finalhistos, finalresults);
		AnalyzeEffCorr(	Good1Track_Agl , finalhistos, finalresults);

		AnalyzeEffCorr(	GoodTime_TOF , finalhistos, finalresults);
		AnalyzeEffCorr(	RICHEffCorr_NaF , finalhistos, finalresults);
		AnalyzeEffCorr(	RICHEffCorr_Agl , finalhistos, finalresults);
		AnalyzeEffCorr(	RICHQualEffCorr_NaF , finalhistos, finalresults);
		AnalyzeEffCorr(	RICHQualEffCorr_Agl , finalhistos, finalresults);
	}

	return ;
}


void AnalyzeEffCorr(EffCorr * Correction, FileSaver  finalhistos, FileSaver  finalresults){
	Correction   -> Eval_Efficiencies();
	Correction   -> Eval_Corrections();
	Correction   -> SaveResults(finalresults);
}	


