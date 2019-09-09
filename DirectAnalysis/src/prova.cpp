#include "Analyzer.h"
#include "AllRangesEfficiency.h"
#include "TGraphErrors.h"
#include "filesaver.h"


void Analyzer::BookEfficiencyAnalysis(FileSaver finalhistos, FileSaver finalresults, bool refill)
{

	cout<<"****************************** BINS ***************************************"<<endl;
    	SetUpTOIBinning();

	bool checkfile = finalhistos.CheckFile();
	check_file = checkfile;
	cout<<"****************************** Efficiecny ANALYIS ******************************************"<<endl;

       Efficiency * For_Acceptance_P = new Efficiency(finalhistos,"For_Acceptance_P","For_Acceptance_P",ForAcceptance,"IsProtonMC","IsProtonMC");
       Efficiency * For_Acceptance_D = new Efficiency(finalhistos,"For_Acceptance_D","For_Acceptance_D",ForAcceptance,"IsDeutonMC","IsDeutonMC");


       Efficiency * CountsTestEff[10];
       Efficiency * CountsTestEff_rig[10];


	for(int i=0;i<10;i++) {
		CountsTestEff[i] = new Efficiency(finalhistos,("CountsTestEff_v" + to_string(i)).c_str(),("CountsTestEff_v" + to_string(i)).c_str(),PRB,"IsProtonMC",("IsPositive&IsProtonMC&IsStandardSel_v"+to_string(i)).c_str());
	}		

	for(int i=0;i<10;i++) {
		CountsTestEff_rig[i] = new Efficiency(finalhistos,("CountsTestEff_rig_v" + to_string(i)).c_str(),("CountsTestEff_rig_v" + to_string(i)).c_str(),PRB,"IsProtonMC",("IsPositive&IsProtonMC&IsStandardSel_v"+to_string(i)).c_str());
	}		



	Efficiency * RigBinBaselineEff_Trig = new Efficiency(finalhistos,"RigBinBaselineEff_Trig","RigBinBaselineEff_Trig",PRB,"IsProtonMC","IsProtonMC&IsPhysTrig&IsPositive&IsBaseline");
	Efficiency * RigBinBaselineL1Eff_Trig = new Efficiency(finalhistos,"RigBinBaselineL1Eff_Trig","RigBinBaselineL1Eff_Trig",PRB,"IsProtonMC","IsProtonMC&IsPhysTrig&IsPositive&IsBaseline&L1LooseCharge1");
	Efficiency * RigBinQualEff     = new Efficiency(finalhistos,"RigBinQualEff"    ,"RigBinQualEff"    ,PRB,"IsProtonMC","IsProtonMC&IsPhysTrig&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning");



     	AllRangesEfficiency * FullSetTOT_P = new AllRangesEfficiency(finalhistos,"FullSetTOT_P","FullSetTOT",
	"IsProtonMC",			      
	"IsProtonMC",			   	      
	"IsProtonMC",
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime",			      
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut",			   	      
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut",
	ToFPB,NaFPB,AglPB,
	refill);

	AllRangesEfficiency * FullSetTOT_D = new AllRangesEfficiency(finalhistos,"FullSetTOT_D","FullSetTOT",
	"IsProtonMC",			      
	"IsProtonMC",			   	      
	"IsProtonMC",
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime",			      
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut",			   	      
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut",
	ToFDB,NaFDB,AglDB,
	refill);


	AllRangesEfficiency * Baseline_P = new AllRangesEfficiency(finalhistos,"Baseline_P","Baseline",
	"IsProtonMC",			      
	"IsProtonMC",			   	      
	"IsProtonMC",
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1",			      
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1",			   	      
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1",
		ToFDB,NaFDB,AglDB,
	refill);


	AllRangesEfficiency * Baseline_D = new AllRangesEfficiency(finalhistos,"Baseline_D","Baseline",
	"IsDeutonMC",			      
	"IsDeutonMC",			   	      
	"IsDeutonMC",
	"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1",			      
	"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1",			   	      
	"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1",
		ToFDB,NaFDB,AglDB,
	refill);


	Efficiency * Cascade0 = new Efficiency(finalhistos,"Cascade1","Cascade1",PRB,"IsProtonMC","IsProtonMC&IsPositive&IsPhysTrig");
	Efficiency * Cascade1 = new Efficiency(finalhistos,"Cascade1","Cascade1",PRB,"IsProtonMC","IsProtonMC&IsPositive&IsBaseline");
	Efficiency * Cascade2 = new Efficiency(finalhistos,"Cascade2","Cascade2",PRB,"IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1");
	Efficiency * Cascade3 = new Efficiency(finalhistos,"Cascade3","Cascade3",PRB,"IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning");
	Efficiency * Cascade4 = new Efficiency(finalhistos,"Cascade4","Cascade4",PRB,"IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime");
	Efficiency * Cascade5 = new Efficiency(finalhistos,"Cascade5","Cascade5",PRB,"IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF");
	Efficiency * Cascade6 = new Efficiency(finalhistos,"Cascade6","Cascade6",PRB,"IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl");
	Efficiency * Cascade7 = new Efficiency(finalhistos,"Cascade7","Cascade7",PRB,"IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut");
	Efficiency * Cascade8 = new Efficiency(finalhistos,"Cascade8","Cascade8",PRB,"IsProtonMC","IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut");

	AllRangesEfficiency * Trigger_P_PID = new AllRangesEfficiency(finalhistos,"Trigger_P_PID","Trigger",
	"IsPurePMC&IsPositive&IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge1Track&L1LooseCharge1",			      
	"IsPurePMC&IsPositive&IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge1Track&L1LooseCharge1",			   	      
	"IsPurePMC&IsPositive&IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge1Track&L1LooseCharge1",
	"IsPurePMC&IsPositive&IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge1Track&L1LooseCharge1&IsPhysTrig",			      
	"IsPurePMC&IsPositive&IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge1Track&L1LooseCharge1&IsPhysTrig",			   	      
	"IsPurePMC&IsPositive&IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge1Track&L1LooseCharge1&IsPhysTrig",
		ToFDB,NaFDB,AglDB,
	refill);


	AllRangesEfficiency * Trigger_D_PID = new AllRangesEfficiency(finalhistos,"Trigger_D_PID","Trigger",
	"IsPureDMC&IsPositive&IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge1Track&L1LooseCharge1",			      
	"IsPureDMC&IsPositive&IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge1Track&L1LooseCharge1",			   	      
	"IsPureDMC&IsPositive&IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge1Track&L1LooseCharge1",
	"IsPureDMC&IsPositive&IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge1Track&L1LooseCharge1&IsPhysTrig",			      
	"IsPureDMC&IsPositive&IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge1Track&L1LooseCharge1&IsPhysTrig",			   	      
	"IsPureDMC&IsPositive&IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge1Track&L1LooseCharge1&IsPhysTrig",
		ToFDB,NaFDB,AglDB,
	refill);


	AllRangesEfficiency * MinimBias_P_PID = new AllRangesEfficiency(finalhistos,"MinimBias_P_PID","MinimBias",
	"IsPurePMC&IsPhysTrig",			      
	"IsPurePMC&IsPhysTrig",			   	      
	"IsPurePMC&IsPhysTrig",
	"IsPurePMC&IsPositive&IsBaseline&L1LooseCharge1",			      
	"IsPurePMC&IsPositive&IsBaseline&L1LooseCharge1",			   	      
	"IsPurePMC&IsPositive&IsBaseline&L1LooseCharge1",
		ToFDB,NaFDB,AglDB,
	refill);


	AllRangesEfficiency * MinimBias_D_PID = new AllRangesEfficiency(finalhistos,"MinimBias_D_PID","MinimBias",
	"IsPureDMC&IsPhysTrig",			      
	"IsPureDMC&IsPhysTrig",			   	      
	"IsPureDMC&IsPhysTrig",
	"IsPureDMC&IsPositive&IsBaseline&L1LooseCharge1",			      
	"IsPureDMC&IsPositive&IsBaseline&L1LooseCharge1",			   	      
	"IsPureDMC&IsPositive&IsBaseline&L1LooseCharge1",
		ToFDB,NaFDB,AglDB,
	refill);


	AllRangesEfficiency * Cleaning_P_PID = new AllRangesEfficiency(finalhistos,"Cleaning_P_PID","Cleaning",
	"IsPurePMC&IsPositive&IsBaseline&L1LooseCharge1",			      
	"IsPurePMC&IsPositive&IsBaseline&L1LooseCharge1",			   	      
	"IsPurePMC&IsPositive&IsBaseline&L1LooseCharge1",
	"IsPurePMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",			      
	"IsPurePMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",			   	      
	"IsPurePMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",
		ToFDB,NaFDB,AglDB,
	refill);


	AllRangesEfficiency * Cleaning_D_PID = new AllRangesEfficiency(finalhistos,"Cleaning_D_PID","Cleaning",
	"IsPureDMC&IsPositive&IsBaseline&L1LooseCharge1",			      
	"IsPureDMC&IsPositive&IsBaseline&L1LooseCharge1",			   	      
	"IsPureDMC&IsPositive&IsBaseline&L1LooseCharge1",
	"IsPureDMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",			      
	"IsPureDMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",			   	      
	"IsPureDMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",
		ToFDB,NaFDB,AglDB,
	refill);


	AllRangesEfficiency * RICH_P_PID = new AllRangesEfficiency(finalhistos,"RICH_P_PID","RICH",
	"IsPurePMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",			      
	"IsPurePMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",			   	      
	"IsPurePMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",
	"IsPurePMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",			      
	"IsPurePMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF",			   	      
	"IsPurePMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl",
		ToFDB,NaFDB,AglDB,
	refill);


	AllRangesEfficiency * RICH_D_PID = new AllRangesEfficiency(finalhistos,"RICH_D_PID","RICH",
	"IsPureDMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",			      
	"IsPureDMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",			   	      
	"IsPureDMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",
	"IsPureDMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",			      
	"IsPureDMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF",			   	      
	"IsPureDMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl",
		ToFDB,NaFDB,AglDB,
	refill);


	AllRangesEfficiency * RICHQual_P_PID = new AllRangesEfficiency(finalhistos,"RICHQual_P_PID","RICHQual",
	"IsPurePMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",			      
	"IsPurePMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF",			   	      
	"IsPurePMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl",
	"IsPurePMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",			      
	"IsPurePMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut",			   	      
	"IsPurePMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut",
		ToFDB,NaFDB,AglDB,
	refill);


	AllRangesEfficiency * RICHQual_D_PID = new AllRangesEfficiency(finalhistos,"RICHQual_D_PID","RICHQual",
	"IsPureDMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",			      
	"IsPureDMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF",			   	      
	"IsPureDMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl",
	"IsPureDMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",			      
	"IsPureDMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut",			   	      
	"IsPureDMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut",
		ToFDB,NaFDB,AglDB,
	refill);


	AllRangesEfficiency * GoldenTOF_P_PID = new AllRangesEfficiency(finalhistos,"GoldenTOF_P_PID","GoldenTOF",
	"IsPurePMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",			      
	"IsPurePMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF",			   	      
	"IsPurePMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl",
	"IsPurePMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime",			      
	"IsPurePMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut",			   	      
	"IsPurePMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut",
		ToFDB,NaFDB,AglDB,
	refill);
	//GoldenTOF_P_PID->SetBins(ToFPB,NaFPB,AglPB);


	AllRangesEfficiency * GoldenTOF_D_PID = new AllRangesEfficiency(finalhistos,"GoldenTOF_D_PID","GoldenTOF",
	"IsPureDMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",			      
	"IsPureDMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF",			   	      
	"IsPureDMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl",
	"IsPureDMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime",			      
	"IsPureDMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut",			   	      
	"IsPureDMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut",
		ToFDB,NaFDB,AglDB,
	refill);


	AllRangesEfficiency * Trigger_P = new AllRangesEfficiency(finalhistos,"Trigger_P","Trigger",
	"IsProtonMC&IsPositive&IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge1Track&L1LooseCharge1",			      
	"IsProtonMC&IsPositive&IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge1Track&L1LooseCharge1",			   	      
	"IsProtonMC&IsPositive&IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge1Track&L1LooseCharge1",
	"IsProtonMC&IsPositive&IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge1Track&L1LooseCharge1&IsPhysTrig",			      
	"IsProtonMC&IsPositive&IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge1Track&L1LooseCharge1&IsPhysTrig",			   	      
	"IsProtonMC&IsPositive&IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge1Track&L1LooseCharge1&IsPhysTrig",
		ToFDB,NaFDB,AglDB,
	refill);


	AllRangesEfficiency * Trigger_D = new AllRangesEfficiency(finalhistos,"Trigger_D","Trigger",
	"IsDeutonMC&IsPositive&IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge1Track&L1LooseCharge1",			      
	"IsDeutonMC&IsPositive&IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge1Track&L1LooseCharge1",			   	      
	"IsDeutonMC&IsPositive&IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge1Track&L1LooseCharge1",
	"IsDeutonMC&IsPositive&IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge1Track&L1LooseCharge1&IsPhysTrig",			      
	"IsDeutonMC&IsPositive&IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge1Track&L1LooseCharge1&IsPhysTrig",			   	      
	"IsDeutonMC&IsPositive&IsDownGoing&IsGoodTrack&IsGoodChi2&IsCharge1Track&L1LooseCharge1&IsPhysTrig",
		ToFDB,NaFDB,AglDB,
	refill);


	AllRangesEfficiency * MinimBias_P = new AllRangesEfficiency(finalhistos,"MinimBias_P","MinimBias",
	"IsProtonMC&IsPhysTrig",			      
	"IsProtonMC&IsPhysTrig",			   	      
	"IsProtonMC&IsPhysTrig",
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1",			      
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1",			   	      
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1",
		ToFDB,NaFDB,AglDB,
	refill);


	AllRangesEfficiency * MinimBias_D = new AllRangesEfficiency(finalhistos,"MinimBias_D","MinimBias",
	"IsDeutonMC&IsPhysTrig",			      
	"IsDeutonMC&IsPhysTrig",			   	      
	"IsDeutonMC&IsPhysTrig",
	"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1",			      
	"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1",			   	      
	"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1",
		ToFDB,NaFDB,AglDB,
	refill);


	AllRangesEfficiency * Cleaning_P = new AllRangesEfficiency(finalhistos,"Cleaning_P","Cleaning",
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1",			      
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1",			   	      
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1",
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",			      
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",			   	      
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",
		ToFDB,NaFDB,AglDB,
	refill);


	AllRangesEfficiency * Cleaning_D = new AllRangesEfficiency(finalhistos,"Cleaning_D","Cleaning",
	"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1",			      
	"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1",			   	      
	"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1",
	"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",			      
	"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",			   	      
	"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",
		ToFDB,NaFDB,AglDB,
	refill);


	AllRangesEfficiency * RICH_P = new AllRangesEfficiency(finalhistos,"RICH_P","RICH",
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",			      
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",			   	      
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",			      
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF",			   	      
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl",
		ToFDB,NaFDB,AglDB,
	refill);


	AllRangesEfficiency * RICH_D = new AllRangesEfficiency(finalhistos,"RICH_D","RICH",
	"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",			      
	"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",			   	      
	"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",
	"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",			      
	"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF",			   	      
	"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl",
		ToFDB,NaFDB,AglDB,
	refill);


	AllRangesEfficiency * RICHQual_P = new AllRangesEfficiency(finalhistos,"RICHQual_P","RICHQual",
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",			      
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF",			   	      
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl",
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",			      
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut",			   	      
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut",
		ToFDB,NaFDB,AglDB,
	refill);


	AllRangesEfficiency * RICHQual_D = new AllRangesEfficiency(finalhistos,"RICHQual_D","RICHQual",
	"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",			      
	"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF",			   	      
	"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl",
	"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",			      
	"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut",			   	      
	"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut",
		ToFDB,NaFDB,AglDB,
	refill);


	AllRangesEfficiency * GoldenTOF_P = new AllRangesEfficiency(finalhistos,"GoldenTOF_P","GoldenTOF",
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",			      
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF",			   	      
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl",
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime",			      
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut",			   	      
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut",
		ToFDB,NaFDB,AglDB,
	refill);


	AllRangesEfficiency * GoldenTOF_D = new AllRangesEfficiency(finalhistos,"GoldenTOF_D","GoldenTOF",
	"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning",			      
	"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF",			   	      
	"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl",
	"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsGoodTime",			      
	"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut",			   	      
	"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut",
		ToFDB,NaFDB,AglDB,
	refill);


	AllRangesEfficiency * Fragmentation_P = new AllRangesEfficiency(finalhistos,"Fragmentation_P","Fragmentation",
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1",			      
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1",			   	      
	"IsProtonMC&IsPositive&IsBaseline&L1LooseCharge1",
	"IsPurePMC&IsPositive&IsBaseline&L1LooseCharge1",			      
	"IsPurePMC&IsPositive&IsBaseline&L1LooseCharge1",			   	      
	"IsPurePMC&IsPositive&IsBaseline&L1LooseCharge1",
		ToFDB,NaFDB,AglDB,
	refill);


	AllRangesEfficiency * Fragmentation_D = new AllRangesEfficiency(finalhistos,"Fragmentation_D","Fragmentation",
	"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1",			      
	"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1",			   	      
	"IsDeutonMC&IsPositive&IsBaseline&L1LooseCharge1",
	"IsPureDMC&IsPositive&IsBaseline&L1LooseCharge1",			      
	"IsPureDMC&IsPositive&IsBaseline&L1LooseCharge1",			   	      
	"IsPureDMC&IsPositive&IsBaseline&L1LooseCharge1",
		ToFDB,NaFDB,AglDB,
	refill);


	For_Acceptance_P->SetNotWeightedMC();
	For_Acceptance_D->SetNotWeightedMC();
//	RigBinBaselineEff_Trig->SetNotWeightedMC();
//	RigBinBaselineL1Eff_Trig->SetNotWeightedMC();
//	RigBinQualEff->SetNotWeightedMC();
	for(int i=0;i<10;i++) CountsTestEff[i]->SetNotWeightedMC();	
//	for(int i=0;i<10;i++) CountsTestEff_rig[i]->SetNotWeightedMC();	
//	FullSetTOT_P->SetNotWeightedMC();	
//	FullSetTOT_D->SetNotWeightedMC();	
//	Baseline_P->SetNotWeightedMC();	
//	Baseline_D->SetNotWeightedMC();	


	For_Acceptance_P->SetDefaultOutFile(finalhistos);
	For_Acceptance_D->SetDefaultOutFile(finalhistos);
	RigBinBaselineEff_Trig->SetDefaultOutFile(finalhistos);
	RigBinBaselineL1Eff_Trig->SetDefaultOutFile(finalhistos);
	RigBinQualEff    ->SetDefaultOutFile(finalhistos);
	for(int i=0;i<10;i++) CountsTestEff[i]->SetDefaultOutFile(finalhistos);	
	for(int i=0;i<10;i++) CountsTestEff_rig[i]->SetDefaultOutFile(finalhistos);	


	FullSetTOT_P 	 ->SetDefaultOutFile(finalhistos);
	FullSetTOT_D 	 ->SetDefaultOutFile(finalhistos);
	Baseline_P 	 ->SetDefaultOutFile(finalhistos);
	Baseline_P 	 ->SetDefaultOutFile(finalhistos);
	Baseline_D 	 ->SetDefaultOutFile(finalhistos);
	Cascade0         ->SetDefaultOutFile(finalhistos);
	Cascade1  	 ->SetDefaultOutFile(finalhistos);
	Cascade2  	 ->SetDefaultOutFile(finalhistos);
	Cascade3 	 ->SetDefaultOutFile(finalhistos);
	Cascade4 	 ->SetDefaultOutFile(finalhistos);
	Cascade5 	 ->SetDefaultOutFile(finalhistos);
	Cascade6 	 ->SetDefaultOutFile(finalhistos);
	Cascade7 	 ->SetDefaultOutFile(finalhistos);
	Cascade8 	 ->SetDefaultOutFile(finalhistos);
	Trigger_P  	 ->SetDefaultOutFile(finalhistos);
	Trigger_D 	 ->SetDefaultOutFile(finalhistos);
	MinimBias_P 	 ->SetDefaultOutFile(finalhistos);
	MinimBias_D 	 ->SetDefaultOutFile(finalhistos);
	Cleaning_P 	 ->SetDefaultOutFile(finalhistos);
	Cleaning_D 	 ->SetDefaultOutFile(finalhistos);
	RICH_P 		 ->SetDefaultOutFile(finalhistos);
	RICH_D 		 ->SetDefaultOutFile(finalhistos);
	RICHQual_P 	 ->SetDefaultOutFile(finalhistos);
	RICHQual_D 	 ->SetDefaultOutFile(finalhistos);
	GoldenTOF_P	 ->SetDefaultOutFile(finalhistos);
	GoldenTOF_D	 ->SetDefaultOutFile(finalhistos);
	Trigger_P_PID 	 ->SetDefaultOutFile(finalhistos);
	Trigger_D_PID 	 ->SetDefaultOutFile(finalhistos);
	MinimBias_P_PID  ->SetDefaultOutFile(finalhistos);
	MinimBias_D_PID  ->SetDefaultOutFile(finalhistos);
	Cleaning_P_PID 	 ->SetDefaultOutFile(finalhistos);
	Cleaning_D_PID 	 ->SetDefaultOutFile(finalhistos);
	RICH_P_PID 	 ->SetDefaultOutFile(finalhistos);
	RICH_D_PID 	 ->SetDefaultOutFile(finalhistos);
	RICHQual_P_PID 	 ->SetDefaultOutFile(finalhistos);
	RICHQual_D_PID 	 ->SetDefaultOutFile(finalhistos);
	GoldenTOF_P_PID	 ->SetDefaultOutFile(finalhistos);
	GoldenTOF_D_PID	 ->SetDefaultOutFile(finalhistos);
	Fragmentation_P	 ->SetDefaultOutFile(finalhistos);
	Fragmentation_D  ->SetDefaultOutFile(finalhistos);


	Filler.AddObject2beFilled(For_Acceptance_P      ,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(For_Acceptance_D      ,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(RigBinBaselineEff_Trig,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(RigBinBaselineL1Eff_Trig,GetRigidity,GetRigidity);
	
	Filler.AddObject2beFilled(RigBinQualEff         ,GetRigidity,GetRigidity);
	for(int i=0;i<10;i++) Filler.AddObject2beFilled(CountsTestEff[i]        ,GetGenMomentum,GetGenMomentum);
	for(int i=0;i<10;i++) Filler.AddObject2beFilled(CountsTestEff_rig[i]        ,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(Cascade0,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(Cascade1,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(Cascade2,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(Cascade3,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(Cascade4,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(Cascade5,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(Cascade6,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(Cascade7,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(Cascade8,GetGenMomentum,GetGenMomentum);


	Filler.AddObject2beFilled(FullSetTOT_P,GetBetaMeas,GetBetaMeas);
	Filler.AddObject2beFilled(FullSetTOT_D,GetBetaMeas,GetBetaMeas);
	Filler.AddObject2beFilled(Baseline_P,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(Baseline_D,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(Trigger_P,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(Trigger_D,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(MinimBias_P,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(MinimBias_D,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(Cleaning_P,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(Cleaning_D,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(RICH_P,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(RICH_D,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(RICHQual_P,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(RICHQual_D,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(GoldenTOF_P,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(GoldenTOF_D,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(Trigger_P_PID,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(Trigger_D_PID,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(MinimBias_P_PID,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(MinimBias_D_PID,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(Cleaning_P_PID,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(Cleaning_D_PID,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(RICH_P_PID,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(RICH_D_PID,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(RICHQual_P_PID,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(RICHQual_D_PID,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(GoldenTOF_P_PID,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(GoldenTOF_D_PID,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(Fragmentation_P,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(Fragmentation_D,GetBetaGen,GetBetaGen);

	Filler.ReinitializeAll(refill);

	if(!refill&&checkfile) {	
		For_Acceptance_P->Eval_Efficiency();
		For_Acceptance_D->Eval_Efficiency();
	
		RigBinBaselineEff_Trig->Eval_Efficiency();
		RigBinBaselineL1Eff_Trig->Eval_Efficiency();
		RigBinQualEff->Eval_Efficiency();
		for(int i=0;i<10;i++) CountsTestEff[i]->Eval_Efficiency();
		for(int i=0;i<10;i++) CountsTestEff_rig[i]->Eval_Efficiency();
		FullSetTOT_P 	->Eval_Efficiency();
		FullSetTOT_D 	->Eval_Efficiency();
		Baseline_P 	->Eval_Efficiency();
		Baseline_D 	->Eval_Efficiency();
		Cascade0  	->Eval_Efficiency();
		Cascade1  	->Eval_Efficiency();
		Cascade2  	->Eval_Efficiency();
		Cascade3 	->Eval_Efficiency();
		Cascade4 	->Eval_Efficiency();
		Cascade5 	->Eval_Efficiency();
		Cascade6 	->Eval_Efficiency();
		Cascade7 	->Eval_Efficiency();
		Cascade8 	->Eval_Efficiency();
		Trigger_P 	->Eval_Efficiency();
		Trigger_D 	->Eval_Efficiency();
		MinimBias_P 	->Eval_Efficiency();
		MinimBias_D 	->Eval_Efficiency();
		Cleaning_P 	->Eval_Efficiency();
		Cleaning_D 	->Eval_Efficiency();
		RICH_P 	->Eval_Efficiency();
		RICH_D 	->Eval_Efficiency();
		RICHQual_P 	->Eval_Efficiency();
		RICHQual_D 	->Eval_Efficiency();
		GoldenTOF_P	->Eval_Efficiency();
		GoldenTOF_D     ->Eval_Efficiency();
		Trigger_P_PID 	->Eval_Efficiency();
		Trigger_D_PID 	->Eval_Efficiency();
		MinimBias_P_PID 	->Eval_Efficiency();
		MinimBias_D_PID 	->Eval_Efficiency();
		Cleaning_P_PID 	->Eval_Efficiency();
		Cleaning_D_PID 	->Eval_Efficiency();
		RICH_P_PID 	->Eval_Efficiency();
		RICH_D_PID 	->Eval_Efficiency();
		RICHQual_P_PID 	->Eval_Efficiency();
		RICHQual_D_PID 	->Eval_Efficiency();
		GoldenTOF_P_PID	->Eval_Efficiency();
		GoldenTOF_D_PID ->Eval_Efficiency();
		Fragmentation_P->Eval_Efficiency();
		Fragmentation_D->Eval_Efficiency();


		/*
		   FullSetTOT_P      ->Eval_StatError();
		   FullSetTOT_D      ->Eval_StatError();
		   FullSetTOT_P      ->Eval_SystError(FullSet_P,FullSet_D);
		   FullSetTOT_D      ->Eval_SystError(FullSet_P,FullSet_D);
		   */	

		RigBinBaselineEff_Trig->SaveResults(finalresults);
		RigBinBaselineL1Eff_Trig->SaveResults(finalresults);
		RigBinQualEff->SaveResults(finalresults);
		for(int i=0;i<10;i++) CountsTestEff[i]->Eval_Efficiency();
		for(int i=0;i<10;i++) CountsTestEff_rig[i]->Eval_Efficiency();


		For_Acceptance_P->SaveResults(finalresults);
		For_Acceptance_D->SaveResults(finalresults);
		FullSetTOT_P 	->SaveResults(finalresults);
		FullSetTOT_D 	->SaveResults(finalresults);
		Baseline_P 	->SaveResults(finalresults);
		Baseline_D 	->SaveResults(finalresults);
		for(int i=0;i<10;i++) CountsTestEff[i]->SaveResults(finalresults);
		for(int i=0;i<10;i++) CountsTestEff_rig[i]->SaveResults(finalresults);
		Cascade1        ->SaveResults(finalresults);
		Cascade1  	->SaveResults(finalresults);
		Cascade2  	->SaveResults(finalresults);
		Cascade3 	->SaveResults(finalresults);
		Cascade4 	->SaveResults(finalresults);
		Cascade5 	->SaveResults(finalresults);
		Cascade6 	->SaveResults(finalresults);
		Cascade7 	->SaveResults(finalresults);
		Cascade8 	->SaveResults(finalresults);
		Trigger_P 	->SaveResults(finalresults);
		Trigger_D 	->SaveResults(finalresults);
		MinimBias_P 	->SaveResults(finalresults);
		MinimBias_D 	->SaveResults(finalresults);
		Cleaning_P 	->SaveResults(finalresults);
		Cleaning_D 	->SaveResults(finalresults);
		RICH_P 	->SaveResults(finalresults);
		RICH_D 	->SaveResults(finalresults);
		RICHQual_P 	->SaveResults(finalresults);
		RICHQual_D 	->SaveResults(finalresults);
		GoldenTOF_P	->SaveResults(finalresults);
		GoldenTOF_D     ->SaveResults(finalresults);
		Trigger_P_PID 	->SaveResults(finalresults);
		Trigger_D_PID 	->SaveResults(finalresults);
		MinimBias_P_PID 	->SaveResults(finalresults);
		MinimBias_D_PID 	->SaveResults(finalresults);
		Cleaning_P_PID 	->SaveResults(finalresults);
		Cleaning_D_PID 	->SaveResults(finalresults);
		RICH_P_PID 	->SaveResults(finalresults);
		RICH_D_PID 	->SaveResults(finalresults);
		RICHQual_P_PID 	->SaveResults(finalresults);
		RICHQual_D_PID 	->SaveResults(finalresults);
		GoldenTOF_P_PID	->SaveResults(finalresults);
		GoldenTOF_D_PID ->SaveResults(finalresults);
		Fragmentation_P->SaveResults(finalresults);
		Fragmentation_D->SaveResults(finalresults);
	}

	return ;
}




