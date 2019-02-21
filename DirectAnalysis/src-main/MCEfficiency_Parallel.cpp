#include <bitset>
#include "TROOT.h"
#include "TNtuple.h"
#include <TSpline.h>
#include "../../include/binning.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include <TVector3.h>
#include "TMath.h"
#include <TFile.h>
#include "TFile.h"
#include "TH2.h"
#include "TF2.h"
#include <TVector3.h>
#include "TMath.h"
#include "TGraphErrors.h"
#include "TFractionFitter.h"
#include "TRandom3.h"
#include "TChain.h"
#include "../include/InputFileReader.h"
#include "../include/Globals.h"

#include "../include/Variables.hpp"
#include "../include/Cuts.h"
#include "../include/ParallelFiller.h"

#include "../include/filesaver.h"

#include "../include/Tool.h"
#include "../include/Efficiency.h"
#include "../include/AllRangesEfficiency.h"

int main(int argc, char * argv[])
{

	TH1::SetDefaultSumw2();
        cout<<"****************************** FILES OPENING ***************************************"<<endl;
	string INPUT1 = "";
	string INPUT2 = "";
	string OUTPUT = "";
	
	if(argc<=2) { 
		OUTPUT = argv[1];
	}	
	
	else {
	INPUT1 = argv[1];
	INPUT2 = argv[2];
	OUTPUT = argv[3];
	}
	string refill="";
	if(argc > 4 ) 	refill = argv[4];	
	
	bool Refill = false;
	if(refill!="") Refill=true;


        FileSaver finalHistos;
        finalHistos.setName(OUTPUT.c_str());
        bool checkfile = finalHistos.CheckFile();

	FileSaver finalResults;
        finalResults.setName((OUTPUT+"_Results").c_str());

	TChain * chainDT = InputFileReader(INPUT1.c_str(),"Event");
    	TChain * chainMC = InputFileReader(INPUT2.c_str(),"Event");


	cout<<"****************************** BINS ***************************************"<<endl;

        SetBins();

        PRB.Print();

        cout<<"**TOF**"<<endl;
        ToFDB.Print();

        cout<<"**NaF**"<<endl;
        NaFPB.Print();

        cout<<"**Agl**"<<endl;
        AglDB.Print();

        ToFDB.UseBetaEdges();
        NaFDB.UseBetaEdges();
        AglDB.UseBetaEdges();

        PRB.UseREdges();


        cout<<endl;


	cout<<"****************************** VARIABLES ***************************************"<<endl;

        Variables * vars = new Variables();

	cout<<"****************************** ANALYIS ******************************************"<<endl;


	Efficiency * RigBinBaselineEff = new Efficiency(finalHistos,"RigBinBaselineEff","RigBinBaselineEff",PRB,"IsProtonMC","IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1");
	Efficiency * RigBinQualEff     = new Efficiency(finalHistos,"RigBinQualEff"    ,"RigBinQualEff"    ,PRB,"IsProtonMC","IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning");


     	AllRangesEfficiency * FullSetTOT_P = new AllRangesEfficiency(finalHistos,"FullSetTOT_P","FullSetTOT",
	"IsProtonMC",			      
	"IsProtonMC",			   	      
	"IsProtonMC",
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsGoodTime",			      
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut",			   	      
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut",
	Refill);

	AllRangesEfficiency * FullSetTOT_D = new AllRangesEfficiency(finalHistos,"FullSetTOT_D","FullSetTOT",
	"IsDeutonMC",			      
	"IsDeutonMC",			   	      
	"IsDeutonMC",
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsGoodTime",			      
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut",			   	      
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut",
	Refill);

	AllRangesEfficiency * Baseline_P = new AllRangesEfficiency(finalHistos,"Baseline_P","Baseline",
	"IsProtonMC",			      
	"IsProtonMC",			   	      
	"IsProtonMC",
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1",			      
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1",			   	      
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1",
	Refill);

	AllRangesEfficiency * Baseline_D = new AllRangesEfficiency(finalHistos,"Baseline_D","Baseline",
	"IsDeutonMC",			      
	"IsDeutonMC",			   	      
	"IsDeutonMC",
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1",			      
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1",			   	      
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1",
	Refill);

	Efficiency * Cascade0 = new Efficiency(finalHistos,"Cascade1","Cascade1",PRB,"IsProtonMC","IsProtonMC&IsPositive&IsPhysTrig");
	Efficiency * Cascade1 = new Efficiency(finalHistos,"Cascade1","Cascade1",PRB,"IsProtonMC","IsProtonMC&IsPositive&IsMinimumBias");
	Efficiency * Cascade2 = new Efficiency(finalHistos,"Cascade2","Cascade2",PRB,"IsProtonMC","IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1");
	Efficiency * Cascade3 = new Efficiency(finalHistos,"Cascade3","Cascade3",PRB,"IsProtonMC","IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning");
	Efficiency * Cascade4 = new Efficiency(finalHistos,"Cascade4","Cascade4",PRB,"IsProtonMC","IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsGoodTime");
	Efficiency * Cascade5 = new Efficiency(finalHistos,"Cascade5","Cascade5",PRB,"IsProtonMC","IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF");
	Efficiency * Cascade6 = new Efficiency(finalHistos,"Cascade6","Cascade6",PRB,"IsProtonMC","IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl");
	Efficiency * Cascade7 = new Efficiency(finalHistos,"Cascade7","Cascade7",PRB,"IsProtonMC","IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut");
	Efficiency * Cascade8 = new Efficiency(finalHistos,"Cascade8","Cascade8",PRB,"IsProtonMC","IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut");

	AllRangesEfficiency * Trigger_P_PID = new AllRangesEfficiency(finalHistos,"Trigger_P_PID","Trigger",
	"IsPurePMC&IsPositive&IsMinimumBias_notrigg&IsLooseCharge1",			      
	"IsPurePMC&IsPositive&IsMinimumBias_notrigg&IsLooseCharge1",			   	      
	"IsPurePMC&IsPositive&IsMinimumBias_notrigg&IsLooseCharge1",
	"IsPurePMC&IsPositive&IsMinimumBias&IsLooseCharge1",			      
	"IsPurePMC&IsPositive&IsMinimumBias&IsLooseCharge1",			   	      
	"IsPurePMC&IsPositive&IsMinimumBias&IsLooseCharge1",
	Refill);

	AllRangesEfficiency * Trigger_D_PID = new AllRangesEfficiency(finalHistos,"Trigger_D_PID","Trigger",
	"IsPureDMC&IsPositive&IsMinimumBias_notrigg&IsLooseCharge1",			      
	"IsPureDMC&IsPositive&IsMinimumBias_notrigg&IsLooseCharge1",			   	      
	"IsPureDMC&IsPositive&IsMinimumBias_notrigg&IsLooseCharge1",
	"IsPureDMC&IsPositive&IsMinimumBias&IsLooseCharge1",			      
	"IsPureDMC&IsPositive&IsMinimumBias&IsLooseCharge1",			   	      
	"IsPureDMC&IsPositive&IsMinimumBias&IsLooseCharge1",
	Refill);

	AllRangesEfficiency * MinimBias_P_PID = new AllRangesEfficiency(finalHistos,"MinimBias_P_PID","MinimBias",
	"IsPurePMC&IsPhysTrig",			      
	"IsPurePMC&IsPhysTrig",			   	      
	"IsPurePMC&IsPhysTrig",
	"IsPurePMC&IsPositive&IsMinimumBias&IsLooseCharge1",			      
	"IsPurePMC&IsPositive&IsMinimumBias&IsLooseCharge1",			   	      
	"IsPurePMC&IsPositive&IsMinimumBias&IsLooseCharge1",
	Refill);

	AllRangesEfficiency * MinimBias_D_PID = new AllRangesEfficiency(finalHistos,"MinimBias_D_PID","MinimBias",
	"IsPureDMC&IsPhysTrig",			      
	"IsPureDMC&IsPhysTrig",			   	      
	"IsPureDMC&IsPhysTrig",
	"IsPureDMC&IsPositive&IsMinimumBias&IsLooseCharge1",			      
	"IsPureDMC&IsPositive&IsMinimumBias&IsLooseCharge1",			   	      
	"IsPureDMC&IsPositive&IsMinimumBias&IsLooseCharge1",
	Refill);

	AllRangesEfficiency * Cleaning_P_PID = new AllRangesEfficiency(finalHistos,"Cleaning_P_PID","Cleaning",
	"IsPurePMC&IsPositive&IsMinimumBias&IsLooseCharge1",			      
	"IsPurePMC&IsPositive&IsMinimumBias&IsLooseCharge1",			   	      
	"IsPurePMC&IsPositive&IsMinimumBias&IsLooseCharge1",
	"IsPurePMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",			      
	"IsPurePMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",			   	      
	"IsPurePMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",
	Refill);

	AllRangesEfficiency * Cleaning_D_PID = new AllRangesEfficiency(finalHistos,"Cleaning_D_PID","Cleaning",
	"IsPureDMC&IsPositive&IsMinimumBias&IsLooseCharge1",			      
	"IsPureDMC&IsPositive&IsMinimumBias&IsLooseCharge1",			   	      
	"IsPureDMC&IsPositive&IsMinimumBias&IsLooseCharge1",
	"IsPureDMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",			      
	"IsPureDMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",			   	      
	"IsPureDMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",
	Refill);

	AllRangesEfficiency * RICH_P_PID = new AllRangesEfficiency(finalHistos,"RICH_P_PID","RICH",
	"IsPurePMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",			      
	"IsPurePMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",			   	      
	"IsPurePMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",
	"IsPurePMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",			      
	"IsPurePMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF",			   	      
	"IsPurePMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl",
	Refill);

	AllRangesEfficiency * RICH_D_PID = new AllRangesEfficiency(finalHistos,"RICH_D_PID","RICH",
	"IsPureDMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",			      
	"IsPureDMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",			   	      
	"IsPureDMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",
	"IsPureDMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",			      
	"IsPureDMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF",			   	      
	"IsPureDMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl",
	Refill);

	AllRangesEfficiency * RICHQual_P_PID = new AllRangesEfficiency(finalHistos,"RICHQual_P_PID","RICHQual",
	"IsPurePMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",			      
	"IsPurePMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF",			   	      
	"IsPurePMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl",
	"IsPurePMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",			      
	"IsPurePMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut",			   	      
	"IsPurePMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut",
	Refill);

	AllRangesEfficiency * RICHQual_D_PID = new AllRangesEfficiency(finalHistos,"RICHQual_D_PID","RICHQual",
	"IsPureDMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",			      
	"IsPureDMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF",			   	      
	"IsPureDMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl",
	"IsPureDMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",			      
	"IsPureDMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut",			   	      
	"IsPureDMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut",
	Refill);

	AllRangesEfficiency * GoldenTOF_P_PID = new AllRangesEfficiency(finalHistos,"GoldenTOF_P_PID","GoldenTOF",
	"IsPurePMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",			      
	"IsPurePMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF",			   	      
	"IsPurePMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl",
	"IsPurePMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsGoodTime",			      
	"IsPurePMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut",			   	      
	"IsPurePMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut",
	Refill);

	AllRangesEfficiency * GoldenTOF_D_PID = new AllRangesEfficiency(finalHistos,"GoldenTOF_D_PID","GoldenTOF",
	"IsPureDMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",			      
	"IsPureDMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF",			   	      
	"IsPureDMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl",
	"IsPureDMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsGoodTime",			      
	"IsPureDMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut",			   	      
	"IsPureDMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut",
	Refill);


	AllRangesEfficiency * Trigger_P = new AllRangesEfficiency(finalHistos,"Trigger_P","Trigger",
	"IsProtonMC&IsPositive&IsMinimumBias_notrigg&IsLooseCharge1",			      
	"IsProtonMC&IsPositive&IsMinimumBias_notrigg&IsLooseCharge1",			   	      
	"IsProtonMC&IsPositive&IsMinimumBias_notrigg&IsLooseCharge1",
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1",			      
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1",			   	      
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1",
	Refill);

	AllRangesEfficiency * Trigger_D = new AllRangesEfficiency(finalHistos,"Trigger_D","Trigger",
	"IsDeutonMC&IsPositive&IsMinimumBias_notrigg&IsLooseCharge1",			      
	"IsDeutonMC&IsPositive&IsMinimumBias_notrigg&IsLooseCharge1",			   	      
	"IsDeutonMC&IsPositive&IsMinimumBias_notrigg&IsLooseCharge1",
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1",			      
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1",			   	      
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1",
	Refill);

	AllRangesEfficiency * MinimBias_P = new AllRangesEfficiency(finalHistos,"MinimBias_P","MinimBias",
	"IsProtonMC&IsPhysTrig",			      
	"IsProtonMC&IsPhysTrig",			   	      
	"IsProtonMC&IsPhysTrig",
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1",			      
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1",			   	      
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1",
	Refill);

	AllRangesEfficiency * MinimBias_D = new AllRangesEfficiency(finalHistos,"MinimBias_D","MinimBias",
	"IsDeutonMC&IsPhysTrig",			      
	"IsDeutonMC&IsPhysTrig",			   	      
	"IsDeutonMC&IsPhysTrig",
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1",			      
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1",			   	      
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1",
	Refill);

	AllRangesEfficiency * Cleaning_P = new AllRangesEfficiency(finalHistos,"Cleaning_P","Cleaning",
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1",			      
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1",			   	      
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1",
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",			      
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",			   	      
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",
	Refill);

	AllRangesEfficiency * Cleaning_D = new AllRangesEfficiency(finalHistos,"Cleaning_D","Cleaning",
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1",			      
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1",			   	      
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1",
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",			      
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",			   	      
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",
	Refill);

	AllRangesEfficiency * RICH_P = new AllRangesEfficiency(finalHistos,"RICH_P","RICH",
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",			      
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",			   	      
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",			      
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF",			   	      
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl",
	Refill);

	AllRangesEfficiency * RICH_D = new AllRangesEfficiency(finalHistos,"RICH_D","RICH",
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",			      
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",			   	      
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",			      
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF",			   	      
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl",
	Refill);

	AllRangesEfficiency * RICHQual_P = new AllRangesEfficiency(finalHistos,"RICHQual_P","RICHQual",
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",			      
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF",			   	      
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl",
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",			      
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut",			   	      
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut",
	Refill);

	AllRangesEfficiency * RICHQual_D = new AllRangesEfficiency(finalHistos,"RICHQual_D","RICHQual",
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",			      
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF",			   	      
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl",
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",			      
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut",			   	      
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut",
	Refill);

	AllRangesEfficiency * GoldenTOF_P = new AllRangesEfficiency(finalHistos,"GoldenTOF_P","GoldenTOF",
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",			      
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF",			   	      
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl",
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsGoodTime",			      
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut",			   	      
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut",
	Refill);

	AllRangesEfficiency * GoldenTOF_D = new AllRangesEfficiency(finalHistos,"GoldenTOF_D","GoldenTOF",
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning",			      
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF",			   	      
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl",
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsGoodTime",			      
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromNaF&RICHBDTCut",			   	      
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsFromAgl&RICHBDTCut",
	Refill);

	AllRangesEfficiency * Fragmentation_P = new AllRangesEfficiency(finalHistos,"Fragmentation_P","Fragmentation",
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1",			      
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1",			   	      
	"IsProtonMC&IsPositive&IsMinimumBias&IsLooseCharge1",
	"IsPurePMC&IsPositive&IsMinimumBias&IsLooseCharge1",			      
	"IsPurePMC&IsPositive&IsMinimumBias&IsLooseCharge1",			   	      
	"IsPurePMC&IsPositive&IsMinimumBias&IsLooseCharge1",
	Refill);

	AllRangesEfficiency * Fragmentation_D = new AllRangesEfficiency(finalHistos,"Fragmentation_D","Fragmentation",
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1",			      
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1",			   	      
	"IsDeutonMC&IsPositive&IsMinimumBias&IsLooseCharge1",
	"IsPureDMC&IsPositive&IsMinimumBias&IsLooseCharge1",			      
	"IsPureDMC&IsPositive&IsMinimumBias&IsLooseCharge1",			   	      
	"IsPureDMC&IsPositive&IsMinimumBias&IsLooseCharge1",
	Refill);

	RigBinBaselineEff->SetNotWeightedMC();
	RigBinQualEff->SetNotWeightedMC();
	FullSetTOT_P->SetNotWeightedMC();	
	FullSetTOT_D->SetNotWeightedMC();	
	Baseline_P->SetNotWeightedMC();	
	Baseline_D->SetNotWeightedMC();	


	ParallelFiller<Tool *> Filler;

	Filler.AddObject2beFilled(RigBinBaselineEff,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(RigBinQualEff,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(Cascade0,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(Cascade1,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(Cascade2,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(Cascade3,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(Cascade4,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(Cascade5,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(Cascade6,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(Cascade7,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(Cascade8,GetGenMomentum,GetGenMomentum);
	

	Filler.AddObject2beFilled(FullSetTOT_P,GetBetaGen,GetBetaGen);
	Filler.AddObject2beFilled(FullSetTOT_D,GetBetaGen,GetBetaGen);
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
	Filler.ReinitializeAll(Refill);

	//main loop
	Filler.LoopOnMC(DBarReader(chainMC, true ),vars);

	RigBinBaselineEff->Save(finalHistos);
	RigBinQualEff->Save(finalHistos);
	FullSetTOT_P 	->Save(finalHistos);
	FullSetTOT_D 	->Save(finalHistos);
	Baseline_P 	->Save(finalHistos);
	Baseline_D 	->Save(finalHistos);
	Cascade0        ->Save(finalHistos);
	Cascade1  	->Save(finalHistos);
	Cascade2  	->Save(finalHistos);
	Cascade3 	->Save(finalHistos);
	Cascade4 	->Save(finalHistos);
	Cascade5 	->Save(finalHistos);
	Cascade6 	->Save(finalHistos);
	Cascade7 	->Save(finalHistos);
	Cascade8 	->Save(finalHistos);
	Trigger_P 	->Save(finalHistos);
	Trigger_D 	->Save(finalHistos);
	MinimBias_P 	->Save(finalHistos);
	MinimBias_D 	->Save(finalHistos);
	Cleaning_P 	->Save(finalHistos);
	Cleaning_D 	->Save(finalHistos);
	RICH_P 	->Save(finalHistos);
	RICH_D 	->Save(finalHistos);
	RICHQual_P 	->Save(finalHistos);
	RICHQual_D 	->Save(finalHistos);
	GoldenTOF_P	->Save(finalHistos);
	GoldenTOF_D	->Save(finalHistos);
	Trigger_P_PID 	->Save(finalHistos);
	Trigger_D_PID 	->Save(finalHistos);
	MinimBias_P_PID 	->Save(finalHistos);
	MinimBias_D_PID 	->Save(finalHistos);
	Cleaning_P_PID 	->Save(finalHistos);
	Cleaning_D_PID 	->Save(finalHistos);
	RICH_P_PID 	->Save(finalHistos);
	RICH_D_PID 	->Save(finalHistos);
	RICHQual_P_PID 	->Save(finalHistos);
	RICHQual_D_PID 	->Save(finalHistos);
	GoldenTOF_P_PID	->Save(finalHistos);
	GoldenTOF_D_PID	->Save(finalHistos);
	Fragmentation_P	->Save(finalHistos);
	Fragmentation_D ->Save(finalHistos);


	RigBinBaselineEff->Eval_Efficiency();
	RigBinQualEff->Eval_Efficiency();
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

	RigBinBaselineEff->SaveResults(finalResults);
	RigBinQualEff->SaveResults(finalResults);
	FullSetTOT_P 	->SaveResults(finalResults);
	FullSetTOT_D 	->SaveResults(finalResults);
	Baseline_P 	->SaveResults(finalResults);
	Baseline_D 	->SaveResults(finalResults);
	Cascade1        ->SaveResults(finalResults);
	Cascade1  	->SaveResults(finalResults);
	Cascade2  	->SaveResults(finalResults);
	Cascade3 	->SaveResults(finalResults);
	Cascade4 	->SaveResults(finalResults);
	Cascade5 	->SaveResults(finalResults);
	Cascade6 	->SaveResults(finalResults);
	Cascade7 	->SaveResults(finalResults);
	Cascade8 	->SaveResults(finalResults);
	Trigger_P 	->SaveResults(finalResults);
	Trigger_D 	->SaveResults(finalResults);
	MinimBias_P 	->SaveResults(finalResults);
	MinimBias_D 	->SaveResults(finalResults);
	Cleaning_P 	->SaveResults(finalResults);
	Cleaning_D 	->SaveResults(finalResults);
	RICH_P 	->SaveResults(finalResults);
	RICH_D 	->SaveResults(finalResults);
	RICHQual_P 	->SaveResults(finalResults);
	RICHQual_D 	->SaveResults(finalResults);
	GoldenTOF_P	->SaveResults(finalResults);
	GoldenTOF_D     ->SaveResults(finalResults);
	Trigger_P_PID 	->SaveResults(finalResults);
	Trigger_D_PID 	->SaveResults(finalResults);
	MinimBias_P_PID 	->SaveResults(finalResults);
	MinimBias_D_PID 	->SaveResults(finalResults);
	Cleaning_P_PID 	->SaveResults(finalResults);
	Cleaning_D_PID 	->SaveResults(finalResults);
	RICH_P_PID 	->SaveResults(finalResults);
	RICH_D_PID 	->SaveResults(finalResults);
	RICHQual_P_PID 	->SaveResults(finalResults);
	RICHQual_D_PID 	->SaveResults(finalResults);
	GoldenTOF_P_PID	->SaveResults(finalResults);
	GoldenTOF_D_PID ->SaveResults(finalResults);
	Fragmentation_P->SaveResults(finalResults);
	Fragmentation_D->SaveResults(finalResults);


	return 0;
}


	

