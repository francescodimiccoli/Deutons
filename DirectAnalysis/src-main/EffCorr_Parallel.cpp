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
#include "../include/Globals.h"
#include "../include/InputFileReader.h"

#include "../include/Variables.hpp"
#include "../include/Cuts.h"
#include "../include/ParallelFiller.h"
#include "../include/LatReweighter.h"

#include "../include/filesaver.h"

#include "../include/Efficiency.h"
#include "../include/EffCorr.h"
#include "../include/EffCorrTemplate.h"

void AnalyzeEffCorr(EffCorr * Correction, FileSaver  finalHistos, FileSaver  finalResults);

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
	cout<<"ecco"<<endl;
	
	bool Refill = false;
	if(refill!="") Refill=true;

	TChain * chain_RTI = InputFileReader(INPUT1.c_str(),"RTI");

	TChain * chainDT = InputFileReader(INPUT1.c_str(),"Event");
	TChain * chainMC = InputFileReader(INPUT2.c_str(),"Event");
	//TChain * chainDT = InputFileReader(INPUT1.c_str(),"template_stuff");
	//TChain * chainMC = InputFileReader(INPUT2.c_str(),"template_stuffMC");
	//TChain * chainDT = InputFileReader(INPUT1.c_str(),"parametri_geo");
	//TChain * chainMC = InputFileReader(INPUT2.c_str(),"parametri_MC");

	FileSaver LatWeights;
	LatWeights.setName("/afs/cern.ch/work/f/fdimicco/private/Deutons/DirectAnalysis/LatWeights/Weights.root");
	LatReweighter * weighter = new LatReweighter(LatWeights,"LatWeights");	


	FileSaver finalHistos;



	finalHistos.setName(OUTPUT.c_str());

	FileSaver finalResults;
	finalResults.setName((OUTPUT+"_Results").c_str());


	bool checkfile = finalHistos.CheckFile();

	cout<<"****************************** BINS ***************************************"<<endl;
	SetUpEffCorrBinning();
	

	cout<<"****************************** VARIABLES ***************************************"<<endl;

        Variables * vars = new Variables;

	cout<<"****************************** ANALYIS ******************************************"<<endl;

	BadEventSimulator * NaFBadEvSimulator= new BadEventSimulator("IsFromNaF",30,0.82,1);
        BadEventSimulator * AglBadEvSimulator= new BadEventSimulator("IsFromAgl",120,0.96,1);



	//Baseline efficiency corrections
	std::string before;
        std::string after;
        before = "IsPositive&IsMinimumBias_notrigg&IsLooseCharge1";
        after  = "IsPositive&IsMinimumBias&IsLooseCharge1";
	EffCorr * TriggerEffCorr_HE  = new EffCorr(finalHistos,"TriggerEffCorr_HE" ,"Trigger Eff. Corr",PRB  ,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC"); 
	EffCorr * TriggerEffCorr_TOF = new EffCorr(finalHistos,"TriggerEffCorr_TOF","Trigger Eff. Corr",ToFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * TriggerEffCorr_NaF = new EffCorr(finalHistos,"TriggerEffCorr_NaF","Trigger Eff. Corr",NaFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * TriggerEffCorr_Agl = new EffCorr(finalHistos,"TriggerEffCorr_Agl","Trigger Eff. Corr",AglPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");

	before = "IsPositive&IsPhysTrig&IsGoodTOFStandaloneQ1&IsL1HitNearExtrapol&IsCleanL1Hit";
        after  = "IsPositive&IsPhysTrig&IsGoodTOFStandaloneQ1&IsL1HitNearExtrapol&IsCleanL1Hit&L1LooseCharge1";
	EffCorr * L1PickUpEffCorr_HE = new EffCorr(finalHistos,"L1PickUpEffCorr_HE","L1PickUp Eff. Corr",PRB,before,after,"IsPrimary",    "IsPurePMC","IsPureDMC","IsDeutonMC"); 
	EffCorr * L1PickUpEffCorr_TOF = new EffCorr(finalHistos,"L1PickUpEffCorr_TOF","L1PickUp Eff. Corr",ToFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * L1PickUpEffCorr_NaF = new EffCorr(finalHistos,"L1PickUpEffCorr_NaF","L1PickUp Eff. Corr",NaFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * L1PickUpEffCorr_Agl = new EffCorr(finalHistos,"L1PickUpEffCorr_Agl","L1PickUp Eff. Corr",AglPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");

	before = "IsPositive&IsMinimumBias_notof&IsLooseCharge1";
        after  = "IsPositive&IsMinimumBias&IsLooseCharge1";
	EffCorr * MinTOFEffCorr_HE  = new EffCorr(finalHistos,"MinTOFEffCorr_HE" ,"Min TOF Eff. Corr",PRB,before,after,  "IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC"); 
	EffCorr * MinTOFEffCorr_TOF = new EffCorr(finalHistos,"MinTOFEffCorr_TOF","Min TOF Eff. Corr",ToFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * MinTOFEffCorr_NaF = new EffCorr(finalHistos,"MinTOFEffCorr_NaF","Min TOF Eff. Corr",NaFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * MinTOFEffCorr_Agl = new EffCorr(finalHistos,"MinTOFEffCorr_Agl","Min TOF Eff. Corr",AglPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");

        before = "IsMinimumBias_notrack&IsGoodTOFStandaloneQ1&IsExtrapolInsideL8";
        after  = "IsMinimumBias&IsGoodTOFStandaloneQ1&IsExtrapolInsideL8"; 
	EffCorr * TrackerEffCorr_HE = new EffCorr(finalHistos,"TrackerEffCorr_HE","Tracker Eff. Corr",PRB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");

	// Good Z=1 and Golden Efficiency corrections
	before = "IsPositive&IsMinimumBias&IsLooseCharge1";
        after  = "IsPositive&IsMinimumBias&IsLooseCharge1&IsGoodChi2"; 
	EffCorr * GoodChi_HE  = new EffCorr(finalHistos,"GoodChiEffCorr_HE","GoodChi Eff. Corr",PRB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodChi_TOF = new EffCorr(finalHistos,"GoodChiEffCorr_TOF","GoodChi Eff. Corr",ToFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodChi_NaF = new EffCorr(finalHistos,"GoodChiEffCorr_NaF","GoodChi Eff. Corr",NaFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodChi_Agl = new EffCorr(finalHistos,"GoodChiEffCorr_Agl","GoodChi Eff. Corr",AglPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");

	before = "IsPositive&IsMinimumBias&IsLooseCharge1";
        after  = "IsPositive&IsMinimumBias&IsLooseCharge1&IsCharge1UTOF"; 
	EffCorr * GoodUtof_HE  = new EffCorr(finalHistos,"GoodUtofEffCorr_HE","GoodUtof Eff. Corr",PRB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodUtof_TOF = new EffCorr(finalHistos,"GoodUtofEffCorr_TOF","GoodUtof Eff. Corr",ToFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodUtof_NaF = new EffCorr(finalHistos,"GoodUtofEffCorr_NaF","GoodUtof Eff. Corr",NaFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodUtof_Agl = new EffCorr(finalHistos,"GoodUtofEffCorr_Agl","GoodUtof Eff. Corr",AglPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");

	before = "IsPositive&IsMinimumBias&IsLooseCharge1";
        after  = "IsPositive&IsMinimumBias&IsLooseCharge1&IsCharge1LTOF";
        EffCorr * GoodLtof_HE  = new EffCorr(finalHistos,"GoodLTOFEffCorr_HE","GoodQTrack Eff. Corr",ToFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
        EffCorr * GoodLtof_TOF = new EffCorr(finalHistos,"GoodLTOFEffCorr_TOF","GoodQTrack Eff. Corr",ToFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
        EffCorr * GoodLtof_NaF = new EffCorr(finalHistos,"GoodLTOFEffCorr_NaF","GoodQTrack Eff. Corr",NaFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
        EffCorr * GoodLtof_Agl = new EffCorr(finalHistos,"GoodLTOFEffCorr_Agl","GoodQTrack Eff. Corr",AglPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");

	before = "IsPositive&IsMinimumBias&IsLooseCharge1";
        after  = "IsPositive&IsMinimumBias&IsLooseCharge1&Is1TrTrack"; 
	EffCorr * Good1Track_HE  = new EffCorr(finalHistos,"Good1TrackEffCorr_HE","Good1Track Eff. Corr",ToFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * Good1Track_TOF = new EffCorr(finalHistos,"Good1TrackEffCorr_TOF","Good1Track Eff. Corr",ToFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * Good1Track_NaF = new EffCorr(finalHistos,"Good1TackEffCorr_NaF", "Good1Track Eff. Corr",NaFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * Good1Track_Agl = new EffCorr(finalHistos,"Good1TrackEffCorr_Agl","Good1Track Eff. Corr",AglPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");

	before = "IsPositive&IsMinimumBias&IsLooseCharge1";
        after  = "IsPositive&IsMinimumBias&IsLooseCharge1&IsCharge1Track"; 
	EffCorr * GoodQTrack_HE  = new EffCorr(finalHistos,"GoodQTrackEffCorr_HE","GoodQTrack Eff. Corr",ToFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodQTrack_TOF = new EffCorr(finalHistos,"GoodQTrackEffCorr_TOF","GoodQTrack Eff. Corr",ToFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodQTrack_NaF = new EffCorr(finalHistos,"GoodQTrackEffCorr_NaF","GoodQTrack Eff. Corr",NaFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodQTrack_Agl = new EffCorr(finalHistos,"GoodQTrackEffCorr_Agl","GoodQTrack Eff. Corr",AglPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");


	before = "IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning";
        after  = "IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning&IsGoodTime"; 
	EffCorr * GoodTime_TOF = new EffCorr(finalHistos,"GoodTimeEffCorr_TOF","GoodTime Eff. Corr",ToFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	
	before = "IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning";
        after  = "IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning";
	EffCorr * RICHEffCorr_NaF = new EffCorr(finalHistos,"RICHCorrection_NaF","RICH Eff. Corr",NaFPB,before,(after+"&IsFromNaF").c_str(),"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * RICHEffCorr_Agl = new EffCorr(finalHistos,"RICHCorrection_Agl","RICH Eff. Corr",AglPB,before,(after+"&IsFromAgl").c_str(),"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");


	before = "IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning";
        after  = "IsPositive&IsMinimumBias&IsLooseCharge1&IsCleaning";
	EffCorr * RICHQualEffCorr_NaF = new EffCorr(finalHistos,"RICHQualCorrection_NaF","RICH Qual Eff. Corr",NaFPB,(before+"&IsFromNaF").c_str(),(after+"&IsFromNaF&RICHBDTCut").c_str(),"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * RICHQualEffCorr_Agl = new EffCorr(finalHistos,"RICHqualCorrection_Agl","RICH Qual. Eff. Corr",AglPB,(before+"&IsFromAgl").c_str(),(after+"&IsFromAgl&RICHBDTCut").c_str(),"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");


	ParallelFiller<EffCorr*> Filler2;

	Filler2.AddObject2beFilled(TriggerEffCorr_HE,GetRigidity,GetRigidity);
	Filler2.AddObject2beFilled(TriggerEffCorr_TOF,GetRigidity,GetRigidity);
	Filler2.AddObject2beFilled(TriggerEffCorr_NaF,GetRigidity,GetRigidity);	
	Filler2.AddObject2beFilled(TriggerEffCorr_Agl,GetRigidity,GetRigidity);

	Filler2.AddObject2beFilled(L1PickUpEffCorr_HE,GetRigidity,GetRigidity);
	Filler2.AddObject2beFilled(L1PickUpEffCorr_TOF,GetRigidity,GetRigidity);
	Filler2.AddObject2beFilled(L1PickUpEffCorr_NaF,GetRigidity,GetRigidity);	
	Filler2.AddObject2beFilled(L1PickUpEffCorr_Agl,GetRigidity,GetRigidity);

	Filler2.AddObject2beFilled(MinTOFEffCorr_HE , GetRigidity,GetRigidity);
	Filler2.AddObject2beFilled(MinTOFEffCorr_TOF,GetRigidity,GetRigidity);
	Filler2.AddObject2beFilled(MinTOFEffCorr_NaF,GetRigidity,GetRigidity);	
	Filler2.AddObject2beFilled(MinTOFEffCorr_Agl,GetRigidity,GetRigidity);
	
	Filler2.AddObject2beFilled(TrackerEffCorr_HE,GetMomentumProxy,GetMomentumProxy);

	Filler2.AddObject2beFilled(GoodChi_HE,GetRigidity,GetRigidity);	
	Filler2.AddObject2beFilled(GoodChi_TOF,GetRigidity,GetRigidity);	
	Filler2.AddObject2beFilled(GoodChi_NaF,GetRigidity,GetRigidity);
	Filler2.AddObject2beFilled(GoodChi_Agl,GetRigidity,GetRigidity);
	
	Filler2.AddObject2beFilled(GoodUtof_HE,GetRigidity,GetRigidity);
	Filler2.AddObject2beFilled(GoodUtof_TOF,GetRigidity,GetRigidity);
	Filler2.AddObject2beFilled(GoodUtof_NaF,GetRigidity,GetRigidity);	
	Filler2.AddObject2beFilled(GoodUtof_Agl,GetRigidity,GetRigidity);

	Filler2.AddObject2beFilled(GoodLtof_HE,GetRigidity,GetRigidity);	
	Filler2.AddObject2beFilled(GoodLtof_TOF,GetRigidity,GetRigidity);	
	Filler2.AddObject2beFilled(GoodLtof_NaF,GetRigidity,GetRigidity);
	Filler2.AddObject2beFilled(GoodLtof_Agl,GetRigidity,GetRigidity);	

	Filler2.AddObject2beFilled(GoodQTrack_HE,GetRigidity,GetRigidity);
	Filler2.AddObject2beFilled(GoodQTrack_TOF,GetRigidity,GetRigidity);
	Filler2.AddObject2beFilled(GoodQTrack_NaF,GetRigidity,GetRigidity);	
	Filler2.AddObject2beFilled(GoodQTrack_Agl,GetRigidity,GetRigidity);

	Filler2.AddObject2beFilled(Good1Track_HE,GetRigidity,GetRigidity);
	Filler2.AddObject2beFilled(Good1Track_TOF,GetRigidity,GetRigidity);
	Filler2.AddObject2beFilled(Good1Track_NaF,GetRigidity,GetRigidity);	
	Filler2.AddObject2beFilled(Good1Track_Agl,GetRigidity,GetRigidity);

	Filler2.AddObject2beFilled(GoodTime_TOF,GetRigidity,GetRigidity);	
	Filler2.AddObject2beFilled(RICHEffCorr_NaF,GetRigidity,GetRigidity);
	Filler2.AddObject2beFilled(RICHEffCorr_Agl,GetRigidity,GetRigidity);	
	Filler2.AddObject2beFilled(RICHQualEffCorr_NaF,GetRigidity,GetRigidity);
	Filler2.AddObject2beFilled(RICHQualEffCorr_Agl,GetRigidity,GetRigidity);	

	Filler2.ReinitializeAll(Refill);
	//main loops 2
	Filler2.LoopOnMC(DBarReader(chainMC, true ),vars);
	Filler2.LoopOnData(DBarReader(chainDT, false,chain_RTI),vars);

	//saving
	AnalyzeEffCorr(	TriggerEffCorr_HE  , finalHistos, finalResults);
	AnalyzeEffCorr(	TriggerEffCorr_TOF , finalHistos, finalResults);
	AnalyzeEffCorr(	TriggerEffCorr_NaF , finalHistos, finalResults);
	AnalyzeEffCorr(	TriggerEffCorr_Agl , finalHistos, finalResults);

	AnalyzeEffCorr(	L1PickUpEffCorr_HE , finalHistos, finalResults);
	AnalyzeEffCorr(	L1PickUpEffCorr_TOF, finalHistos, finalResults);
	AnalyzeEffCorr(	L1PickUpEffCorr_NaF, finalHistos, finalResults);
	AnalyzeEffCorr(	L1PickUpEffCorr_Agl, finalHistos, finalResults);

	AnalyzeEffCorr(	MinTOFEffCorr_HE  , finalHistos, finalResults);
	AnalyzeEffCorr(	MinTOFEffCorr_TOF , finalHistos, finalResults);
	AnalyzeEffCorr(	MinTOFEffCorr_NaF , finalHistos, finalResults);
	AnalyzeEffCorr(	MinTOFEffCorr_Agl , finalHistos, finalResults);

	AnalyzeEffCorr(	TrackerEffCorr_HE , finalHistos, finalResults);

	AnalyzeEffCorr(	GoodChi_HE  , finalHistos, finalResults);  
	AnalyzeEffCorr(	GoodChi_TOF , finalHistos, finalResults);  
	AnalyzeEffCorr(	GoodChi_NaF , finalHistos, finalResults);    
	AnalyzeEffCorr(	GoodChi_Agl , finalHistos, finalResults);

	AnalyzeEffCorr(	GoodUtof_HE , finalHistos, finalResults);
	AnalyzeEffCorr(	GoodUtof_TOF , finalHistos, finalResults);
	AnalyzeEffCorr(	GoodUtof_NaF , finalHistos, finalResults);
	AnalyzeEffCorr(	GoodUtof_Agl , finalHistos, finalResults);

	AnalyzeEffCorr(	GoodLtof_HE  , finalHistos, finalResults);
	AnalyzeEffCorr(	GoodLtof_TOF , finalHistos, finalResults);
	AnalyzeEffCorr(	GoodLtof_NaF , finalHistos, finalResults);
	AnalyzeEffCorr(	GoodLtof_Agl , finalHistos, finalResults);

	AnalyzeEffCorr(	GoodQTrack_HE , finalHistos, finalResults);
	AnalyzeEffCorr(	GoodQTrack_TOF , finalHistos, finalResults);
	AnalyzeEffCorr(	GoodQTrack_NaF , finalHistos, finalResults);
	AnalyzeEffCorr(	GoodQTrack_Agl , finalHistos, finalResults);

	AnalyzeEffCorr(	Good1Track_HE , finalHistos, finalResults);
	AnalyzeEffCorr(	Good1Track_TOF , finalHistos, finalResults);
	AnalyzeEffCorr(	Good1Track_NaF , finalHistos, finalResults);
	AnalyzeEffCorr(	Good1Track_Agl , finalHistos, finalResults);

	AnalyzeEffCorr(	GoodTime_TOF , finalHistos, finalResults);
	AnalyzeEffCorr(	RICHEffCorr_NaF , finalHistos, finalResults);
        AnalyzeEffCorr(	RICHEffCorr_Agl , finalHistos, finalResults);
        AnalyzeEffCorr(	RICHQualEffCorr_NaF , finalHistos, finalResults);
        AnalyzeEffCorr(	RICHQualEffCorr_Agl , finalHistos, finalResults);


	return 0;
}


void AnalyzeEffCorr(EffCorr * Correction, FileSaver  finalHistos, FileSaver  finalResults){
	Correction   -> Save(finalHistos);
	Correction   -> Eval_Efficiencies();
	Correction   -> Eval_Corrections();
	Correction   -> SaveResults(finalResults);
}	


