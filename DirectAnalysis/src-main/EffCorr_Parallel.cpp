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


int main(int argc, char * argv[])
{

	TH1::SetDefaultSumw2();
        cout<<"****************************** FILES OPENING ***************************************"<<endl;

	string INPUT1(argv[1]);
	string INPUT2(argv[2]);
	string OUTPUT(argv[3]);

	string refill="";
	if(argc > 4 ) 	refill = argv[4];	

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

	EffCorr * HEPPresEffCorr = new EffCorr(finalHistos,"HEPPresEffCorr","HEPPresEffCorr",PRB,"IsPositive&IsMinimumBias&IsLooseCharge1","IsPositive&IsMinimumBias&IsLooseCharge1&IsGolden","","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * HEPQualEffCorr = new EffCorr(finalHistos,"HEPQualEffCorr","HEPQualEffCorr",PRB,"IsPositive&IsMinimumBias&IsLooseCharge1","IsPositive&IsMinimumBias&IsLooseCharge1&IsGolden","","IsPurePMC","IsPureDMC","IsDeutonMC");


	std::string before;
        std::string after;
        before = "IsPositive&IsMinimumBias_notrigg&IsLooseCharge1";
        after  = "IsPositive&IsMinimumBias&IsLooseCharge1";
	EffCorr * TriggerEffCorr_HE = new EffCorr(finalHistos,"TriggerEffCorr_HE","Trigger Eff. Corr",PRB,before,after,"IsPrimary",    "IsProtonMC"    ,"IsPureDMC","IsDeutonMC"); 
	EffCorr * TriggerEffCorr_TOF = new EffCorr(finalHistos,"TriggerEffCorr_TOF","Trigger Eff. Corr",ToFPB,before,after,"IsPrimary","IsProtonMC","IsPureDMC","IsDeutonMC");
	EffCorr * TriggerEffCorr_NaF = new EffCorr(finalHistos,"TriggerEffCorr_NaF","Trigger Eff. Corr",NaFPB,before,after,"IsPrimary","IsProtonMC","IsPureDMC","IsDeutonMC");
	EffCorr * TriggerEffCorr_Agl = new EffCorr(finalHistos,"TriggerEffCorr_Agl","Trigger Eff. Corr",AglPB,before,after,"IsPrimary","IsProtonMC","IsPureDMC","IsDeutonMC");

	before;
        after;
        before = "IsPositive&IsMinimumBias_notrack&IsLooseCharge1&QualChargeCut_notrack";
        after  = "IsPositive&IsMinimumBias&IsLooseCharge1&QualChargeCut_notrack"; 
	EffCorr * TrackerEffCorr_TOF = new EffCorr(finalHistos,"TrackerEffCorr_TOF","Tracker Eff. Corr",ToFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * TrackerEffCorr_NaF = new EffCorr(finalHistos,"TrackerEffCorr_NaF","Tracker Eff. Corr",NaFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * TrackerEffCorr_Agl = new EffCorr(finalHistos,"TrackerEffCorr_Agl","Tracker Eff. Corr",AglPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");

	before = "IsPositive&IsMinimumBias&IsLooseCharge1";
        after  = "IsPositive&IsMinimumBias&IsLooseCharge1&IsGoodChi2"; 
	EffCorr * GoodChi_TOF = new EffCorr(finalHistos,"GoodChiEffCorr_TOF","GoodChi Eff. Corr",ToFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodChi_NaF = new EffCorr(finalHistos,"GoodChiEffCorr_NaF","GoodChi Eff. Corr",NaFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodChi_Agl = new EffCorr(finalHistos,"GoodChiEffCorr_Agl","GoodChi Eff. Corr",AglPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");

	before = "IsPositive&IsMinimumBias&IsLooseCharge1";
        after  = "IsPositive&IsMinimumBias&IsLooseCharge1&IsCharge1UTOF"; 
	EffCorr * GoodUtof_TOF = new EffCorr(finalHistos,"GoodUtofEffCorr_TOF","GoodUtof Eff. Corr",ToFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodUtof_NaF = new EffCorr(finalHistos,"GoodUtofEffCorr_NaF","GoodUtof Eff. Corr",NaFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodUtof_Agl = new EffCorr(finalHistos,"GoodUtofEffCorr_Agl","GoodUtof Eff. Corr",AglPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");

	before = "IsPositive&IsMinimumBias&IsLooseCharge1";
        after  = "IsPositive&IsMinimumBias&IsLooseCharge1&IsCharge1LTOF";
        EffCorr * GoodLtof_TOF = new EffCorr(finalHistos,"GoodLTOFEffCorr_TOF","GoodQTrack Eff. Corr",ToFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
        EffCorr * GoodLtof_NaF = new EffCorr(finalHistos,"GoodLTOFEffCorr_NaF","GoodQTrack Eff. Corr",NaFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
        EffCorr * GoodLtof_Agl = new EffCorr(finalHistos,"GoodLTOFEffCorr_Agl","GoodQTrack Eff. Corr",AglPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");


	before = "IsPositive&IsMinimumBias&IsLooseCharge1";
        after  = "IsPositive&IsMinimumBias&IsLooseCharge1&Is1TrTrack"; 
	EffCorr * Good1Track_TOF = new EffCorr(finalHistos,"Good1TrackEffCorr_TOF","Good1Track Eff. Corr",ToFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * Good1Track_NaF = new EffCorr(finalHistos,"Good1TackEffCorr_NaF", "Good1Track Eff. Corr",NaFPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * Good1Track_Agl = new EffCorr(finalHistos,"Good1TrackEffCorr_Agl","Good1Track Eff. Corr",AglPB,before,after,"IsPrimary","IsPurePMC","IsPureDMC","IsDeutonMC");

	before = "IsPositive&IsMinimumBias&IsLooseCharge1";
        after  = "IsPositive&IsMinimumBias&IsLooseCharge1&IsCharge1Track"; 
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
//	Filler2.AddObject2beFilled(HEPPresEffCorr,GetRigidity,GetRigidity);
//	Filler2.AddObject2beFilled(HEPQualEffCorr,GetRigidity,GetRigidity);

	Filler2.AddObject2beFilled(TriggerEffCorr_HE,GetRigidity,GetRigidity);
	Filler2.AddObject2beFilled(TriggerEffCorr_TOF,GetRigidity,GetRigidity);
	Filler2.AddObject2beFilled(TriggerEffCorr_NaF,GetRigidity,GetRigidity);	
	Filler2.AddObject2beFilled(TriggerEffCorr_Agl,GetRigidity,GetRigidity);
/*	
	Filler2.AddObject2beFilled(TrackerEffCorr_TOF,GetRigidity,GetRigidity);
	Filler2.AddObject2beFilled(TrackerEffCorr_NaF,GetRigidity,GetRigidity);	
	Filler2.AddObject2beFilled(TrackerEffCorr_Agl,GetRigidity,GetRigidity);
	Filler2.AddObject2beFilled(GoodChi_TOF,GetRigidity,GetRigidity);	
	Filler2.AddObject2beFilled(GoodChi_NaF,GetRigidity,GetRigidity);
	Filler2.AddObject2beFilled(GoodChi_Agl,GetRigidity,GetRigidity);	
	Filler2.AddObject2beFilled(GoodUtof_TOF,GetRigidity,GetRigidity);
	Filler2.AddObject2beFilled(GoodUtof_NaF,GetRigidity,GetRigidity);	
	Filler2.AddObject2beFilled(GoodUtof_Agl,GetRigidity,GetRigidity);
	Filler2.AddObject2beFilled(GoodLtof_TOF,GetRigidity,GetRigidity);	
	Filler2.AddObject2beFilled(GoodLtof_NaF,GetRigidity,GetRigidity);
	Filler2.AddObject2beFilled(GoodLtof_Agl,GetRigidity,GetRigidity);	
	Filler2.AddObject2beFilled(GoodQTrack_TOF,GetRigidity,GetRigidity);
	Filler2.AddObject2beFilled(GoodQTrack_NaF,GetRigidity,GetRigidity);	
	Filler2.AddObject2beFilled(GoodQTrack_Agl,GetRigidity,GetRigidity);
	Filler2.AddObject2beFilled(Good1Track_TOF,GetRigidity,GetRigidity);
	Filler2.AddObject2beFilled(Good1Track_NaF,GetRigidity,GetRigidity);	
	Filler2.AddObject2beFilled(Good1Track_Agl,GetRigidity,GetRigidity);
	Filler2.AddObject2beFilled(GoodTime_TOF,GetRigidity,GetRigidity);	
	Filler2.AddObject2beFilled(RICHEffCorr_NaF,GetRigidity,GetRigidity);
	Filler2.AddObject2beFilled(RICHEffCorr_Agl,GetRigidity,GetRigidity);	
	Filler2.AddObject2beFilled(RICHQualEffCorr_NaF,GetRigidity,GetRigidity);
	Filler2.AddObject2beFilled(RICHQualEffCorr_Agl,GetRigidity,GetRigidity);	
*/
	Filler2.ReinitializeAll(Refill);
	//main loops 2
	Filler2.LoopOnMC(DBarReader(chainMC, true ),vars);
	Filler2.LoopOnData(DBarReader(chainDT, false),vars);

	//saving
	HEPPresEffCorr -> Save(finalHistos);
	HEPQualEffCorr -> Save(finalHistos);
	TriggerEffCorr_HE  -> Save(finalHistos);
	TriggerEffCorr_TOF -> Save(finalHistos);
	TriggerEffCorr_NaF -> Save(finalHistos);
	TriggerEffCorr_Agl -> Save(finalHistos);
	TrackerEffCorr_TOF -> Save(finalHistos);
	TrackerEffCorr_NaF -> Save(finalHistos);
	TrackerEffCorr_Agl -> Save(finalHistos);
	GoodChi_TOF -> Save(finalHistos);
	GoodChi_NaF -> Save(finalHistos);
	GoodChi_Agl -> Save(finalHistos);
	GoodUtof_TOF -> Save(finalHistos);
	GoodUtof_NaF -> Save(finalHistos);
	GoodUtof_Agl -> Save(finalHistos);
	GoodLtof_TOF -> Save(finalHistos);
	GoodLtof_NaF -> Save(finalHistos);
	GoodLtof_Agl -> Save(finalHistos);
	GoodQTrack_TOF -> Save(finalHistos);
	GoodQTrack_NaF -> Save(finalHistos);
	GoodQTrack_Agl -> Save(finalHistos);
	Good1Track_TOF -> Save(finalHistos);
	Good1Track_NaF -> Save(finalHistos);
	Good1Track_Agl -> Save(finalHistos);
	GoodTime_TOF -> Save(finalHistos);
	RICHEffCorr_NaF -> Save(finalHistos);
        RICHEffCorr_Agl -> Save(finalHistos);
        RICHQualEffCorr_NaF -> Save(finalHistos);
        RICHQualEffCorr_Agl -> Save(finalHistos);

	//analysis
	HEPPresEffCorr -> Eval_Efficiencies();
	HEPQualEffCorr -> Eval_Efficiencies();
	TriggerEffCorr_HE  -> Eval_Efficiencies();
	TriggerEffCorr_TOF -> Eval_Efficiencies();
	TriggerEffCorr_NaF -> Eval_Efficiencies();
	TriggerEffCorr_Agl -> Eval_Efficiencies();
	TrackerEffCorr_TOF -> Eval_Efficiencies();
	TrackerEffCorr_NaF -> Eval_Efficiencies();
	TrackerEffCorr_Agl -> Eval_Efficiencies();
	GoodChi_TOF -> Eval_Efficiencies();
	GoodChi_NaF -> Eval_Efficiencies();
	GoodChi_Agl -> Eval_Efficiencies();
	GoodUtof_TOF -> Eval_Efficiencies();
	GoodUtof_NaF -> Eval_Efficiencies();
	GoodUtof_Agl -> Eval_Efficiencies();
	GoodLtof_TOF -> Eval_Efficiencies();
	GoodLtof_NaF -> Eval_Efficiencies();
	GoodLtof_Agl -> Eval_Efficiencies();
	GoodQTrack_TOF -> Eval_Efficiencies();
	GoodQTrack_NaF -> Eval_Efficiencies();
	GoodQTrack_Agl -> Eval_Efficiencies();
	Good1Track_TOF -> Eval_Efficiencies();
	Good1Track_NaF -> Eval_Efficiencies();
	Good1Track_Agl -> Eval_Efficiencies();
	GoodTime_TOF -> Eval_Efficiencies();
	RICHEffCorr_NaF -> Eval_Efficiencies();
        RICHEffCorr_Agl -> Eval_Efficiencies();
        RICHQualEffCorr_NaF -> Eval_Efficiencies();
        RICHQualEffCorr_Agl -> Eval_Efficiencies();


	HEPPresEffCorr -> Eval_Corrections();
	HEPQualEffCorr -> Eval_Corrections();
	TriggerEffCorr_HE  -> Eval_Corrections();
	TriggerEffCorr_TOF -> Eval_Corrections();
	TriggerEffCorr_NaF -> Eval_Corrections();
	TriggerEffCorr_Agl -> Eval_Corrections();
	TrackerEffCorr_TOF -> Eval_Corrections();
	TrackerEffCorr_NaF -> Eval_Corrections();
	TrackerEffCorr_Agl -> Eval_Corrections();
	GoodChi_TOF -> Eval_Corrections();
	GoodChi_NaF -> Eval_Corrections();
	GoodChi_Agl -> Eval_Corrections();
	GoodUtof_TOF -> Eval_Corrections();
	GoodUtof_NaF -> Eval_Corrections();
	GoodUtof_Agl -> Eval_Corrections();
	GoodLtof_TOF -> Eval_Corrections();
	GoodLtof_NaF -> Eval_Corrections();
	GoodLtof_Agl -> Eval_Corrections();
	GoodQTrack_TOF -> Eval_Corrections();
	GoodQTrack_NaF -> Eval_Corrections();
	GoodQTrack_Agl -> Eval_Corrections();
	Good1Track_TOF -> Eval_Corrections();
	Good1Track_NaF -> Eval_Corrections();
	Good1Track_Agl -> Eval_Corrections();
	GoodTime_TOF -> Eval_Corrections();
	RICHEffCorr_NaF -> Eval_Corrections();
        RICHEffCorr_Agl -> Eval_Corrections();
        RICHQualEffCorr_NaF -> Eval_Corrections();
        RICHQualEffCorr_Agl -> Eval_Corrections();


	HEPPresEffCorr -> SaveResults(finalResults);
	HEPQualEffCorr -> SaveResults(finalResults);
	TriggerEffCorr_HE  -> SaveResults(finalResults);
	TriggerEffCorr_TOF -> SaveResults(finalResults);
	TriggerEffCorr_NaF -> SaveResults(finalResults);
	TriggerEffCorr_Agl -> SaveResults(finalResults);
	TrackerEffCorr_TOF -> SaveResults(finalResults);
	TrackerEffCorr_NaF -> SaveResults(finalResults);
	TrackerEffCorr_Agl -> SaveResults(finalResults);
	GoodChi_TOF -> SaveResults(finalResults);
	GoodChi_NaF -> SaveResults(finalResults);
	GoodChi_Agl -> SaveResults(finalResults);
	GoodUtof_TOF -> SaveResults(finalResults);
	GoodUtof_NaF -> SaveResults(finalResults);
	GoodUtof_Agl -> SaveResults(finalResults);
	GoodLtof_TOF -> SaveResults(finalResults);
	GoodLtof_NaF -> SaveResults(finalResults);
	GoodLtof_Agl -> SaveResults(finalResults);
	GoodQTrack_TOF -> SaveResults(finalResults);
	GoodQTrack_NaF -> SaveResults(finalResults);
	GoodQTrack_Agl -> SaveResults(finalResults);
	Good1Track_TOF -> SaveResults(finalResults);
	Good1Track_NaF -> SaveResults(finalResults);
	Good1Track_Agl -> SaveResults(finalResults);
	GoodTime_TOF -> SaveResults(finalResults);
	RICHEffCorr_NaF -> SaveResults(finalResults);
        RICHEffCorr_Agl -> SaveResults(finalResults);
        RICHQualEffCorr_NaF -> SaveResults(finalResults);
        RICHQualEffCorr_Agl -> SaveResults(finalResults);

	return 0;
}


	

