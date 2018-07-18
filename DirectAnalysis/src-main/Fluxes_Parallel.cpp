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

#include "../include/Efficiency.h"
#include "../include/AllRangesEfficiency.h"
#include "../include/Flux.h"
#include "../include/EffCorr.h"
//#include "../include/EffCorrTemplate.h"


//MC parameters taken from /cvmfs/ams.cern.ch/Offline/AMSDataDir/DataManagement/DataSetsDesc

int main(int argc, char * argv[])
{


        cout<<"****************************** FILES OPENING ***************************************"<<endl;

        string INPUT1(argv[1]);
        string INPUT2(argv[2]);
        string OUTPUT(argv[3]);

	string refill="";
        if(argc > 4 )   refill = argv[4];

        bool Refill = false;
        if(refill!="") Refill=true;

        FileSaver finalHistos;
        finalHistos.setName(OUTPUT.c_str());
        bool checkfile = finalHistos.CheckFile();
	if(!checkfile) Refill=true;

	FileSaver finalResults;
        finalResults.setName((OUTPUT+"_Results").c_str());


//	TChain * chainDT = InputFileReader(INPUT1.c_str(),"parametri_geo");
   //	TChain * chainMC = InputFileReader(INPUT2.c_str(),"parametri_MC");
  	TChain * chainDT = InputFileReader(INPUT1.c_str(),"RTI");
   	TChain * chainMC = InputFileReader(INPUT2.c_str(),"Event");
 
	cout<<"****************************** BINS ***************************************"<<endl;
	SetUpUsualBinning();

	cout<<"****************************** VARIABLES ***************************************"<<endl;

        Variables * vars = new Variables;

	cout<<"****************************** EFFICIENCY CORR. READING  ******************************************"<<endl;
	std::string before;
        std::string after;
        before = "";
        after  = ""; 
	EffCorr * TrackerEffCorr_TOF = new EffCorr(finalResults,"TrackerEffCorr_TOF","Tracker Eff. Corr",ToFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * TrackerEffCorr_NaF = new EffCorr(finalResults,"TrackerEffCorr_NaF","Tracker Eff. Corr",NaFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * TrackerEffCorr_Agl = new EffCorr(finalResults,"TrackerEffCorr_Agl","Tracker Eff. Corr",AglPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodChi_TOF 		= new EffCorr(finalResults,"GoodChiEffCorr_TOF","GoodChi Eff. Corr",ToFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodChi_NaF 		= new EffCorr(finalResults,"GoodChiEffCorr_NaF","GoodChi Eff. Corr",NaFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodChi_Agl 		= new EffCorr(finalResults,"GoodChiEffCorr_Agl","GoodChi Eff. Corr",AglPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodUtof_TOF 		= new EffCorr(finalResults,"GoodUtofEffCorr_TOF","GoodUtof Eff. Corr",ToFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodUtof_NaF 		= new EffCorr(finalResults,"GoodUtofEffCorr_NaF","GoodUtof Eff. Corr",NaFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodUtof_Agl 		= new EffCorr(finalResults,"GoodUtofEffCorr_Agl","GoodUtof Eff. Corr",AglPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
        EffCorr * GoodLtof_TOF 		= new EffCorr(finalResults,"GoodLTOFEffCorr_TOF","GoodQTrack Eff. Corr",ToFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
        EffCorr * GoodLtof_NaF 		= new EffCorr(finalResults,"GoodLTOFEffCorr_NaF","GoodQTrack Eff. Corr",NaFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
        EffCorr * GoodLtof_Agl 		= new EffCorr(finalResults,"GoodLTOFEffCorr_Agl","GoodQTrack Eff. Corr",AglPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * Good1Track_TOF 	= new EffCorr(finalResults,"Good1TrackEffCorr_TOF","Good1Track Eff. Corr",ToFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * Good1Track_NaF 	= new EffCorr(finalResults,"Good1TackEffCorr_NaF", "Good1Track Eff. Corr",NaFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * Good1Track_Agl 	= new EffCorr(finalResults,"Good1TrackEffCorr_Agl","Good1Track Eff. Corr",AglPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodQTrack_TOF 	= new EffCorr(finalResults,"GoodQTrackEffCorr_TOF","GoodQTrack Eff. Corr",ToFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodQTrack_NaF 	= new EffCorr(finalResults,"GoodQTrackEffCorr_NaF","GoodQTrack Eff. Corr",NaFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodQTrack_Agl 	= new EffCorr(finalResults,"GoodQTrackEffCorr_Agl","GoodQTrack Eff. Corr",AglPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * GoodTime_TOF 		= new EffCorr(finalResults,"GoodTimeEffCorr_TOF","GoodTime Eff. Corr",ToFPB,before,after,"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * RICHEffCorr_NaF 	= new EffCorr(finalResults,"RICHCorrection_NaF","RICH Eff. Corr",NaFPB,before,(after+"&IsFromNaF").c_str(),"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * RICHEffCorr_Agl 	= new EffCorr(finalResults,"RICHCorrection_Agl","RICH Eff. Corr",AglPB,before,(after+"&IsFromAgl").c_str(),"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * RICHQualEffCorr_NaF 	= new EffCorr(finalResults,"RICHQualCorrection_NaF","RICH Qual Eff. Corr",NaFPB,(before+"&IsFromNaF").c_str(),(after+"&IsFromNaF&RICHBDTCut").c_str(),"","IsPurePMC","IsPureDMC","IsDeutonMC");
	EffCorr * RICHQualEffCorr_Agl 	= new EffCorr(finalResults,"RICHqualCorrection_Agl","RICH Qual. Eff. Corr",AglPB,(before+"&IsFromNaF").c_str(),(after+"&IsFromAgl&RICHBDTCut").c_str(),"","IsPurePMC","IsPureDMC","IsDeutonMC");

	cout<<"****************************** FLUXES EVALUATION ******************************************"<<endl;

	
	Flux * DFluxTOF = new Flux(finalHistos,finalResults, "DFluxTOF", "FullSetTOT_D_TOF","FullSetTOT","TOFfits/Fit Results/Primary Deuteron Counts","ExposureTOF","Gen. Acceptance",ToFDB);
	Flux * DFluxNaF = new Flux(finalHistos,finalResults, "DFluxNaF", "FullSetTOT_D_NaF","FullSetTOT","NaFfits/Fit Results/Primary Deuteron Counts","ExposureNaF","Gen. Acceptance",NaFDB);
	Flux * DFluxAgl = new Flux(finalHistos,finalResults, "DFluxAgl", "FullSetTOT_D_Agl","FullSetTOT","Aglfits/Fit Results/Primary Deuteron Counts","ExposureAgl","Gen. Acceptance",AglDB);

	Flux * HEPFlux  = new Flux(finalHistos,finalResults,"PFluxHE", "RigBinFullSetEff","RigBinFullSetEff","HEPCounts/HEPCounts/HEPCounts","HEExposure","Gen. Acceptance",PRB);

	Flux * PFluxTOF = new Flux(finalHistos,finalResults, "PFluxTOF", "FullSetTOT_P_TOF","FullSetTOT","TOFfits/Fit Results/Primary Proton Counts","ExposureTOF","Gen. Acceptance",ToFPB);
	Flux * PFluxNaF = new Flux(finalHistos,finalResults, "PFluxNaF", "FullSetTOT_P_NaF","FullSetTOT","NaFfits/Fit Results/Primary Proton Counts","ExposureNaF","Gen. Acceptance",NaFPB);
	Flux * PFluxAgl = new Flux(finalHistos,finalResults, "PFluxAgl", "FullSetTOT_P_Agl","FullSetTOT","Aglfits/Fit Results/Primary Proton Counts","ExposureAgl","Gen. Acceptance",AglPB);

	Flux * DummyDTOF = new Flux(finalHistos,finalResults, "DFluxTOF", "Baseline_D_TOF","Baseline","TOFfits/Fit Results/Primary Deuteron Counts","ExposureTOF","Gen. Acceptance",ToFDB);
	Flux * DummyDNaF = new Flux(finalHistos,finalResults, "DFluxNaF", "Baseline_D_NaF","Baseline","NaFfits/Fit Results/Primary Deuteron Counts","ExposureNaF","Gen. Acceptance",NaFDB);
	Flux * DummyDAgl = new Flux(finalHistos,finalResults, "DFluxAgl", "Baseline_D_Agl","Baseline","Aglfits/Fit Results/Primary Deuteron Counts","ExposureAgl","Gen. Acceptance",AglDB);
                                                                        
	Flux * DummyPTOF = new Flux(finalHistos,finalResults, "PFluxTOF", "Baseline_P_TOF","Baseline","TOFfits/Fit Results/Primary Proton Counts","ExposureTOF","Gen. Acceptance",ToFPB);
	Flux * DummyPNaF = new Flux(finalHistos,finalResults, "PFluxNaF", "Baseline_P_NaF","Baseline","NaFfits/Fit Results/Primary Proton Counts","ExposureNaF","Gen. Acceptance",NaFPB);
	Flux * DummyPAgl = new Flux(finalHistos,finalResults, "PFluxAgl", "Baseline_P_Agl","Baseline","Aglfits/Fit Results/Primary Proton Counts","ExposureAgl","Gen. Acceptance",AglPB);



	cout<<"********** EXPOSURE TIME & GEOM. ACCEPT. ********"<<endl;

	ParallelFiller<Flux *> Filler;
                Filler.AddObject2beFilled(DFluxTOF,GetBetaGen,GetBetaGen,"IsDeutonMC");
		Filler.AddObject2beFilled(DFluxNaF,GetBetaGen,GetBetaGen,"IsDeutonMC");
		Filler.AddObject2beFilled(DFluxAgl,GetBetaGen,GetBetaGen,"IsDeutonMC");
		Filler.AddObject2beFilled(PFluxTOF,GetBetaGen,GetBetaGen,"IsProtonMC");
		Filler.AddObject2beFilled(PFluxNaF,GetBetaGen,GetBetaGen,"IsProtonMC");
		Filler.AddObject2beFilled(PFluxAgl,GetBetaGen,GetBetaGen,"IsProtonMC");
		Filler.AddObject2beFilled(HEPFlux,GetGenMomentum,GetGenMomentum,"IsProtonMC");
		Filler.ReinitializeAll(Refill);
		//main loops
		Filler.ExposureTimeFilling_RTI(DBarReader(chainDT, false ),vars,finalHistos);
		Filler.LoopOnMCForGenAcceptance(DBarReader(chainMC, true ),vars,finalHistos);		

	cout<<"************* D FLUX ************"<<endl;

	DFluxTOF->Set_MCPar(0.5,200,0.009014888);	
	DFluxNaF->Set_MCPar(0.5,200,0.009014888);	
	DFluxAgl->Set_MCPar(0.5,200,0.009014888);	
	DummyDTOF->Set_MCPar(0.5,200,0.009014888);	
	DummyDNaF->Set_MCPar(0.5,200,0.009014888);	
	DummyDAgl->Set_MCPar(0.5,200,0.009014888);	


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
	
	DFluxTOF->SaveResults(finalResults);
	DFluxNaF->SaveResults(finalResults);
	DFluxAgl->SaveResults(finalResults);
	DummyDTOF->SaveResults(finalResults);
	DummyDNaF->SaveResults(finalResults);
	DummyDAgl->SaveResults(finalResults);


	cout<<"************* P FLUX ************"<<endl;

	HEPFlux ->Set_MCPar(0.5,100,0.00856558);	
	PFluxTOF->Set_MCPar(0.5,100,0.00856558);	
	PFluxNaF->Set_MCPar(0.5,100,0.00856558);	
	PFluxAgl->Set_MCPar(0.5,100,0.00856558);	
	DummyPTOF->Set_MCPar(0.5,100,0.00856558);	
	DummyPNaF->Set_MCPar(0.5,100,0.00856558);	
	DummyPAgl->Set_MCPar(0.5,100,0.00856558);	
	
	/*HEPFlux->Set_MCPar(0.5,200,0.009014888);
	PFluxTOF->Set_MCPar(0.5,200,0.009014888);	
	PFluxNaF->Set_MCPar(0.5,200,0.009014888);	
	PFluxAgl->Set_MCPar(0.5,200,0.009014888);	
	DummyPTOF->Set_MCPar(0.5,200,0.009014888);	
	DummyPNaF->Set_MCPar(0.5,200,0.009014888);	
	DummyPAgl->Set_MCPar(0.5,200,0.009014888);	
*/

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
	PFluxTOF-> Eval_Flux();
        PFluxNaF-> Eval_Flux();
        PFluxAgl-> Eval_Flux();
	DummyPTOF-> Eval_Flux();
        DummyPNaF-> Eval_Flux();
        DummyPAgl-> Eval_Flux();

	DummyPTOF->ChangeName("DummyPTOF");
	DummyPNaF->ChangeName("DummyPNaF");
	DummyPAgl->ChangeName("DummyPAgl");
	
	HEPFlux ->SaveResults(finalResults);
	PFluxTOF->SaveResults(finalResults);
	PFluxNaF->SaveResults(finalResults);
	PFluxAgl->SaveResults(finalResults);
	DummyPTOF->SaveResults(finalResults);
	DummyPNaF->SaveResults(finalResults);
	DummyPAgl->SaveResults(finalResults);


	cout<<"************* D/P RATIO ************"<<endl;


	TH1F * DPRatioTOF = DFluxTOF->Eval_FluxRatio(PFluxTOF,"DP ratio TOF");
	TH1F * DPRatioNaF = DFluxNaF->Eval_FluxRatio(PFluxNaF,"DP ratio NaF");
	TH1F * DPRatioAgl = DFluxAgl->Eval_FluxRatio(PFluxAgl,"DP ratio Agl");

	finalResults.Add(DPRatioTOF);
	finalResults.Add(DPRatioNaF);
	finalResults.Add(DPRatioAgl);

        finalResults.writeObjsInFolder("Fluxes/");	

	return 0;
}
