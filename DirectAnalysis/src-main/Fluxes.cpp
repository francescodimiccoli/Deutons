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

#include "../include/GlobalBinning.h"

#include "../include/Commonglobals.cpp"

#include "../include/Variables.hpp"
#include "../include/Cuts.h"


#include "../include/filesaver.h"

#include "../include/Efficiency.h"
#include "../include/AllRangesEfficiency.h"
#include "../include/Flux.h"
#include "../include/EffCorr.h"
//#include "../include/EffCorrTemplate.h"


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

	FileSaver finalResults;
        finalResults.setName((OUTPUT+"_Results").c_str());


        TFile *fileDT =TFile::Open(INPUT1.c_str());
        TFile *fileMC =TFile::Open(INPUT2.c_str());

        TTree *treeMC = (TTree *)fileMC->Get("parametri_geo");
        TTree *treeDT = (TTree *)fileDT->Get("parametri_geo");


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
	ToFPB.UseBetaEdges();
        NaFPB.UseBetaEdges();
        AglPB.UseBetaEdges();

	PRB.UseREdges();


        cout<<endl;


	cout<<"****************************** VARIABLES ***************************************"<<endl;

        Variables * vars = new Variables;

	cout<<"****************************** EFFICIENCY CORR. READING  ******************************************"<<endl;

	EffCorr * HEPPresEffCorr = new EffCorr(finalResults,"HEPPresEffCorr","HEPPresEffCorr",PRB,"PresControlSample","PresControlSample&IsPreselected","","IsProtonMC");
	EffCorr * HEPQualEffCorr = new EffCorr(finalResults,"HEPQualEffCorr","HEPQualEffCorr",PRB,"ControlSample","ControlSample&DistanceCut&LikelihoodCut","","IsProtonMC");
	
	EffCorr * RICHEffCorr_NaF = new EffCorr(finalResults,"RICHCorrection_NaF","RICH Eff. Corr",NaFPB,"ControlSample","ControlSample&IsFromNaF","","IsProtonMC");
	EffCorr * RICHEffCorr_Agl = new EffCorr(finalResults,"RICHCorrection_Agl","RICH Eff. Corr",AglPB,"ControlSample","ControlSample&IsFromAgl","","IsProtonMC");
/*
	EffCorrTemplate* DistCorr_TOF = new EffCorrTemplate(finalResults,"DistanceCorrTOF","Quality Eff. Corr",ToFDB,"ControlSample","ControlSample&DistanceCut","","");	
	EffCorrTemplate* DistCorr_NaF = new EffCorrTemplate(finalResults,"DistanceCorrNaF","Quality Eff. Corr",NaFDB,"ControlSample&IsFromNaF","ControlSample&IsFromNaF&DistanceCut","","",true);	
	EffCorrTemplate* DistCorr_Agl = new EffCorrTemplate(finalResults,"DistanceCorrAgl","Quality Eff. Corr",AglDB,"ControlSample&IsFromAgl","ControlSample&IsFromAgl&DistanceCut","","",true);	
	EffCorrTemplate* LikCorr_TOF = new EffCorrTemplate(finalResults,"LikelihoodCorrTOF","Quality Eff. Corr",ToFDB,"ControlSample&DistanceCut","ControlSample&DistanceCut&LikelihoodCut","","");	
	EffCorrTemplate* LikCorr_NaF = new EffCorrTemplate(finalResults,"LikelihoodCorrNaF","Quality Eff. Corr",NaFDB,"ControlSample&DistanceCut&IsFromNaF","ControlSample&IsFromNaF&DistanceCut&LikelihoodCut","","",true);	
	EffCorrTemplate* LikCorr_Agl = new EffCorrTemplate(finalResults,"LikelihoodCorrAgl","Quality Eff. Corr",AglDB,"ControlSample&DistanceCut&IsFromAgl","ControlSample&IsFromAgl&DistanceCut&LikelihoodCut","","",true);	
*/
	
	cout<<"****************************** FLUXES EVALUATION ******************************************"<<endl;

	cout<<"************* D FLUX ************"<<endl;

	Flux * DFluxTOF = new Flux(finalHistos,finalResults, "DFluxTOF", "FullsetEff_D_TOF","FullsetEfficiency","TOFfits/Fit Results/Primary Deuteron Counts","ExposureTOF","Gen. Acceptance",ToFDB);
	Flux * DFluxNaF = new Flux(finalHistos,finalResults, "DFluxNaF", "FullsetEff_D_NaF","FullsetEfficiency","NaFfits/Fit Results/Primary Deuteron Counts","ExposureNaF","Gen. Acceptance",NaFDB);
	Flux * DFluxAgl = new Flux(finalHistos,finalResults, "DFluxAgl", "FullsetEff_D_Agl","FullsetEfficiency","Aglfits/Fit Results/Primary Deuteron Counts","ExposureAgl","Gen. Acceptance",AglDB);

	DFluxTOF-> Eval_ExposureTime(vars,treeDT,finalHistos,Refill);
        DFluxNaF-> Eval_ExposureTime(vars,treeDT,finalHistos,Refill);
        DFluxAgl-> Eval_ExposureTime(vars,treeDT,finalHistos,Refill);

	DFluxTOF->Set_MCPar(0.5,20,0.0242236931);	
	DFluxNaF->Set_MCPar(0.5,20,0.0242236931);	
	DFluxAgl->Set_MCPar(0.5,20,0.0242236931);	

	DFluxTOF-> Eval_GeomAcceptance(treeMC,finalHistos,"IsDeutonMC",Refill);
        DFluxNaF-> Eval_GeomAcceptance(treeMC,finalHistos,"IsDeutonMC",Refill);
        DFluxAgl-> Eval_GeomAcceptance(treeMC,finalHistos,"IsDeutonMC",Refill);
/*
	DFluxNaF->ApplyEfficCorr(RICHEffCorr_NaF->GetGlobCorrection());
	DFluxAgl->ApplyEfficCorr(RICHEffCorr_Agl->GetGlobCorrection());

	DFluxTOF->ApplyEfficCorr(DistCorr_TOF->GetGlobCorrection_D());
	DFluxNaF->ApplyEfficCorr(DistCorr_NaF->GetGlobCorrection_D());
	DFluxAgl->ApplyEfficCorr(DistCorr_Agl->GetGlobCorrection_D());

	DFluxTOF->ApplyEfficCorr(LikCorr_TOF->GetGlobCorrection_D());
	DFluxNaF->ApplyEfficCorr(LikCorr_NaF->GetGlobCorrection_D());
	DFluxAgl->ApplyEfficCorr(LikCorr_Agl->GetGlobCorrection_D());

*/
	DFluxTOF-> Eval_Flux();
        DFluxNaF-> Eval_Flux();
        DFluxAgl-> Eval_Flux();

	DFluxTOF->SaveResults(finalResults);
	DFluxNaF->SaveResults(finalResults);
	DFluxAgl->SaveResults(finalResults);

	cout<<"************* P FLUX ************"<<endl;

	Flux * HEPFlux  = new Flux(finalHistos,finalResults,"PFluxHE", "RigBinFullSetEff","RigBinFullSetEff","HEPCounts/HEPCounts/HEPCounts","HEExposure","Gen. Acceptance",PRB);

	Flux * PFluxTOF = new Flux(finalHistos,finalResults, "PFluxTOF", "FullsetEff_P_TOF","FullsetEfficiency","TOFfits/Fit Results/Primary Proton Counts","ExposureTOF","Gen. Acceptance",ToFPB);
	Flux * PFluxNaF = new Flux(finalHistos,finalResults, "PFluxNaF", "FullsetEff_P_NaF","FullsetEfficiency","NaFfits/Fit Results/Primary Proton Counts","ExposureNaF","Gen. Acceptance",NaFPB);
	Flux * PFluxAgl = new Flux(finalHistos,finalResults, "PFluxAgl", "FullsetEff_P_Agl","FullsetEfficiency","Aglfits/Fit Results/Primary Proton Counts","ExposureAgl","Gen. Acceptance",AglPB);

	HEPFlux -> Eval_ExposureTime(vars,treeDT,finalHistos,Refill);
	PFluxTOF-> Eval_ExposureTime(vars,treeDT,finalHistos,Refill);
        PFluxNaF-> Eval_ExposureTime(vars,treeDT,finalHistos,Refill);
        PFluxAgl-> Eval_ExposureTime(vars,treeDT,finalHistos,Refill);

	HEPFlux ->Set_MCPar(0.5,100,0.0308232619);	
	PFluxTOF->Set_MCPar(0.5,100,0.0308232619);	
	PFluxNaF->Set_MCPar(0.5,100,0.0308232619);	
	PFluxAgl->Set_MCPar(0.5,100,0.0308232619);	

	HEPFlux -> Eval_GeomAcceptance(treeMC,finalHistos,"IsProtonMC",Refill,true);
	PFluxTOF-> Eval_GeomAcceptance(treeMC,finalHistos,"IsProtonMC",Refill);
        PFluxNaF-> Eval_GeomAcceptance(treeMC,finalHistos,"IsProtonMC",Refill);
        PFluxAgl-> Eval_GeomAcceptance(treeMC,finalHistos,"IsProtonMC",Refill);

//	HEPFlux -> ApplyEfficCorr(HEPPresEffCorr->GetGlobCorrection());
//	HEPFlux -> ApplyEfficCorr(HEPQualEffCorr->GetGlobCorrection());
/*	
	PFluxNaF->ApplyEfficCorr(RICHEffCorr_NaF->GetGlobCorrection());
	PFluxAgl->ApplyEfficCorr(RICHEffCorr_Agl->GetGlobCorrection());
	
	PFluxTOF->ApplyEfficCorr(DistCorr_TOF->GetGlobCorrection_P());
	PFluxNaF->ApplyEfficCorr(DistCorr_NaF->GetGlobCorrection_P());
	PFluxAgl->ApplyEfficCorr(DistCorr_Agl->GetGlobCorrection_P());

	PFluxTOF->ApplyEfficCorr(LikCorr_TOF->GetGlobCorrection_P());
	PFluxNaF->ApplyEfficCorr(LikCorr_NaF->GetGlobCorrection_P());
	PFluxAgl->ApplyEfficCorr(LikCorr_Agl->GetGlobCorrection_P());
*/
	HEPFlux -> Eval_Flux();
	PFluxTOF-> Eval_Flux();
        PFluxNaF-> Eval_Flux();
        PFluxAgl-> Eval_Flux();

	HEPFlux ->SaveResults(finalResults);
	PFluxTOF->SaveResults(finalResults);
	PFluxNaF->SaveResults(finalResults);
	PFluxAgl->SaveResults(finalResults);

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
