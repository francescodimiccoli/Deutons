#include <bitset>
#include "TROOT.h"
#include "TNtuple.h"
#include <TSpline.h>
#include "../include/binning.h"
#include "TFile.h"
#include "TH2.h"
#include "TF2.h"
#include <TVector3.h>
#include "TMath.h"
#include <TFile.h>
#include "TFile.h"
#include "TH2.h"
#include "TF2.h"
#include <TVector3.h>
#include "TMath.h"
#include <TMVA/Reader.h>
#include <TMVA/Tools.h>
#include "TRandom3.h"

#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "../include/GlobalBinning.h"

#include "Commonglobals.cpp"
#include "Variables.hpp"

#include "../include/Cuts.h"
#include "../include/filesaver.h"

#include "../include/FitError.h"
#include "../include/Resolution.h"

#include "Discriminants5D.h"

#include "reweight.h"
#include "histUtils.h"
#include "../include/binning.h"
#include "Functions.hpp"





bool ReadCalibration(string month);
bool ReadPdfForLikelihood();
int AssignMC_type(float Massa_gen);
float mcweight=0;
bool IsMC(TTree * tree);
Reweighter ReweightInitializer();
void CalibrateEdep(Variables *vars);

void ProcessEvent(Variables *vars,bool isMC,Reweighter reweighter);


int main(int argc, char * argv[])
{
	string month = argv[1];
	if(!ReadCalibration(month)) return 0;
	if(!ReadPdfForLikelihood()) return 0;
	string INPUT(argv[2]);
	string OUTPUT(argv[3]);

	TFile *file =TFile::Open(INPUT.c_str());
	TTree *tree = (TTree *)file->Get("parametri_geo");
	TFile * File = new TFile(OUTPUT.c_str(), "RECREATE");

	bool check=tree?true:false;
	if(!check) {cout<<"Corrupted input file"<<endl;return 0;}

	bool isMC = IsMC(tree);

	Variables * vars = new Variables;
	vars->ReadBranches(tree);

	TNtuple * grandezzesepd;
	TNtuple * trig;
	TNtuple * Q;

	if(isMC) cout<<"It is a MC file"<<endl; 
	if(isMC){
		grandezzesepd=new TNtuple("grandezzesepd","grandezzesepd","R:Beta:EdepL1:Massa_gen:Cutmask:PhysBPatt:EdepTOF:EdepTrack:EdepTOFD:Momentogen:BetaRICH_new:LDiscriminant:mcweight:Dist5D:Dist5D_P");
		trig=new TNtuple("trig","trig","Massa_gen:Momento_gen:R_L1:R_pre:Beta_pre:Cutmask:EdepL1:EdepTOFU:EdepTOFD:EdepTrack:BetaRICH:Second:PhysBPatt:mcweight");
		Q = new TNtuple("Q","Q","R:Beta:qL1:Massa_gen:Cutmask:PhysBPatt:qUtof:qInner:qLtof:Momentogen:BetaRICH_new:LDiscriminant:mcweight:Dist5D:Dist5D_P");
	}
	else{
		grandezzesepd = new TNtuple("grandezzesepd","grandezzesepd","R:Beta:EdepL1:Cutmask:Latitude:PhysBPatt:EdepTOFU:EdepTrack:EdepTOFD:IGRFRcutoff:BetaRICH_new:LDiscriminant:Dist5D:Dist5D_P:Rcutoff");
		trig=new TNtuple("trig","trig","U_Time:Latitude:Rcutoff:IGRFRcutoff:R_pre:Beta_pre:Cutmask:EdepL1:EdepTOFU:EdepTOFD:EdepTrack:BetaRICH:Second:PhysBPatt:Livetime");
		Q = new TNtuple("Q","Q","R:Beta:qL1:Cutmask:Latitude:PhysBPatt:qUtof:qInner:qLtof:IGRFRcutoff:BetaRICH_new:LDiscriminant:Dist5D:Dist5D_P:Rcutoff"); 
	}

	Reweighter reweighter;
	if(isMC) reweighter=ReweightInitializer(); 	

	for(int i=0; i<tree->GetEntries(); i++){
		UpdateProgressBar(i, tree->GetEntries());	
		tree->GetEvent(i);
		vars->Update();
		
		if(i==0) Timebeg = vars->U_time;	
		
		ProcessEvent(vars,isMC,reweighter);	

		vars->FillwithAnalysisVariables(grandezzesepd,isMC);
		vars->FillwithRawVariables(trig,isMC);	
		vars->FillwithChargeInfos(Q,isMC); 
	}
	File->Write();
	File->Close();

	return 0;

}


