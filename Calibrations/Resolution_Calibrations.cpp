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
#include "TGraphErrors.h"

#include "../include/GlobalBinning.h"

#include "../Ntuple-making/Commonglobals.cpp"
#include "../Ntuple-making/Variables.hpp"
#include "../include/Cuts.h"

#include "../include/filesaver.h"

#include "../include/FitError.h"
#include "../include/Resolution.h"



float GetInverseRigidity (Variables * vars) {return 1/vars->R-1/vars->Momento_gen;}
float GetGenMomentum     (Variables * vars) {return vars->Momento_gen;}
float GetInverseEdepUToF (Variables * vars) {return 1/vars->EdepTOFU;}
float GetInverseEdepLToF (Variables * vars) {return 1/vars->EdepTOFD;}
float GetBetaTOF	 (Variables * vars) {return vars->Beta;}
float GetInverseEdepTrack(Variables * vars) {return 1/vars->EdepTrack;}
float GetBetaGen	 (Variables * vars) {return pow((pow((vars->Momento_gen/vars->Massa_gen),2)/(pow((vars->Momento_gen/vars->Massa_gen),2)+1)),0.5);}
float GetInverseBetaTOF  (Variables * vars) {return 1/vars->Beta - 1/GetBetaGen(vars);}
float GetInverseBetaRICH (Variables * vars) {return 1/vars->BetaRICH_new - 1/GetBetaGen(vars);}



bool Calculate_Resolution( Resolution * Reso, Variables * vars, bool checkfile, TTree * treeMC,  FileSaver finalHistos, float (*var) (Variables * vars),float (*discr_var) (Variables * vars),  std::vector<float> ExpValues={-1}, bool spline = false){
	if(!checkfile){
                Reso->Fill(treeMC,vars,var,discr_var);
		Reso->Normalize();
        }

        else Reso = new Resolution(finalHistos,Reso->GetName(),Reso->GetBinning());

        if(Reso->CheckHistos()){
                Reso->Eval_Resolution(ExpValues);
		Reso->Save(finalHistos);
        }
   	finalHistos.Add(Reso->Get_Means()	);
   	finalHistos.Add(Reso->Get_Sigmas()	);
	finalHistos.Add(Reso->Get_Resolutions()	);

	if(spline)   finalHistos.Add(Reso->ModelSigmasWithSpline());
	else 	     finalHistos.Add(Reso->ModelSigmasWithPoly());

	if(spline)   finalHistos.Add(Reso->ModelMeansWithSpline());
	else 	     finalHistos.Add(Reso->ModelMeansWithPoly());

	
   	finalHistos.writeObjsInFolder((Reso->GetName()+"/Fit Results").c_str());				
	
	return Reso->CheckHistos();
}

void CalculateMeanRatio( Resolution * Reso1,Resolution * Reso2, FileSaver finalHistos){

	Reso1 = new Resolution(finalHistos,Reso1->GetName(),Reso1->GetBinning());
	Reso2 = new Resolution(finalHistos,Reso2->GetName(),Reso2->GetBinning());

	finalHistos.Add(Reso1->ModelMeansRatio(Reso2));
        finalHistos.writeObjsInFolder((Reso1->GetName()+"/Fit Results").c_str());
	cout<<"Ratio: "<<Reso1->ModelMeansRatio(Reso2)->Eval(0.5)<<endl<<endl;;
	return;	
}


int main(int argc, char * argv[])
{


	cout<<"****************************** FILES OPENING ***************************************"<<endl;

	string INPUT1(argv[1]);
	string INPUT2(argv[2]);
        string OUTPUT(argv[3]);

	FileSaver finalHistos;
        finalHistos.setName(OUTPUT.c_str());
	bool checkfile = finalHistos.CheckFile();	

	TFile *fileDT =TFile::Open(INPUT1.c_str());
        TFile *fileMC =TFile::Open(INPUT2.c_str());
        
	TTree *treeMC = (TTree *)fileMC->Get("parametri_geo");
	TTree *treeDT = (TTree *)fileDT->Get("parametri_geo");


	cout<<"****************************** BINS ***************************************"<<endl;

	SetBins();	

	PResB.Print();

	cout<<"**TOF**"<<endl;
	ToFResB.Print();

	cout<<"**NaF**"<<endl;
	NaFResB.Print();

	cout<<"**Agl**"<<endl;
	AglResB.Print();

	ToFResB.UseBetaEdges();
	NaFResB.UseBetaEdges();
	AglResB.UseBetaEdges();

	PResB.UseREdges();


	cout<<endl;

	cout<<"****************************** VARIABLES ***************************************"<<endl;
	
	Variables * varsMC = new Variables;
        varsMC->ReadBranches(treeMC);

	Variables * varsDT = new Variables;
        varsDT->ReadBranches(treeDT);



	cout<<"****************************** ANALYSIS  ***************************************"<<endl;

	//rigidity resolution vs rigidity 
	Resolution * RigidityResolution_P = new Resolution("RvsR Resolution (P)",PResB,"IsPreselected&IsProtonMC",1000,-0.5,1.5);	
	Calculate_Resolution( RigidityResolution_P, varsMC,checkfile, treeMC, finalHistos, GetInverseRigidity, GetGenMomentum,PResB.RigBinsCent(),true);

	Resolution * RigidityResolution_D = new Resolution("RvsR Resolution (D)",PResB,"IsPreselected&IsDeutonMC",1000,-0.5,1.5);	
	Calculate_Resolution( RigidityResolution_D, varsMC,checkfile, treeMC, finalHistos, GetInverseRigidity, GetGenMomentum,PResB.RigBinsCent(),true);


	//beta resolution vs beta
	Resolution * BetaTOFResolution_P = new Resolution("BetaTOFvsBeta Resolution (P)",ToFResB,"IsPreselected&IsProtonMC",1000,-0.5,1.5);
	Calculate_Resolution( BetaTOFResolution_P, varsMC,checkfile, treeMC, finalHistos, GetInverseBetaTOF, GetBetaGen);

	Resolution * BetaNaFResolution_P = new Resolution("BetaNaFvsBeta Resolution (P)",NaFResB,"IsPreselected&IsProtonMC&IsFromNaF",250,-0.05,0.15);
	Calculate_Resolution( BetaNaFResolution_P, varsMC,checkfile, treeMC, finalHistos, GetInverseBetaRICH, GetBetaGen);

	Resolution * BetaAglResolution_P = new Resolution("BetaAglvsBeta Resolution (P)",AglResB,"IsPreselected&IsProtonMC&IsFromAgl",250,-0.075,0.075);
	Calculate_Resolution( BetaAglResolution_P, varsMC,checkfile, treeMC, finalHistos, GetInverseBetaRICH, GetBetaGen);

	Resolution * BetaTOFResolution_D = new Resolution("BetaTOFvsBeta Resolution (D)",ToFResB,"IsPreselected&IsDeutonMC",1000,-0.5,1.5);
	Calculate_Resolution( BetaTOFResolution_D, varsMC,checkfile, treeMC, finalHistos, GetInverseBetaTOF, GetBetaGen);

	Resolution * BetaNaFResolution_D = new Resolution("BetaNaFvsBeta Resolution (D)",NaFResB,"IsPreselected&IsDeutonMC&IsFromNaF",250,-0.05,0.15);
	Calculate_Resolution( BetaNaFResolution_D, varsMC,checkfile, treeMC, finalHistos, GetInverseBetaRICH, GetBetaGen);

	Resolution * BetaAglResolution_D = new Resolution("BetaAglvsBeta Resolution (D)",AglResB,"IsPreselected&IsDeutonMC&IsFromAgl",250,-0.075,0.075);
	Calculate_Resolution( BetaAglResolution_D, varsMC,checkfile, treeMC, finalHistos, GetInverseBetaRICH, GetBetaGen);


	// E. dep. (U.ToF, L.ToF, Inner Tracker)
	
	Resolution * EdepUTOFResolution_P = new Resolution("EdepUTOFvsBeta Resolution (P)",ToFResB,"IsPreselected&IsProtonMC",1000,-1.5,1.5);
	Calculate_Resolution( EdepUTOFResolution_P, varsMC,checkfile, treeMC, finalHistos, GetInverseEdepUToF, GetBetaGen);

        Resolution * EdepUTOFResolution_D = new Resolution("EdepUTOFvsBeta Resolution (D)",ToFResB,"IsPreselected&IsDeutonMC",1000,-1.5,1.5);
	Calculate_Resolution( EdepUTOFResolution_D, varsMC,checkfile, treeMC, finalHistos, GetInverseEdepUToF, GetBetaGen);
	
	Resolution * EdepLTOFResolution_P = new Resolution("EdepLTOFvsBeta Resolution (P)",ToFResB,"IsPreselected&IsProtonMC",1000,-1.5,1.5);
	Calculate_Resolution( EdepLTOFResolution_P, varsMC,checkfile, treeMC, finalHistos, GetInverseEdepLToF, GetBetaGen);

        Resolution * EdepLTOFResolution_D = new Resolution("EdepUTOFvsBeta Resolution (D)",ToFResB,"IsPreselected&IsDeutonMC",1000,-1.5,1.5);
        Calculate_Resolution( EdepLTOFResolution_D, varsMC,checkfile, treeMC, finalHistos, GetInverseEdepLToF, GetBetaGen);

	Resolution * EdepTrackResolution_P = new Resolution("EdepTrackvsBeta Resolution (P)",ToFResB,"IsPreselected&IsProtonMC",500,-1.5,15);
	Calculate_Resolution( EdepTrackResolution_P, varsMC,checkfile, treeMC, finalHistos, GetInverseEdepTrack, GetBetaGen);
	
	Resolution * EdepTrackResolution_D = new Resolution("EdepTrackvsBeta Resolution (D)",ToFResB,"IsPreselected&IsDeutonMC",500,-1.5,15);
	Calculate_Resolution( EdepTrackResolution_D, varsMC,checkfile, treeMC, finalHistos, GetInverseEdepTrack, GetBetaGen);


	// E.dep. calibration
	
	//Upper ToF 	
	Resolution * EdepUTOFMC_P = new Resolution("EdepUTOFvsBeta Measured MC",ToFResB,"IsPreselected&IsProtonMC&L1LooseCharge1",1000,-1.5,1.5);
	Calculate_Resolution( EdepUTOFMC_P, varsMC,checkfile, treeMC, finalHistos, GetInverseEdepUToF, GetBetaTOF);

	Resolution * EdepUTOFDT_P = new Resolution("EdepUTOFvsBeta Measured DT",ToFResB,"IsPreselected&IsData&L1LooseCharge1",1000,-1.5,1.5);
	Calculate_Resolution( EdepUTOFDT_P, varsDT,checkfile, treeDT, finalHistos, GetInverseEdepUToF, GetBetaTOF);

	CalculateMeanRatio( EdepUTOFDT_P,EdepUTOFMC_P, finalHistos);

	// Inner Tracker
	Resolution * EdepTrackMC_P = new Resolution("EdepTrackvsBeta Measured MC",ToFResB,"IsPreselected&IsProtonsMC&L1LooseCharge1",500,-1.5,15);
        Calculate_Resolution( EdepTrackMC_P, varsMC,checkfile, treeMC, finalHistos, GetInverseEdepTrack, GetBetaTOF);

	Resolution * EdepTrackDT_P = new Resolution("EdepTrackvsBeta Measured DT",ToFResB,"IsPreselected&IsData&L1LooseCharge1",500,-1.5,15);
	Calculate_Resolution( EdepTrackDT_P, varsDT,checkfile, treeDT, finalHistos, GetInverseEdepTrack, GetBetaTOF);

	CalculateMeanRatio( EdepTrackDT_P,EdepTrackMC_P, finalHistos);


	//Lower ToF
	Resolution * EdepLTOFMC_P = new Resolution("EdepLTOFvsBeta Measured MC",ToFResB,"IsPreselected&IsProtonMC&L1LooseCharge1",1000,-1.5,1.5);
	Calculate_Resolution( EdepLTOFMC_P, varsMC,checkfile, treeMC, finalHistos, GetInverseEdepLToF, GetBetaTOF);

	Resolution * EdepLTOFDT_P = new Resolution("EdepLTOFvsBeta Measured DT",ToFResB,"IsPreselected&IsData&L1LooseCharge1",1000,-1.5,1.5);
	Calculate_Resolution( EdepLTOFDT_P, varsDT,checkfile, treeDT, finalHistos, GetInverseEdepLToF, GetBetaTOF);

	CalculateMeanRatio( EdepLTOFDT_P,EdepLTOFMC_P, finalHistos);

	return 0;
}


