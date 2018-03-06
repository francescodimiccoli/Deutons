#include <bitset>
#include "TROOT.h"
#include "TNtuple.h"
#include <TSpline.h>
#include "../include/binning.h"
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
#include "Globals.h"
#include "../include/InputFileReader.h"

#include "../include/filesaver.h"
#include "../include/Variables.hpp"
#include "../include/FitError.h"
#include "../include/Resolution.h"





bool Calculate_Resolution( Resolution * Reso, Variables * vars, bool checkfile, DBarReader reader,  FileSaver finalHistos, float (*var) (Variables * vars),float (*discr_var) (Variables * vars),  std::vector<float> ExpValues={-1}, bool spline = false,float low_limit=0.05,float high_limit=0.85,bool fixedwindow=false){
	if(!checkfile){
                Reso->Fill(reader,vars,var,discr_var);
		Reso->Normalize();
        }

        else Reso = new Resolution(finalHistos,Reso->GetName(),Reso->GetBinning());
	if(fixedwindow) Reso->SetFixedFitWindow(low_limit,high_limit);
        if(Reso->CheckHistos()){
                Reso->Eval_Resolution(ExpValues,low_limit,high_limit);
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

	string refill="";
	if(argc > 4 ) 	refill = argv[4];	

	bool Refill = false;
	if(refill!="") Refill=true;

	TChain * chainDT = InputFileReader(INPUT1.c_str(),"Event");
	TChain * chainMC = InputFileReader(INPUT2.c_str(),"Event");

	FileSaver finalHistos;
        finalHistos.setName(OUTPUT.c_str());
        bool checkfile = finalHistos.CheckFile();
	
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

	Variables * varsDT = new Variables;

	std::vector<float> meanedep;
        for(int i=0;i<ToFResB.size();i++) meanedep.push_back(1.);

	cout<<"****************************** ANALYSIS  ***************************************"<<endl;
	
	//rigidity resolution vs rigidity 
	Resolution * RigidityResolution_P = new Resolution("RvsR Resolution (P)",PResB,"IsPositive&IsPreselected&IsProtonMC",1000,-0.5,1.5);	
	Calculate_Resolution( RigidityResolution_P, varsMC,checkfile, DBarReader(chainMC, true), finalHistos, GetInverseRigidity, GetGenMomentum,PResB.RigBinsCent(),true);

	Resolution * RigidityResolution_D = new Resolution("RvsR Resolution (D)",PResB,"IsPositive&IsPreselected&IsDeutonMC",1000,-0.5,1.5);	
	Calculate_Resolution( RigidityResolution_D, varsMC,checkfile, DBarReader(chainMC, true), finalHistos, GetInverseRigidity, GetGenMomentum,PResB.RigBinsCent(),true);


	//beta resolution vs beta
	Resolution * BetaTOFResolution_P = new Resolution("BetaTOFvsBeta Resolution (P)",ToFResB,"IsPreselected&IsProtonMC",1000,-0.5,1.5);
	Calculate_Resolution( BetaTOFResolution_P, varsMC,checkfile, DBarReader(chainMC, true ), finalHistos, GetInverseBetaTOF, GetBetaGen);

	Resolution * BetaNaFResolution_P = new Resolution("BetaNaFvsBeta Resolution (P)",NaFResB,"IsPreselected&IsProtonMC&IsFromNaF",250,-0.05,0.15);
	Calculate_Resolution( BetaNaFResolution_P, varsMC,checkfile, DBarReader(chainMC, true ), finalHistos, GetInverseBetaRICH, GetBetaGen);

	Resolution * BetaAglResolution_P = new Resolution("BetaAglvsBeta Resolution (P)",AglResB,"IsPreselected&IsProtonMC&IsFromAgl",250,-0.02,0.075);
	Calculate_Resolution( BetaAglResolution_P, varsMC,checkfile, DBarReader(chainMC, true ), finalHistos, GetInverseBetaRICH, GetBetaGen,AglResB.BetaBinsCent(),false,0.315,0.93);

	Resolution * BetaTOFResolution_D = new Resolution("BetaTOFvsBeta Resolution (D)",ToFResB,"IsPreselected&IsDeutonMC",1000,-0.5,1.5);
	Calculate_Resolution( BetaTOFResolution_D, varsMC,checkfile, DBarReader(chainMC, true ), finalHistos, GetInverseBetaTOF, GetBetaGen);

	Resolution * BetaNaFResolution_D = new Resolution("BetaNaFvsBeta Resolution (D)",NaFResB,"IsPreselected&IsDeutonMC&IsFromNaF",250,-0.05,0.15);
	Calculate_Resolution( BetaNaFResolution_D, varsMC,checkfile, DBarReader(chainMC, true ), finalHistos, GetInverseBetaRICH, GetBetaGen);

	Resolution * BetaAglResolution_D = new Resolution("BetaAglvsBeta Resolution (D)",AglResB,"IsPreselected&IsDeutonMC&IsFromAgl",500,-0.02,0.075);
	Calculate_Resolution( BetaAglResolution_D, varsMC,checkfile, DBarReader(chainMC, true ), finalHistos, GetInverseBetaRICH, GetBetaGen,AglResB.BetaBinsCent(),false,0.315,0.93);

	// E. dep. (U.ToF, L.ToF, Inner Tracker)
	
	Resolution * EdepUTOFResolution_P = new Resolution("EdepUTOFvsBeta Resolution (P)",ToFResB,"IsPreselected&IsProtonMC",1000,-1.5,1.5);
	Calculate_Resolution( EdepUTOFResolution_P, varsMC,checkfile, DBarReader(chainMC, true ), finalHistos, GetInverseEdepUToF, GetBetaGen,meanedep,true,0.35,0.95);

        Resolution * EdepUTOFResolution_D = new Resolution("EdepUTOFvsBeta Resolution (D)",ToFResB,"IsPreselected&IsDeutonMC",1000,-1.5,1.5);
	Calculate_Resolution( EdepUTOFResolution_D, varsMC,checkfile, DBarReader(chainMC, true ), finalHistos, GetInverseEdepUToF, GetBetaGen,meanedep,true,0.35,0.95);
	
	Resolution * EdepLTOFResolution_P = new Resolution("EdepLTOFvsBeta Resolution (P)",ToFResB,"IsPreselected&IsProtonMC",1000,-1.5,1.5);
	Calculate_Resolution( EdepLTOFResolution_P, varsMC,checkfile, DBarReader(chainMC, true ), finalHistos, GetInverseEdepLToF, GetBetaGen,meanedep,true,0.35,0.95);

        Resolution * EdepLTOFResolution_D = new Resolution("EdepLTOFvsBeta Resolution (D)",ToFResB,"IsPreselected&IsDeutonMC",1000,-1.5,1.5);
        Calculate_Resolution( EdepLTOFResolution_D, varsMC,checkfile, DBarReader(chainMC, true ), finalHistos, GetInverseEdepLToF, GetBetaGen,meanedep,true,0.35,0.95);

	Resolution * EdepTrackResolution_P = new Resolution("EdepTrackvsBeta Resolution (P)",ToFResB,"IsPreselected&IsProtonMC",500,-1.5,15);
	Calculate_Resolution( EdepTrackResolution_P, varsMC,checkfile, DBarReader(chainMC, true ), finalHistos, GetInverseEdepTrack, GetBetaGen,meanedep,true,0.1,0.75);
	
	Resolution * EdepTrackResolution_D = new Resolution("EdepTrackvsBeta Resolution (D)",ToFResB,"IsPreselected&IsDeutonMC",500,-1.5,15);
	Calculate_Resolution( EdepTrackResolution_D, varsMC,checkfile, DBarReader(chainMC, true ), finalHistos, GetInverseEdepTrack, GetBetaGen,meanedep,true,0.1,0.75);


	// E.dep. calibration
	
	//Upper ToF 	
	Resolution * EdepUTOFMC_P = new Resolution("EdepUTOFvsBeta Measured MC",ToFResB,"IsPreselected&IsProtonMC&L1LooseCharge1",1000,-1.5,1.5);
	Calculate_Resolution( EdepUTOFMC_P, varsMC,checkfile, DBarReader(chainMC, true ), finalHistos, GetInverseEdepUToF, GetBetaTOF,meanedep,true,0.35,0.95);

	Resolution * EdepUTOFDT_P = new Resolution("EdepUTOFvsBeta Measured DT",ToFResB,"IsPreselected&IsData&L1LooseCharge1",1000,-1.5,1.5);
	Calculate_Resolution( EdepUTOFDT_P, varsDT,checkfile, DBarReader(chainDT, false ), finalHistos, GetInverseEdepUToF, GetBetaTOF,meanedep,true,0.35,0.95);

	CalculateMeanRatio( EdepUTOFDT_P,EdepUTOFMC_P, finalHistos);

	// Inner Tracker
	Resolution * EdepTrackMC_P = new Resolution("EdepTrackvsBeta Measured MC",ToFResB,"IsPreselected&IsProtonsMC&L1LooseCharge1",500,-1.5,15);
        Calculate_Resolution( EdepTrackMC_P, varsMC,checkfile, DBarReader(chainMC, true ), finalHistos, GetInverseEdepTrack, GetBetaTOF,meanedep,true,0.1,0.75);

	Resolution * EdepTrackDT_P = new Resolution("EdepTrackvsBeta Measured DT",ToFResB,"IsPreselected&IsData&L1LooseCharge1",500,-1.5,15);
	Calculate_Resolution( EdepTrackDT_P, varsDT,checkfile, DBarReader(chainDT, false ), finalHistos, GetInverseEdepTrack, GetBetaTOF,meanedep,true,0.1,0.75);

	CalculateMeanRatio( EdepTrackDT_P,EdepTrackMC_P, finalHistos);


	//Lower ToF
	Resolution * EdepLTOFMC_P = new Resolution("EdepLTOFvsBeta Measured MC",ToFResB,"IsPreselected&IsProtonMC&L1LooseCharge1",1000,-1.5,1.5);
	Calculate_Resolution( EdepLTOFMC_P, varsMC,checkfile, DBarReader(chainMC, true ), finalHistos, GetInverseEdepLToF, GetBetaTOF,meanedep,true,0.35,0.95);

	Resolution * EdepLTOFDT_P = new Resolution("EdepLTOFvsBeta Measured DT",ToFResB,"IsPreselected&IsData&L1LooseCharge1",1000,-1.5,1.5);
	Calculate_Resolution( EdepLTOFDT_P, varsDT,checkfile, DBarReader(chainDT, false ), finalHistos, GetInverseEdepLToF, GetBetaTOF,meanedep,true,0.35,0.95);

	CalculateMeanRatio( EdepLTOFDT_P,EdepLTOFMC_P, finalHistos);

	//TRD
	Resolution * EdepTRDMC_P = new Resolution("EdepTRDvsBeta Measured MC",ToFResB,"IsPreselected&IsProtonMC&L1LooseCharge1",1000,-0.15,0.15);
	Calculate_Resolution( EdepTRDMC_P, varsMC,checkfile, DBarReader(chainMC, true ), finalHistos, GetInverseEdepTRD, GetBetaTOF,meanedep,true);

	Resolution * EdepTRDDT_P = new Resolution("EdepTRDvsBeta Measured DT",ToFResB,"IsPreselected&IsData&L1LooseCharge1",1000,-0.15,0.15);
	Calculate_Resolution( EdepTRDDT_P, varsDT,checkfile, DBarReader(chainDT, false ), finalHistos, GetInverseEdepTRD, GetBetaTOF,meanedep,true);

	CalculateMeanRatio( EdepTRDDT_P,EdepTRDMC_P, finalHistos);




	// Charge Calibration Check
/*	
	//QUtof 
	Resolution * QUTOFMC_P = new Resolution("QUTOFvsBeta Measured MC",ToFResB,"IsPreselected&IsProtonMC&L1LooseCharge1",300,0,3); 
	Calculate_Resolution( QUTOFMC_P, varsMC,checkfile, DBarReader(chainMC, true ), finalHistos, GetUtofQ, GetBetaTOF,meanedep,false,0.8,1.08,true);

	Resolution * QUTOFDT_P = new Resolution("QUTOFvsBeta Measured DT",ToFResB,"IsPreselected&IsData&L1LooseCharge1",300,0,3);
        Calculate_Resolution( QUTOFDT_P, varsDT,checkfile, DBarReader(chainDT, false ), finalHistos, GetLtofQ, GetBetaTOF,meanedep,false,0.8,1.08,true);	
	
	CalculateMeanRatio( QUTOFDT_P,QUTOFMC_P, finalHistos);

	//QLtof
	Resolution * QLTOFMC_P = new Resolution("QLTOFvsBeta Measured MC",ToFResB,"IsPreselected&IsProtonMC&L1LooseCharge1",300,0,3); 
	Calculate_Resolution( QLTOFMC_P, varsMC,checkfile, DBarReader(chainMC, true ), finalHistos, GetLtofQ, GetBetaTOF,meanedep,false,0.8,1.08,true);

	Resolution * QLTOFDT_P = new Resolution("QLTOFvsBeta Measured DT",ToFResB,"IsPreselected&IsData&L1LooseCharge1",300,0,3);
	Calculate_Resolution( QLTOFDT_P, varsDT,checkfile, DBarReader(chainDT, false ), finalHistos, GetLtofQ, GetBetaTOF,meanedep,false,0.8,1.08,true);	

	CalculateMeanRatio( QLTOFDT_P,QLTOFMC_P, finalHistos);

	//QInner
	Resolution * QInnerMC_P = new Resolution("QInnervsBeta Measured MC",ToFResB,"IsPreselected&IsProtonMC&L1LooseCharge1",300,0,3); 
	Calculate_Resolution( QInnerMC_P, varsMC,checkfile, DBarReader(chainMC, true ), finalHistos, GetInnerQ, GetBetaTOF,meanedep,false,0.8,1.08,true);

	Resolution * QInnerDT_P = new Resolution("QInnervsBeta Measured DT",ToFResB,"IsPreselected&IsData&L1LooseCharge1",300,0,3);
        Calculate_Resolution( QInnerDT_P, varsDT,checkfile, DBarReader(chainDT, false ), finalHistos, GetInnerQ, GetBetaTOF,meanedep,false,0.8,1.08,true);	

	CalculateMeanRatio( QInnerDT_P,QInnerMC_P, finalHistos);
*/

	

	return 0;
}


