#include <bitset>
#include "TROOT.h"
#include "TNtuple.h"
#include <TSpline.h>
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include <TVector3.h>
#include "TMath.h"
#include <TFile.h>
#include <string>
#include "TFile.h"
#include "TH2.h"
#include "TF2.h"
#include <TVector3.h>
#include "TChain.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TFractionFitter.h"
#include "TRandom3.h"

#include "InputFileReader.h"
#include "DBarReader.h"

#include "filesaver.h"
#include "binning.h"
#include "Globals.h"
#include "Variables.hpp"
#include "ParallelFiller.h"

#include "HistoBooker.h"

int main(int argc, char * argv[])
{


	cout<<"****************************** FILES OPENING ***************************************"<<endl;

	TH1::SetDefaultSumw2();     	
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
	SetUpUsualBinning();

	cout<<"****************************** VARIABLES ***************************************"<<endl;
	Variables * vars = new Variables;

	cout<<"****************************** ANALYIS ******************************************"<<endl;
	HistoBooker Booker;

//	Booker.BookSingleScatter("RigvsBetaTOF",300,0,4,300,0.4,1 ,"IsPreselected&LikelihoodCut&DistanceCut",GetRigidity,GetBetaTOF,GetBetaTOF);
//	Booker.BookSingleScatter("RigvsBetaNaF",300,1,5,300,0.7,1 ,"IsPreselected&LikelihoodCut&DistanceCut&IsFromNaF",GetRigidity,GetBetaRICH,GetBetaRICH);
//	Booker.BookSingleScatter("RigvsBetaAgl",300,2,10,300,0.9,1,"IsPreselected&LikelihoodCut&DistanceCut&IsFromAgl",GetRigidity,GetBetaRICH,GetBetaRICH);
/*
	Booker.BookSingleScatter("NegMassvsTRDLIkTOF",300,0,4,300,-20,20 ,"IsPreselected&LikelihoodCut&DistanceCut&TofBetaSafetyCut",GetNegRecMassTOF,GetTRDePLikRatio,GetBetaTOF);
	Booker.BookSingleScatter("NegMassvsTRDLIkNaF",300,0,4,300,-20,20 ,"IsPreselected&LikelihoodCut&DistanceCut&IsFromNaF&NafBetaSafetyCut",GetNegRecMassRICH,GetTRDePLikRatio,GetBetaRICH);
	Booker.BookSingleScatter("NegMassvsTRDLIkAgl",300,0,4,300,-20,20 ,"IsPreselected&LikelihoodCut&DistanceCut&IsFromAgl&AglBetaSafetyCut",GetNegRecMassRICH,GetTRDePLikRatio,GetBetaRICH);

	Booker.BookSingleHisto("MassNaF_noSel",100,0,5,"IsPositive&IsPreselected&DistanceCut&NaFBetaSafetyCut&IsFromNaF_nosel",GetRecMassRICH,GetBetaRICH);
	Booker.BookSingleHisto("MassNaF_Sel",100,0,5,"IsPositive&IsPreselected&DistanceCut&NaFBetaSafetyCut&IsFromNaF",GetRecMassRICH,GetBetaRICH);
	Booker.BookSingleHisto("MassNaF_Qual",100,0,5,"IsPositive&IsPreselected&DistanceCut&Likelihood&NaFBetaSafetyCut&IsFromNaF",GetRecMassRICH,GetBetaRICH);
	Booker.BookSingleHisto("MassAgl_noSel",100,0,5,"IsPositive&IsPreselected&DistanceCut&AglBetaSafetyCut&IsFromAgl_nosel",GetRecMassRICH,GetBetaRICH);
	Booker.BookSingleHisto("MassAgl_Sel",100,0,5,"IsPositive&IsPreselected&DistanceCut&AglBetaSafetyCut&IsFromAgl",GetRecMassRICH,GetBetaRICH);
	Booker.BookSingleHisto("MassAgl_Qual",100,0,5,"IsPositive&IsPreselected&DistanceCut&Likelihood&AglBetaSafetyCut&IsFromAgl",GetRecMassRICH,GetBetaRICH);

	Booker.BookSingleHisto("AntiDMassTOF",100,0,4,"IsPreselected&LikelihoodCut&DistanceCut&TofBetaSafetyCut",GetNegRecMassTOF);
	Booker.BookSingleHisto("AntiDMassNaF",100,0,4,"IsPreselected&LikelihoodCut&DistanceCut&IsFromNaF&NafBetaSafetyCut",GetNegRecMassRICH);
	Booker.BookSingleHisto("AntiDMassAgl",100,0,4,"IsPreselected&LikelihoodCut&DistanceCut&IsFromAgl&AglBetaSafetyCut",GetNegRecMassRICH);

	
	Booker.FillEverything(DBarReader(chainDT, false));
	Booker.SaveEverything(finalHistos);
*/
	HistoBooker BookerMC;

	BookerMC.BookSingleScatter("BetagenvsBetaMeasTOFP",300,0.3,1,300,0.3,1,"IsPreselected&LikelihoodCut&DistanceCut&IsProtonMC",GetBetaGen,GetBetaTOF,GetBetaTOF);
	BookerMC.BookSingleScatter("BetagenvsBetaMeasNaFP",300,0.7,1,300,0.7,1,"IsPreselected&LikelihoodCut&DistanceCut&IsProtonMC&IsFromNaF",GetBetaGen,GetBetaRICH,GetBetaRICH);
	BookerMC.BookSingleScatter("BetagenvsBetaMeasAglP",300,0.9,1,300,0.9,1,"IsPreselected&LikelihoodCut&DistanceCut&IsProtonMC&IsFromAgl",GetBetaGen,GetBetaRICH,GetBetaRICH);

	BookerMC.BookSingleScatter("BetagenvsBetaMeasTOFHe",300,0.3,1,300,0.3,1,"IsPreselectedHe&LikelihoodCut&&IsHeliumMC"	    ,GetBetaGen,GetBetaTOF,GetBetaTOF);
	BookerMC.BookSingleScatter("BetagenvsBetaMeasNaFHe",300,0.7,1,300,0.7,1,"IsPreselectedHe&LikelihoodCut&IsHeliumMC&IsFromNaF",GetBetaGen,GetBetaRICH,GetBetaRICH);
	BookerMC.BookSingleScatter("BetagenvsBetaMeasAglHe",300,0.9,1,300,0.9,1,"IsPreselectedHe&LikelihoodCut&IsHeliumMC&IsFromAgl",GetBetaGen,GetBetaRICH,GetBetaRICH);

	BookerMC.BookSingleScatter("BetagenvsBetaMeasTOFD",300,0.3,1,300,0.3,1,"IsPreselectedHe&LikelihoodCut&&IsDeutonMC"	    ,GetBetaGen,GetBetaTOF,GetBetaTOF);
	BookerMC.BookSingleScatter("BetagenvsBetaMeasNaFD",300,0.7,1,300,0.7,1,"IsPreselectedHe&LikelihoodCut&IsDeutonMC&IsFromNaF",GetBetaGen,GetBetaRICH,GetBetaRICH);
	BookerMC.BookSingleScatter("BetagenvsBetaMeasAglD",300,0.9,1,300,0.9,1,"IsPreselectedHe&LikelihoodCut&IsDeutonMC&IsFromAgl",GetBetaGen,GetBetaRICH,GetBetaRICH);

	BookerMC.BookSingleScatter("RICHBDTvsMassNaFP",100,0,4.5,100,-0.5,0.6,"IsPositive&IsPreselected&LikelihoodCut&DistanceCut&IsFromNaF&NaFBetaSafetyCut&IsProtonMC",GetRecMassRICH,GetRICHBDT,GetRICHBDT);
	BookerMC.BookSingleScatter("RICHBDTvsMassAglP",100,0,4.5,100,-0.5,0.6,"IsPositive&IsPreselected&LikelihoodCut&DistanceCut&IsFromAgl&AglBetaSafetyCut&IsProtonMC",GetRecMassRICH,GetRICHBDT,GetRICHBDT);

	BookerMC.BookSingleScatter("RICHBDTvsMassNaFD",100,0,4.5,100,-0.5,0.6,"IsPositive&IsPreselected&LikelihoodCut&DistanceCut&IsFromNaF&NaFBetaSafetyCut&IsDeutonMC",GetRecMassRICH,GetRICHBDT,GetRICHBDT);
	BookerMC.BookSingleScatter("RICHBDTvsMassAglD",100,0,4.5,100,-0.5,0.6,"IsPositive&IsPreselected&LikelihoodCut&DistanceCut&IsFromAgl&AglBetaSafetyCut&IsDeutonMC",GetRecMassRICH,GetRICHBDT,GetRICHBDT);


	BookerMC.FillEverything(DBarReader(chainMC, true));
        BookerMC.SaveEverything(finalHistos);
	
	return 0;
}




