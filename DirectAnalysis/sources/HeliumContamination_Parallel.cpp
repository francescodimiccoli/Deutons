#include <bitset>
#include "TROOT.h"
#include "TTree.h"
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
#include "../include/BetaSmearing.h"
#include "../include/Cuts.h"
#include "../include/ParallelFiller.h"

#include "../include/filesaver.h"

#include "../include/Efficiency.h"


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

        FileSaver finalHistos;
        finalHistos.setName(OUTPUT.c_str());
        bool checkfile = finalHistos.CheckFile();
	if(!checkfile) Refill=true;

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

        PRB.UseREdges();


        cout<<endl;


	cout<<"****************************** VARIABLES ***************************************"<<endl;

        Variables * vars = new Variables;

	cout<<"****************************** ANALYIS ******************************************"<<endl;

	BadEventSimulator * NaFBadEvSimulator= new BadEventSimulator("IsFromNaF",22,0.8,1);
        BadEventSimulator * AglBadEvSimulator= new BadEventSimulator("IsFromAgl",250,0.95,1);

	// TOF
	
	Efficiency * HeliumFragmTOF = new Efficiency(finalHistos,"HeliumFrragmTOF","HeliumFragmentation",ToFDB,"IsPreselected&IsHeliumMC&IsGoodHe" ,"IsPreselected&IsHeliumMC&LikelihoodCut&DistanceCut");
	Efficiency * HeliumFragmIntoDTOF = new Efficiency(finalHistos,"HeliumFragmIntoDTOF","HeliumFragmentation",ToFDB,"IsPreselected&IsHeliumMC&LikelihoodCut&DistanceCut" ,"IsPreselected&IsHeliumMC&LikelihoodCut&DistanceCut&DeutonsMassCut");
	Efficiency * HeContaminationTOF = new Efficiency(finalHistos,"HeContTOF","HeliumFragmentation",ToFDB,"IsPreselected&LikelihoodCut&DistanceCut","IsPreselected&IsGoodHe");

	// NaF
	
	Efficiency * HeliumFragmNaF = new Efficiency(finalHistos,"HeliumFrragmNaF","HeliumFragmentation",NaFDB,"IsPreselected&IsHeliumMC&IsFromNaF&IsGoodHe" ,"IsPreselected&IsHeliumMC&LikelihoodCut&DistanceCut&IsFromNaF");
	Efficiency * HeliumFragmIntoDNaF = new Efficiency(finalHistos,"HeliumFragmIntoDNaF","HeliumFragmentation",NaFDB,"IsPreselected&IsHeliumMC&LikelihoodCut&DistanceCut&IsFromNaF" ,"IsPreselected&IsHeliumMC&LikelihoodCut&DistanceCut&IsFromNaF&DeutonsMassCut");
	Efficiency * HeContaminationNaF = new Efficiency(finalHistos,"HeContNaF","HeliumFragmentation",NaFDB,"IsPreselected&LikelihoodCut&DistanceCut&IsFromNaF"  ,"IsPreselected&IsFromNaF&IsGoodHe"); 

	HeliumFragmNaF->SetUpBadEventSimulator(NaFBadEvSimulator);
	HeliumFragmIntoDNaF->SetUpBadEventSimulator(NaFBadEvSimulator);
	HeContaminationNaF->SetUpBadEventSimulator(NaFBadEvSimulator);

	//Agl
	
	Efficiency * HeliumFragmAgl = new Efficiency(finalHistos,"HeliumFrragmAgl","HeliumFragmentation",AglDB,"IsPreselected&IsHeliumMC&IsFromAgl&IsGoodHe" ,"IsPreselected&IsHeliumMC&LikelihoodCut&DistanceCut&IsFromAgl");
	Efficiency * HeliumFragmIntoDAgl = new Efficiency(finalHistos,"HeliumFragmIntoDAgl","HeliumFragmentation",AglDB,"IsPreselected&IsHeliumMC&LikelihoodCut&DistanceCut&IsFromNaF" ,"IsPreselected&IsHeliumMC&LikelihoodCut&DistanceCut&IsFromAgl&DeutonsMassCut");
	Efficiency * HeContaminationAgl = new Efficiency(finalHistos,"HeContAgl","HeliumFragmentation",AglDB,"IsPreselected&LikelihoodCut&DistanceCut&IsFromAgl" ,"IsPreselected&IsFromAgl&IsGoodHe");
	
	HeliumFragmAgl->SetUpBadEventSimulator(AglBadEvSimulator);
	HeliumFragmIntoDAgl->SetUpBadEventSimulator(AglBadEvSimulator);
	HeContaminationAgl->SetUpBadEventSimulator(AglBadEvSimulator);

		ParallelFiller<Efficiency *> Filler;
		Filler.AddObject2beFilled(HeliumFragmTOF,GetSmearedBetaTOF,GetSmearedBetaTOF);
		Filler.AddObject2beFilled(HeliumFragmIntoDTOF,GetSmearedBetaTOF,GetSmearedBetaTOF);
		Filler.AddObject2beFilled(HeContaminationTOF,GetSmearedBetaTOF,GetSmearedBetaTOF);
		Filler.AddObject2beFilled(HeliumFragmNaF,GetSmearedBetaRICH,GetSmearedBetaRICH);
		Filler.AddObject2beFilled(HeliumFragmIntoDNaF,GetSmearedBetaRICH,GetSmearedBetaRICH);
		Filler.AddObject2beFilled(HeContaminationNaF,GetSmearedBetaRICH,GetSmearedBetaRICH);
		Filler.AddObject2beFilled(HeliumFragmAgl,GetSmearedBetaRICH,GetSmearedBetaRICH);
		Filler.AddObject2beFilled(HeliumFragmIntoDAgl,GetSmearedBetaRICH,GetSmearedBetaRICH);
		Filler.AddObject2beFilled(HeContaminationAgl,GetSmearedBetaRICH,GetSmearedBetaRICH);
		Filler.ReinitializeAll(Refill);
		//main loop
		Filler.LoopOnMC(treeMC,vars);


	HeliumFragmTOF->Save(finalHistos);
	HeliumFragmTOF->Eval_Efficiency();
	HeliumFragmTOF->SaveResults(finalResults);


	HeliumFragmIntoDTOF->Save(finalHistos);
	HeliumFragmIntoDTOF->Eval_Efficiency();
	HeliumFragmIntoDTOF->SaveResults(finalResults);

	HeContaminationTOF->Save(finalHistos);
	HeContaminationTOF->Eval_Efficiency();

	HeContaminationTOF->ComposeEfficiency(HeliumFragmTOF);
//	HeContaminationTOF->ComposeEfficiency(HeliumFragmIntoDTOF);
	HeContaminationTOF->SaveResults(finalResults);

	

	HeliumFragmNaF->Save(finalHistos);
	HeliumFragmNaF->Eval_Efficiency();
	HeliumFragmNaF->SaveResults(finalResults);

	HeliumFragmIntoDNaF->Save(finalHistos);
	HeliumFragmIntoDNaF->Eval_Efficiency();
	HeliumFragmIntoDNaF->SaveResults(finalResults);
	
	HeContaminationNaF->Save(finalHistos);
	HeContaminationNaF->Eval_Efficiency();
	HeContaminationNaF->ComposeEfficiency(HeliumFragmNaF);
//	HeContaminationNaF->ComposeEfficiency(HeliumFragmIntoDNaF);
	HeContaminationNaF->SaveResults(finalResults);
	


	HeliumFragmAgl->Save(finalHistos);
	HeliumFragmAgl->Eval_Efficiency();
	HeliumFragmAgl->SaveResults(finalResults);
	
	HeliumFragmIntoDAgl->Save(finalHistos);
	HeliumFragmIntoDAgl->Eval_Efficiency();
	HeliumFragmIntoDAgl->SaveResults(finalResults);
	
	HeContaminationAgl->Save(finalHistos);
	HeContaminationAgl->Eval_Efficiency();
	HeContaminationAgl->ComposeEfficiency(HeliumFragmAgl);
//	HeContaminationAgl->ComposeEfficiency(HeliumFragmIntoDAgl);
	HeContaminationAgl->SaveResults(finalResults);
	


	
	
	return 0;
}


	

