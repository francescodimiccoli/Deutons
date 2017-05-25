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

#include "../include/GlobalBinning.h"

#include "../Ntuple-making/Commonglobals.cpp"

#include "../include/Variables.hpp"
#include "../include/Cuts.h"


#include "../include/filesaver.h"

#include "../include/Efficiency.h"


int main(int argc, char * argv[])
{


        cout<<"****************************** FILES OPENING ***************************************"<<endl;

        string INPUT1(argv[1]);
        string INPUT2(argv[2]);
        string OUTPUT(argv[3]);

        FileSaver finalHistos;
        finalHistos.setName(OUTPUT.c_str());
        bool checkfile = finalHistos.CheckFile();

	FileSaver finalResults;
        finalResults.setName((OUTPUT+"_Results").c_str());


        TFile *fileDT =TFile::Open(INPUT1.c_str());
        TFile *fileMC =TFile::Open(INPUT2.c_str());

        TNtuple *treeMC = (TNtuple *)fileMC->Get("grandezzesepd");
        TNtuple *treeDT = (TNtuple *)fileDT->Get("grandezzesepd");


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



	// TOF
	
	Efficiency * HeliumFragmTOF = new Efficiency(finalHistos,"HeliumFrragmTOF","HeliumFragmentation",ToFDB,"IsPreselected&IsHeliumMC&IsGoodHe" ,"IsPreselected&IsHeliumMC&LikelihoodCut&DistanceCut");
	
	Efficiency * HeliumFragmIntoDTOF = new Efficiency(finalHistos,"HeliumFragmIntoDTOF","HeliumFragmentation",ToFDB,"IsPreselected&IsHeliumMC&LikelihoodCut&DistanceCut" ,"IsPreselected&IsHeliumMC&LikelihoodCut&DistanceCut&DeutonsMassCut");

	Efficiency * HeContaminationTOF = new Efficiency(finalHistos,"HeContTOF","HeliumFragmentation",ToFDB,"IsPreselected&LikelihoodCut&DistanceCut"		,"IsPreselected&IsGoodHe");

	HeliumFragmTOF->Fill(treeMC, vars,GetBetaTOF);
	
	HeliumFragmIntoDTOF->Fill(treeMC, vars,GetBetaTOF);

	HeContaminationTOF->Fill(treeDT, vars,GetBetaTOF);

	HeliumFragmTOF->Save(finalHistos);
	HeliumFragmTOF->Eval_Efficiency();
	HeliumFragmTOF->SaveResults(finalResults);


	HeliumFragmIntoDTOF->Save(finalHistos);
	HeliumFragmIntoDTOF->Eval_Efficiency();
	HeliumFragmIntoDTOF->SaveResults(finalResults);

	HeContaminationTOF->Save(finalHistos);
	HeContaminationTOF->Eval_Efficiency();
	HeContaminationTOF->ComposeEfficiency(HeliumFragmTOF);
	HeContaminationTOF->ComposeEfficiency(HeliumFragmIntoDTOF);
	HeContaminationTOF->SaveResults(finalResults);

	
	// NaF
	
	Efficiency * HeliumFragmNaF = new Efficiency(finalHistos,"HeliumFrragmNaF","HeliumFragmentation",NaFDB,"IsPreselected&IsHeliumMC&IsFromNaF&IsGoodHe" ,"IsPreselected&IsHeliumMC&LikelihoodCut&DistanceCut&IsFromNaF");
	
	Efficiency * HeliumFragmIntoDNaF = new Efficiency(finalHistos,"HeliumFragmIntoDNaF","HeliumFragmentation",NaFDB,"IsPreselected&IsHeliumMC&LikelihoodCut&DistanceCut&IsFromNaF" ,"IsPreselected&IsHeliumMC&LikelihoodCut&DistanceCut&IsFromNaF&DeutonsMassCut");

	Efficiency * HeContaminationNaF = new Efficiency(finalHistos,"HeContNaF","HeliumFragmentation",NaFDB,"IsPreselected&LikelihoodCut&DistanceCut&IsFromNaF"  ,"IsPreselected&IsFromNaF&IsGoodHe"); 
	

	HeliumFragmNaF->Fill(treeMC, vars,GetBetaRICH);
	
	HeliumFragmIntoDNaF->Fill(treeMC, vars,GetBetaRICH);

	HeContaminationNaF->Fill(treeDT, vars,GetBetaRICH);
			

	HeliumFragmNaF->Save(finalHistos);
	HeliumFragmNaF->Eval_Efficiency();
	HeliumFragmNaF->SaveResults(finalResults);
	

	HeliumFragmIntoDNaF->Save(finalHistos);
	HeliumFragmIntoDNaF->Eval_Efficiency();
	HeliumFragmIntoDNaF->SaveResults(finalResults);
	
	HeContaminationNaF->Save(finalHistos);
	HeContaminationNaF->Eval_Efficiency();
	HeContaminationNaF->ComposeEfficiency(HeliumFragmNaF);
	HeContaminationNaF->ComposeEfficiency(HeliumFragmIntoDNaF);
	HeContaminationNaF->SaveResults(finalResults);
	

	//Agl
	
	Efficiency * HeliumFragmAgl = new Efficiency(finalHistos,"HeliumFrragmAgl","HeliumFragmentation",AglDB,"IsPreselected&IsHeliumMC&IsFromAgl&IsGoodHe" ,"IsPreselected&IsHeliumMC&LikelihoodCut&DistanceCut&IsFromAgl");
	
	Efficiency * HeliumFragmIntoDAgl = new Efficiency(finalHistos,"HeliumFragmIntoDAgl","HeliumFragmentation",AglDB,"IsPreselected&IsHeliumMC&LikelihoodCut&DistanceCut&IsFromNaF" ,"IsPreselected&IsHeliumMC&LikelihoodCut&DistanceCut&IsFromAgl&DeutonsMassCut");

	Efficiency * HeContaminationAgl = new Efficiency(finalHistos,"HeContAgl","HeliumFragmentation",AglDB,"IsPreselected&LikelihoodCut&DistanceCut&IsFromAgl" ,"IsPreselected&IsFromAgl&IsGoodHe");
	


	HeliumFragmAgl->Fill(treeMC, vars,GetBetaRICH);
	
	HeliumFragmIntoDAgl->Fill(treeMC, vars,GetBetaRICH);

	HeContaminationAgl->Fill(treeDT, vars,GetBetaRICH);

	HeliumFragmAgl->Save(finalHistos);
	HeliumFragmAgl->Eval_Efficiency();
	HeliumFragmAgl->SaveResults(finalResults);

	
	HeliumFragmIntoDAgl->Save(finalHistos);
	HeliumFragmIntoDAgl->Eval_Efficiency();
	HeliumFragmIntoDAgl->SaveResults(finalResults);
	
	HeContaminationAgl->Save(finalHistos);
	HeContaminationAgl->Eval_Efficiency();
	HeContaminationAgl->ComposeEfficiency(HeliumFragmAgl);
	HeContaminationAgl->ComposeEfficiency(HeliumFragmIntoDAgl);
	HeContaminationAgl->SaveResults(finalResults);
	


	
	
	return 0;
}


	

