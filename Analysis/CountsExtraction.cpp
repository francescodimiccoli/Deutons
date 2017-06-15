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
#include "../include/TemplateFITbetasmear.h"




int main(int argc, char * argv[])
{


        cout<<"****************************** FILES OPENING ***************************************"<<endl;

        string INPUT1(argv[1]);
        string INPUT2(argv[2]);
        string OUTPUT(argv[3]);

        FileSaver finalHistos;
        finalHistos.setName(OUTPUT.c_str());

	FileSaver finalResults;
        finalResults.setName((OUTPUT+"_Results").c_str());


        bool checkfile = finalHistos.CheckFile();

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
	
	TemplateFIT * TOFfits= new TemplateFIT("TOFfits","HeContTOF",ToFDB,"IsPreselected&LikelihoodCut&DistanceCut",100,0.1,4);
	if(!checkfile){
		TOFfits->Fill(treeMC,treeDT,vars,GetRecMassTOF,GetBetaTOF);
		TOFfits->DisableFit();
		TOFfits->Save(finalHistos);
	}
	else { TOFfits= new TemplateFIT(finalHistos,"TOFfits","HeContTOF",ToFDB);
	
	//	TOFfits->Save(finalHistos);
		TOFfits->ExtractCounts(finalHistos,finalResults);	
		TOFfits->SaveFitResults(finalResults);
	}

	TemplateFIT * NaFfits= new TemplateFIT("NaFfits","HeContNaF",NaFDB,"IsPreselected&LikelihoodCut&DistanceCut&IsFromNaF",100,0.1,4,true,11,400);
	if(!checkfile){
		NaFfits->Fill(treeMC,treeDT,vars,GetRecMassRICH,GetBetaRICH);
		NaFfits->DisableFit();
		NaFfits->Save(finalHistos);
	}
	else {NaFfits= new TemplateFIT(finalHistos,"NaFfits","HeContNaF",NaFDB,true,11,400,200);
	
	//	NaFfits->Save(finalHistos);
		NaFfits->SetFitRange(0.6,3);
		NaFfits->ExtractCounts(finalHistos,finalResults);
		NaFfits->SaveFitResults(finalResults);
	}
	
	
	TemplateFIT * Aglfits= new TemplateFIT("Aglfits","HeContAgl",AglDB,"IsPreselected&LikelihoodCut&DistanceCut&IsFromAgl",100,0.1,4,true,11,110);
	if(!checkfile){
		Aglfits->Fill(treeMC,treeDT,vars,GetRecMassRICH,GetBetaRICH);
		Aglfits->DisableFit();
		Aglfits->Save(finalHistos);
	}
	else {Aglfits= new TemplateFIT(finalHistos,"Aglfits","HeContAgl",AglDB,true,11,110,80);
	
	//	Aglfits->Save(finalHistos);
		Aglfits->SetFitRange(0.6,3);
		Aglfits->ExtractCounts(finalHistos,finalResults);
		Aglfits->SaveFitResults(finalResults);
	}

	return 0;
}



	


