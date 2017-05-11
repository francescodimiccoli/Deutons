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
        NaFDB.Print();

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

	TemplateFIT * TOFfits= new TemplateFIT("TOFfits",ToFDB,"IsPreselected&LikelihoodCut&DistanceCut",100,0.1,4);
	if(!checkfile){
		TOFfits->Fill(treeMC,treeDT,vars,GetRecMassTOF,GetBetaTOF);
		TOFfits->DisableFit();
		TOFfits->Save(finalHistos);
	}
	else { TOFfits= new TemplateFIT(finalHistos,"TOFfits",ToFDB);
	
		TOFfits->Save(finalHistos);
		TOFfits->ExtractCounts(finalHistos);	
		TOFfits->SaveFitResults(finalHistos);
	}

/*	TemplateFIT * NaFfits= new TemplateFIT("NaFfits",NaFDB,"IsPreselected&LikelihoodCut&DistanceCut&IsFromNaF",100,0.1,4);
	if(!checkfile){
		NaFfits->Fill(treeMC,treeDT,vars,GetRecMassRICH,GetBetaRICH);
		NaFfits->DisableFit();
		NaFfits->Save(finalHistos);
	}
	else {NaFfits= new TemplateFIT(finalHistos,"NaFfits",NaFDB);
	
		NaFfits->Save(finalHistos);
		NaFfits->SetSystematicParameters(3,0.04,0.04);
		NaFfits->SetFitRange(0.7,2.5);
		NaFfits->ExtractCounts(finalHistos);
		NaFfits->SaveFitResults(finalHistos);
	}

	TemplateFIT * Aglfits= new TemplateFIT("Aglfits",AglDB,"IsPreselected&LikelihoodCut&DistanceCut&IsFromAgl",100,0.1,4);
	if(!checkfile){
		Aglfits->Fill(treeMC,treeDT,vars,GetRecMassRICH,GetBetaRICH);
		Aglfits->DisableFit();
		Aglfits->Save(finalHistos);
	}
	else {Aglfits= new TemplateFIT(finalHistos,"Aglfits",AglDB);
	
		Aglfits->Save(finalHistos);
		Aglfits->SetSystematicParameters(3,0.04,0.04);
		Aglfits->SetFitRange(0.7,2.5);
		Aglfits->ExtractCounts(finalHistos);
		Aglfits->SaveFitResults(finalHistos);
	}
*/
	return 0;
}



	


