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
#include "../include/AllRangesEfficiency.h"
#include "../include/Flux.h"


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

        TNtuple *treeMC = (TNtuple *)fileMC->Get("grandezzesepd");
        TNtuple *treeDT = (TNtuple *)fileDT->Get("grandezzesepd");
	TNtuple *RawDT = (TNtuple *)fileDT->Get("trig");



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

	cout<<"****************************** FLUXES EVALUATION ******************************************"<<endl;

	Flux * DFluxTOF = new Flux(finalHistos,finalResults, "DFluxTOF", "FullsetEff_D_TOF","FullsetEfficiency","TOFfits/Fit Results/Primary Deuteron Counts","ExposureTOF","Gen. Acceptance",ToFDB);
	Flux * DFluxNaF = new Flux(finalHistos,finalResults, "DFluxNaF", "FullsetEff_D_NaF","FullsetEfficiency","NaFfits/Fit Results/Primary Deuteron Counts","ExposureNaF","Gen. Acceptance",NaFDB);
	Flux * DFluxAgl = new Flux(finalHistos,finalResults, "DFluxAgl", "FullsetEff_D_Agl","FullsetEfficiency","Aglfits/Fit Results/Primary Deuteron Counts","ExposureAgl","Gen. Acceptance",AglDB);

	DFluxTOF-> Eval_ExposureTime(RawDT,finalHistos);
        DFluxNaF-> Eval_ExposureTime(RawDT,finalHistos);
        DFluxAgl-> Eval_ExposureTime(RawDT,finalHistos);

	DFluxTOF->Set_MCPar(0.5,20,0.0242236931);	
	DFluxNaF->Set_MCPar(0.5,20,0.0242236931);	
	DFluxAgl->Set_MCPar(0.5,20,0.0242236931);	

	DFluxTOF-> Eval_GeomAcceptance(treeMC,finalHistos,"IsDeutonMC");
        DFluxNaF-> Eval_GeomAcceptance(treeMC,finalHistos,"IsDeutonMC");
        DFluxAgl-> Eval_GeomAcceptance(treeMC,finalHistos,"IsDeutonMC");

	DFluxTOF-> Eval_Flux();
        DFluxNaF-> Eval_Flux();
        DFluxAgl-> Eval_Flux();

	DFluxTOF->SaveResults(finalResults);
	DFluxNaF->SaveResults(finalResults);
	DFluxAgl->SaveResults(finalResults);



	return 0;
}
