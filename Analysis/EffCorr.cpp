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
#include "../include/EffCorr.h"

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

        TNtuple *treeMC = (TNtuple *)fileMC->Get("Q");
        TNtuple *treeDT = (TNtuple *)fileDT->Get("Q");


	cout<<"****************************** BINS ***************************************"<<endl;

        SetBins();

        PRB.Print();

        cout<<"**TOF**"<<endl;
        ToFDB.Print();

        cout<<"**NaF**"<<endl;
        NaFPB.Print();

        cout<<"**Agl**"<<endl;
        AglDB.Print();

        ToFPB.UseREdges();
        NaFPB.UseREdges();
        AglPB.UseREdges();

        PRB.UseREdges();


        cout<<endl;


	cout<<"****************************** VARIABLES ***************************************"<<endl;

        Variables * vars = new Variables;

	cout<<"****************************** ANALYIS ******************************************"<<endl;

	EffCorr * RICHEffCorr_NaF = new EffCorr(finalHistos,"RICHCorrection_NaF","RICH Eff. Corr",NaFPB,"ControlSample","ControlSample&IsFromNaF","","IsProtonMC");
	EffCorr * RICHEffCorr_Agl = new EffCorr(finalHistos,"RICHCorrection_Agl","RICH Eff. Corr",AglPB,"ControlSample","ControlSample&IsFromNaF","","IsProtonMC");


	RICHEffCorr_NaF -> Fill(treeMC,treeDT,vars,GetRigidity,Refill);	
	RICHEffCorr_Agl -> Fill(treeMC,treeDT,vars,GetRigidity,Refill);	


	RICHEffCorr_NaF -> Save(finalHistos);
	RICHEffCorr_Agl -> Save(finalHistos);


	RICHEffCorr_NaF -> Eval_Efficiencies();
	RICHEffCorr_Agl -> Eval_Efficiencies();
	
	RICHEffCorr_NaF -> SaveResults(finalResults);
 	RICHEffCorr_Agl -> SaveResults(finalResults);
 


	return 0;
}


	

