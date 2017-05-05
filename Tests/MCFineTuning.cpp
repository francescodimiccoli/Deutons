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

#include "../include/MCTuning.h"



void UpdateZoneLivetime (float Livetime, float Rcutoff, TH1F * esposizionegeo){

        for(int i=0;i<esposizionegeo->GetNbinsX();i++)
                        if(esposizionegeo->GetBinLowEdge(i+1)>=1.0*Rcutoff){
                                esposizionegeo -> SetBinContent(i+1, esposizionegeo -> GetBinContent(i+1) + Livetime) ;
        }
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

        TNtuple *treeMC = (TNtuple *)fileMC->Get("grandezzesepd");
        TNtuple *treeDT = (TNtuple *)fileDT->Get("grandezzesepd");
	TNtuple *RawDT  = (TNtuple *)fileDT->Get("trig");

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

	cout<<"************ Exposure TIME **************"<<endl;
	
	TH1F * ExposureTime;

	if(!checkfile){

		float Livetime,U_time,Rcutoff;
		RawDT->SetBranchAddress("Livetime"           ,&Livetime);
		RawDT->SetBranchAddress("U_Time"             ,&U_time);
		RawDT->SetBranchAddress("Rcutoff"            ,&Rcutoff);


		ExposureTime = new TH1F("Exposure Time","Exposure Time",500,0,10);

		vars->ReadAnalysisBranches(treeDT);

		int ActualTime=0;
		for(int i=0;i<RawDT->GetEntries();i++){
			vars->AnalysisVariablseReset();
			UpdateProgressBar(i, treeDT->GetEntries());
			RawDT->GetEvent(i);
			if((int)U_time!=ActualTime) {
				UpdateZoneLivetime(Livetime,Rcutoff,ExposureTime);
				ActualTime=U_time;
			}
		}

		finalHistos.Add(ExposureTime);
		finalHistos.Add(ExtractCutoffWeight(ExposureTime));	
		finalHistos.writeObjsInFolder("");	
	}
	
	cout<<"****************************** ANALYIS ******************************************"<<endl;

	Tuning * MCTuning;
	MCTuning = new Tuning(ToFDB,10,50,40,ExposureTime);
	if(!checkfile){
		MCTuning->UseCutoffFilterMode();
		MCTuning->Fill(treeMC,treeDT,vars);
		MCTuning->Save(finalHistos);	
	}
	else { 
		MCTuning = new Tuning(finalHistos,ToFDB,10,50,40,ExposureTime); 
		MCTuning->UseCutoffFilterMode();
	}

	MCTuning->Normalize();
	MCTuning->EvalResiduals();	
	MCTuning->SaveResults(finalHistos);

	return 0;
}



	


