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
#include "../include/TemplateFITbetasmear.h"
#include "../include/ChargeFitter.h"

void ExtractContaminationWeight(TemplateFIT * HeContTemplate, FileSaver finalResults);

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
/*
	BadEventSimulator * NaFBadEvSimulator= new BadEventSimulator("IsFromNaF",22,0.8,1);
        BadEventSimulator * AglBadEvSimulator= new BadEventSimulator("IsFromAgl",250,0.95,1);


	TemplateFIT * HeContTemplate= new TemplateFIT("HeContTemplate",ToFDB,"IsPreselected&LikelihoodCut&DistanceCut",100,0.1,4);
        if((!checkfile)||Refill){
                HeContTemplate->Fill(treeMC,treeDT,vars,GetRecMassTOF,GetBetaTOF);
                HeContTemplate->DisableFit();
                HeContTemplate->Save(finalHistos);
        }
        else { HeContTemplate= new TemplateFIT(finalHistos,"HeContTemplate",ToFDB);
                HeContTemplate->ExtractCounts(finalHistos);
                HeContTemplate->SaveFitResults(finalResults);
        	ExtractContaminationWeight(HeContTemplate,finalResults);
	}
*/	
	
	Efficiency * TestReweigthingP  = new Efficiency(finalHistos,"TestReweigthingP","TestReweigthingP",PRB,"IsProtonMC" ,"IsProtonMC");
	Efficiency * TestReweigthingHe = new Efficiency(finalHistos,"TestReweigthingHe","TestReweigthingHe",PRB,"IsHeliumMC" ,"IsHeliumMC");


	ParallelFiller<Efficiency *> Filler1;
	Filler1.AddObject2beFilled(TestReweigthingP,GetRigidity,GetRigidity);
	Filler1.AddObject2beFilled(TestReweigthingHe,GetRigidity,GetRigidity);
	Filler1.ReinitializeAll(false);
	//main loop
	Filler1.LoopOnMC(treeMC,vars);
	
	TestReweigthingP ->Save(finalHistos);
	TestReweigthingHe->Save(finalHistos);
	
	TestReweigthingP ->Save(finalResults);
	TestReweigthingHe->Save(finalResults);




	ChargeFitter * L1Charge = new ChargeFitter(finalHistos,"L1Charge","HeliumContamination",PRB,"IsPreselectedInner","InnerAndL1Charge2");
	L1Charge->Fill(treeMC,treeDT,vars,GetL1Q,GetRigidity,Refill);
	L1Charge->Save(finalHistos);
	L1Charge->ModelBefores();
	L1Charge->Eval_Efficiency();
	L1Charge->SaveResults(finalResults);

	return 0;
}


void ExtractContaminationWeight(TemplateFIT * HeContTemplate, FileSaver finalResults){
	
	TH1F * ContWeights = new TH1F("ContWeights","ContWeights",HeContTemplate->GetBinning().size(),0,HeContTemplate->GetBinning().size());
	TGraphErrors * ContWeightsFit = new TGraphErrors();
	TF1 * Model = new TF1("Model","pol3",0,1.1);

	for(int i=0;i<ContWeights->GetNbinsX();i++)
		if(HeContTemplate->GetHeContaminationWeight(i)<0.024){
			ContWeights->SetBinContent(i+1,HeContTemplate->GetHeContaminationWeight(i));
			ContWeights->SetBinError(i+1,HeContTemplate->GetHeContaminationErr(i));
			ContWeightsFit->SetPoint(i,HeContTemplate->GetBinning().BetaBinCent(i),HeContTemplate->GetHeContaminationWeight(i));
			ContWeightsFit->SetPointError(i,0,HeContTemplate->GetHeContaminationErr(i));
		}

	ContWeightsFit->SetMarkerColor(3);
	ContWeightsFit->Fit("Model","","",0.52,0.79);
	

	finalResults.Add(ContWeights);
	finalResults.Add(ContWeightsFit);
	finalResults.Add(Model);

	finalResults.writeObjsInFolder(HeContTemplate->GetName().c_str());	
	return;
}

	

