#include <bitset>
#include "TROOT.h"
#include "TNtuple.h"
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

#include "../../include/GlobalBinning.h"

#include "../../Ntuple-making/Commonglobals.cpp"

#include "../../include/Variables.hpp"
#include "../../include/Cuts.h"


#include "../../include/filesaver.h"

#include "../../include/Efficiency.h"
#include "../../include/EffCorr.h"
#include "../../include/EffCorrTemplate.h"


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

	ToFDB.UseBetaEdges();
        NaFDB.UseBetaEdges();
        AglDB.UseBetaEdges();


        PRB.UseREdges();


        cout<<endl;


	cout<<"****************************** VARIABLES ***************************************"<<endl;

        Variables * vars = new Variables;

	cout<<"****************************** ANALYIS ******************************************"<<endl;

	BadEventSimulator * NaFBadEvSimulator= new BadEventSimulator("IsFromNaF",30,0.82,1);
        BadEventSimulator * AglBadEvSimulator= new BadEventSimulator("IsFromAgl",120,0.96,1);

	EffCorr * HEPPresEffCorr = new EffCorr(finalHistos,"HEPPresEffCorr","HEPPresEffCorr",PRB,"PresControlSample","PresControlSample&IsPreselected","","IsProtonMC");
	EffCorr * HEPQualEffCorr = new EffCorr(finalHistos,"HEPQualEffCorr","HEPQualEffCorr",PRB,"ControlSample","ControlSample&DistanceCut&LikelihoodCut","","IsProtonMC");
	
	EffCorr * RICHEffCorr_NaF = new EffCorr(finalHistos,"RICHCorrection_NaF","RICH Eff. Corr",NaFPB,"ControlSample","ControlSample&IsFromNaF","","IsProtonMC");
	EffCorr * RICHEffCorr_Agl = new EffCorr(finalHistos,"RICHCorrection_Agl","RICH Eff. Corr",AglPB,"ControlSample","ControlSample&IsFromAgl","","IsProtonMC");

	EffCorrTemplate* DistCorr_TOF = new EffCorrTemplate(finalHistos,"DistanceCorrTOF","Quality Eff. Corr",ToFDB,"ControlSample","ControlSample&DistanceCut","","");	
	EffCorrTemplate* DistCorr_NaF = new EffCorrTemplate(finalHistos,"DistanceCorrNaF","Quality Eff. Corr",NaFDB,"ControlSample&IsFromNaF","ControlSample&IsFromNaF&DistanceCut","","",true);	
	EffCorrTemplate* DistCorr_Agl = new EffCorrTemplate(finalHistos,"DistanceCorrAgl","Quality Eff. Corr",AglDB,"ControlSample&IsFromAgl","ControlSample&IsFromAgl&DistanceCut","","",true);	
	EffCorrTemplate* LikCorr_TOF = new EffCorrTemplate(finalHistos,"LikelihoodCorrTOF","Quality Eff. Corr",ToFDB,"ControlSample&DistanceCut","ControlSample&DistanceCut&LikelihoodCut","","");	
	EffCorrTemplate* LikCorr_NaF = new EffCorrTemplate(finalHistos,"LikelihoodCorrNaF","Quality Eff. Corr",NaFDB,"ControlSample&DistanceCut&IsFromNaF","ControlSample&IsFromNaF&DistanceCut&LikelihoodCut","","",true);	
	EffCorrTemplate* LikCorr_Agl = new EffCorrTemplate(finalHistos,"LikelihoodCorrAgl","Quality Eff. Corr",AglDB,"ControlSample&DistanceCut&IsFromAgl","ControlSample&IsFromAgl&DistanceCut&LikelihoodCut","","",true);	

	DistCorr_NaF->SetUpBadEventSimulator(NaFBadEvSimulator);
	DistCorr_Agl->SetUpBadEventSimulator(AglBadEvSimulator);
	LikCorr_NaF->SetUpBadEventSimulator(NaFBadEvSimulator);
	LikCorr_Agl->SetUpBadEventSimulator(AglBadEvSimulator);


	HEPPresEffCorr -> Fill(treeMC,treeDT,vars,GetRigidity,Refill);
	HEPQualEffCorr -> Fill(treeMC,treeDT,vars,GetRigidity,Refill);	
	RICHEffCorr_NaF -> Fill(treeMC,treeDT,vars,GetRigidity,Refill);	
	RICHEffCorr_Agl -> Fill(treeMC,treeDT,vars,GetRigidity,Refill);	
	DistCorr_TOF -> Fill(treeMC,treeDT,vars,GetRecMassTOF,GetBetaTOF,Refill);	
	DistCorr_NaF -> Fill(treeMC,treeDT,vars,GetRecMassRICH,GetBetaRICH,Refill);	
	DistCorr_Agl -> Fill(treeMC,treeDT,vars,GetRecMassRICH,GetBetaRICH,Refill);	
	LikCorr_TOF -> Fill(treeMC,treeDT,vars,GetRecMassTOF,GetBetaTOF,Refill);	
	LikCorr_NaF -> Fill(treeMC,treeDT,vars,GetRecMassRICH,GetBetaRICH,Refill);	
	LikCorr_Agl -> Fill(treeMC,treeDT,vars,GetRecMassRICH,GetBetaRICH,Refill);	

	HEPPresEffCorr -> Save(finalHistos);
	HEPQualEffCorr -> Save(finalHistos);
	RICHEffCorr_NaF -> Save(finalHistos);
	RICHEffCorr_Agl -> Save(finalHistos);
	DistCorr_TOF -> Save(finalHistos); 
	DistCorr_NaF -> Save(finalHistos); 
	DistCorr_Agl -> Save(finalHistos); 
	LikCorr_TOF -> Save(finalHistos); 
	LikCorr_NaF -> Save(finalHistos); 
	LikCorr_Agl -> Save(finalHistos); 


	HEPPresEffCorr -> Eval_Efficiencies();
	HEPQualEffCorr -> Eval_Efficiencies();
	RICHEffCorr_NaF -> Eval_Efficiencies();
	RICHEffCorr_Agl -> Eval_Efficiencies();
	DistCorr_TOF -> Eval_Efficiencies();
	DistCorr_NaF -> Eval_Efficiencies();
	DistCorr_Agl -> Eval_Efficiencies();
	LikCorr_TOF -> Eval_Efficiencies();
	LikCorr_NaF -> Eval_Efficiencies();
	LikCorr_Agl -> Eval_Efficiencies();


	HEPPresEffCorr -> Eval_Corrections();
	HEPQualEffCorr -> Eval_Corrections();
	RICHEffCorr_NaF -> Eval_Corrections();
	RICHEffCorr_Agl -> Eval_Corrections();
	DistCorr_TOF -> Eval_Corrections();
	DistCorr_NaF -> Eval_Corrections();
	DistCorr_Agl -> Eval_Corrections();
	LikCorr_TOF -> Eval_Corrections();
	LikCorr_NaF -> Eval_Corrections();
	LikCorr_Agl -> Eval_Corrections();

	
	HEPPresEffCorr -> SaveResults(finalResults);
	HEPQualEffCorr -> SaveResults(finalResults);
	RICHEffCorr_NaF -> SaveResults(finalResults);
 	RICHEffCorr_Agl -> SaveResults(finalResults);
	DistCorr_TOF -> SaveResults(finalResults);
        DistCorr_NaF -> SaveResults(finalResults);
 	DistCorr_Agl -> SaveResults(finalResults);
	LikCorr_TOF -> SaveResults(finalResults);
        LikCorr_NaF -> SaveResults(finalResults);
 	LikCorr_Agl -> SaveResults(finalResults);

	return 0;
}


	

