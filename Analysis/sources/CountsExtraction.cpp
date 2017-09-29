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
#include "../../include/TemplateFITbetasmear.h"


void ExtractSimpleCountNr(FileSaver finalhistos, FileSaver finalResults, TNtuple* tree,Binning Bins,float (*discr_var) (Variables * vars),std::string name,std::string cut, bool refill);

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

	FileSaver finalResults;
        finalResults.setName((OUTPUT+"_Results").c_str());


        bool checkfile = finalHistos.CheckFile();

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

        ToFDB.UseBetaEdges();
        NaFDB.UseBetaEdges();
        AglDB.UseBetaEdges();

        PRB.UseREdges();


        cout<<endl;


	cout<<"****************************** VARIABLES ***************************************"<<endl;

        Variables * vars = new Variables;
	TH1F * HeContTOF=0x0;
	TH1F * HeContNaF=0x0;
	TH1F * HeContAgl=0x0;



	cout<<"****************************** ANALYIS ******************************************"<<endl;
	if(finalResults.GetFile()){
	      HeContTOF = (TH1F *) finalResults.Get("HeliumFragmentation/HeContTOF/HeContTOF_Eff");
	      HeContNaF = (TH1F *) finalResults.Get("HeliumFragmentation/HeContNaF/HeContNaF_Eff");
	      HeContAgl = (TH1F *) finalResults.Get("HeliumFragmentation/HeContAgl/HeContAgl_Eff");
	}	

	BadEventSimulator * NaFBadEvSimulator= new BadEventSimulator("IsFromNaF",22,0.8,1); 
	BadEventSimulator * AglBadEvSimulator= new BadEventSimulator("IsFromAgl",250,0.95,1); 
	

	TemplateFIT * TOFfits= new TemplateFIT("TOFfits",ToFDB,"IsPreselected&LikelihoodCut&DistanceCut",100,0.1,4);
	if((!checkfile)||Refill){
		TOFfits->Fill(treeMC,treeDT,vars,GetRecMassTOF,GetBetaTOF);
		TOFfits->DisableFit();
		TOFfits->Save(finalHistos);
	}
	else { TOFfits= new TemplateFIT(finalHistos,"TOFfits",ToFDB);
	
//		TOFfits->DisableFit();
		TOFfits->SetHeliumContamination(HeContTOF);				
		TOFfits->ExtractCounts(finalHistos);	
		TOFfits->SaveFitResults(finalResults);
	}

	TemplateFIT * NaFfits= new TemplateFIT("NaFfits",NaFDB,"IsPreselected&LikelihoodCut&DistanceCut&IsFromNaF",100,0.1,4,true,11,400);
	if((!checkfile)||Refill){
		NaFfits->SetUpBadEventSimulator(NaFBadEvSimulator);
		NaFfits->Fill(treeMC,treeDT,vars,GetRecMassRICH,GetBetaRICH);
		NaFfits->DisableFit();
		NaFfits->Save(finalHistos);
	}
	else {NaFfits= new TemplateFIT(finalHistos,"NaFfits",NaFDB,true,11,400,200);
	
		NaFfits->SetFitRange(0.6,3);
	//	NaFfits->DisableFit();
		NaFfits->SetHeliumContamination(HeContNaF);			
		NaFfits->ExtractCounts(finalHistos);
		NaFfits->SaveFitResults(finalResults);
	}
	
	
	TemplateFIT * Aglfits= new TemplateFIT("Aglfits",AglDB,"IsPreselected&LikelihoodCut&DistanceCut&IsFromAgl",100,0.1,4,true,11,110);
	if((!checkfile)||Refill){
		Aglfits->SetUpBadEventSimulator(AglBadEvSimulator);
		Aglfits->Fill(treeMC,treeDT,vars,GetRecMassRICH,GetBetaRICH);
		Aglfits->DisableFit();
		Aglfits->Save(finalHistos);
	}
	else {Aglfits= new TemplateFIT(finalHistos,"Aglfits",AglDB,true,11,110,80);
	
		Aglfits->SetFitRange(0.6,3);
	//	Aglfits->DisableFit();
		Aglfits->SetHeliumContamination(HeContAgl);
		Aglfits->ExtractCounts(finalHistos);
		Aglfits->SaveFitResults(finalResults);
	}

	ExtractSimpleCountNr(finalHistos,finalResults,treeDT,PRB,GetRigidity,"HEPCounts","IsPreselected&LikelihoodCut&DistanceCut&IsPrimary",Refill);

	return 0;
}


void ExtractSimpleCountNr(FileSaver finalhistos, FileSaver finalResults, TNtuple* tree,Binning Bins,float (*discr_var) (Variables * vars),std::string name,std::string cut, bool refill){

	TH1F * Counts;
	if(refill){
		Counts = new TH1F(name.c_str(),name.c_str(),Bins.size(),0,Bins.size());
		cout<<name.c_str()<<" Filling ... "<< endl;
		Variables * vars = new Variables;
		vars->ReadAnalysisBranches(tree);
	
		for(int i=0;i<tree->GetEntries()/FRAC;i++){
			vars->AnalysisVariablseReset();		
			UpdateProgressBar(i, tree->GetEntries()/FRAC);
			tree->GetEvent(i);
			int kbin;
			kbin = 	Bins.GetBin(discr_var(vars));
			if(ApplyCuts(cut.c_str(),vars)&&kbin>0)
				Counts->Fill(kbin);
		}
	}
	else Counts = (TH1F*) finalhistos.Get((name + "/" + name+ "/" + name).c_str());
	
      finalhistos.Add(Counts);
      finalhistos.writeObjsInFolder((name + "/" + name).c_str());

      finalResults.Add(Counts);
      finalResults.writeObjsInFolder((name + "/" + name).c_str());
      return;	
}
	


