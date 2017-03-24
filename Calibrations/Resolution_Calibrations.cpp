#include <bitset>
#include "TROOT.h"
#include "TNtuple.h"
#include <TSpline.h>
#include "../include/binning.h"
#include "TFile.h"
#include "TH2.h"
#include "TF2.h"
#include <TVector3.h>
#include "TMath.h"
#include <TFile.h>
#include "TFile.h"
#include "TH2.h"
#include "TF2.h"
#include <TVector3.h>
#include "TMath.h"
#include "TGraphErrors.h"

#include "../include/GlobalBinning.h"

#include "../Ntuple-making/Commonglobals.cpp"
#include "../Ntuple-making/Variables.hpp"
#include "../include/filesaver.h"

#include "../include/FitError.h"
#include "../include/Resolution.h"



bool Calculate_Resolution( Resolution * Reso, bool checkfile, TTree * treeMC,  FileSaver finalHistos , bool spline = false){
	
	if(!checkfile){
                Reso->Fill(treeMC);
		Reso->Normalize();
        }

        else Reso = new Resolution(finalHistos,Reso->GetName(),Reso->GetBinning());

        if(Reso->CheckHistos()){
                Reso->Eval_Resolution();
		Reso->Save(finalHistos);
        }
	 
	finalHistos.Add(Reso->Get_Resolutions() );
   	finalHistos.Add(Reso-> Get_Means()	);
   	finalHistos.Add(Reso->Get_Sigmas()	);

	if(spline)  finalHistos.Add(Reso->ModelSigmaWithSpline());
	else 	    finalHistos.Add(Reso->ModelSigmaWithPoly());
	
   	finalHistos.writeObjsInFolder((Reso->GetName()+"/Fit Results").c_str());				
	
	return Reso->CheckHistos();
}


int main(int argc, char * argv[])
{


	cout<<"****************************** FILES OPENING ***************************************"<<endl;

	string INPUT(argv[1]);
        string OUTPUT(argv[2]);

	FileSaver finalHistos;
        finalHistos.setName(OUTPUT.c_str());
	bool checkfile = finalHistos.CheckFile();	


        TFile *fileMC =TFile::Open(INPUT.c_str());
        TTree *treeMC = (TTree *)fileMC->Get("parametri_geo");


	cout<<"****************************** BINS ***************************************"<<endl;

	SetBins();	

	PRB.Print();
	DRB.Print();

	cout<<"**TOF**"<<endl;
	ToFDB.Print();
	ToFPB.Print();

	cout<<"**NaF**"<<endl;
	NaFDB.Print();
	NaFPB.Print();

	cout<<"**Agl**"<<endl;
	AglDB.Print();
	AglPB.Print();

	ToFDB.UseBetaEdges();
	ToFPB.UseBetaEdges();
	NaFDB.UseBetaEdges();
	NaFPB.UseBetaEdges();
	AglDB.UseBetaEdges();
	AglPB.UseBetaEdges();

	DRB.UseREdges();
	PRB.UseREdges();


	cout<<endl;

	cout<<"****************************** VARIABLES ***************************************"<<endl;
	
/*	Variables * varsMC = new Variables;
        varsMC->ReadBranches(treeMC);
*/
	std::string IsPreselected = "(CUTMASK&187)==187"  ;
	std::string IsFromNaF 	  = "(RICHmask&1023)==512";
	std::string IsFromAgl     = "(RICHmask&1023)==0"  ;
	std::string Beta_gen	  = "(GenMomentum^2/(GenMomentum^2+1))^0.5";
	std::string InverseBeta_gen = "((GenMomentum^2 + 1)/GenMomentum^2)^0.5";

	cout<<"****************************** ANALYSIS ***************************************"<<endl;
	
	Resolution * RigidityResolution = new Resolution("Rigidity Resolution",PRB,(IsPreselected).c_str(),"1/R - 1/GenMomentum","GenMomentum",1000,-0.5,1.5);
	Calculate_Resolution( RigidityResolution, checkfile, treeMC, finalHistos,true);


	Resolution * BetaTOFResolution_P = new Resolution("BetaTOF Resolution (P)",ToFPB,(IsPreselected).c_str(),("1/BetaHR -" + InverseBeta_gen).c_str(),Beta_gen.c_str(),1000,-0.5,1.5);
	Calculate_Resolution( BetaTOFResolution_P, checkfile, treeMC, finalHistos);


	Resolution * BetaNaFResolution_P = new Resolution("BetaNaF Resolution (P)",NaFPB,(IsPreselected + "&&" + IsFromNaF ).c_str(),("1/BetaRICH -" + InverseBeta_gen).c_str(),Beta_gen.c_str(),250,-0.05,0.15);
	Calculate_Resolution( BetaNaFResolution_P, checkfile, treeMC, finalHistos);


	Resolution * BetaAglResolution_P = new Resolution("BetaAgl Resolution (P)",AglPB,(IsPreselected + "&&" + IsFromAgl ).c_str(),("1/BetaRICH -" + InverseBeta_gen).c_str(),Beta_gen.c_str(),250,-0.075,0.075);
	Calculate_Resolution( BetaAglResolution_P, checkfile, treeMC, finalHistos);




	return 0;
}


