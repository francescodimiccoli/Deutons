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



int main(int argc, char * argv[]){

	cout<<"****************************** FILES OPENING ***************************************"<<endl;

	string INPUT(argv[1]);

	FileSaver finalHistos;
        finalHistos.setName(INPUT.c_str());
	bool checkfile = finalHistos.CheckFile();	


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
	
	cout<<"****************************** PLOTTING ***************************************"<<endl;





	return 0;
	
}
