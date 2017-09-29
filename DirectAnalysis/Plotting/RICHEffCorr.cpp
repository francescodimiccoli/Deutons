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
#include "TCanvas.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "../include/GlobalBinning.h"
#include "TKey.h"
#include "TFractionFitter.h"

#include "../../Ntuple-making/Commonglobals.cpp"
#include "../include/Variables.hpp"
#include "../include/Cuts.h"
#include "../include/filesaver.h"

#include "../include/FitError.h"
#include "../include/Resolution.h"
#include "../include/Efficiency.h"

#include "../include/PlottingFunctions.h"

#include "../include/EffCorr.h"


int main(int argc, char * argv[]){


cout<<"****************************** FILES OPENING ***************************************"<<endl;

        string INPUT(argv[1]);
        string OUTPUT(argv[2]);

        FileSaver finalHistos;
        FileSaver Plots;

        finalHistos.setName(INPUT.c_str());
        Plots.setName(OUTPUT.c_str());


        bool checkfile = finalHistos.CheckFile();


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

	cout<<"**************************** PLOTTING ***************************************"<<endl;

	EffCorr * RICHEffCorr_NaF = new EffCorr(finalHistos,"RICHCorrection_NaF","RICH Eff. Corr",NaFPB,"ControlSample","ControlSample&IsFromNaF","","IsProtonMC");
	EffCorr * RICHEffCorr_Agl = new EffCorr(finalHistos,"RICHCorrection_Agl","RICH Eff. Corr",AglPB,"ControlSample","ControlSample&IsFromAgl","","IsProtonMC");

	TCanvas * c3 = new TCanvas("Latitude Corrections"); 
        c3->SetCanvasSize(2000,1500);

	for(int lat=0;lat<10;lat++)
		PlotTH1FintoGraph(gPad,NaFDB, (TH1F*)RICHEffCorr_NaF->GetCorrectionLat(lat),"Kinetic Energy [GeV/nucl.]", "Efficiency Correction Factor",55+5*lat,true,"Psame",0.2,4.5,0.1,1.45,("Lat zone " + to_string(lat)).c_str());

	TCanvas * c4 = new TCanvas("Global Corrections"); 
        c4->SetCanvasSize(2000,1500);

	PlotTH1FintoGraph(gPad,NaFDB, (TH1F*)RICHEffCorr_NaF->GetGlobCorrection(),"Kinetic Energy [GeV/nucl.]", "Efficiency Correction Factor",2,true,"Psame",0.2,4.5,0.1,1.45,"RICH NaF Efficiency Correction");


	
	Plots.Add(c3);
	Plots.Add(c4);
	Plots.writeObjsInFolder("Correction for NaF");
	
	TCanvas * c5 = new TCanvas("Latitude Corrections"); 
        c5->SetCanvasSize(2000,1500);

	for(int lat=0;lat<10;lat++)
		PlotTH1FintoGraph(gPad,AglDB, (TH1F*)RICHEffCorr_Agl->GetCorrectionLat(lat),"Kinetic Energy [GeV/nucl.]", "Efficiency Correction Factor",55+5*lat,true,"Psame",2,10,0.1,1.45,("Lat zone " + to_string(lat)).c_str());

	TCanvas * c6 = new TCanvas("Global Corrections"); 
        c6->SetCanvasSize(2000,1500);

	PlotTH1FintoGraph(gPad,AglDB, (TH1F*)RICHEffCorr_Agl->GetGlobCorrection(),"Kinetic Energy [GeV/nucl.]", "Efficiency Correction Factor",2,true,"Psame",2,10,0.1,1.45,"RICH Agl Efficiency Correction");


	
	Plots.Add(c5);
	Plots.Add(c6);
	Plots.writeObjsInFolder("Correction for Agl");


	
	return 0;
}

