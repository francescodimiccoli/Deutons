#ifndef RANGEMERGER_H
#define RANGEMERGER_H

#include "binning.h"
#include <fstream>
#include <vector>
#include "string.h"
#include "TH1F.h"
#include "TF1.h"
#include <iostream>
#include <fstream>
#include <sstream>


class RangeMerger{

	private:
	Binning Global_P;
	Binning Global_D;

	Binning ToFD;
	Binning NaFD;
	Binning AglD;
	Binning ToFP;
	Binning NaFP;
        Binning AglP;

	public:
	RangeMerger(){
		
		Particle proton(0.9382720813, 1, 1);  // proton mass 938 MeV
		Particle deuton(1.8756129   , 1, 2);  // deuterium mass 1876 MeV, Z=1, A=2

		 Binning       global_P(proton);
                 Binning       global_D(deuton);
                 Binning       toFD(deuton);
                 Binning       naFD(deuton);
                 Binning       aglD(deuton);
                 Binning       toFP(proton);
                 Binning       naFP(proton);
                 Binning       aglP(proton);

		Global_P = global_P;
                Global_D = global_D;
		ToFD =     toFD;
		NaFD =     naFD;
		AglD =     aglD;
		ToFP =     toFP;
		NaFP = 	   naFP;
		AglP =     aglP;
	}

	Binning GetToFDBins() {return ToFD;}
	Binning GetNaFDBins() {return NaFD;}
	Binning GetAglDBins() {return AglD;}
	Binning GetToFPBins() {return ToFP;}
	Binning GetNaFPBins() {return NaFP;}
	Binning GetAglPBins() {return AglP;}
	Binning GetGlobalPBins() {return Global_P;}
	Binning GetGlobalDBins() {return Global_D;}


	void Reset();
	void setBinsFromRDatacard(std::string datacard, TF1 * ResponseTOF, TF1 * ResponseNaF, TF1 * ResponseAgl);
	
	void UseBetaEdges();
	void UseREdges();
	void UseRTOIEdges();
	void UseBetaTOIEdges();
	
	int GetGlobalBinRTOICenter(int n);
	int GetGlobalBinEkinTOICenter(int n);
	int GetGlobalBinBetaTOICenter(int n);
	int GetGlobalBinRCenter(int n);
	int GetGlobalBinEkinCenter(int n);
	int GetGlobalBinBetaCenter(int n);

	int GetToFBinP(int globalbin);
	int GetNaFBinP(int globalbin);
	int GetAglBinP(int globalbin);
	int GetToFBinD(int globalbin);
	int GetNaFBinD(int globalbin);
	int GetAglBinD(int globalbin);


	TH1F * MergeSubDResult_P(TH1F * ResultTOF, TH1F * ResultNaF, TH1F * ResultAgl);
	TH1F * MergeSubDResult_D(TH1F * ResultTOF, TH1F * ResultNaF, TH1F * ResultAgl);

	TH1F * MergedRatio(TH1F * Result_D, TH1F * Result_P);	
	TH1F * MergedRatio_Ekin(TH1F * Result_D, TH1F * Result_P);	


};

#endif

