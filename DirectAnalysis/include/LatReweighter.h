#ifndef LATREW_H
#define LATREW_H

#include <TF1.h>
#include "Globals.h"
#include "BadEventSimulator.h"
#include "binning.h"
#include "Livetime.h"
#include "filesaver.h"
#include "DBarReader.h"
#include "filesaver.h"



class LatReweighter{

	private:
	TH1F * AllOrbit;
	TH1F * HighLat;
	TH1F * Weights;
	std::string cut;
	std::string basename;


	public: 
	LatReweighter(std::string Basename, std::string Cut,int Nbins=500,int xmin=0,int xmax=150){

		AllOrbit = new TH1F((Basename+"_AllOrbit").c_str(),(Basename+"_AllOrbit").c_str(),Nbins,xmin,xmax);
		HighLat = new TH1F((Basename+"_HighLat").c_str(),(Basename+"_HighLat").c_str(),Nbins,xmin,xmax);
		Weights = new TH1F((Basename+"_Weights").c_str(),(Basename+"_Weights").c_str(),Nbins,xmin,xmax);
		basename = Basename;
		cut = Cut;
	};

	LatReweighter(FileSaver FinalHisto,std::string Basename);
	void LoopOnData(TTree * treeDT, Variables * vars, bool Refill);
	void LoopOnData(DBarReader readerDT, Variables * vars, bool Refill);
	void Save(FileSaver finalhisto);	
	void CalculateWeights();
	float GetWeight( float R);
};

#endif
