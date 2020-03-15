#ifndef LATREW_H
#define LATREW_H

#include <TF1.h>
#include "TH1.h"
#include "Globals.h"
#include "GlobalPaths.h"
#include "binning.h"
#include "Livetime.h"
#include "filesaver.h"
//#include "DBarReader.h"
#include "filesaver.h"
//#include "Variables.hpp"


class LatReweighter{

	private:
	TH1F * AllOrbit;
	TH1F * HighLat;
	TH1F * Weights;
	TH1F * ExposureTime;
	TH1F * Reference;
	TF1  * WeightModel;	
	std::string cut;
	std::string basename;
	Binning Bins;


	public: 
	LatReweighter(std::string Basename, std::string Cut,int Nbins=500,int xmin=0,int xmax=150){

		Binning bins(proton);
		Bins = bins;
		Bins.Reset();
		Bins.setBinsFromRigidity(100,0.5,50,ResponseTOF,0.00347548,5.8474);
		Bins.UseREdges();
		Bins.Print();

		AllOrbit = new TH1F((Basename+"_AllOrbit").c_str(),(Basename+"_AllOrbit").c_str(),Nbins,xmin,xmax);
		HighLat = new TH1F((Basename+"_HighLat").c_str(),(Basename+"_HighLat").c_str(),Nbins,xmin,xmax);
		Weights = new TH1F((Basename+"_Weights").c_str(),(Basename+"_Weights").c_str(),Nbins,xmin,xmax);
		ExposureTime = new TH1F ((Basename+"_Weights").c_str(),(Basename+"_Weights").c_str(),Bins.size(),0,Bins.size());
		basename = Basename;
		cut = Cut;
	};

	LatReweighter(FileSaver FinalHisto,std::string Basename);
/*	void LoopOnData(TTree * treeDT, Variables * vars, bool Refill);
	void LoopOnData(DBarReader readerDT, Variables * vars, bool Refill);
	void LoopOnRTI(DBarReader readerDT, Variables * vars,bool Refill);
*/	void Save(FileSaver finalhisto);	
	void CalculateWeights();
	float GetCleaningWeight(float Rmeas, float Rgen, float SF);
	float GetTimeDepWeight( float R);
	void ModelWithSpline();
	void SaveResults(FileSaver finalhisto);	

};

#endif
