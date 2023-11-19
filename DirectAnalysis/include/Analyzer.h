#ifndef ANAL_H
#define ANAL_H

#include "Tool.h"
#include "ParallelFiller.h"
#include "filesaver.h"
#include "TChain.h"
#include "InputFileReader.h"
#include "DBarReader.h"
#include "Flux.h"
#include "Globals.h"


class Analyzer{

	typedef void (*ComputingFunction) (FileSaver finalhistos, FileSaver finalresults);
	
	private:
	TChain * chain_RTI;
        TChain * chainDT;  
        TChain * chainMC;  
	TChain * chainDT_Cpct;  
        TChain * chainMC1_Cpct;  
        TChain * chainMC2_Cpct;  
        TChain * chainMC3_Cpct;  
        TChain * chainMC4_Cpct;  





	int timeindex=0;
	std::string filelistMC1;	
	std::string filelistMC2;	
	std::string filelistMC3;	
	std::string filelistMC4;	
	std::string filelistDT;	
	ParallelFiller<Acceptance *> FillerAcc;
	ParallelFiller<Tool *> Filler;
	ParallelFiller<Flux *> Filler_RTI;
	std::vector<ComputingFunction> ComputingFunctions;	
	bool check_file = false;

	public:
	Analyzer(std::string INPUT1, std::string INPUT2){
		chain_RTI  	= InputFileReader(INPUT1.c_str(),"RTI");
		chainDT    	= InputFileReader(INPUT1.c_str(),"Event");
		chainMC    	= InputFileReader((INPUT2+"_P").c_str(),"Event");
		chainDT_Cpct    = InputFileReader(INPUT1.c_str(),"Compact");
		chainMC1_Cpct    = InputFileReader((INPUT2+"_P").c_str(),"Compact");
		chainMC2_Cpct    = InputFileReader((INPUT2+"_D").c_str(),"Compact");
		chainMC3_Cpct    = InputFileReader((INPUT2+"_He").c_str(),"Compact");
		chainMC4_Cpct    = InputFileReader((INPUT2+"_He3").c_str(),"Compact");
		timeindex = FindTimeIndex(INPUT1);
		filelistMC1 = (INPUT2+"_P").c_str();
		filelistMC2 = (INPUT2+"_D").c_str();
		filelistMC3 = (INPUT2+"_He").c_str();
		filelistMC4 = (INPUT2+"_He3").c_str();
		filelistDT = INPUT1.c_str();

	};
	bool CheckFile() {return check_file;}
	void FillAll();
	void SaveAll();

	void BookCountsAnalysis(FileSaver finalhistos, FileSaver finalresults, bool refill);
	void BookAcceptanceAnalysis(FileSaver finalhistos, FileSaver finalresults, bool refill);
	void BookFluxAnalysis(FileSaver finalhistos, FileSaver finalresults, bool refill);
	void BookTestCascade(FileSaver finalhistos, FileSaver finalresults, bool refill);
};


#endif
