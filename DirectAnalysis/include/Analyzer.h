#ifndef ANAL_H
#define ANAL_H

#include "Tool.h"
#include "ParallelFiller.h"
#include "filesaver.h"
#include "TChain.h"
#include "InputFileReader.h"
#include "DBarReader.h"
#include "Flux.h"


class Analyzer{

	typedef void (*ComputingFunction) (FileSaver finalhistos, FileSaver finalresults);
	
	private:
	TChain * chain_RTI;
        TChain * chainDT;  
        TChain * chainMC;  
	ParallelFiller<Tool *> Filler;
	ParallelFiller<Flux *> Filler_RTI;
	std::vector<ComputingFunction> ComputingFunctions;	
	bool check_file = false;

	public:
	Analyzer(std::string INPUT1, std::string INPUT2){
		chain_RTI  = InputFileReader(INPUT1.c_str(),"RTI");
		chainDT    = InputFileReader(INPUT1.c_str(),"Event");
		chainMC    = InputFileReader(INPUT2.c_str(),"Event");
	};
	bool CheckFile() {return check_file;}
	void FillAll();
	void SaveAll();

	void BookCountsAnalysis(FileSaver finalhistos, FileSaver finalresults, bool refill);
	void BookEfficiencyAnalysis(FileSaver finalhistos, FileSaver finalresults, bool refill);
	void BookEffCorrAnalysis(FileSaver finalhistos, FileSaver finalresults, bool refill);
	void BookFluxAnalysis(FileSaver finalhistos, FileSaver finalresults, bool refill);

};

#endif
