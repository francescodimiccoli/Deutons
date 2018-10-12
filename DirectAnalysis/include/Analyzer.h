#ifndef ANAL_H
#define ANAL_H

#include "Tool.h"
#include "ParallelFiller.h"
#include "filesaver.h"
#include "TChain.h"
#include "InputFileReader.h"
#include "DBarReader.h"
#include "Variables.hpp"

class Analyzer{

	typedef void (*ComputingFunction) (FileSaver finalhistos, FileSaver finalresults);
	
	private:
	TChain * chain_RTI;
        TChain * chainDT;  
        TChain * chainMC;  
	ParallelFiller<Tool *> Filler;
	std::vector<ComputingFunction> ComputingFunctions;	

	public:
	Analyzer(std::string INPUT1, std::string INPUT2){
		chain_RTI  = InputFileReader(INPUT1.c_str(),"RTI");
		chainDT    = InputFileReader(INPUT1.c_str(),"Event");
		chainMC    = InputFileReader(INPUT2.c_str(),"Event");
	};
	void FillAll();
	void BookCountsAnalysis(FileSaver finalhistos, FileSaver finalresults, bool refill);

};

#endif
