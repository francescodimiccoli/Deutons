#ifndef ANAL2_H
#define ANAL2_H


#include "Analyzer.h"

void Analyzer::FillAll(){
        Variables * vars = new Variables();
	Filler.LoopOnMC  (DBarReader(chainMC, true ),vars);
        Filler.LoopOnData(DBarReader(chainDT, false,chain_RTI),vars);
        
}

#endif
