#ifndef ANAL2_H
#define ANAL2_H


#include "Analyzer.h"

void Analyzer::FillAll(){
	Variables * vars = new Variables();
	cout<<Filler.GetObjectNr()<<" Event Objects to be filled... "<< endl;
	cout<<Filler_RTI.GetObjectNr()<<" RTI Objects to be filled... "<< endl;
	if(Filler.GetRefillFlag()&&Filler.GetObjectNr()) {
		Filler.LoopOnMC  (DBarReader(chainMC, true ),vars);
		Filler.LoopOnData(DBarReader(chainDT, false,chain_RTI),vars);
	} 
	if(Filler_RTI.GetRefillFlag()&&Filler_RTI.GetObjectNr()) 
		Filler_RTI.ExposureTimeFilling_RTI(DBarReader(chain_RTI, false ),vars);

}

void Analyzer::SaveAll(){

	for(int i=0;i<Filler.GetObjectNr();i++) Filler.GetObject(i)->Save();
 
}

#endif