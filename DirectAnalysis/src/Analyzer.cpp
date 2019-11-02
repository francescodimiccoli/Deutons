#ifndef ANAL2_H
#define ANAL2_H


#include "Analyzer.h"

void Analyzer::FillAll(){
	Variables * vars = new Variables();
	cout<<Filler.GetObjectNr()<<" Event Objects to be filled... "<< endl;
	cout<<Filler_RTI.GetObjectNr()<<" RTI Objects to be filled... "<< endl;
	if(Filler.GetRefillFlag()&&Filler.GetObjectNr()) {
		Filler.LoopOnData(DBarReader(chainDT, false,chain_RTI,chainDT_Cpct),vars);
	//	Filler.LoopOnMC  (DBarReader(chainMC, true ,chain_RTI,chainMC_Cpct),vars);
	} 
	if(Filler_RTI.GetRefillFlag()&&Filler_RTI.GetObjectNr()) 
		Filler_RTI.ExposureTimeFilling_RTI(DBarReader(chain_RTI, false ),vars);
		//fill using data
	//	Filler_RTI.ExposureTimeFilling(DBarReader(chainDT, false,chain_RTI),vars);

}

void Analyzer::SaveAll(){

	for(int i=0;i<Filler.GetObjectNr();i++) Filler.GetObject(i)->Save();
 
}

#endif
