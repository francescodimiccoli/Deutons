#ifndef ANAL2_H
#define ANAL2_H


#include "Analyzer.h"

void Analyzer::FillAll(){
	Variables * vars = new Variables(timeindex);
	cout<<Filler.GetObjectNr()<<" Event Objects to be filled... "<< endl;
	cout<<Filler_RTI.GetObjectNr()<<" RTI Objects to be filled... "<< endl;
	cout<<FillerAcc.GetObjectNr()<<" Acceptance Objects to be filled... "<< endl;
	if(Filler.GetRefillFlag()&&Filler.GetObjectNr()) {
		//Filler.LoopOnData(DBarReader(chainDT, false,chain_RTI,chainDT_Cpct),vars);
		Filler.LoopOnMC  (DBarReader(chainMC, true ,chain_RTI,chainMC1_Cpct),DBarReader(chainMC, true ,chain_RTI,chainMC2_Cpct),DBarReader(chainMC, true ,chain_RTI,chainMC3_Cpct),vars);
	}
	 
	if(Filler_RTI.GetRefillFlag()&&Filler_RTI.GetObjectNr()) 
		Filler_RTI.ExposureTimeFilling_RTI(DBarReader(chain_RTI, false ),vars);
		
//	if(FillerAcc.GetRefillFlag()&&FillerAcc.GetObjectNr()) 
//		FillerAcc.LoopOnMCForGenAcceptance(DBarReader(chainMC, true ,chain_RTI,chainMC_Cpct),vars);
	
}

void Analyzer::SaveAll(){

	for(int i=0;i<Filler.GetObjectNr();i++) Filler.GetObject(i)->Save();
 
}

#endif
