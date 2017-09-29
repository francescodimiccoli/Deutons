#ifndef BADEVENT_H
#define BADEVENT_H

#include "Cuts.h"
#include "Variables.hpp"

class BadEventSimulator{
	private:
	std::string simulatorcut;
	int BadEventsFrequency;
	float badrangemin;
	float badrangemax;
	int EvNum=0;
	TF1 * BadEvModel;

	public:
	BadEventSimulator(std::string cut, int Freq, float min, float max) {
			simulatorcut = cut; 
			BadEventsFrequency=Freq; 
			badrangemin=min; 
			badrangemax=max;
			BadEvModel=new TF1("Linear Model","[0]*x+[1]",badrangemin,badrangemax);
			BadEvModel->SetParameter(0,1/(badrangemax-badrangemin));
			BadEvModel->SetParameter(1,-badrangemin/(badrangemax-badrangemin));
	};
	
	void   LoadEvent(Variables * vars) {if(ApplyCuts(simulatorcut,vars)) EvNum++; return;};
	float SimulateBadEvents(float Beta) {
		if(EvNum%BadEventsFrequency==0) { float betabad=BadEvModel->GetRandom(); return betabad;}
	        return Beta;

	}
};


#endif
