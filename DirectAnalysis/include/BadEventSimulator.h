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
	TF1 * BadEvPDF;
	TRandom3 *R;

	public:
	BadEventSimulator(std::string cut, int Freq, float min, float max) {
			simulatorcut = cut; 
			BadEventsFrequency=Freq; 
			badrangemin=min; 
			badrangemax=max;
			BadEvPDF=new TF1("Noise Model","[0]*x+[1]",badrangemin,badrangemax);
			R=new TRandom3();
	};

	void SetFrequency(int freq){BadEventsFrequency = freq;}	
	void  LoadEvent(Variables * vars) {if(ApplyCuts(simulatorcut,vars)) EvNum++; return;};
	float SimulateBadEvents(float Beta) {
		if(EvNum%BadEventsFrequency==0) { 
			BadEvPDF->SetParameter(0,1/(Beta*(Beta-badrangemin)));
			BadEvPDF->SetParameter(1,-badrangemin/(Beta*(Beta-badrangemin)));
			float x=Beta;
			float y=2;
			int i=0;
			while(BadEvPDF->Eval(x)<y){
				x=R->Uniform((double)badrangemin,(double)Beta);
				float y=R->Rndm();
				if(BadEvPDF->Eval(x)>y) { return x;}
			}
			return x;	
		}
		return Beta;
	}
};



#endif
