#ifndef TOOL_H
#define TOOL_H

#include "filesaver.h"
#include "ParallelFiller.h"

class Tool{

	protected:
	FileSaver finalhistos;

	public:
	virtual void SetDefaultOutFile(FileSaver FinalHistos){finalhistos = FinalHistos; return;}
	virtual void SetUpBadEventSimulator(BadEventSimulator * Sim) {return;};
	virtual void LoadEventIntoBadEvSim(Variables * vars) { return;}
	virtual bool ReinitializeHistos(bool refill){ return refill;}
	virtual void FillEventByEventMC(Variables * vars, float (*var) (Variables * vars), float (*discr_var) (Variables * vars)){ return; }
	virtual void FillEventByEventData(Variables * vars, float (*var) (Variables * vars), float (*discr_var) (Variables * vars)){ return;}
	virtual void Save() {return;}
	virtual void SaveResults(FileSaver finalhistos){return;}

};
#endif
