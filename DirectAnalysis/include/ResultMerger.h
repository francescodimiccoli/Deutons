#ifndef RESULTMERGER_H
#define RESULTMERGER_H


#include "binning.h"
#include "Globals.h"
#include "particle.h"
#include "PlottingFunctions.h"
#include "Flux.h"
#include "Acceptance.h"


class ResultMerger{

	public:
	
	TH1F* flux;
	TH1F* flux_stat;
	TH1F* flux_unf;
	TH1F* flux_unf_stat;

	TH1F* flux27;
	TH1F* flux_stat27;
	TH1F* flux_unf27;
	TH1F* flux_unf_stat27;

	TH1F* flux_ekin;
	TH1F* flux_ekin_unf;


	TH1F* effAcc;
	TH1F* effAccMC;
	TH1F* effAccMC_raw;
	TH1F* counts;
	TH1F* unfolding;
	TH1F* roounfolding;


	TH1F * statErr;
	TH1F * systErr;
	TH1F * accErr;
	TH1F * unfErr;			
	TH1F * roounfErr;			


	std::string Name;

	ResultMerger(FileSaver finalhostos, std::string name, RangeMerger Global, Flux *FluxTOF,  Flux *FluxNaF, Flux *FluxAgl, Particle particle);

	void SaveResults(FileSaver finalhistos);
	

};



#endif
