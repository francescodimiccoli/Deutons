#ifndef RATIO_H
#define RATIO_H


#include "Analyzer.h"
#include "Efficiency.h"
#include "Flux.h"
#include "ResultMerger.h"

TH1F * DoRatio(TH1F* numerator, TH1F * denominator,std::string name);

class RatioMaker{

	public:
	TH1F * ratio		; 	 
        TH1F * ratio_stat	; 
        TH1F * ratio_unf 	; 
        TH1F * ratio_unf_stat   ;
        TH1F * ratio_errstat    ;
        TH1F * ratio_errsyst    ;
	TH1F * ratiounfolding   ;
	TH1F * ratioacceptance  ;
	TH1F * ratiocounts      ;
	ResultMerger * Num;
	 ResultMerger * Den;
	std::string Name;

	RatioMaker(ResultMerger * num, ResultMerger * den, std::string name);
	void DoRatioErrStat(float C);
	void DoRatioErrSyst(float C);
	
	void SaveResults(FileSaver finalhistos);
};


#endif
