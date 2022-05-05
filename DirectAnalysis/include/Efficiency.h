#ifndef EFFICIENCY_H
#define EFFICIENCY_H

#include "BetaSmearing.h"
#include "Cuts.h"
#include "DBarReader.h"
#include "TGraphErrors.h"
#include "filesaver.h"
#include "Tool.h"
#include "GlobalPaths.h"
#include "PlottingFunctions.h"


class Tool;

class Efficiency : public Tool{

	private:

	bool notweighted_phtr = false;
	bool notweighted = false;
	TH1 * before;
	TH1 * after;
	TH1 * Eff;

	TGraphErrors * FittedEff;
	FileSaver  file;
	Binning bins;
	std::string cut_before;
	std::string cut_after;
	std::string basename;
	std::string directory;
	bool fitrequested=false;
	bool IsExtern = false;
	BadEventSimulator * BadEvSim=0x0;

	TH1F * Stat_Error=0x0;
	TH1F * Syst_Error=0x0;

	bool Refill = true;
	public:

	Efficiency(FileSaver  File, std::string Basename,std::string Directory, Binning Bins, std::string Cut_before,std::string Cut_after,bool external = false){
		file=File;
		bins=Bins;
		basename=Basename;
		cut_before=Cut_before;
		cut_after=Cut_after;
		directory=Directory;
		IsExtern = external;

		ReadFile(external);		
		before = new TH1F((basename+"_before").c_str(),(basename+"_before").c_str(),bins.size(),0,bins.size());
		after  = new TH1F((basename+"_after" ).c_str(),(basename+"_after" ).c_str(),bins.size(),0,bins.size());
	};

	Efficiency(FileSaver  File, std::string Basename,std::string Directory, Binning Bins, std::string Cut_before,std::string Cut_after,std::vector<float> LatZones, bool external = false){
		file=File;
		bins=Bins;
		basename=Basename;
		cut_before=Cut_before;
		cut_after=Cut_after;
		directory=Directory;
		IsExtern = external;

	
		ReadFile(external);
		before = new TH2F((basename+"_before").c_str(),(basename+"_before").c_str(),bins.size(),0,bins.size(),LatZones.size(),0,LatZones.size());
		after  = new TH2F((basename+"_after" ).c_str(),(basename+"_after" ).c_str(),bins.size(),0,bins.size(),LatZones.size(),0,LatZones.size());
	};

	Efficiency(FileSaver  File, std::string Basename,std::string Directory, Binning Bins,bool external = false){
		file=File;
		bins=Bins;
		basename=Basename;
		directory=Directory;
		IsExtern = external;
		ReadFile(external);
	}

	void ReadFile(bool external = false){
		TFile * ff;
		if(external) { 
			ff = TFile::Open((outdir+"/ExternalMCEff.root").c_str());
			cout<<"*********************** EXTERNAL MC FILE FOUND: "<<ff<<"***********************"<<endl;
			if(!ff) ff = file.GetFile();
		}

		else ff = file.GetFile();
		cout<<ff<<endl;
		if(ff){
			before =(TH1 *) ff->Get((directory+"/"+basename+"/"+basename+"_before").c_str());
                	after  =(TH1 *) ff->Get((directory+"/"+basename+"/"+basename+"_after").c_str());
                	Eff    =(TH1 *) ff->Get((directory+"/"+basename+"/"+basename+"_Eff").c_str());
        		Stat_Error = (TH1F *) ff->Get((directory+"/"+basename+"/"+basename+"_Stat_Error").c_str());
			Syst_Error = (TH1F *) ff->Get((directory+"/"+basename+"/"+basename+"_Syst_Error").c_str());
		}
		
	}
	
	virtual bool ReinitializeHistos(bool refill);
	virtual void Fill( TTree * treeMC, Variables * vars, float (*discr_var) (Variables * vars),bool refill=false);
	virtual void Fill(DBarReader readerMC, Variables * vars, float (*discr_var) (Variables * vars),bool refill);
	virtual void FillEventByEventMC(Variables * vars, float (*var) (Variables * vars), float (*discr_var) (Variables * vars));
	virtual void FillEventByEventData(Variables * vars, float (*var) (Variables * vars), float (*discr_var) (Variables * vars));
	virtual void SetNotWeightedMC(){notweighted = true; return;}
	virtual void SetNotWeightedMC_phtr(){notweighted_phtr = true; return;}


	virtual void SetUpBadEventSimulator(BadEventSimulator * Sim) {BadEvSim = Sim; return; };
	BadEventSimulator * GetBadEventSimulator() {return BadEvSim;};
	virtual void LoadEventIntoBadEvSim(Variables * vars) {if(BadEvSim) BadEvSim->LoadEvent(vars);}	
	
	virtual void Save();
	virtual void SaveResults(FileSaver finalhistos);
	virtual void Eval_Efficiency();
	virtual void Eval_TrigEfficiency();
	
	virtual void Eval_FittedEfficiency();

	virtual void CloneEfficiency(Efficiency * Second);
	virtual void ComposeEfficiency(Efficiency * Second);
	virtual void Eval_StatError();
	virtual void Eval_SystError(Efficiency * First,  Efficiency * Second);
	
	Binning GetBins() {return bins;}	
	TH1 * GetEfficiency() {return Eff;}
	TGraphErrors * GetFittedEfficiency() {return FittedEff;}
	TH1 * GetBefore()     {return before;}
	TH1 * GetAfter()      {return after;}

	TH1F * GetSyst_Error() {return Syst_Error;}
	TH1F * GetStat_Error() {return Stat_Error;}

	std::string GetCut_Before()	{return cut_before;}	
	std::string GetCut_After()	{return cut_after;}
};

#endif

