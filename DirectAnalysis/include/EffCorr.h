#ifndef EFFCORR_H
#define EFFCORR_H

#include "Tool.h"
#include "Efficiency.h"


class EffCorr : public Tool{

	private:

	Efficiency * EffMC;
	Efficiency * EffMC2;
	Efficiency * EffMCpid;

	Efficiency * EffData;
	Efficiency * EffData_glob;
	
	TH1F * ExposureZones;

	TH1F * LatCorrections[10];
	TH1F * LatEfficiencies[10];
	
	TH1F * GlobalEfficiency;

	TH1F * GlobalCorrection;
	TH1F * GlobalCorrection2;
	TH1F * GlobalCorrectionpid;
	
	TH1F * Syst_Err;
        TH1F * Stat_Err;

	bool IsTrigEffCorr = false;

	std::string basename;
	std::string directory;

	BadEventSimulator * BadEvSim;	
	public:

	EffCorr(FileSaver  File, std::string Basename,std::string Directory, Binning Bins, std::string Cut_before,std::string Cut_after,std::string Cut_Data,std::string Cut_MC,std::string Cut_MC2,std::string Cut_MCpid){
	
		EffMC   = new Efficiency(File, (Basename+"_MC" ).c_str(),Directory,Bins, (Cut_before+"&"+Cut_MC  ).c_str(),(Cut_after+"&"+Cut_MC  ).c_str());
		EffMC2 = new Efficiency(File, (Basename+"_MC2").c_str(),Directory,Bins, (Cut_before+"&"+Cut_MC2  ).c_str(),(Cut_after+"&"+Cut_MC2  ).c_str());
		EffMCpid = new Efficiency(File, (Basename+"_MCpid").c_str(),Directory,Bins, (Cut_before+"&"+Cut_MCpid  ).c_str(),(Cut_after+"&"+Cut_MCpid  ).c_str());

		EffData = new Efficiency(File, (Basename+"_lat").c_str(),Directory,Bins, (Cut_before+"&"+Cut_Data).c_str(),(Cut_after+"&"+Cut_Data).c_str(),LatEdges);
		EffData_glob = new Efficiency(File, (Basename+"_glob").c_str(),Directory,Bins,(Cut_before+"&"+Cut_Data).c_str(),(Cut_after+"&"+Cut_Data).c_str());
		
		TFile * file = File.GetFile();
		if(file){
			for(int lat=0;lat<10;lat++) LatCorrections[lat]=(TH1F*) file->Get((Directory+"/"+Basename+"/" + Basename + "_Corr_lat"+to_string(lat)).c_str());
			for(int lat=0;lat<10;lat++) LatEfficiencies[lat]=(TH1F*) file->Get((Directory+"/"+Basename+"/" + Basename + "_Eff_lat"+to_string(lat)).c_str());
	
			GlobalEfficiency = (TH1F*) file->Get((Directory+"/"+Basename+"/" + Basename + "_glob_Eff").c_str());
			GlobalCorrection = (TH1F*) file->Get((Directory+"/"+Basename+"/" + Basename + "_Corr_glob").c_str());
			GlobalCorrection2 = (TH1F*) file->Get((Directory+"/"+Basename+"/" + Basename + "_Corr_glob2").c_str());
			GlobalCorrectionpid = (TH1F*) file->Get((Directory+"/"+Basename+"/" + Basename + "_Corr_globpid").c_str());
			Syst_Err    = (TH1F*) file->Get((Directory+"/"+Basename+"/" + Basename + "Syst_Err").c_str());
			Stat_Err = (TH1F*) file->Get((Directory+"/"+Basename+"/" + Basename + "Stat_Err").c_str());




		}
		basename = Basename;	
		directory = Directory;
	}	

	virtual void SetUpBadEventSimulator(BadEventSimulator * Sim) {
		BadEvSim = Sim;
	};
	
	virtual void LoadEventIntoBadEvSim(Variables * vars) {
		EffMC->LoadEventIntoBadEvSim(vars);
		EffMC2->LoadEventIntoBadEvSim(vars);
		EffMCpid->LoadEventIntoBadEvSim(vars);
	}
	virtual bool ReinitializeHistos(bool refill){
		bool checkifsomeismissing=false;
		bool allfound=true;
		if(!(EffMC -> ReinitializeHistos(refill))) checkifsomeismissing   = true;
         	if(!(EffMC2 -> ReinitializeHistos(refill))) checkifsomeismissing   = true;
		if(!(EffMCpid -> ReinitializeHistos(refill))) checkifsomeismissing   = true;
	        if(!(EffData -> ReinitializeHistos(refill))) checkifsomeismissing = true;
		if(!(EffData_glob -> ReinitializeHistos(refill))) checkifsomeismissing = true;
	 	 
		if(checkifsomeismissing||refill) allfound=false;
		return allfound;
	}
	virtual void FillEventByEventMC(Variables * vars, float (*var) (Variables * vars), float (*discr_var) (Variables * vars)){
		EffMC -> FillEventByEventMC(vars,var,discr_var);
		EffMC2 -> FillEventByEventMC(vars,var,discr_var);
		EffMCpid -> FillEventByEventMC(vars,var,discr_var);
	}
	
	virtual void FillEventByEventData(Variables * vars, float (*var) (Variables * vars), float (*discr_var) (Variables * vars)){
		EffData -> FillEventByEventData(vars,var,discr_var);
		EffData_glob -> FillEventByEventData(vars,var,discr_var);
	
	}
	virtual void Fill(TTree * treeMC,TTree * treeDT, Variables * vars, float (*discr_var) (Variables * vars),bool refill=false);
	virtual void Save();
	virtual void Eval_Efficiencies();
	virtual void SaveResults(FileSaver finalhistos);
	virtual void Eval_Corrections();
	void SetToConstantValue(float value);
	void SetDefaultOutFile(FileSaver FinalHistos);
	void SetAsTrigEffCorr() {IsTrigEffCorr = true;};
	
	void Eval_Errors();

	TH1F * GetMCEfficiency()	{return (TH1F*)EffMC  -> GetEfficiency();}
	TH1F * GetMCEfficiency2()	{return (TH1F*)EffMC2  -> GetEfficiency();}
	TH1F * GetMCEfficiency_noPID()	{return (TH1F*)EffMCpid  -> GetEfficiency();}
	
	TH1F * GetCorrectionLat(int lat)  {return LatCorrections[lat];}
	TH1F * GetEfficiencyLat(int lat)  {return LatEfficiencies[lat];}
	
	TH1F * GetGlobEfficiency()	  {return GlobalEfficiency;}
	TH1F * GetGlobCorrection()	  {return GlobalCorrection;}
	TH1F * GetGlobCorrection2()	  {return GlobalCorrection2;}
	TH1F * GetGlobCorrection_noPID()	  {return GlobalCorrectionpid;}

	TH1F * GetStat_Err() 		  {return Stat_Err;}
	TH1F * GetSyst_Err() 		  {return Syst_Err;}


};

#endif
