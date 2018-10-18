#ifndef EFFCORR_H
#define EFFCORR_H

#include "Tool.h"
#include "Efficiency.h"


class EffCorr : public Tool{

	private:

	Efficiency * EffMC;
	Efficiency * EffMC2;
	Efficiency * EffMCnopid;

	Efficiency * EffData;
	
	TH1F * ExposureZones;

	TH1F * LatCorrections[10];
	TH1F * LatEfficiencies[10];
	
	TH1F * GlobalEfficiency;

	TH1F * GlobalCorrection;
	TH1F * GlobalCorrection2;
	TH1F * GlobalCorrectionnopid;



	std::string basename;
	std::string directory;

	BadEventSimulator * BadEvSim;	
	public:

	EffCorr(FileSaver  File, std::string Basename,std::string Directory, Binning Bins, std::string Cut_before,std::string Cut_after,std::string Cut_Data,std::string Cut_MC,std::string Cut_MC2,std::string Cut_MCnopid){
	
		EffMC   = new Efficiency(File, (Basename+"_MC" ).c_str(),Directory,Bins, (Cut_before+"&"+Cut_MC  ).c_str(),(Cut_after+"&"+Cut_MC  ).c_str());
		EffMC2 = new Efficiency(File, (Basename+"_MC2").c_str(),Directory,Bins, (Cut_before+"&"+Cut_MC2  ).c_str(),(Cut_after+"&"+Cut_MC2  ).c_str());
		EffMCnopid = new Efficiency(File, (Basename+"_MCnopid").c_str(),Directory,Bins, (Cut_before+"&"+Cut_MCnopid  ).c_str(),(Cut_after+"&"+Cut_MCnopid  ).c_str());

		EffData = new Efficiency(File, (Basename+"_lat").c_str(),Directory,Bins, (Cut_before+"&"+Cut_Data).c_str(),(Cut_after+"&"+Cut_Data).c_str(),LatEdges);
		TFile * file = File.GetFile();
		if(file){
			for(int lat=0;lat<10;lat++) LatCorrections[lat]=(TH1F*) file->Get((Directory+"/"+Basename+"/" + Basename + "_Corr_lat"+to_string(lat)).c_str());
			for(int lat=0;lat<10;lat++) LatEfficiencies[lat]=(TH1F*) file->Get((Directory+"/"+Basename+"/" + Basename + "_Eff_lat"+to_string(lat)).c_str());
	
			GlobalEfficiency = (TH1F*) file->Get((Directory+"/"+Basename+"/" + Basename + "_Eff_glob").c_str());
			GlobalCorrection = (TH1F*) file->Get((Directory+"/"+Basename+"/" + Basename + "_Corr_glob").c_str());
			GlobalCorrection2 = (TH1F*) file->Get((Directory+"/"+Basename+"/" + Basename + "_Corr_glob2").c_str());
			GlobalCorrectionnopid = (TH1F*) file->Get((Directory+"/"+Basename+"/" + Basename + "_Corr_globnopid").c_str());
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
		EffMCnopid->LoadEventIntoBadEvSim(vars);
	}
	virtual bool ReinitializeHistos(bool refill){
		bool checkifsomeismissing=false;
		bool allfound=true;
		if(!(EffMC -> ReinitializeHistos(refill))) checkifsomeismissing   = true;
         	if(!(EffMC2 -> ReinitializeHistos(refill))) checkifsomeismissing   = true;
		if(!(EffMCnopid -> ReinitializeHistos(refill))) checkifsomeismissing   = true;
	        if(!(EffData -> ReinitializeHistos(refill))) checkifsomeismissing = true;
	 	if(checkifsomeismissing||refill) allfound=false;
		return allfound;
	}
	virtual void FillEventByEventMC(Variables * vars, float (*var) (Variables * vars), float (*discr_var) (Variables * vars)){
		EffMC -> FillEventByEventMC(vars,var,discr_var);
		EffMC2 -> FillEventByEventMC(vars,var,discr_var);
		EffMCnopid -> FillEventByEventMC(vars,var,discr_var);
	}
	
	virtual void FillEventByEventData(Variables * vars, float (*var) (Variables * vars), float (*discr_var) (Variables * vars)){
		EffData -> FillEventByEventData(vars,var,discr_var);
	}
	virtual void Fill(TTree * treeMC,TTree * treeDT, Variables * vars, float (*discr_var) (Variables * vars),bool refill=false);
	virtual void Save();
	virtual void Eval_Efficiencies();
	virtual void SaveResults(FileSaver finalhistos);
	virtual void Eval_Corrections();
	void SetToConstantValue(float value);
	void SetDefaultOutFile(FileSaver FinalHistos);

	TH1F * GetMCEfficiency()	{return (TH1F*)EffMC  -> GetEfficiency();}
	TH1F * GetMCEfficiency2()	{return (TH1F*)EffMC2  -> GetEfficiency();}
	TH1F * GetMCEfficiency_noPID()	{return (TH1F*)EffMCnopid  -> GetEfficiency();}
	
	TH1F * GetCorrectionLat(int lat)  {return LatCorrections[lat];}
	TH1F * GetEfficiencyLat(int lat)  {return LatEfficiencies[lat];}
	
	TH1F * GetGlobEfficiency()	  {return GlobalEfficiency;}
	TH1F * GetGlobCorrection()	  {return GlobalCorrection;}
	TH1F * GetGlobCorrection2()	  {return GlobalCorrection2;}
	TH1F * GetGlobCorrection_noPID()	  {return GlobalCorrectionnopid;}



};

#endif
