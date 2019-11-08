#ifndef EFFCORR_H
#define EFFCORR_H

#include "Tool.h"
#include "Efficiency.h"
#include "Globals.h"

class EffCorr : public Tool{

	private:

	Efficiency * EffMC;
	Efficiency * EffMCpid2;
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
	bool IsEkinCorrection = true;

	std::string basename;
	std::string directory;

	BadEventSimulator * BadEvSim;	
	
	TF1 * CorrectionModel;
	TSpline3 * CorrectionModel_Spline;


	protected:
	Binning Bins;
	Binning Bins_D;

	public:

	EffCorr(FileSaver  File, std::string Basename,std::string Directory, bool ekin, std::string Cut_before,std::string Cut_after,std::string Cut_Data,std::string Cut_MC,std::string Cut_MC2,std::string Cut_MCpid){
	

		
		IsEkinCorrection = ekin;
		SetBinsForCorrections(IsEkinCorrection);
		Bins = ForEffCorr;
		Bins_D = ForEffCorr_D;	

		EffMC   = new Efficiency(File, (Basename+"_MC" ).c_str(),Directory,Bins, (Cut_before+"&"+Cut_MC  ).c_str(),(Cut_after+"&"+Cut_MC  ).c_str(),true);
		EffMCpid2 = new Efficiency(File, (Basename+"_MC2").c_str(),Directory,Bins_D, (Cut_before+"&"+Cut_MC2  ).c_str(),(Cut_after+"&"+Cut_MC2  ).c_str(),true);
		EffMCpid = new Efficiency(File, (Basename+"_MCpid").c_str(),Directory,Bins, (Cut_before+"&"+Cut_MCpid  ).c_str(),(Cut_after+"&"+Cut_MCpid  ).c_str(),true);

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
			CorrectionModel_Spline = (TSpline3 *) file->Get((Directory+"/"+Basename+"/" + Basename + "_CorrSpline").c_str());

		}
		basename = Basename;	
		directory = Directory;
	}

	Binning GetBins(){ return Bins;}	
	bool IsEkin(){return IsEkinCorrection;}
	TSpline3 * GetCorrectionModel(){return CorrectionModel_Spline;}

	void SetBinsForCorrections(bool ekin){
		if(ekin) {
				ForEffCorr.Reset();
				ForEffCorr_D.Reset();
				ForEffCorr.setBinsFromEkPerMass(nbinsr,0.15,15,ResponseTOF,0.00347548,5.8474);
				ForEffCorr_D.setBinsFromEkPerMass(nbinsr,0.15,15,ResponseTOF,0.00347548,5.8474);
				ForEffCorr.UseREdges();
				ForEffCorr_D.UseREdges();
				ForEffCorr.Print();
				ForEffCorr_D.Print();
		}

		else {
				ForEffCorr.Reset();
				ForEffCorr_D.Reset();
				ForEffCorr.setBinsFromRigidity(nbinsr,0.5,30,ResponseTOF,0.00347548,5.8474);
				ForEffCorr_D.setBinsFromRigidity(nbinsr,0.5,30,ResponseTOF,0.00347548,5.8474);
				ForEffCorr.UseREdges();
				ForEffCorr_D.UseREdges();
				ForEffCorr.Print();
				ForEffCorr_D.Print();
	
		}

	}

	virtual void SetUpBadEventSimulator(BadEventSimulator * Sim) {
		BadEvSim = Sim;
	};
	
	virtual void LoadEventIntoBadEvSim(Variables * vars) {
		EffMC->LoadEventIntoBadEvSim(vars);
		EffMCpid2->LoadEventIntoBadEvSim(vars);
		EffMCpid->LoadEventIntoBadEvSim(vars);
	}
	virtual bool ReinitializeHistos(bool refill){
		bool checkifsomeismissing=false;
		bool allfound=true;
		if(!(EffMC -> ReinitializeHistos(refill))) checkifsomeismissing   = true;
         	if(!(EffMCpid2 -> ReinitializeHistos(refill))) checkifsomeismissing   = true;
		if(!(EffMCpid -> ReinitializeHistos(refill))) checkifsomeismissing   = true;
	        if(!(EffData -> ReinitializeHistos(refill))) checkifsomeismissing = true;
		if(!(EffData_glob -> ReinitializeHistos(refill))) checkifsomeismissing = true;
	 	 
		if(checkifsomeismissing||refill) allfound=false;
		return allfound;
	}
	virtual void FillEventByEventMC(Variables * vars, float (*var) (Variables * vars), float (*discr_var) (Variables * vars)){
		EffMC -> FillEventByEventMC(vars,var,discr_var);
		EffMCpid2 -> FillEventByEventMC(vars,GetGenMomentum,GetGenMomentum);
		EffMCpid -> FillEventByEventMC(vars,GetGenMomentum,GetGenMomentum);
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
	void ModelWithSpline();

	void Eval_Errors();

	TH1F * GetMCEfficiency()	{return (TH1F*)EffMC  -> GetEfficiency();}
	TH1F * GetMCEfficiency2()	{return (TH1F*)EffMCpid2  -> GetEfficiency();}
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
