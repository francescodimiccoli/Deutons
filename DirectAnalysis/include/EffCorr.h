#ifndef EFFCORR_H
#define EFFCORR_H

#include "Tool.h"
#include "Efficiency.h"
#include "Globals.h"
#include "PlottingFunctions.h"


class EffCorr : public Tool{

	private:

	Efficiency * EffMC;
	Efficiency * EffMC2;
	Efficiency * EffMCpid2;
	Efficiency * EffMCpid;

	Efficiency * EffData;
	Efficiency * EffData_glob;
	
	TH1F * ExposureZones;

	TH1F * LatCorrections[10];
	TH1F * LatEfficiencies[10];
	
	TH1F * GlobalEfficiency;
	TH1F * GlobalEfficiency_plot;

	TH1F * GlobalCorrection;
	TH1F * GlobalCorrectionpid;

	TGraphErrors * GlobalCorrection_timeavg;
	
	
	TH1F * Syst_Err;
        TH1F * Stat_Err;
	TH1F * Syst_Stat;
        
	bool IsTrigEffCorr = false;
	bool IsEkinCorrection = true;

	std::string basename;
	std::string directory;

	BadEventSimulator * BadEvSim;	
	
	TF1 * CorrectionModel;
	TF1 * DataEffModel;
	TSpline3 * CorrectionModel_Spline;
	TSpline3 * DataEffModel_Spline;

	TH1F * syst_stat=0x0;
	TH2F * avg_time=0x0;
	int time=1429054000;


	bool splinemodel=false;

	protected:
	Binning Bins;
	Binning Bins_D;

	public:

	EffCorr(FileSaver  File, std::string Basename,std::string Directory, bool ekin, int nucl, std::string Cut_before,std::string Cut_after,std::string Cut_Data,bool notweighted = false,std::string cut_MC="IsProtonMC"){
	

		
		IsEkinCorrection = ekin;
		SetBinsForCorrections(IsEkinCorrection);
		if(nucl<2) Bins = ForEffCorr;
		else Bins = ForEffCorr_D;
		std::string Cut_MC = cut_MC;
		std::string Cut_MCpid = "IsPurePMC";

		EffMC     = new Efficiency(File, (Basename+"_MC" ).c_str(),Directory,Bins, (Cut_before+"&"+Cut_MC  ).c_str(),(Cut_after+"&"+Cut_MC  ).c_str(),false);
		EffMCpid  = new Efficiency(File, (Basename+"_MCpid").c_str(),Directory,Bins, (Cut_before+"&"+Cut_MCpid  ).c_str(),(Cut_after+"&"+Cut_MCpid  ).c_str(),false);
	
		EffData = new Efficiency(File, (Basename+"_lat").c_str(),Directory,Bins, (Cut_before+"&"+Cut_Data).c_str(),(Cut_after+"&"+Cut_Data).c_str(),LatEdges);
		EffData_glob = new Efficiency(File, (Basename+"_glob").c_str(),Directory,Bins,(Cut_before+"&"+Cut_Data).c_str(),(Cut_after+"&"+Cut_Data).c_str());
	
		TFile * file = File.GetFile();	
		if(file){
			for(int lat=0;lat<10;lat++) LatCorrections[lat]=(TH1F*) file->Get((Directory+"/"+Basename+"/" + Basename + "_Corr_lat"+to_string(lat)).c_str());
			for(int lat=0;lat<10;lat++) LatEfficiencies[lat]=(TH1F*) file->Get((Directory+"/"+Basename+"/" + Basename + "_Eff_lat"+to_string(lat)).c_str());
	
			GlobalEfficiency = (TH1F*) file->Get((Directory+"/"+Basename+"/" + Basename + "_glob_Eff").c_str());
			GlobalCorrection = (TH1F*) file->Get((Directory+"/"+Basename+"/" + Basename + "_Corr_glob").c_str());
			GlobalCorrection_timeavg = (TGraphErrors*) file->Get((Directory+"/"+Basename+"/" + Basename + "_Corr_avgtime").c_str());
			GlobalCorrectionpid = (TH1F*) file->Get((Directory+"/"+Basename+"/" + Basename + "_Corr_globpid").c_str());
			Syst_Err    = (TH1F*) file->Get((Directory+"/"+Basename+"/" + Basename + "Syst_Err").c_str());
			Stat_Err = (TH1F*) file->Get((Directory+"/"+Basename+"/" + Basename + "Stat_Err").c_str());
			Syst_Stat = (TH1F*) file->Get((Directory+"/"+Basename+"/" + Basename + "Syst_Stat").c_str());
			CorrectionModel_Spline = (TSpline3 *) file->Get((Directory+"/"+Basename+"/" + Basename + "_CorrSpline").c_str());
			DataEffModel_Spline = (TSpline3 *) file->Get((Directory+"/"+Basename+"/" + Basename + "_EffSpline").c_str());


		}
		basename = Basename;	
		directory = Directory;
	
		if(notweighted){
			EffMC     ->SetNotWeightedMC_phtr();	
			EffMCpid  ->SetNotWeightedMC_phtr();	
		}
	}

	EffCorr(FileSaver File,std::string Basename, std::string Directory, Binning bins,TH1F * DataBefore, TH1F * DataAfter, TH1F * MCBefore, TH1F * MCAfter){

		Bins = bins;

		basename = Basename;	
		directory = Directory;
	
		EffMC     = new Efficiency(File, (Basename+"_MC" ).c_str(),directory,Bins, "","",false);
		EffMCpid  = new Efficiency(File, (Basename+"_MCpid").c_str(),directory,Bins, "","",false);
	
		EffData = new Efficiency(File, (Basename+"_lat").c_str(),directory,Bins, "","",LatEdges);
		EffData_glob = new Efficiency(File, (Basename+"_glob").c_str(),directory,Bins,"","");
	
		EffMC->SetBefore((TH1F*)MCBefore->Clone((Basename+"_MC_before" ).c_str()));	
		EffMC->SetAfter( (TH1F*)MCAfter->Clone((Basename+"_MC_after" ).c_str()));	
		EffMCpid->SetBefore( (TH1F*)MCBefore->Clone((Basename+"_MCpid_before" ).c_str()));	
		EffMCpid->SetAfter( (TH1F*)MCAfter->Clone((Basename+"_MCpid_after" ).c_str()));	

		EffData_glob->SetBefore( (TH1F*)DataBefore->Clone((Basename+"_glob_before" ).c_str()));	
		EffData_glob->SetAfter( (TH1F*)DataAfter->Clone((Basename+"_glob_after" ).c_str()));	

	}

	Binning GetBins(){ return Bins;}	
	bool IsEkin(){return IsEkinCorrection;}
	TF1 * GetCorrectionModel(){return CorrectionModel;}
	TF1 * GetDataEffModel(){return DataEffModel;}
	std::string  GetName(){return basename;}

	void SetBinsForCorrections(bool ekin){
		if(ekin) {
				ForEffCorr.Reset();
				ForEffCorr_D.Reset();
		/*		ForEffCorr.setBinsFromRDatacard ((workdir+"/bindatacard_PMIT.data").c_str(), 0.1, 0.9999999 ,ResponseTOF,0.00347548,5.8474); ;
				ForEffCorr_D.setBinsFromRDatacard ((workdir+"/bindatacard_PMIT.data").c_str(), 0.1, 0.9999999 ,ResponseTOF,0.00347548,5.8474); 
		*/
				float ekmin=0.1, ekmax=25;
				ForEffCorr.setBinsFromEkPerMass (50, ekmin, ekmax,ResponseTOF,0.00347548,5.8474);
				ForEffCorr_D.setBinsFromEkPerMass(50, ekmin, ekmax,ResponseTOF,0.00347548,5.8474);
				ForEffCorr.UseREdges();
				ForEffCorr_D.UseREdges();
				cout<<"FOREFFCORRBINS: "<<basename<<endl;
				ForEffCorr.Print();
				ForEffCorr_D.Print();
		}

		else {
				ForEffCorr.Reset();
				ForEffCorr_D.Reset();
				ForEffCorr.setBinsFromRDatacard ((workdir+"/bindatacard_PMIT.data").c_str(), 0.1, 0.9999999 ,ResponseTOF,0.00347548,5.8474); ;
				ForEffCorr_D.setBinsFromRDatacard ((workdir+"/bindatacard_PMIT.data").c_str(), 0.1, 0.9999999 ,ResponseTOF,0.00347548,5.8474); 
				ForEffCorr.UseREdges();
				ForEffCorr_D.UseREdges();
				cout<<"FOREFFCORRBINS: "<<basename<<endl;
				ForEffCorr.Print();
				ForEffCorr_D.Print();
	
		}

	}

	virtual void SetUpBadEventSimulator(BadEventSimulator * Sim) {
		BadEvSim = Sim;
	};
	
	virtual void LoadEventIntoBadEvSim(Variables * vars) {
		EffMC->LoadEventIntoBadEvSim(vars);
		EffMCpid->LoadEventIntoBadEvSim(vars);
	}
	virtual bool ReinitializeHistos(bool refill){
		bool checkifsomeismissing=false;
		bool allfound=true;
		if(!(EffMC -> ReinitializeHistos(refill))) checkifsomeismissing   = true;
		if(!(EffMCpid -> ReinitializeHistos(refill))) checkifsomeismissing   = true;
	        if(!(EffData -> ReinitializeHistos(refill))) checkifsomeismissing = true;
		if(!(EffData_glob -> ReinitializeHistos(refill))) checkifsomeismissing = true;
	 	 
		if(checkifsomeismissing||refill) allfound=false;
		return allfound;
	}
	virtual void FillEventByEventMC(Variables * vars, float (*var) (Variables * vars), float (*discr_var) (Variables * vars)){
	//	EffMC ->  FillEventByEventMC(vars,GetGenRigidity,GetGenRigidity);
	//	EffMCpid -> FillEventByEventMC(vars ,GetGenRigidity,GetGenRigidity);
		EffMC ->  FillEventByEventMC(vars,var,discr_var);
		EffMCpid -> FillEventByEventMC(vars ,var,discr_var);
	}
	
	virtual void FillEventByEventData(Variables * vars, float (*var) (Variables * vars), float (*discr_var) (Variables * vars)){
	//	EffData -> FillEventByEventData(vars,var,discr_var);
		EffData_glob -> FillEventByEventData(vars,var,discr_var);
	
	}
	virtual void Fill(TTree * treeMC,TTree * treeDT, Variables * vars, float (*discr_var) (Variables * vars),bool refill=false);
	virtual void Save();
	virtual void Eval_Efficiencies();
	virtual void SaveResults(FileSaver finalhistos);
	virtual void Eval_Corrections(float shift=0);
	void SetToConstantValue(float value);
	void SetDefaultOutFile(FileSaver FinalHistos);
	void SetAsTrigEffCorr() {IsTrigEffCorr = true;};
	void ModelWithSpline(float shift);
	void ModelWithSimple(float shift);
	void SetSplineModel(){splinemodel=true;}

	void Eval_Errors();
	void Set_SystStat(TH1F * syst, TH2F * avg,std::string timename) { 
		if(syst) syst_stat = (TH1F*) syst->Clone();  
		if(avg)  avg_time=(TH2F*)avg->Clone();
		time=std::atoi(timename.substr(timename.find("-")+1,10).c_str())-4665600;
		std::cout<<"Set_SystPar: "<<time<<" "<<syst_stat<<" "<<avg_time<<std::endl;
	}

	TH1F * GetMCEfficiency()	{return (TH1F*)EffMC     -> GetEfficiency();}
	TH1F * GetMCEfficiency_noPID()	{return (TH1F*)EffMCpid  -> GetEfficiency();}


	TH1F * GetMCBefore()		{return (TH1F*)EffMC  -> GetBefore();}
	TH1F * GetMCBefore_noPID()	{return (TH1F*)EffMCpid  -> GetBefore();}
	
	TH1F * GetMCAfter()		{return (TH1F*)EffMC  -> GetAfter();}
	TH1F * GetMCAfter_noPID()	{return (TH1F*)EffMCpid  -> GetAfter();}
	
	TH1F * GetCorrectionLat(int lat)  {return LatCorrections[lat];}
	TH1F * GetEfficiencyLat(int lat)  {return LatEfficiencies[lat];}
	
	TH1F * GetGlobEfficiency()	  {return GlobalEfficiency;}
	TH1F * GetGlobEfficiency_plot()	  {return GlobalEfficiency_plot;}
	TH1F * GetGlobCorrection()	  {return GlobalCorrection;}
	TH1F * GetGlobCorrection_noPID()  {return GlobalCorrectionpid;}

	TH1F * GetStat_Err() 		  {return Stat_Err;}
	TH1F * GetSyst_Err() 		  {return Syst_Err;}
	TH1F * GetSyst_Stat() 		  {return Syst_Stat;}

	TH2F * Get_AvgTime() 		{return avg_time;}
	TGraphErrors * Get_TimeAvgCorrection() {return GlobalCorrection_timeavg;}
};

class Model{

	private:
		int nodes;
		std::vector<float> xs;

	public:
		Model(std::vector<float> X){
			nodes = X.size();
			for(int i=0; i<nodes;i++) xs.push_back(X[i]);
		}

	double Function(double *x, double *p); 
};


#endif
