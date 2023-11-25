#ifndef FLUX_H
#define FLUX_H

#include "Livetime.h"
#include "binning.h"
#include "Cuts.h"
#include <fstream>
#include "rundb.h"
#include "Efficiency.h"
#include "EffCorr.h"
#include "Acceptance.h"
#include "PlottingFunctions.h"
#include "Unfolding.h"
#include "TCanvas.h"

class Flux{

	private:
	Acceptance * EffAcceptance;
	TH1F * Counts=0x0;
	TH1F * Counts_density=0x0;
	TH1F * Counts_density27=0x0;
	TH1F * Counts_density_unfolded=0x0;
	TH1F * Counts_density_roounf=0x0;
	TH1F * Prior=0x0;
	TH1F * Counts_density_hRecoroounf=0x0;
	TH1F * Counts_density_Ekin=0x0;
	TH1F * Counts_density_Ekin_unf=0x0;
	TH1F * Unfolding_factor=0x0;
	TH1F * Unfolding_factor_raw=0x0;

	TH1F * Counts_ForAverage=0x0;
	TH1F * Exposure_ForAverage=0x0;

	float systfragm = 0.02;

	bool rooUnfolding=false;


	TGraphErrors * Unfolding_factor_timeavg=0x0;
	TGraphErrors * Unfolding_factor_timeavg_err=0x0;
	TH2F * avg_time=0x0;
	int time=1429054000;

	TH1F * RooUnfolding_factor=new TH1F();;
	TGraph * RooUnfolding_factor_g=new TGraphErrors;
	TH1F * RooUnfolding_factor_raw=new TH1F();
	
	TH1F * Eff_Acceptance=0x0;
	TH1F * MC_Acceptance=0x0;
	TH1F * MC_Acceptance_raw=0x0;
	TH1F * ExposureTime=0x0;
	TH1F * ForAcceptance=0x0;

	TH1F * Eff_Acceptance_plot_rig=0x0;
	TH1F * ExposureTime_plot_rig=0x0;
	TH1F * Eff_Acceptance_plot=0x0;
	TH1F * ExposureTime_plot=0x0;

	TH2F * migr_matr=0x0;
		
	Binning bins;
	std::string basename;	
	std::string exposurename;
	TF1 * AcceptanceModel;	

	TH1F * FluxEstim=0x0;
	TH1F * FluxEstim_ekin_unf=0x0;
	TH1F * FluxEstim_rig=0x0;
	TH1F * FluxEstim_rig_stat=0x0;
	TH1F * FluxEstim_unf=0x0;
	TH1F * FluxEstim_unf_stat=0x0;
	TH1F * DeltaE=0x0;

	std::vector<EffCorr*> EfficiencyCorrections;
	std::vector<EffCorr*> EfficiencyFromData;

	FileSaver finalhistos;

	TCanvas * test=0x0;

	TH1F * Counts_statErr=0x0;
	TH1F * TOT_statErr=0x0;
	TH1F * Counts_systErr=0x0;
	TH1F * Counts_systErr_unb=0x0;
	TH1F * Acc_Err=0x0;
	TH1F * Acc_ErrStat=0x0;
	TH1F * Unfolding_Err=0x0;
	TH1F * RooUnfolding_Err=0x0;


	bool forcefolded = false;


	TH2F * migr_matr_rootunfold=0x0;
	TH1F* roounfold_true=0x0;
	TH1F* roounfold_meas=0x0;


	bool Usepdf=true;

	public:

	Flux(FileSaver File, FileSaver FileRes, std::string Basename,std::string Accname, std::string AccDir,std::string CountsName,std::string CountsStatName,std::string CountsSystName,std::string ExposureName, Binning Bins){
	
		EffAcceptance = new Acceptance(FileRes,Accname,AccDir,Bins,Bins);
		Eff_Acceptance = (TH1F *) EffAcceptance->GetEffAcc_gen();	
		MC_Acceptance = (TH1F *) EffAcceptance->GetEffAccMC_gen();	
		MC_Acceptance_raw = (TH1F *) EffAcceptance->GetEffAccMC_gen_raw();	



		if(EffAcceptance) {
		cout<<"************ EFFECTIVE ACCEPTANCE FOUND: ********************"<<endl;
		cout<<Eff_Acceptance<<" "<<EffAcceptance->GetEffAcc()<<" "<<EffAcceptance->GetEffAcc_gen()<<endl; }

		if(FileRes.CheckFile()) Counts = (TH1F *) FileRes.Get((CountsName).c_str());	 
		if(FileRes.CheckFile()) Counts_statErr = (TH1F *) FileRes.Get((CountsStatName).c_str());	 
		if(FileRes.CheckFile()) Counts_systErr = (TH1F *) FileRes.Get((CountsSystName).c_str());	 
		if(FileRes.CheckFile()) Counts_systErr_unb = (TH1F *) FileRes.Get((CountsSystName+"_unb").c_str());	 


		//Counts_statErr->Scale(1.5);

		ExposureTime = (TH1F *) File.Get(("Fluxes/"+Basename+"/"+ExposureName).c_str());
		cout<<("Fluxes/"+Basename+"/"+ExposureName).c_str()<<endl;
		cout<<ExposureTime<<"; File Name: "<<File.GetName()<<endl;
		bins = Bins;		
		basename = Basename;
		exposurename = ExposureName;
		EfficiencyCorrections.clear();
		EfficiencyFromData.clear();

		if(migr_matr) cout<<"**************** MIGRATION MATRIX FOUND: ********************"<<endl;
		
		cout<<"Matrix: "<<migr_matr<<endl;


		if(bins.IsUsingBetaEdges()) {
			migr_matr_rootunfold = (TH2F*) EffAcceptance->GetMigr_Unf();
			roounfold_true =  (TH1F*) EffAcceptance->GetMigr_True();
			roounfold_meas =  (TH1F*) EffAcceptance->GetMigr_Meas();
		}


	}


	bool ReinitializeHistos(bool refill) {
		if(!finalhistos.CheckFile()||refill) { 

			cout<<"BINS SIZE: "<<bins.size()<<endl;
			ExposureTime = new TH1F(exposurename.c_str(),exposurename.c_str(),bins.size(),0,bins.size());
		
			Counts=0x0;
			return false;
		}
		else {
		return true;	}
	}
	void Eval_ExposureTime(Variables * vars, TTree * treeDT,FileSaver finalhistos,bool refill);

        void AverageCountsWithOther(FileSaver FileRes, std::string namecounts, std::string nameexposure,float toll=0.05);
	
	void Eval_Flux(float corr_acc=1, float fit_min=0, float fit_max=0,int knots=10,float offset=0.0,bool regularize=false);
	void Eval_Errors();
	void SaveResults(FileSaver finalhistos,std::string bn="");
	void ChangeName (std::string newname) {basename = newname; return;}
	Binning GetBins(){return bins;}
	std::string GetName(){return basename;}

	void SetDefaultOutFile(FileSaver FinalHistos) {finalhistos = FinalHistos; return;}
	FileSaver GetOutFileSaver() {return finalhistos;}
	void Unfold_Counts(float fit_min, float fit_max,int knots,float offset, bool regularize=false);
	void ModelAcceptanceWithSpline(float shift);

	void ActivateRooUnfolding(bool usepdf=true){ rooUnfolding=true; Usepdf=usepdf;}

	void Set_UnfoldingTime(TGraphErrors * avg,TGraphErrors * avg_err) {
		Unfolding_factor_timeavg = (TGraphErrors *) avg->Clone(); 
		Unfolding_factor_timeavg_err = (TGraphErrors *) avg_err->Clone(); 
	}

	void Roounfold(int iterations,bool usepdf);
	
	
	void SetForceFolded(){forcefolded=true;}
	TH1F * GetExposureTime(){return ExposureTime;}
	
	TH1F * GetFlux(){return FluxEstim;}
	TH1F * GetFlux_ekin_unf(){return FluxEstim_ekin_unf;}
	TH1F * GetFlux_rig(){return FluxEstim_rig;}
	TH1F * GetFlux_rig_stat(){return FluxEstim_rig_stat;}
	TH1F * GetFlux_unf(){return FluxEstim_unf;}
	TH1F * GetFlux_unf_stat(){return FluxEstim_unf_stat;}


	TH1F * GetEffAcceptance(){return Eff_Acceptance;}
	TH1F * GetMCAcceptance(){return MC_Acceptance;}
	TH1F * GetMCAcceptance_raw(){return MC_Acceptance_raw;}
	TH1F * GetCounts() {return Counts;}
	TH1F * GetStatError() {return TOT_statErr;}
	TH1F * GetAccError() {return Acc_Err;}
	TH1F * GetSystError() {return Counts_systErr;}
	TH1F * GetSystError_unb() {return Counts_systErr_unb;}
	TH1F * GetUnfError() {return Unfolding_Err;}
	TH1F * GetUnfoldingFactor() {return Unfolding_factor;}
	TH1F * GetRooUnfoldingFactor() {return RooUnfolding_factor;}
	TH1F * GetRooUnfoldingFactor_Raw() {return RooUnfolding_factor_raw;}
	TH1F * GetRooUnfError() {return RooUnfolding_Err;}
	
	void Average_with_another(TH1F* ratio);
	TH1F * Eval_FluxRatio(Flux * Denominator,std::string name);
	void SetSystFragm(float sfr) {systfragm=sfr;}
};


#endif
