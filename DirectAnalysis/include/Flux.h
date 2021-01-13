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
	TH1F * Counts_density_unfolded=0x0;
	TH1F * Counts_density_Ekin=0x0;
	TH1F * Unfolding_factor=0x0;
	
	TH1F * Eff_Acceptance=0x0;
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

	TH1F * FluxEstim=0x0;
	TH1F * FluxEstim_rig=0x0;
	TH1F * FluxEstim_unf=0x0;


	std::vector<EffCorr*> EfficiencyCorrections;
	std::vector<EffCorr*> EfficiencyFromData;

	FileSaver finalhistos;

	TCanvas * test;


	TH1F * Counts_statErr=0x0;
	TH1F * Counts_systErr=0x0;
	TH1F * Acc_Err=0x0;
	TH1F * Unfolding_Err=0x0;


	public:

	Flux(FileSaver File, FileSaver FileRes, std::string Basename,std::string Accname, std::string AccDir,std::string CountsName,std::string CountsStatName,std::string ExposureName, Binning Bins){
	
		EffAcceptance = new Acceptance(FileRes,Accname,AccDir,Bins,Bins);
		Eff_Acceptance = (TH1F *) EffAcceptance->GetEffAcc();	

		if(EffAcceptance) {
		cout<<"************ EFFECTIVE ACCEPTANCE FOUND: ********************"<<endl;
		cout<<Eff_Acceptance<<endl; }

		if(FileRes.CheckFile()) Counts = (TH1F *) FileRes.Get((CountsName).c_str());	 
		if(FileRes.CheckFile()) Counts_statErr = (TH1F *) FileRes.Get((CountsStatName).c_str());	 

		ExposureTime = (TH1F *) File.Get(("Fluxes/"+Basename+"/"+ExposureName).c_str());
		cout<<("Fluxes/"+Basename+"/"+ExposureName).c_str()<<" "<<ExposureTime<<"; File Name: "<<File.GetName()<<endl;
		bins = Bins;		
		basename = Basename;
		exposurename = ExposureName;
		EfficiencyCorrections.clear();
		EfficiencyFromData.clear();

		if(bins.IsUsingBetaEdges()) migr_matr = (TH2F*) EffAcceptance->GetMigr_beta();
		else migr_matr = (TH2F*) EffAcceptance->GetMigr_rig();
		if(migr_matr) cout<<"**************** MIGRATION MATRIX FOUND: ********************"<<endl;
		
		cout<<"Matrix: "<<migr_matr<<endl;

		test = new TCanvas();

	}


	Flux(FileSaver FileRes, std::string Basename, std::string Accname, std::string AccDir,std::string CountsName,std::string ExposureName, Binning Bins){
		
		EffAcceptance = new Acceptance(FileRes,Accname,AccDir,Bins,Bins);
		Eff_Acceptance = (TH1F *) EffAcceptance->GetEffAcc();	

		if(EffAcceptance) {
			cout<<"************ EFFECTIVE ACCEPTANCE FOUND: ********************"<<endl;
			cout<<Eff_Acceptance<<endl; }

		TFile * fileres = FileRes.GetFile();

			if(FileRes.CheckFile()) {
			Counts = (TH1F *) fileres->Get((CountsName).c_str());
			ExposureTime = (TH1F *) fileres->Get(("Fluxes/"+Basename+"/"+ExposureName).c_str());

			FluxEstim = (TH1F *) fileres->Get(("Fluxes/"+Basename+"/"+Basename+"_Flux").c_str());
			FluxEstim_rig = (TH1F *) fileres->Get(("Fluxes/"+Basename+"/"+Basename+"_Flux_rig").c_str());
			FluxEstim_unf = (TH1F *) fileres->Get(("Fluxes/"+Basename+"/"+Basename+"_Flux_unf").c_str());
			cout<<("Fluxes/"+Basename+"/"+Basename+"_Flux").c_str()<<" "<<FluxEstim<<"; File Name: "<<FileRes.GetName()<<endl;
			cout<<("Fluxes/"+Basename+"/"+Basename+"_Flux_rig").c_str()<<" "<<FluxEstim_rig<<"; File Name: "<<FileRes.GetName()<<endl;
			EfficiencyCorrections.clear();
			EfficiencyFromData.clear();
		
		}
		bins = Bins;		
		basename = Basename;
	}
	bool ReinitializeHistos(bool refill) {
		if(!finalhistos.CheckFile()||refill) { 
			ExposureTime = new TH1F(exposurename.c_str(),exposurename.c_str(),bins.size(),0,bins.size());
		
			Counts=0x0;
			return false;
		}
		else {
		return true;	}
	}
	void Eval_ExposureTime(Variables * vars, TTree * treeDT,FileSaver finalhistos,bool refill);
	
	void Eval_Flux(float corr_acc=1, float fit_min=0, float fit_max=0,int knots=10,float offset=0.0,int uptobin=99);
	void Eval_Errors();
	void SaveResults(FileSaver finalhistos);
	void ChangeName (std::string newname) {basename = newname; return;}
	Binning GetBins(){return bins;}
	std::string GetName(){return basename;}

	void SetDefaultOutFile(FileSaver FinalHistos) {finalhistos = FinalHistos; return;}
	FileSaver GetOutFileSaver() {return finalhistos;}
	void Unfold_Counts(float fit_min, float fit_max,int knots,float offset, int uptobin);
	
	TH1F * GetExposureTime(){return ExposureTime;}
	
	TH1F * GetFlux(){return FluxEstim;}
	TH1F * GetFlux_rig(){return FluxEstim_rig;}
	TH1F * GetFlux_unf(){return FluxEstim_unf;}

	TH1F * GetEffAcceptance(){return Eff_Acceptance;}
	TH1F * GetCounts() {return Counts;}
	TH1F * GetStatError() {return Counts_statErr;}

	TH1F * Eval_FluxRatio(Flux * Denominator,std::string name);

};


#endif
