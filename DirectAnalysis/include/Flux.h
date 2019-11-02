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
	
class Flux{

	private:
	Acceptance * EffAcceptance;
	TH1F * Counts=0x0;
	TH1F * Eff_Acceptance=0x0;
	TH1F * ExposureTime=0x0;
	TH1F * ForAcceptance=0x0;
	
	Binning bins;
	std::string basename;	
	std::string exposurename;

	TH1F * FluxEstim=0x0;
	TH1F * FluxEstim_rig=0x0;
	

	TH1F * Acc_StatErr=0x0;
	TH1F * Acc_SystErr=0x0;
	TH1F * Counts_Err=0x0;

	std::vector<EffCorr*> EfficiencyCorrections;
	std::vector<EffCorr*> EfficiencyFromData;

	FileSaver finalhistos;

	public:


	Flux(FileSaver File, FileSaver FileRes, std::string Basename,std::string Accname, std::string AccDir,std::string CountsName,std::string ExposureName, Binning Bins){
	
		EffAcceptance = new Acceptance(FileRes,Accname,AccDir,Bins);
		Eff_Acceptance = (TH1F *) EffAcceptance->GetEffAcc();	
		Acc_StatErr= (TH1F *) EffAcceptance->GetEffAcc_staterr();
		Acc_SystErr= (TH1F *) EffAcceptance->GetEffAcc_systerr();

		if(EffAcceptance) {
		cout<<"************ EFFECTIVE ACCEPTANCE FOUND: ********************"<<endl;
		cout<<Eff_Acceptance<<" "<<Acc_StatErr<<" "<<Acc_SystErr<<endl; }

		if(FileRes.CheckFile()) Counts = (TH1F *) FileRes.Get((CountsName).c_str());	 
		ExposureTime = (TH1F *) File.Get(("Fluxes/"+Basename+"/"+ExposureName).c_str());
		cout<<("Fluxes/"+Basename+"/"+ExposureName).c_str()<<" "<<ExposureTime<<"; File Name: "<<File.GetName()<<endl;
		bins = Bins;		
		basename = Basename;
		exposurename = ExposureName;
		EfficiencyCorrections.clear();
		EfficiencyFromData.clear();


	}


	Flux(FileSaver FileRes, std::string Basename, std::string Accname, std::string AccDir,std::string CountsName,std::string ExposureName, Binning Bins){
		
		EffAcceptance = new Acceptance(FileRes,Accname,AccDir,Bins);
		Eff_Acceptance = (TH1F *) EffAcceptance->GetEffAcc();	
		Acc_StatErr= (TH1F *) EffAcceptance->GetEffAcc_staterr();
		Acc_SystErr= (TH1F *) EffAcceptance->GetEffAcc_systerr();

		if(EffAcceptance) {
			cout<<"************ EFFECTIVE ACCEPTANCE FOUND: ********************"<<endl;
			cout<<Eff_Acceptance<<" "<<Acc_StatErr<<" "<<Acc_SystErr<<endl; }

		TFile * fileres = FileRes.GetFile();

		Counts_Err= (TH1F *) fileres->Get(("Fluxes/"+Basename+"/Counts Error").c_str());
			if(FileRes.CheckFile()) {
			Counts = (TH1F *) fileres->Get((CountsName).c_str());
			ExposureTime = (TH1F *) fileres->Get(("Fluxes/"+Basename+"/"+ExposureName).c_str());

			FluxEstim = (TH1F *) fileres->Get(("Fluxes/"+Basename+"/"+Basename+"_Flux").c_str());
			FluxEstim_rig = (TH1F *) fileres->Get(("Fluxes/"+Basename+"/"+Basename+"_Flux_rig").c_str());
			cout<<("Fluxes/"+Basename+"/"+ExposureName).c_str()<<" "<<ExposureTime<<"; File Name: "<<FileRes.GetName()<<endl;
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
	
	void Eval_Flux();
	void SaveResults(FileSaver finalhistos);
	void ChangeName (std::string newname) {basename = newname; return;}
	Binning GetBins(){return bins;}
	std::string GetName(){return basename;}

	void SetDefaultOutFile(FileSaver FinalHistos) {finalhistos = FinalHistos; return;}
	FileSaver GetOutFileSaver() {return finalhistos;}
	
	TH1F * GetExposureTime(){return ExposureTime;}
	
	TH1F * GetFlux(){return FluxEstim;}
	TH1F * GetFlux_rig(){return FluxEstim_rig;}

	TH1F * GetEffAcceptance(){return Eff_Acceptance;}
	TH1F * GetAcc_StatErr() { return  Acc_StatErr;}
	TH1F * GetAcc_SystErr()	{return Acc_SystErr;}
	TH1F * GetCounts_Err() {return Counts_Err;}

	TH1F * GetCounts() {return Counts;}


	TH1F * Eval_FluxRatio(Flux * Denominator,std::string name);

};


#endif
