#ifndef ACCEPTANCE_H
#define ACCEPTANCE_H

#include "Tool.h"
#include "Efficiency.h"
#include "Globals.h"
#include "EffCorr.h"
#include "rundb.h"

struct MCPar{
	float Rmin,Rmax,Trigrate,gen_factor,art_ratio;
	long int tot_ev,tot_trig;
	std::string filename;
	void Eval_trigrate();
};



class Acceptance : public Tool{

	private:
	Efficiency * FullSetEff;
	Efficiency * For_Acceptance;
	std::string directory;
	std::string basename;
	Binning bins;
	MCPar param;

	TH1F * Acc_StatErr=0x0;
	TH1F * Acc_SystErr=0x0;

	std::vector<EffCorr*> EfficiencyCorrections;
	std::vector<EffCorr*> EfficiencyFromData;

	FileSaver finalhistos;

	TH1F * EffAcceptance  =0x0;
	TH1F * EffAcceptanceMC=0x0;

	public:
	//standard constructor
	Acceptance(FileSaver File, std::string Basename, std::string Directory,std::string Cut_before,std::string Cut_after,Binning Bins){
		FullSetEff     = new Efficiency(File, (Basename+"_FullSetMC" ).c_str(),Directory,Bins, Cut_before.c_str(),Cut_after.c_str());
	        For_Acceptance = new Efficiency(File,(Basename +"_For_Acceptance").c_str(),Directory,ForAcceptance,Cut_before.c_str(),Cut_before.c_str());
      		For_Acceptance->SetNotWeightedMC();
 		FullSetEff->SetNotWeightedMC(); // da rimuovere quando avrai vera MC di DEUTONI
		directory = Directory;
		basename = Basename;

		EffAcceptance = new TH1F((Basename +"_Eff_Acceptance").c_str(),(Basename +"_Eff_Acceptance").c_str(),Bins.size(),0,Bins.size());
		EffAcceptanceMC = new TH1F((Basename +"_Eff_AcceptanceMC").c_str(),(Basename +"_Eff_AcceptanceMC").c_str(),Bins.size(),0,Bins.size());
		Acc_StatErr = new TH1F((Basename +"_Acc_StatErr").c_str(),(Basename +"_Acc_StatErr").c_str(),Bins.size(),0,Bins.size());
		Acc_SystErr = new TH1F((Basename +"_Acc_SystErr").c_str(),(Basename +"_Acc_SystErr").c_str(),Bins.size(),0,Bins.size());
		bins = Bins;
	}
	//reading constructor
	Acceptance(FileSaver FileRes, std::string Basename, std::string Directory,Binning Bins){
		FullSetEff     = new Efficiency(FileRes, (Basename+"_FullSetMC" ).c_str(),Directory,Bins);
	        For_Acceptance = new Efficiency(FileRes,(Basename +"_For_Acceptance").c_str(),Directory,ForAcceptance);
      		For_Acceptance->SetNotWeightedMC();
		FullSetEff->SetNotWeightedMC(); // da rimuovere quando avrai vera MC di DEUTONI
 		directory = Directory;
		basename = Basename;

		EffAcceptance = (TH1F*) FileRes.Get((directory + "/" + basename + "/"+ Basename +"_Eff_Acceptance").c_str());
		EffAcceptanceMC = (TH1F*) FileRes.Get((directory + "/" + basename + "/"+ Basename +"_Eff_AcceptanceMC").c_str());
		Acc_StatErr = (TH1F*) FileRes.Get((directory + "/" + basename + "/"+ Basename +"_Eff_Acceptance").c_str());
		Acc_SystErr = (TH1F*) FileRes.Get((directory + "/" + basename + "/"+ Basename +"_Eff_Acceptance").c_str());
		bins = Bins;
	}



	void Set_MCPar(float rmin, float rmax, float Gen_factor, std::string Filename, float Art_ratio=1);
	void ApplyEfficCorr(EffCorr * Correction);
	void ApplyEfficFromData(EffCorr * Correction);

	bool ReinitializeHistos(bool refill){
		bool checkifsomeismissing=false;
		bool allfound=true;
		if(!(FullSetEff -> ReinitializeHistos(refill))) checkifsomeismissing   = true;
		if(!(For_Acceptance -> ReinitializeHistos(refill))) checkifsomeismissing   = true;
		if(checkifsomeismissing||refill) allfound=false;
		return allfound;
	}
	void FillEventByEventMC(Variables * vars, float (*var) (Variables * vars), float (*discr_var) (Variables * vars)){
		FullSetEff 	-> FillEventByEventMC(vars,var,discr_var);
		For_Acceptance  ->FillEventByEventMC(vars,GetGenMomentum,GetGenMomentum);
	}
	void Save();
	void SaveResults(FileSaver finalhistos);
	void SetDefaultOutFile(FileSaver FinalHistos);
	void EvalEffAcc();
	TH1F * GetEffAcc() {return EffAcceptance;}
	TH1F * GetEffAccMC() {return EffAcceptanceMC;}
	TH1F * GetEffAcc_staterr() {return Acc_StatErr;}
	TH1F * GetEffAcc_systerr() {return Acc_SystErr;}
		
};

#endif
