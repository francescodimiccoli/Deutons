#ifndef ACCEPTANCE_H
#define ACCEPTANCE_H

#include "Tool.h"
#include "Efficiency.h"
#include "Globals.h"
#include "EffCorr.h"
#include "Unfolding.h"
#include "rundb.h"

struct MCPar{
	float Rmin,Rmax,Trigrate,gen_factor,art_ratio;
	long int tot_ev,tot_trig;
	std::string filename;
	std::string runlist;
	void Eval_trigrate();
};



class Acceptance : public Tool{

	private:
	Efficiency * FullSetEff;
	Efficiency * FullSetEff_gen;
	Efficiency * For_Acceptance;
	std::string directory;
	std::string basename;
	std::string cut;
	Binning bins;
	Binning binsfu;
	MCPar param;

	TH2F * Migr_rig;
	TH2F * Migr_beta;
	TH2F * Migr_R;
	TH2F * Migr_B;



	std::vector<EffCorr*> EfficiencyCorrections;
	std::vector<EffCorr*> EfficiencyFromData;

	FileSaver finalhistos;

	TH1F * EffAcceptance  =0x0;
	TH1F * EffAcceptanceMC=0x0;

	public:
	//standard constructor
	Acceptance(FileSaver File, std::string Basename, std::string Directory,std::string Cut_before,std::string Cut_after,Binning Bins,Binning BinsForUnfolding){
		FullSetEff     = new Efficiency(File, (Basename+"_FullSetMC" ).c_str(),Directory,Bins, Cut_before.c_str(),Cut_after.c_str(),false);
		FullSetEff_gen     = new Efficiency(File, (Basename+"_FullSetMCgen" ).c_str(),Directory,Bins, Cut_before.c_str(),Cut_after.c_str(),false);
	        For_Acceptance = new Efficiency(File,(Basename +"_For_Acceptance").c_str(),Directory,ForAcceptance,Cut_before.c_str(),Cut_before.c_str(),false);
      		For_Acceptance->SetNotWeightedMC();
 		//FullSetEff->SetNotWeightedMC(); // da rimuovere quando avrai vera MC di DEUTONI
		FullSetEff_gen->SetNotWeightedMC(); 

		directory = Directory;
		basename = Basename;
		cut = Cut_after;

		EffAcceptance = new TH1F((Basename +"_Eff_Acceptance").c_str(),(Basename +"_Eff_Acceptance").c_str(),Bins.size(),0,Bins.size());
		EffAcceptanceMC = new TH1F((Basename +"_Eff_AcceptanceMC").c_str(),(Basename +"_Eff_AcceptanceMC").c_str(),Bins.size(),0,Bins.size());
		bins = Bins;

		binsfu=BinsForUnfolding;

		float rig_gen[binsfu.size()+1];
		float rig_meas[binsfu.size()+1];
		float beta_meas[binsfu.size()+1];


		for(int i=0;i<binsfu.size()+1;i++) { 
			rig_gen[i]=binsfu.RigTOIBins()[i];	
			rig_meas[i]=binsfu.RigBins()[i];	
			beta_meas[i]=binsfu.BetaBins()[i];
		}	
		finalhistos = File;
		TFile * file = finalhistos.GetFile();
		if(file) 
		{
			Migr_rig = (TH2F*) file->Get((directory + "/" + basename + "/"+ Basename +"_UnfoldingRig").c_str());
                       	Migr_beta= (TH2F*) file->Get((directory + "/" + basename + "/"+ Basename +"_Unfolding").c_str());
		        Migr_R= (TH2F*) file->Get((directory + "/" + basename + "/"+ Basename +"_R").c_str());
			Migr_B= (TH2F*) file->Get((directory + "/" + basename + "/"+ Basename +"_B").c_str());
			int i = slicenormalizex(Migr_rig);
			int j = slicenormalizex(Migr_beta);
		}
		else{
			Migr_rig = new TH2F((Basename +"_UnfoldingRig").c_str(),(Basename+"rig;R_{GEN}[GV];R_{meas}[GV]").c_str(),binsfu.size(),rig_gen,binsfu.size(),rig_meas);
			Migr_beta= new TH2F((Basename +"_Unfolding").c_str(),(Basename+"beta;R_{GEN}[GV];#beta_{meas}").c_str(),binsfu.size(),rig_gen,binsfu.size(),beta_meas);
			Migr_R= new TH2F((Basename +"_R").c_str(),(Basename+"mass;R_{GEN}[GV];R_{meas}").c_str(),200,0,30,200,0,30);
			Migr_B= new TH2F((Basename +"_B").c_str(),(Basename+"mass;R_{GEN}[GV];#beta_{meas}").c_str(),200,0,30,200,0.8,1.2);
		}
	}


	//reading constructor

	Acceptance(FileSaver FileRes, std::string Basename, std::string Directory,Binning Bins,Binning BinsForUnfolding){
		FullSetEff     = new Efficiency(FileRes, (Basename+"_FullSetMC" ).c_str(),Directory,Bins,false);
		FullSetEff_gen     = new Efficiency(FileRes, (Basename+"_FullSetMCgen" ).c_str(),Directory,Bins,false);
	        For_Acceptance = new Efficiency(FileRes,(Basename +"_For_Acceptance").c_str(),Directory,ForAcceptance,false);
      		For_Acceptance->SetNotWeightedMC();
		//FullSetEff->SetNotWeightedMC(); // da rimuovere quando avrai vera MC di DEUTONI
 		FullSetEff_gen->SetNotWeightedMC();
 
 		directory = Directory;
		basename = Basename;
		TFile * fileres = FileRes.GetFile();

		if(fileres){
			EffAcceptance   = (TH1F*) fileres->Get((directory + "/" + basename + "/"+ Basename +"_Eff_Acceptance").c_str());
			EffAcceptanceMC = (TH1F*) fileres->Get((directory + "/" + basename + "/"+ Basename +"_Eff_AcceptanceMC").c_str());
	
			Migr_rig = (TH2F*) fileres->Get((directory + "/" + basename + "/"+ Basename +"_UnfoldingRig").c_str());
                       	Migr_beta= (TH2F*) fileres->Get((directory + "/" + basename + "/"+ Basename +"_Unfolding").c_str());
		        Migr_R= (TH2F*) fileres->Get((directory + "/" + basename + "/"+ Basename +"_R").c_str());
			Migr_B= (TH2F*) fileres->Get((directory + "/" + basename + "/"+ Basename +"_B").c_str());


		}
		bins = Bins;
		binsfu= BinsForUnfolding;

	}



	void Set_MCPar(float rmin, float rmax, float Gen_factor, std::string Filename, std::string Runlist, float Art_ratio=1);
	void ApplyEfficCorr(EffCorr * Correction);
	void ApplyEfficFromData(EffCorr * Correction);

	bool ReinitializeHistos(bool refill){
		bool checkifsomeismissing=false;
		bool allfound=true;
		if(!(FullSetEff -> ReinitializeHistos(refill))) checkifsomeismissing   = true;
		if(!(FullSetEff_gen -> ReinitializeHistos(refill))) checkifsomeismissing   = true;
		if(!(For_Acceptance -> ReinitializeHistos(refill))) checkifsomeismissing   = true;

		if(!Migr_rig||!Migr_beta){
			float rig_gen[binsfu.size()+1];
			float rig_meas[binsfu.size()+1];
			float beta_meas[binsfu.size()+1];


			for(int i=0;i<binsfu.size()+1;i++) { 
				rig_gen[i]=binsfu.RigTOIBins()[i];	
				rig_meas[i]=binsfu.RigBins()[i];	
				beta_meas[i]=binsfu.BetaBins()[i];
			}	
			Migr_rig = new TH2F((basename +"_UnfoldingRig").c_str(),(basename+"rig;R_{GEN}[GV];R_{meas}[GV]").c_str(),binsfu.size(),rig_gen,binsfu.size(),rig_meas);
			Migr_beta= new TH2F((basename +"_Unfolding").c_str(),(basename+"beta;R_{GEN}[GV];#beta_{meas}").c_str(),binsfu.size(),rig_gen,binsfu.size(),beta_meas);
			Migr_R= new TH2F((basename +"_R").c_str(),(basename+"mass;R_{GEN}[GV];R_{meas}").c_str(),200,0,30,200,0,30);
			Migr_B= new TH2F((basename +"_B").c_str(),(basename+"mass;R_{GEN}[GV];#beta_{meas}").c_str(),200,0,30,200,0.8,1.2);
			} 		
		if(checkifsomeismissing||refill) allfound=false;
		return allfound;
	}
	void FillEventByEventMC(Variables * vars, float (*var) (Variables * vars), float (*discr_var) (Variables * vars)){
		FullSetEff 	-> FillEventByEventMC(vars,var,discr_var);

		if(FullSetEff_gen->GetBins().IsUsingBetaEdges()) FullSetEff_gen 	-> FillEventByEventMC(vars,GetBetaGen,GetBetaGen);
		else FullSetEff_gen        -> FillEventByEventMC(vars,GetGenMomentum,GetGenMomentum);

	//	For_Acceptance  -> FillEventByEventMC(vars,GetGenMomentum,GetGenMomentum);

		if(ApplyCuts(cut,vars)){
			Migr_rig->Fill(GetGenMomentum(vars),GetRigidity(vars),vars->mcweight);
			Migr_beta->Fill(GetGenMomentum(vars),var(vars),vars->mcweight);
			Migr_R->Fill(GetGenMomentum(vars),GetRigidity(vars),vars->mcweight);
			Migr_B->Fill(GetGenMomentum(vars),var(vars),vars->mcweight);
		}		

	}
	void FillTotalMCEvents(Variables * vars, float (*var) (Variables * vars), float (*discr_var) (Variables * vars)){
		For_Acceptance  -> FillEventByEventMC(vars,GetGenMomentum,GetGenMomentum);
	}

	void Save();
	void SaveResults(FileSaver finalhistos);
	void SetDefaultOutFile(FileSaver FinalHistos);
	void EvalEffAcc(int timeindex,float SF);
	TH1F * GetEffAcc() {return EffAcceptance;}
	TH1F * GetEffAccMC() {return EffAcceptanceMC;}
	TH2F * GetMigr_rig() {return Migr_rig;}	
	TH2F * GetMigr_beta() {return Migr_beta;}	

};

#endif
