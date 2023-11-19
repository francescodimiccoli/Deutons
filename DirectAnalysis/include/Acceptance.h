#ifndef ACCEPTANCE_H
#define ACCEPTANCE_H

#include "Tool.h"
#include "Efficiency.h"
#include "Globals.h"
#include "EffCorr.h"
#include "Unfolding.h"
#include "rundb.h"

struct MCPar{
	float Rmin,Rmax,Trigrate,gen_factor,art_ratio,compact_event_ratio;
	long int tot_ev,tot_trig;
	long int TOT_EV,TOT_TRIG;
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
	std::string dataeff_cut;
	
	Binning bins;
	Binning binsfu;
	MCPar param;

	TH2F * EnLoss;
	TH2F * Migr_Unf;
	TH1F * Migr_True;
	TH1F * Migr_Meas;

	bool IsModeled=false;

	std::vector<EffCorr*> EfficiencyCorrections;
	std::vector<EffCorr*> EfficiencyFromData;
	std::vector<float> EffCorrshift;
	std::vector<float> EffDatashift;


	FileSaver finalhistos;

	TF1 * EffAcceptanceTime  	  =0x0;
	
	TH1F * EffAcceptance  	  =0x0;
	TH1F * EffAcceptanceMC	  =0x0;
	TH1F * EffAcceptance_gen  =0x0;
	TH1F * EffAcceptanceMC_gen=0x0;
	TH1F * EffAcceptanceMC_gen_raw=0x0;


	TH1F * EffAcceptance_StatErr  	  =0x0;
	TH1F * EffAcceptance_SystErr  	  =0x0;
	

	//STEP correction
	bool Applystep=false;
	float rigstep=0;
	float stepsize=0;


	public:
	//standard constructor
	Acceptance(FileSaver File, std::string Basename, std::string Directory,std::string Cut_before,std::string Cut_after,Binning Bins,Binning BinsForUnfolding,std::string dataeffcut=""){
		FullSetEff     = new Efficiency(File, (Basename+"_FullSetMC" ).c_str(),Directory,Bins, Cut_before.c_str(),Cut_after.c_str(),false);
		FullSetEff_gen     = new Efficiency(File, (Basename+"_FullSetMCgen" ).c_str(),Directory,Bins, Cut_before.c_str(),Cut_after.c_str(),false);
	        For_Acceptance = new Efficiency(File,(Basename +"_For_Acceptance").c_str(),Directory,ForAcceptance,Cut_before.c_str(),Cut_before.c_str(),false);
      		For_Acceptance->SetNotWeightedMC();
 		//FullSetEff->SetNotWeightedMC(); // da rimuovere quando avrai vera MC di DEUTONI
		FullSetEff_gen->SetNotWeightedMC(); 

		directory = Directory;
		basename = Basename;
		cut = Cut_after;
		dataeff_cut=dataeffcut;
		EffAcceptance = new TH1F((Basename +"_Eff_Acceptance").c_str(),(Basename +"_Eff_Acceptance").c_str(),Bins.size(),0,Bins.size());
		EffAcceptanceMC = new TH1F((Basename +"_Eff_AcceptanceMC").c_str(),(Basename +"_Eff_AcceptanceMC").c_str(),Bins.size(),0,Bins.size());
		EffAcceptance_gen = new TH1F((Basename +"_Eff_Acceptance_gen").c_str(),(Basename +"_Eff_Acceptance_gen").c_str(),Bins.size(),0,Bins.size());
		EffAcceptanceMC_gen = new TH1F((Basename +"_Eff_AcceptanceMC_gen").c_str(),(Basename +"_Eff_AcceptanceMC_gen").c_str(),Bins.size(),0,Bins.size());
		EffAcceptanceMC_gen_raw = new TH1F((Basename +"_Eff_AcceptanceMC_gen_raw").c_str(),(Basename +"_Eff_AcceptanceMC_gen_raw").c_str(),Bins.size(),0,Bins.size());
		bins = Bins;

		binsfu=BinsForUnfolding;

		float rig_gen[binsfu.size()+1];
		float rig_meas[bins.size()+1];
		float beta_meas[bins.size()+1];
		float beta_gen[binsfu.size()+1];
		float ekin_meas[bins.size()+1];
		float ekin_true[binsfu.size()+1]; 

		for(int i=0;i<binsfu.size()+1;i++) { 
			rig_gen[i]	=0;
			rig_meas[i]	=0;;
			beta_meas[i]	=0;
			beta_gen[i]	=0;
			ekin_true[i]	=0;;
			ekin_meas[i]	=0;
			rig_gen[i]=binsfu.RigBins()[i];	
			rig_meas[i]=bins.RigBins()[i];	
			beta_meas[i]=bins.BetaBins()[i];
			beta_gen[i]=bins.BetaTOIBins()[i];
			ekin_true[i]=binsfu.EkPerMasTOIBins()[i];
			ekin_meas[i]=bins.EkPerMasTOIBins()[i];
	}

//		for(int i=0;i<bins.size()+1;i++) ekin_meas[i]=bins.EkPerMasBins()[i];
			
		finalhistos = File;
		TFile * file = finalhistos.GetFile();
		if(file) 
		{
			cout<<basename<<endl;
			EnLoss= (TH2F*) file->Get((directory + "/" + basename + "/"+ Basename +"_EnLoss").c_str());
			Migr_Unf= (TH2F*) file->Get((directory + "/" + basename + "/"+ Basename +"_RooUnfold").c_str());
			Migr_True= (TH1F*) file->Get((directory + "/" + basename + "/"+ Basename +"_RooTrue").c_str());
			Migr_Meas= (TH1F*) file->Get((directory + "/" + basename + "/"+ Basename +"_RooMeas").c_str());
		//	int i = slicenormalizex(Migr_rig);
		//	int j = slicenormalizex(Migr_beta);
		}
		else{
			EnLoss= new TH2F((Basename +"_EnLoss").c_str(),(Basename+"mass;1/#beta_{GEN}[GV];1/#beta_{meas}").c_str(),400,0.8,2.,400,0.8,2.);
			if(Bins.IsUsingBetaEdges()){
				Migr_Unf= new TH2F((Basename +"_RooUnfold").c_str(),(Basename+"beta;Ekin_{meas}[GeV/n];Ekin_{GEN}").c_str(),binsfu.size(),ekin_true,binsfu.size(),ekin_true);
				Migr_True= new TH1F((Basename +"_RooTrue").c_str(),(Basename+"beta;Ekin_{GEN}[GeV/nV];").c_str(),binsfu.size(),ekin_true);
				Migr_Meas= new TH1F((Basename +"_RooMeas").c_str(),(Basename+"beta;Ekin_{meas}[GeV/n];").c_str(),binsfu.size(),ekin_true);
			}
			else{
					for(int i=0;i<binsfu.size();i++) {
						rig_gen[i]=binsfu.RigBins()[i];	
						cout<<basename<<" RGEN "<<rig_gen[i]<<" "<<binsfu.RigBins()[i]<<endl;
					}	
				Migr_Unf= new TH2F((Basename +"_RooUnfold").c_str(),(Basename+"beta;R_{meas}[GV];R_{GEN}").c_str(),binsfu.size(),rig_gen,binsfu.size(),rig_gen);
				Migr_True= new TH1F((Basename +"_RooTrue").c_str(),(Basename+"R;R_{GEN}[GV];").c_str(),binsfu.size(),rig_gen);
				Migr_Meas= new TH1F((Basename +"_RooMeas").c_str(),(Basename+"R;R_{meas}[GV];").c_str(),binsfu.size(),rig_gen);
			}

		}
	
	}


	//reading constructor

	Acceptance(FileSaver FileRes, std::string Basename, std::string Directory,Binning Bins,Binning BinsForUnfolding){
		FullSetEff     = new Efficiency(FileRes, (Basename+"_FullSetMC" ).c_str(),Directory,Bins);
		FullSetEff_gen     = new Efficiency(FileRes, (Basename+"_FullSetMCgen" ).c_str(),Directory,Bins);
	        For_Acceptance = new Efficiency(FileRes,(Basename +"_For_Acceptance").c_str(),Directory,ForAcceptance);
      		For_Acceptance->SetNotWeightedMC();
		//FullSetEff->SetNotWeightedMC(); // da rimuovere quando avrai vera MC di DEUTONI
 		FullSetEff_gen->SetNotWeightedMC();
 
 		directory = Directory;
		basename = Basename;
		TFile * fileres = FileRes.GetFile();

		if(fileres){
			EffAcceptance       = (TH1F*) fileres->Get((directory + "/" + basename + "/"+ Basename +"_Eff_Acceptance").c_str());
			EffAcceptanceMC     = (TH1F*) fileres->Get((directory + "/" + basename + "/"+ Basename +"_Eff_AcceptanceMC").c_str());
			EffAcceptance_gen   = (TH1F*) fileres->Get((directory + "/" + basename + "/"+ Basename +"_Eff_Acceptance_gen").c_str());
			EffAcceptanceMC_gen = (TH1F*) fileres->Get((directory + "/" + basename + "/"+ Basename +"_Eff_AcceptanceMC_gen").c_str());
			EffAcceptanceMC_gen_raw = (TH1F*) fileres->Get((directory + "/" + basename + "/"+ Basename +"_Eff_AcceptanceMC_gen_raw").c_str());
	

			EffAcceptance_StatErr       = (TH1F*) fileres->Get((directory + "/" + basename + "/"+ Basename +"_Eff_Acceptance_StatErr").c_str());
			EffAcceptance_SystErr       = (TH1F*) fileres->Get((directory + "/" + basename + "/"+ Basename +"_Eff_Acceptance_SystErr").c_str());
		

		    	EnLoss= (TH2F*) fileres->Get((directory + "/" + basename + "/"+ Basename +"_EnLoss").c_str());
			Migr_Unf= (TH2F*) fileres->Get((directory + "/" + basename + "/"+ Basename +"_RooUnfold").c_str());
			Migr_True= (TH1F*) fileres->Get((directory + "/" + basename + "/"+ Basename +"_RooTrue").c_str());
			Migr_Meas= (TH1F*) fileres->Get((directory + "/" + basename + "/"+ Basename +"_RooMeas").c_str());
		   

		}
		bins = Bins;
		binsfu= BinsForUnfolding;

	}


	void ApplyStepCorr(float R,float step){ Applystep=true; rigstep=R; stepsize=step;}

	void Set_MCPar(float rmin, float rmax, float Gen_factor, std::string Filename, std::string Runlist, float Compact_event_ratio, float Art_ratio=1);
	void ApplyEfficCorr(EffCorr * Correction,float shift=0);
	void ApplyEfficFromData(EffCorr * Correction,float shift=0);

	bool ReinitializeHistos(bool refill){
		bool checkifsomeismissing=false;
		bool allfound=true;
		if(!(FullSetEff -> ReinitializeHistos(refill))) checkifsomeismissing   = true;
		if(!(FullSetEff_gen -> ReinitializeHistos(refill))) checkifsomeismissing   = true;
		if(!(For_Acceptance -> ReinitializeHistos(refill))) checkifsomeismissing   = true;

		if(!Migr_Unf||!Migr_True||!Migr_Meas){

		float rig_gen[binsfu.size()+1];
		float rig_meas[binsfu.size()+1];
		float beta_meas[bins.size()+1];
		float beta_gen[binsfu.size()+1];
		float ekin_meas[bins.size()+1];
		float ekin_true[binsfu.size()+1]; 

		for(int i=0;i<binsfu.size()+1;i++) { 
			rig_gen[i]	=0;
			rig_meas[i]	=0;;
			beta_meas[i]	=0;
			beta_gen[i]	=0;
			ekin_true[i]	=0;;
			ekin_meas[i]	=0;
			rig_gen[i]=binsfu.RigBins()[i];	
			rig_meas[i]=binsfu.RigBins()[i];	
			beta_meas[i]=binsfu.BetaBins()[i];
			beta_gen[i]=binsfu.BetaTOIBins()[i];
			ekin_true[i]=binsfu.EkPerMasTOIBins()[i];

		}	

		for(int i=0;i<bins.size()+1;i++) ekin_meas[i]=bins.EkPerMasBins()[i];
			
			EnLoss= new TH2F((basename +"_EnLoss").c_str(),(basename+"mass;1/#beta_{GEN}[GV];1/#beta_{meas}").c_str(),400,0.8,2.,400,0.8,2.);
			if(bins.IsUsingBetaEdges()){
				Migr_Unf= new TH2F((basename +"_RooUnfold").c_str(),(basename+"beta;Ekin_{meas}[GeV/n];#Ekin_{GEN}").c_str(),binsfu.size(),ekin_true,binsfu.size(),ekin_true);
				Migr_True= new TH1F((basename +"_RooTrue").c_str(),(basename+"beta;Ekin_{GEN}[GeV/nV];").c_str(),binsfu.size(),ekin_true);
				Migr_Meas= new TH1F((basename +"_RooMeas").c_str(),(basename+"beta;Ekin_{meas}[GeV/n];").c_str(),binsfu.size(),ekin_true);
			}
			else{
				for(int i=0;i<binsfu.size();i++) {
						rig_gen[i]=binsfu.RigBins()[i];	
						cout<<basename<<" RGEN "<<rig_gen[i]<<" "<<binsfu.RigBins()[i]<<endl;
					}	
		
			Migr_Unf= new TH2F((basename +"_RooUnfold").c_str(),(basename+"R;R_{meas}[GV];R_{GEN}").c_str(),binsfu.size(),rig_gen,binsfu.size(),rig_gen);
				Migr_True= new TH1F((basename +"_RooTrue").c_str(),(basename+"R;R_{GEN}[GV];").c_str(),binsfu.size(),rig_gen);
				Migr_Meas= new TH1F((basename +"_RooMeas").c_str(),(basename+"R;R_{meas}[GV];").c_str(),binsfu.size(),rig_gen);
				cout<<Migr_Unf->GetNbinsX()<<endl;
		}

	} 		
		if(checkifsomeismissing||refill) allfound=false;
		return allfound;
	}
	
	void FillEventByEventMC(Variables * vars, float (*var) (Variables * vars), float (*discr_var) (Variables * vars)){
		FullSetEff 	-> FillEventByEventMC(vars,var,discr_var);

		if(FullSetEff_gen->GetBins().IsUsingBetaEdges()) FullSetEff_gen 	-> FillEventByEventMC(vars,GetBetaGen_cpct,GetBetaGen_cpct);
		else FullSetEff_gen        -> FillEventByEventMC(vars,GetGenRigidity,GetGenRigidity);

		For_Acceptance  -> FillEventByEventMC(vars,GetGenMomentum_10,GetGenMomentum_10);

		if(ApplyCuts((cut+"&"+dataeff_cut).c_str(),vars)){
			float Ekin_n = 0.938*((1/sqrt(1-pow(var(vars),2)))-1); 
			float betagen = GetBetaFromR(vars->Massa_gen,GetGenRigidity(vars),vars->Charge_gen);
			float Ekin_gen = 0.938*((1/sqrt(1-pow(betagen,2)))-1); 
			if(FullSetEff_gen->GetBins().IsUsingBetaEdges()){
				Migr_True->Fill(Ekin_gen,vars->mcweight);
				Migr_Unf->Fill(Ekin_n,Ekin_gen,vars->mcweight);	
				Migr_Meas->Fill(Ekin_n,vars->mcweight);
			}
			else{
				Migr_True->Fill(GetGenRigidity(vars),vars->mcweight);
				Migr_Unf->Fill(vars->R,GetGenRigidity(vars),vars->mcweight);	
				Migr_Meas->Fill(vars->R,vars->mcweight);
			}
			EnLoss->Fill(1/betagen,1/var(vars),vars->mcweight);
				
		}		

	}
	void FillTotalMCEvents(Variables * vars, float (*var) (Variables * vars), float (*discr_var) (Variables * vars)){
		For_Acceptance  -> FillEventByEventMC(vars,GetGenMomentum,GetGenMomentum);
	}

	void Save();
	void SaveResults(FileSaver finalhistos);
	void SetDefaultOutFile(FileSaver FinalHistos);
	void EvalEffAcc(int timeindex,float SF, bool IsHe=false);
	TH1F * GetEffAcc() {return EffAcceptance;}
	TH1F * GetEffAccMC() {return EffAcceptanceMC;}
	TH1F * GetEffAcc_gen() {return EffAcceptance_gen;}
	TH1F * GetEffAccMC_gen() {return EffAcceptanceMC_gen;}
	TH1F * GetEffAccMC_gen_raw() {return EffAcceptanceMC_gen_raw;}
	TH2F * GetMigr_Unf() {return Migr_Unf;}	
	TH1F * GetMigr_True() {return Migr_True;}	
	TH1F * GetMigr_Meas() {return Migr_Meas;}	


	TH1F * GetStat_Err(){return EffAcceptance_StatErr;}
	TH1F * GetSyst_Err(){return EffAcceptance_SystErr;}

	void Set_AcceptanceTime(TF1 * avg) {
		EffAcceptanceTime = (TF1 *) avg->Clone(); 
	}

	void SetModeled() {IsModeled=true;}


};

#endif
