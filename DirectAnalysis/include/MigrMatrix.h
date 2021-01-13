#ifndef MIGRMATRIX_H
#define MIGRMATRIX_H


#include "Tool.h"
#include "Globals.h"
#include "filesaver.h"

class MigrMatrix: public Tool{

	private:

	TH2F * Migr_rig;
	TH2F * Migr_beta;
	Binning bins;
	FileSaver file;

	public:
	MigrMatrix(FileSaver File,std::string name,  std::string cut, Binning Bins){ ;
		bins=Bins;
		file = File;

		float rig_gen[bins.size()+1];
		float rig_meas[bins.size()+1];
		float beta_meas[bins.size()+1];


		for(int i=0;i<bins.size()+1;i++) { 
			rig_gen[i]=Bins.RigTOIBins()[i];	
			rig_meas[i]=Bins.RigBins()[i];	
			beta_meas[i]=Bins.BetaBins()[i];
		}	


		Migr_rig = new TH2F((name+"rig" ).c_str(),(name+"rig;R_{GEN}[GV];R_{meas}[GV]").c_str(),bins.size(),rig_gen,bins.size(),rig_meas);
		Migr_beta= new TH2F((name+"beta").c_str(),(name+"beta;R_{GEN}[GV];#beta_{meas}").c_str(),bins.size(),rig_gen,bins.size(),rig_meas);
	}

	bool ReinitializeHistos(bool refill){
		bool checkifsomeismissing=false;
		bool allfound=true;
		TFile * f = TFile::Open
		checkifsomeismissing   = true;
		checkifsomeismissing   = true;
		if(checkifsomeismissing||refill) allfound=false;
		return allfound;
	}

	void FillEventByEventMC(Variables * vars, float (*var) (Variables * vars), float (*discr_var) (Variables * vars)){
		FullSetEff 	-> FillEventByEventMC(vars,var,discr_var);

		if(FullSetEff_gen->GetBins().IsUsingBetaEdges()) FullSetEff_gen 	-> FillEventByEventMC(vars,GetBetaGen,GetBetaGen);
		else FullSetEff_gen        -> FillEventByEventMC(vars,GetGenMomentum,GetGenMomentum);

		For_Acceptance  -> FillEventByEventMC(vars,GetGenMomentum,GetGenMomentum);
	
	}
	
	void Save();
	

};




#endif
