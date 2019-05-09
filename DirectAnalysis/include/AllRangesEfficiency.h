#include "Efficiency.h"

struct AllRangesEfficiency : public Efficiency{

	private:
	bool refill=false;
	BadEventSimulator * BadEvSimTOF=0x0;
	BadEventSimulator * BadEvSimNaF=0x0;
	BadEventSimulator * BadEvSimAgl=0x0;

	public:
	Efficiency * EffTOF;
	Efficiency * EffNaF;
	Efficiency * EffAgl;


	AllRangesEfficiency(FileSaver  File, std::string Basename,std::string Directory,std::string Cut_before,std::string Cut_after,Binning binsTOF,Binning binsNaF,Binning binsAgl,bool Refill=false):Efficiency(File,(Basename+"_TOF").c_str(),Directory,binsTOF,Cut_before,Cut_after){
	
		EffTOF = new Efficiency(File,(Basename+"_TOF").c_str(),Directory,binsTOF,(Cut_before+"&IsOnlyFromToF").c_str(),(Cut_after+"&IsOnlyFromToF").c_str());	
		EffNaF = new Efficiency(File,(Basename+"_NaF").c_str(),Directory,binsNaF,(Cut_before+"&IsFromNaF"    ).c_str(),(Cut_after+"&IsFromNaF"    ).c_str());
		EffAgl = new Efficiency(File,(Basename+"_Agl").c_str(),Directory,binsAgl,(Cut_before+"&IsFromAgl"    ).c_str(),(Cut_after+"&IsFromAgl"    ).c_str());
		refill=Refill;
	}

	AllRangesEfficiency(FileSaver  File, std::string Basename,std::string Directory,std::string Cut_beforeTOF,std::string Cut_beforeNaF,std::string Cut_beforeAgl,std::string Cut_afterTOF,std::string Cut_afterNaF,std::string Cut_afterAgl,Binning binsTOF,Binning binsNaF,Binning binsAgl,bool Refill=false):Efficiency(File,(Basename+"_TOF").c_str(),Directory,binsTOF,Cut_beforeTOF,Cut_afterTOF){
	
		EffTOF = new Efficiency(File,(Basename+"_TOF").c_str(),Directory,binsTOF,Cut_beforeTOF,Cut_afterTOF);	
		EffNaF = new Efficiency(File,(Basename+"_NaF").c_str(),Directory,binsNaF,Cut_beforeNaF,Cut_afterNaF);
		EffAgl = new Efficiency(File,(Basename+"_Agl").c_str(),Directory,binsAgl,Cut_beforeAgl,Cut_afterAgl);
		refill=Refill;
	}

	AllRangesEfficiency(FileSaver  File, std::string Basename,std::string Directory,Binning binsTOF,Binning binsNaF,Binning binsAgl):Efficiency(File,(Basename+"_TOF").c_str(),Directory,binsTOF){
	
		EffTOF = new Efficiency(File,(Basename+"_TOF").c_str(),Directory,binsTOF);	
		EffNaF = new Efficiency(File,(Basename+"_NaF").c_str(),Directory,binsNaF);
		EffAgl = new Efficiency(File,(Basename+"_Agl").c_str(),Directory,binsAgl);
	}

	virtual void SetUpBadEventSimulator(BadEventSimulator * SimTOF, BadEventSimulator * SimNaF, BadEventSimulator * SimAgl) {
		BadEvSimTOF = SimTOF;
 		BadEvSimNaF = SimNaF;
	        BadEvSimAgl = SimAgl;				
	};
	
	virtual void LoadEventIntoBadEvSim(Variables * vars) {
		EffTOF->LoadEventIntoBadEvSim(vars);
		EffNaF->LoadEventIntoBadEvSim(vars);
		EffAgl->LoadEventIntoBadEvSim(vars);
	}
	
	virtual bool ReinitializeHistos(bool refill){
		bool checkifsomeismissing=false;
		bool allfound=true;
		if(!(EffTOF -> ReinitializeHistos(refill))) checkifsomeismissing = true;
                if(!(EffNaF -> ReinitializeHistos(refill))) checkifsomeismissing = true;
	        if(!(EffAgl -> ReinitializeHistos(refill))) checkifsomeismissing = true;
		if(checkifsomeismissing||refill) allfound=false;
		return allfound;
	}

	virtual void Fill(TTree * tree, Variables * vars){
		EffTOF -> Fill(tree,vars,GetBetaGen,refill);
                EffNaF -> Fill(tree,vars,GetBetaGen,refill);
	        EffAgl -> Fill(tree,vars,GetBetaGen,refill);
	}

	virtual void FillEventByEventMC(Variables * vars, float (*var) (Variables * vars), float (*discr_var) (Variables * vars)){
		EffTOF -> FillEventByEventMC(vars,var,discr_var);
                EffNaF -> FillEventByEventMC(vars,var,discr_var);
	        EffAgl -> FillEventByEventMC(vars,var,discr_var);
	}
	
	virtual void FillEventByEventData(Variables * vars, float (*var) (Variables * vars), float (*discr_var) (Variables * vars)){
		EffTOF -> FillEventByEventData(vars,var,discr_var);
                EffNaF -> FillEventByEventData(vars,var,discr_var);
	        EffAgl -> FillEventByEventData(vars,var,discr_var);
	}
       
	virtual void Save(){
                EffTOF -> Save();
                EffNaF -> Save();
                EffAgl -> Save();
        }

	 virtual void Eval_Efficiency(){
                EffTOF -> Eval_Efficiency();
                EffNaF -> Eval_Efficiency();
                EffAgl -> Eval_Efficiency();
        }
	
	virtual void Eval_FittedEfficiency(){
                EffTOF -> Eval_FittedEfficiency();
                EffNaF -> Eval_FittedEfficiency();
                EffAgl -> Eval_FittedEfficiency();
        }
	virtual void SaveResults(FileSaver finalResults){
                EffTOF -> SaveResults(finalResults);
		EffNaF -> SaveResults(finalResults);
                EffAgl -> SaveResults(finalResults);
        }

	virtual void ComposeEfficiency( AllRangesEfficiency * Second){
		EffTOF->ComposeEfficiency(Second->EffTOF);
		EffNaF->ComposeEfficiency(Second->EffNaF);
		EffAgl->ComposeEfficiency(Second->EffAgl);
	}

	virtual void CloneEfficiency( AllRangesEfficiency * Second){
		EffTOF->CloneEfficiency(Second->EffTOF);
		EffNaF->CloneEfficiency(Second->EffNaF);
		EffAgl->CloneEfficiency(Second->EffAgl);
	}


	virtual void Eval_StatError(){
		EffTOF->Eval_StatError();
		EffNaF->Eval_StatError();
		EffAgl->Eval_StatError();
	} 
	virtual void Eval_SystError( AllRangesEfficiency * First, AllRangesEfficiency * Second){
		EffTOF->Eval_SystError(First->EffTOF,Second->EffTOF);
		EffNaF->Eval_SystError(First->EffNaF,Second->EffNaF);
		EffAgl->Eval_SystError(First->EffAgl,Second->EffAgl);
	}

	 virtual void SetNotWeightedMC() {
		EffTOF->SetNotWeightedMC(); 
		EffNaF->SetNotWeightedMC();
        	EffAgl->SetNotWeightedMC(); 
	}
	void SetDefaultOutFile(FileSaver FinalHistos){
		 EffTOF->SetDefaultOutFile(FinalHistos);
		 EffNaF->SetDefaultOutFile(FinalHistos);
		 EffAgl->SetDefaultOutFile(FinalHistos);
	 }

};


