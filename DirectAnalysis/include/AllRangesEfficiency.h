struct AllRangesEfficiency{

	private:
	bool refill=false;
	BadEventSimulator * BadEvSimTOF=0x0;
	BadEventSimulator * BadEvSimNaF=0x0;
	BadEventSimulator * BadEvSimAgl=0x0;

	public:
	Efficiency * EffTOF;
	Efficiency * EffNaF;
	Efficiency * EffAgl;
	

	AllRangesEfficiency(FileSaver  File, std::string Basename,std::string Directory,std::string Cut_before,std::string Cut_after,bool Refill=false){
	
		EffTOF = new Efficiency(File,(Basename+"_TOF").c_str(),Directory,ToFDB,Cut_before,Cut_after);	
		EffNaF = new Efficiency(File,(Basename+"_NaF").c_str(),Directory,NaFDB,(Cut_before+"&IsFromNaF").c_str(),(Cut_after+"&IsFromNaF").c_str());
		EffAgl = new Efficiency(File,(Basename+"_Agl").c_str(),Directory,AglDB,(Cut_before+"&IsFromAgl").c_str(),(Cut_after+"&IsFromAgl").c_str());
		refill=Refill;
	}

	AllRangesEfficiency(FileSaver  File, std::string Basename,std::string Directory,std::string Cut_beforeTOF,std::string Cut_beforeNaF,std::string Cut_beforeAgl,std::string Cut_afterTOF,std::string Cut_afterNaF,std::string Cut_afterAgl,bool Refill=false){
	
		EffTOF = new Efficiency(File,(Basename+"_TOF").c_str(),Directory,ToFDB,Cut_beforeTOF,Cut_afterTOF);	
		EffNaF = new Efficiency(File,(Basename+"_NaF").c_str(),Directory,NaFDB,Cut_beforeNaF,Cut_afterNaF);
		EffAgl = new Efficiency(File,(Basename+"_Agl").c_str(),Directory,AglDB,Cut_beforeAgl,Cut_afterAgl);
		refill=Refill;
	}

	AllRangesEfficiency(FileSaver  File, std::string Basename,std::string Directory){
	
		EffTOF = new Efficiency(File,(Basename+"_TOF").c_str(),Directory,ToFDB);	
		EffNaF = new Efficiency(File,(Basename+"_NaF").c_str(),Directory,NaFDB);
		EffAgl = new Efficiency(File,(Basename+"_Agl").c_str(),Directory,AglDB);
	}

	void SetUpBadEventSimulator(BadEventSimulator * SimTOF, BadEventSimulator * SimNaF, BadEventSimulator * SimAgl) {
		BadEvSimTOF = SimTOF;
 		BadEvSimNaF = SimNaF;
	        BadEvSimAgl = SimAgl;				
	};
	
	void LoadEventIntoBadEvSim(Variables * vars) {
		EffTOF->LoadEventIntoBadEvSim(vars);
		EffNaF->LoadEventIntoBadEvSim(vars);
		EffAgl->LoadEventIntoBadEvSim(vars);
	}
	
	bool ReinitializeHistos(bool refill){
		bool checkifsomeismissing=false;
		bool allfound=true;
		if(!(EffTOF -> ReinitializeHistos(refill))) checkifsomeismissing = true;
                if(!(EffNaF -> ReinitializeHistos(refill))) checkifsomeismissing = true;
	        if(!(EffAgl -> ReinitializeHistos(refill))) checkifsomeismissing = true;
		if(checkifsomeismissing||refill) allfound=false;
		return allfound;
	}

	void Fill(TTree * tree, Variables * vars){
		EffTOF -> Fill(tree,vars,GetBetaGen,refill);
                EffNaF -> Fill(tree,vars,GetBetaGen,refill);
	        EffAgl -> Fill(tree,vars,GetBetaGen,refill);
	}

	void FillEventByEventMC(Variables * vars, float (*var) (Variables * vars), float (*discr_var) (Variables * vars)){
		EffTOF -> FillEventByEventMC(vars,var,discr_var);
                EffNaF -> FillEventByEventMC(vars,var,discr_var);
	        EffAgl -> FillEventByEventMC(vars,var,discr_var);
	}
	
	void FillEventByEventData(Variables * vars, float (*var) (Variables * vars), float (*discr_var) (Variables * vars)){
		EffTOF -> FillEventByEventData(vars,var,discr_var);
                EffNaF -> FillEventByEventData(vars,var,discr_var);
	        EffAgl -> FillEventByEventData(vars,var,discr_var);
	}
       
	 void Save(FileSaver finalHistos){
                EffTOF -> Save(finalHistos);
                EffNaF -> Save(finalHistos);
                EffAgl -> Save(finalHistos);
        }

	void Eval_Efficiency(){
                EffTOF -> Eval_Efficiency();
                EffNaF -> Eval_Efficiency();
                EffAgl -> Eval_Efficiency();
        }
	
	void Eval_FittedEfficiency(){
                EffTOF -> Eval_FittedEfficiency();
                EffNaF -> Eval_FittedEfficiency();
                EffAgl -> Eval_FittedEfficiency();
        }
	void SaveResults(FileSaver finalResults){
                EffTOF -> SaveResults(finalResults);
		EffNaF -> SaveResults(finalResults);
                EffAgl -> SaveResults(finalResults);
        }

	void ComposeEfficiency( AllRangesEfficiency * Second){
		EffTOF->ComposeEfficiency(Second->EffTOF);
		EffNaF->ComposeEfficiency(Second->EffNaF);
		EffAgl->ComposeEfficiency(Second->EffAgl);
	}

	void CloneEfficiency( AllRangesEfficiency * Second){
		EffTOF->CloneEfficiency(Second->EffTOF);
		EffNaF->CloneEfficiency(Second->EffNaF);
		EffAgl->CloneEfficiency(Second->EffAgl);
	}


};


