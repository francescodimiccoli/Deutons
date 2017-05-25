

class Efficiency{

	private:

	TH1F * before;
	TH1F * after;
	TH1F * Eff;
	FileSaver  file;
	Binning bins;
	std::string cut_before;
	std::string cut_after;
	std::string basename;
	std::string directory;
	public:

	Efficiency(FileSaver  File, std::string Basename,std::string Directory, Binning Bins, std::string Cut_before,std::string Cut_after){
		file=File;
		bins=Bins;
		basename=Basename;
		cut_before=Cut_before;
		cut_after=Cut_after;
		directory=Directory;
		before = new TH1F((basename+"_before").c_str(),(basename+"_before").c_str(),bins.size(),0,bins.size());
		after  = new TH1F((basename+"_after" ).c_str(),(basename+"_after" ).c_str(),bins.size(),0,bins.size());
	};

	Efficiency(FileSaver  File, std::string Basename,std::string Directory, Binning Bins ){
	file=File;
                bins=Bins;
                basename=Basename;
		before =(TH1F *) file.Get((directory+"/"+basename+"/"+basename+"_before").c_str());
		after  =(TH1F *) file.Get((directory+"/"+basename+"/"+basename+"_after").c_str());
		Eff    =(TH1F *) file.Get((directory+"/"+basename+"/"+basename+"_Eff").c_str());

	}

	void Fill( TNtuple * treeMC, Variables * vars, float (*discr_var) (Variables * vars),bool refill=false);
	void FillEventByEvent(float discr_var, bool CUT_BEFORE, bool CUT_AFTER, float weight);
	void Save(FileSaver finalhistos);
	void SaveResults(FileSaver finalhistos);
	void Eval_Efficiency();

	void ComposeEfficiency(Efficiency * Second);
	
	TH1F * GetEfficiency() {return Eff;}
	TH1F * GetBefore()     {return before;}
	TH1F * GetAfter()      {return after;}
	
};

void Efficiency::ComposeEfficiency(Efficiency * Second){
	if(Eff){
		Eff->Sumw2();
		if(Second->GetEfficiency()) Eff->Multiply(Second->GetEfficiency());
	}
	return;
}

void Efficiency::Fill(TNtuple * tree, Variables * vars, float (*discr_var) (Variables * vars),bool refill){

	cout<<file.Get((directory+"/"+basename+"/"+basename+"_before").c_str())<<" "<<file.Get((directory+"/"+basename+"/"+basename+"_after").c_str())<<endl;

	if(( (TH1F*)  file.Get((directory+"/"+basename+"/"+basename+"_before").c_str()) &&
				(TH1F*)  file.Get((directory+"/"+basename+"/"+basename+"_after").c_str()) ) && !refill) {

		before =(TH1F *) file.Get((directory+"/"+basename+"/"+basename+"_before").c_str());
		after  =(TH1F *) file.Get((directory+"/"+basename+"/"+basename+"_after").c_str());

	}	 
	else {

		cout<<basename.c_str()<<" Filling ..."<< endl;
		vars->ReadAnalysisBranches(tree);
		for(int i=0;i<tree->GetEntries();i++){

			vars->AnalysisVariablseReset();		
			UpdateProgressBar(i, tree->GetEntries());
			tree->GetEvent(i);
			FillEventByEvent(discr_var(vars),ApplyCuts(cut_before,vars),ApplyCuts(cut_after,vars),vars->mcweight);
		}

	}
	return;
}



void Efficiency::FillEventByEvent(float discr_var, bool CUT_BEFORE, bool CUT_AFTER, float weight){

	int kbin;
        kbin =  bins.GetBin(discr_var);
	if(kbin>0){
		if(CUT_BEFORE) before->Fill(kbin,weight);
		if(CUT_AFTER ) after->Fill(kbin,weight);
	}
        return;
}


void Efficiency::Save(FileSaver finalhistos){
	finalhistos.Add(before);
	finalhistos.Add(after); 	
        finalhistos.writeObjsInFolder((directory+"/"+basename).c_str());
}


void Efficiency::Eval_Efficiency(){
	Eff = (TH1F *)after->Clone();
	Eff -> Sumw2();
	Eff -> Divide(before);
	Eff ->SetName((basename+"_Eff").c_str());
	Eff ->SetTitle((basename+" Efficiency").c_str());
	return;
}

void Efficiency::SaveResults(FileSaver finalhistos){
	finalhistos.Add(Eff); 	
        finalhistos.writeObjsInFolder((directory+"/"+basename).c_str());
}



