

class Efficiency{

	private:

	TH1 * before;
	TH1 * after;
	TH1 * Eff;

	TGraphErrors * FittedEff;
	FileSaver  file;
	Binning bins;
	std::string cut_before;
	std::string cut_after;
	std::string basename;
	std::string directory;
	bool fitrequested=false;

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

	Efficiency(FileSaver  File, std::string Basename,std::string Directory, Binning Bins, std::string Cut_before,std::string Cut_after,std::vector<float> LatZones){
		file=File;
		bins=Bins;
		basename=Basename;
		cut_before=Cut_before;
		cut_after=Cut_after;
		directory=Directory;
		before = new TH2F((basename+"_before").c_str(),(basename+"_before").c_str(),bins.size(),0,bins.size(),LatZones.size()-1,0,LatZones.size()-1);
		after  = new TH2F((basename+"_after" ).c_str(),(basename+"_after" ).c_str(),bins.size(),0,bins.size(),LatZones.size()-1,0,LatZones.size()-1);
	};



	Efficiency(FileSaver  File, std::string Basename,std::string Directory, Binning Bins ){
	file=File;
                bins=Bins;
                basename=Basename;
		directory=Directory;
		before =(TH1 *) file.Get((directory+"/"+basename+"/"+basename+"_before").c_str());
		after  =(TH1 *) file.Get((directory+"/"+basename+"/"+basename+"_after").c_str());
		Eff    =(TH1 *) file.Get((directory+"/"+basename+"/"+basename+"_Eff").c_str());
	}

	void Fill( TNtuple * treeMC, Variables * vars, float (*discr_var) (Variables * vars),bool refill=false,bool weight=true);
	void FillEventByEvent(float discr_var, bool CUT_BEFORE, bool CUT_AFTER, float weight);
	void FillEventByEventLatitude(float discr_var, bool CUT_BEFORE, bool CUT_AFTER,int latzone);	
	
	void Save(FileSaver finalhistos);
	void SaveResults(FileSaver finalhistos);
	void Eval_Efficiency();
	void Eval_FittedEfficiency();

	void CloneEfficiency(Efficiency * Second);
	void ComposeEfficiency(Efficiency * Second);
	
	TH1 * GetEfficiency() {return Eff;}
	TGraphErrors * GetFittedEfficiency() {return FittedEff;}
	TH1 * GetBefore()     {return before;}
	TH1 * GetAfter()      {return after;}

	std::string GetCut_Before()	{return cut_before;}	
	std::string GetCut_After()	{return cut_after;}
};

void Efficiency::ComposeEfficiency(Efficiency * Second){
	if(Eff){
		Eff->Sumw2();
		if(Second->GetEfficiency()) Eff->Multiply(Second->GetEfficiency());
	}
	return;
}

void Efficiency::CloneEfficiency(Efficiency * Second){

	before = (TH1F*) Second->GetBefore()->Clone();	
	after  = (TH1F*) Second->GetAfter()->Clone();

	before->SetName((basename+"_before").c_str());
	after ->SetName((basename+"_after").c_str());

	before->SetTitle((basename+"_before").c_str());
	after ->SetTitle((basename+"_after").c_str());
	
	return;
}

void Efficiency::Fill(TNtuple * tree, Variables * vars, float (*discr_var) (Variables * vars),bool refill,bool weight){

	cout<<file.Get((directory+"/"+basename+"/"+basename+"_before").c_str())<<" "<<file.Get((directory+"/"+basename+"/"+basename+"_after").c_str())<<endl;
	if(( (TH1F*)  file.Get((directory+"/"+basename+"/"+basename+"_before").c_str()) &&
				(TH1F*)  file.Get((directory+"/"+basename+"/"+basename+"_after").c_str()) ) && !refill) {

		before =(TH1F *) file.Get((directory+"/"+basename+"/"+basename+"_before").c_str());
		after  =(TH1F *) file.Get((directory+"/"+basename+"/"+basename+"_after").c_str());

	}	 
	else {

		cout<<basename.c_str()<<" Filling ..."<< endl;
		vars->ReadAnalysisBranches(tree);
		for(int i=0;i<tree->GetEntries()/FRAC;i++){

			vars->AnalysisVariablseReset();		
			UpdateProgressBar(i, tree->GetEntries()/FRAC);
			tree->GetEvent(i);
			if(weight) vars->mcweight=1;
			if(after->GetNbinsY()==1) FillEventByEvent(discr_var(vars),ApplyCuts(cut_before,vars),ApplyCuts(cut_after,vars),vars->mcweight);
			else FillEventByEventLatitude(discr_var(vars),ApplyCuts(cut_before,vars),ApplyCuts(cut_after,vars),GetLatitude(vars));
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

void Efficiency::FillEventByEventLatitude(float discr_var, bool CUT_BEFORE, bool CUT_AFTER,int latzone){

	int kbin;
        kbin =  bins.GetBin(discr_var);
		if(kbin>0){
			if(CUT_BEFORE) before->Fill(kbin,latzone);
			if(CUT_AFTER ) after->Fill(kbin,latzone);
		}
        return;
}


void Efficiency::Save(FileSaver finalhistos){
	finalhistos.Add(before);
	finalhistos.Add(after); 	
        finalhistos.writeObjsInFolder((directory+"/"+basename).c_str());
}


void Efficiency::Eval_Efficiency(){
	Eff = (TH1 *)after->Clone();
	Eff -> Sumw2();
	Eff -> Divide(before);
	Eff ->SetName((basename+"_Eff").c_str());
	Eff ->SetTitle((basename+" Efficiency").c_str());
	return;
}

void Efficiency::Eval_FittedEfficiency(){
	fitrequested=true;
	FittedEff= new TGraphErrors();
	FittedEff->SetName((basename+"_FitEff").c_str());
	FittedEff->SetName((basename+"_FitEff").c_str());

	for(int i=0;i<bins.size();i++){
		cout<<bins.EkPerMassBinCent(i)<<" "<<Eff->GetBinContent(i+1)<<endl;
		FittedEff->SetPoint(i,bins.EkPerMassBinCent(i), Eff->GetBinContent(i+1));
		FittedEff->SetPointError(i,0, Eff->GetBinError(i+1));	
	}

	FittedEff->Fit("pol3");	
	return;

}

void Efficiency::SaveResults(FileSaver finalhistos){
	finalhistos.Add(Eff); 	
	finalhistos.writeObjsInFolder((directory+"/"+basename).c_str());
	if(fitrequested){
		cout<<FittedEff<<endl;
		finalhistos.Add(FittedEff);
		finalhistos.writeObjsInFolder((directory+"/"+basename).c_str());
	}
}



