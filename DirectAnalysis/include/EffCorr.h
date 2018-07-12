class EffCorr{

	private:

	Efficiency * EffMC;
	Efficiency * EffMC2;
	Efficiency * EffMCnopid;

	Efficiency * EffData;
	
	TH1F * ExposureZones;

	TH1F * LatCorrections[10];
	TH1F * GlobalCorrection;
	TH1F * GlobalCorrection2;
	TH1F * GlobalCorrectionnopid;



	std::string basename;
	std::string directory;

	BadEventSimulator * BadEvSim;	
	public:

	EffCorr(FileSaver  File, std::string Basename,std::string Directory, Binning Bins, std::string Cut_before,std::string Cut_after,std::string Cut_Data,std::string Cut_MC,std::string Cut_MC2,std::string Cut_MCnopid){
		cout<<"cia"<<endl;
	
		EffMC   = new Efficiency(File, (Basename+"_MC" ).c_str(),Directory,Bins, (Cut_before+"&"+Cut_MC  ).c_str(),(Cut_after+"&"+Cut_MC  ).c_str());
		EffMC2 = new Efficiency(File, (Basename+"_MC2").c_str(),Directory,Bins, (Cut_before+"&"+Cut_MC2  ).c_str(),(Cut_after+"&"+Cut_MC2  ).c_str());
		EffMCnopid = new Efficiency(File, (Basename+"_MCnopid").c_str(),Directory,Bins, (Cut_before+"&"+Cut_MCnopid  ).c_str(),(Cut_after+"&"+Cut_MCnopid  ).c_str());

		EffData = new Efficiency(File, (Basename+"_lat").c_str(),Directory,Bins, (Cut_before+"&"+Cut_Data).c_str(),(Cut_after+"&"+Cut_Data).c_str(),LatEdges);
		TFile * file = File.GetFile();
		if(file){
			for(int lat=0;lat<10;lat++) LatCorrections[lat]=(TH1F*) file->Get((Directory+"/"+Basename+"/" + Basename + "_Corr_lat"+to_string(lat)).c_str());
			GlobalCorrection = (TH1F*) file->Get((Directory+"/"+Basename+"/" + Basename + "_Corr_glob").c_str());
			GlobalCorrection2 = (TH1F*) file->Get((Directory+"/"+Basename+"/" + Basename + "_Corr_glob2").c_str());
			GlobalCorrectionnopid = (TH1F*) file->Get((Directory+"/"+Basename+"/" + Basename + "_Corr_globnopid").c_str());
		}
		basename = Basename;	
		directory = Directory;
	}	

	virtual void SetUpBadEventSimulator(BadEventSimulator * Sim) {
		BadEvSim = Sim;
	};
	
	virtual void LoadEventIntoBadEvSim(Variables * vars) {
		EffMC->LoadEventIntoBadEvSim(vars);
		EffMC2->LoadEventIntoBadEvSim(vars);
		EffMCnopid->LoadEventIntoBadEvSim(vars);
	}
	virtual bool ReinitializeHistos(bool refill){
		bool checkifsomeismissing=false;
		bool allfound=true;
		if(!(EffMC -> ReinitializeHistos(refill))) checkifsomeismissing   = true;
         	if(!(EffMC2 -> ReinitializeHistos(refill))) checkifsomeismissing   = true;
		if(!(EffMCnopid -> ReinitializeHistos(refill))) checkifsomeismissing   = true;
	        if(!(EffData -> ReinitializeHistos(refill))) checkifsomeismissing = true;
	 	if(checkifsomeismissing||refill) allfound=false;
		return allfound;
	}
	virtual void FillEventByEventMC(Variables * vars, float (*var) (Variables * vars), float (*discr_var) (Variables * vars)){
		EffMC -> FillEventByEventMC(vars,var,discr_var);
		EffMC2 -> FillEventByEventMC(vars,var,discr_var);
		EffMCnopid -> FillEventByEventMC(vars,var,discr_var);
	}
	
	virtual void FillEventByEventData(Variables * vars, float (*var) (Variables * vars), float (*discr_var) (Variables * vars)){
		EffData -> FillEventByEventData(vars,var,discr_var);
	}
	virtual void Fill(TTree * treeMC,TTree * treeDT, Variables * vars, float (*discr_var) (Variables * vars),bool refill=false);
	virtual void Save(FileSaver finalhistos);
	virtual void Eval_Efficiencies();
	virtual void SaveResults(FileSaver finalhistos);
	virtual void Eval_Corrections();
	TH1F * GetMCEfficiency()	{return (TH1F*)EffMC  -> GetEfficiency();}
	TH1F * GetCorrectionLat(int lat)  {return LatCorrections[lat];}
	TH1F * GetGlobCorrection()	  {return GlobalCorrection;}
	TH1F * GetGlobCorrection2()	  {return GlobalCorrection2;}
	TH1F * GetGlobCorrection_noPID()	  {return GlobalCorrectionnopid;}



};


void EffCorr::Fill(TTree * treeMC,TTree * treeDT, Variables * vars, float (*discr_var) (Variables * vars),bool refill){

	EffMC   -> Fill(treeMC,vars,discr_var,refill);
	EffMC2   -> Fill(treeMC,vars,discr_var,refill);
	EffMCnopid   -> Fill(treeMC,vars,discr_var,refill);
	EffData -> Fill(treeDT,vars,discr_var,refill);
}

void EffCorr::Save(FileSaver finalhistos){
	cout<<"coccole"<<endl;
	EffMC  -> Save(finalhistos);
	EffMC2  -> Save(finalhistos);
	EffMCnopid  -> Save(finalhistos);
	EffData-> Save(finalhistos);
}

void EffCorr::Eval_Efficiencies(){
	EffMC  -> Eval_Efficiency();
	EffMC2  -> Eval_Efficiency();
	EffMCnopid  -> Eval_Efficiency();
	EffData-> Eval_Efficiency();
}

void EffCorr::SaveResults(FileSaver finalhistos){
	EffMC  -> SaveResults(finalhistos);
	EffMC2  -> SaveResults(finalhistos);
	EffMCnopid  -> SaveResults(finalhistos);
	EffData-> SaveResults(finalhistos);

	for(int lat=0;lat<10;lat++) 
		if(LatCorrections[lat]){
			finalhistos.Add(LatCorrections[lat]); 	
			finalhistos.writeObjsInFolder((directory+"/"+basename).c_str());
		}
	finalhistos.Add(GlobalCorrection);
 	finalhistos.Add(GlobalCorrection2);
 	finalhistos.Add(GlobalCorrectionnopid);

       finalhistos.writeObjsInFolder((directory+"/"+basename).c_str());	
}


void EffCorr::Eval_Corrections(){

	for(int lat=0;lat<10;lat++) {
		LatCorrections[lat] = ProjectionXtoTH1F((TH2F*)EffData->GetEfficiency(),(basename + "_Corr_lat"+to_string(lat)).c_str(),lat+1,lat+1);
		LatCorrections[lat]->SetName((basename + "_Corr_lat"+to_string(lat)).c_str());
		LatCorrections[lat]->Sumw2();
		LatCorrections[lat]->Divide(EffMC->GetEfficiency());
	}

	TH1F * Global_Before=ProjectionXtoTH1F((TH2F*)EffData->GetBefore(),(basename + "_Corr_glob").c_str(),0,10);	
	TH1F * Global_After =ProjectionXtoTH1F((TH2F*)EffData->GetAfter() ,(basename + "_Corr_glob").c_str(),0,10);	

	GlobalCorrection = (TH1F *) Global_After->Clone();
	GlobalCorrection -> SetName((basename + "_Corr_glob").c_str());
	GlobalCorrection->Sumw2();
	GlobalCorrection->Divide(Global_Before);	
	
	GlobalCorrection2=(TH1F *) EffMC->GetEfficiency()->Clone();
	GlobalCorrectionnopid=(TH1F *) EffMC->GetEfficiency()->Clone();
	GlobalCorrection2->SetName((basename + "_Corr_glob2").c_str());
	GlobalCorrection2->SetTitle((basename + "_Corr_glob2").c_str());
	GlobalCorrectionnopid->SetName((basename + "_Corr_globnopid").c_str());
	GlobalCorrectionnopid->SetTitle((basename + "_Corr_globnopid").c_str());

	GlobalCorrection->Divide(EffMC->GetEfficiency());	
	GlobalCorrection2->Divide(EffMC2->GetEfficiency());	
	GlobalCorrectionnopid->Divide(EffMCnopid->GetEfficiency());	
	
	return;
}
