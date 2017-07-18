class EffCorr{

	private:

	Efficiency * EffMC;
	Efficiency * EffData;
	
	TH1F * ExposureZones;

	TH1F * LatCorrections[10];
	TH1F * GlobalCorrection;
	std::string basename;
	std::string directory;
	
	public:

	EffCorr(FileSaver  File, std::string Basename,std::string Directory, Binning Bins, std::string Cut_before,std::string Cut_after,std::string Cut_Data,std::string Cut_MC){
		EffMC   = new Efficiency(File, (Basename+"_MC" ).c_str(),Directory,Bins, (Cut_before+"&"+Cut_MC  ).c_str(),(Cut_after+"&"+Cut_MC  ).c_str());
		EffData = new Efficiency(File, (Basename+"_lat").c_str(),Directory,Bins, (Cut_before+"&"+Cut_Data).c_str(),(Cut_after+"&"+Cut_Data).c_str(),LatEdges);

		for(int lat=0;lat<10;lat++) LatCorrections[lat]=(TH1F*) File.Get((Directory+"/"+Basename+"/" + Basename + "_Corr_lat"+to_string(lat)).c_str());
		GlobalCorrection = (TH1F*) File.Get((Directory+"/"+Basename+"/" + Basename + "_Corr_glob").c_str());

			for(int lat=0;lat<10;lat++) cout<<LatCorrections[lat]<<endl;
			cout<<GlobalCorrection<<endl; 
		
		basename = Basename;	
		directory = Directory;
	}	

	void Fill(TNtuple * treeMC,TNtuple * treeDT, Variables * vars, float (*discr_var) (Variables * vars),bool refill=false);
	void Save(FileSaver finalhistos);
	void Eval_Efficiencies();
	void SaveResults(FileSaver finalhistos);
	void Eval_Corrections();
	TH1F * GetCorrectionLat(int lat)  {return LatCorrections[lat];}
	TH1F * GetGlobCorrection()	  {return GlobalCorrection;}

};


void EffCorr::Fill(TNtuple * treeMC,TNtuple * treeDT, Variables * vars, float (*discr_var) (Variables * vars),bool refill){

	EffMC   -> Fill(treeMC,vars,discr_var,refill);
	EffData -> Fill(treeDT,vars,discr_var,refill);
}

void EffCorr::Save(FileSaver finalhistos){
	EffMC  -> Save(finalhistos);
	EffData-> Save(finalhistos);
}

void EffCorr::Eval_Efficiencies(){
	EffMC  -> Eval_Efficiency();
	EffData-> Eval_Efficiency();
}

void EffCorr::SaveResults(FileSaver finalhistos){
	EffMC  -> SaveResults(finalhistos);
	EffData-> SaveResults(finalhistos);

	for(int lat=0;lat<10;lat++) 
		if(LatCorrections[lat]){
			finalhistos.Add(LatCorrections[lat]); 	
			finalhistos.writeObjsInFolder((directory+"/"+basename).c_str());
		}
	finalhistos.Add(GlobalCorrection);
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
	GlobalCorrection->Divide(EffMC->GetEfficiency());	
	
	return;
}
