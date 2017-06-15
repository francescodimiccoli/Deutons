void UpdateZoneLivetime (float Livetime, float Rcutoff, TH1F * esposizionegeo,Binning bins){

        for(int i=0;i<esposizionegeo->GetNbinsX();i++)
                        if(bins.RigBins()[i]>=1.2*Rcutoff){
                                esposizionegeo -> SetBinContent(i+1, esposizionegeo -> GetBinContent(i+1) + Livetime) ;
        }
        return;
}

struct MCPar{
float Rmin,Rmax,Trigrate;
};

class Flux{

	private:
	Efficiency * MCEfficiency;
	TH1F * Counts;
	TH1F * Geom_Acceptance;
	TH1F * ExposureTime;
	Binning bins;
	std::string basename;	
	std::string exposurename;
	std::string geomname;

	TH1F * FluxEstim;
	MCPar param;

	public:


	Flux(FileSaver File, FileSaver FileRes, std::string Basename, std::string Effname, std::string EffDir,std::string CountsName,std::string ExposureName, std::string GeomName,Binning Bins){
		MCEfficiency = new Efficiency(FileRes,Effname,EffDir,Bins);
		Counts = (TH1F *) FileRes.Get((CountsName).c_str());
		ExposureTime = (TH1F *) File.Get(("Fluxes/"+Basename+"/"+ExposureName).c_str());
		Geom_Acceptance = (TH1F *) File.Get(("Fluxes/"+Basename+"/"+GeomName).c_str());

		bins = Bins;		
		basename = Basename;
		exposurename = ExposureName;
		geomname = GeomName;
	}
	Flux(FileSaver File, FileSaver FileRes, std::string Basename, std::string Effname, std::string EffDir,std::string CountsName, Binning Bins){
		MCEfficiency = new Efficiency(FileRes,Effname,EffDir,Bins);
		Counts = (TH1F *) FileRes.Get((CountsName).c_str());

		bins = Bins;		
		basename = Basename;
	}

	void Set_MCPar(float rmin, float rmax, float trigrate);
	void Eval_ExposureTime(TNtuple * treeDT, FileSaver finalhistos);
	void Eval_GeomAcceptance(TNtuple * treeMC, FileSaver finalhistos,std::string cut);
	void Eval_Flux();
	void SaveResults(FileSaver finalhistos);

};

void Flux::Set_MCPar(float rmin, float rmax, float trigrate){
	param.Rmin=rmin;
	param.Rmax=rmax;
	param.Trigrate=trigrate;
}


void Flux::Eval_Flux(){

	if(Counts) FluxEstim = (TH1F *) Counts->Clone();
	else FluxEstim = new TH1F("dummy","dummy",bins.size(),0,bins.size());

	FluxEstim -> SetName((basename+"_Flux").c_str());
	FluxEstim -> SetTitle((basename+"_Flux").c_str());

	FluxEstim -> Divide(MCEfficiency->GetEfficiency());
	if(ExposureTime) FluxEstim -> Divide(ExposureTime);
	if(Geom_Acceptance) FluxEstim -> Divide(Geom_Acceptance);

	for(int i=0;i<bins.size();i++)
		FluxEstim->SetBinContent(i+1,FluxEstim->GetBinContent(i+1)/(bins.EkPerMasBins()[i+1]-bins.EkPerMasBins()[i]));

	return;		
}

void Flux::SaveResults(FileSaver finalhistos){
	finalhistos.Add(FluxEstim); 	
	finalhistos.Add(ExposureTime);
	finalhistos.Add(Geom_Acceptance);
	finalhistos.writeObjsInFolder(("Fluxes/"+basename).c_str());
}

void Flux::Eval_GeomAcceptance(TNtuple * RawMC,FileSaver finalhistos,std::string cut){
	if(1!=0){
		Geom_Acceptance = new TH1F(geomname.c_str(),geomname.c_str(),bins.size(),0,bins.size());

		Variables * vars = new Variables;
		Efficiency * Geom = new Efficiency(finalhistos,"","",bins,cut.c_str() ,cut.c_str());
		Geom->Fill(RawMC, vars,GetBetaGen,true);
	
		Geom_Acceptance = (TH1F *)Geom->GetBefore()->Clone();
		Geom_Acceptance -> SetName(geomname.c_str());
		Geom_Acceptance -> SetTitle(geomname.c_str());
	
		float gen_events = RawMC->GetEntries();
		
		for(int i=0;i<bins.size();i++){
			float gen_bins= bins.RigBins()[i+1]-bins.RigBins()[i];///*gen_events*(pow(param.Trigrate,-1))*/(log(bins.RigBins()[i+1])-log(bins.RigBins()[i]));//(log(param.Rmax)-log(param.Rmin)); 
			Geom_Acceptance -> SetBinContent(i+1,/*Geom_Acceptance -> GetBinContent(i+1)*/gen_bins); 			
			}

		finalhistos.Add(Geom_Acceptance); 	
		finalhistos.writeObjsInFolder(("Fluxes/"+basename).c_str());


	}
	return;
};

void Flux::Eval_ExposureTime(TNtuple * RawDT,FileSaver finalhistos){
		
	if(!ExposureTime){
	
		float Livetime,U_time,Rcutoff;
		RawDT->SetBranchAddress("Livetime"           ,&Livetime);
		RawDT->SetBranchAddress("U_Time"             ,&U_time);
		RawDT->SetBranchAddress("Rcutoff"            ,&Rcutoff);


		ExposureTime = new TH1F(exposurename.c_str(),exposurename.c_str(),bins.size(),0,bins.size());

		int ActualTime=0;
		for(int i=0;i<RawDT->GetEntries();i++){
			UpdateProgressBar(i, RawDT->GetEntries());
			RawDT->GetEvent(i);
			if((int)U_time!=ActualTime) {
				UpdateZoneLivetime(Livetime,Rcutoff,ExposureTime,bins);
				ActualTime=U_time;
			}
		}
		finalhistos.Add(ExposureTime); 	
		finalhistos.writeObjsInFolder(("Fluxes/"+basename).c_str());

	}
	return;
}

