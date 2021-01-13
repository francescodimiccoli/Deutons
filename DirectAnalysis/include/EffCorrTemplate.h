#include "TemplateFITbetasmear.h"
#include "BetaSmearing.h"


class EffCorrTemplate: public EffCorr{

	private:

	Binning bins;
	std::vector<TH1F*> BeforeMC_P;
	std::vector<TH1F*> AfterMC_P;

	std::vector<TH1F*> BeforeMC_D;
	std::vector<TH1F*> AfterMC_D;

	std::vector<std::vector<TH1F*>> BeforeData;
	std::vector<std::vector<TH1F*>> AfterData;
	
	std::vector<TH1F*> BeforeDataGlob;
	std::vector<TH1F*> AfterDataGlob;
	
	std::vector<std::vector<TFit*>> FitBefore;
	std::vector<std::vector<TFit*>> FitAfter;
		
	TH2F * LatEfficiency_P;
	TH2F * LatEfficiency_D;

	TH1F * GlobEfficiency_P;
	TH1F * GlobEfficiency_D;

	TH1F * MCEfficiency_P;
	TH1F * MCEfficiency_D;
	
	TH1F * LatCorrections_P[10];
	TH1F * GlobalCorrection_P;

	TH1F * LatCorrections_D[10];
	TH1F * GlobalCorrection_D;
	
	std::string basename;
	std::string directory;
	std::string cut_before;
	std::string cut_after;
	std::string cut_data;
	std::string cut_MC;

	FileSaver  file;
	
	bool isrich=false;

	BadEventSimulator * BadEvSim=0x0;

	bool Refill = false;
	public:

	EffCorrTemplate(FileSaver  File, std::string Basename,std::string Directory, bool ekin, std::string Cut_before,std::string Cut_after,std::string Cut_Data,std::string Cut_MC,std::string Cut_MC2,std::string Cut_nopid,bool IsRICH=false):
	EffCorr(File,Basename,Directory,ekin,Cut_before,Cut_after,""){
		for(int bin=0;bin<Bins.size();bin++){
			
			std::string name = (Basename+"_MC_Before_"+to_string(bin));
			TH1F * BeforeP = new TH1F(name.c_str(),name.c_str(),30,0,3);	
			BeforeMC_P.push_back(BeforeP);

			name = (Basename+"_MC_After_"+to_string(bin));
                        TH1F * AfterP = new TH1F(name.c_str(),name.c_str(),30,0,3);
			AfterMC_P.push_back(AfterP);

			name = (Basename+"_MCD_Before_"+to_string(bin));
                        TH1F * BeforeD = new TH1F(name.c_str(),name.c_str(),30,0,3);
			BeforeMC_D.push_back(BeforeD);

			name = (Basename+"_MCD_After_"+to_string(bin));
                        TH1F * AfterD = new TH1F(name.c_str(),name.c_str(),30,0,3);
			AfterMC_D.push_back(AfterD);

		}
		for(int lat=0;lat<10;lat++) {
			BeforeData.push_back(std::vector<TH1F*>());
			AfterData.push_back(std::vector<TH1F*>());
			for(int bin=0;bin<Bins.size();bin++){

				string name = (Basename+"_DT_"+to_string(lat)+"_Before_"+to_string(bin));
				TH1F * Before= new TH1F(name.c_str(),name.c_str(),30,0,3);
				BeforeData[lat].push_back(Before);

				name = (Basename+"_DT_"+to_string(lat)+"_After_"+to_string(bin));	
				TH1F * After = new TH1F(name.c_str(),name.c_str(),30,0,3);
				AfterData[lat] .push_back(After);
			}
		}	
		for(int bin=0;bin<Bins.size();bin++){

			string name = (Basename+"_DT_Glob_Before_"+to_string(bin));
			TH1F * BeforeG= new TH1F(name.c_str(),name.c_str(),30,0,3);
			BeforeDataGlob.push_back(BeforeG);

			name = (Basename+"_DT_Glob_After_"+to_string(bin));
			TH1F * AfterG= new TH1F(name.c_str(),name.c_str(),30,0,3);
			AfterDataGlob.push_back(AfterG);

		}	
		if(File.GetFile()){
			for(int lat=0;lat<10;lat++) LatCorrections_P[lat]=(TH1F*) File.Get((Directory+"/"+Basename+"/" + Basename + "_CorrP_lat"+to_string(lat)).c_str());
			GlobalCorrection_P = (TH1F*) File.Get((Directory+"/"+Basename+"/" + Basename + "_CorrP_glob").c_str());
			for(int lat=0;lat<10;lat++) LatCorrections_D[lat]=(TH1F*) File.Get((Directory+"/"+Basename+"/" + Basename + "_CorrD_lat"+to_string(lat)).c_str());
			GlobalCorrection_D = (TH1F*) File.Get((Directory+"/"+Basename+"/" + Basename + "_CorrD_glob").c_str());
		}

		LatEfficiency_P= new TH2F((Basename + "_latP_Eff").c_str(),(Basename + "_latP_Eff").c_str(),Bins.size(),0,Bins.size(),10,0,10);
		LatEfficiency_D= new TH2F((Basename + "_latD_Eff").c_str(),(Basename + "_latD_Eff").c_str(),Bins.size(),0,Bins.size(),10,0,10);
	
		MCEfficiency_P = new TH1F((Basename + "_MCP_Eff").c_str(),(Basename + "_MCP_Eff").c_str(),Bins.size(),0,Bins.size());
		MCEfficiency_D = new TH1F((Basename + "_MCD_Eff").c_str(),(Basename + "_MCD_Eff").c_str(),Bins.size(),0,Bins.size());
		
		GlobEfficiency_P = new TH1F((Basename + "_GlobP_Eff").c_str(),(Basename + "_GlobP_Eff").c_str(),Bins.size(),0,Bins.size());
		GlobEfficiency_D = new TH1F((Basename + "_GlobD_Eff").c_str(),(Basename + "_GlobD_Eff").c_str(),Bins.size(),0,Bins.size());
		

		cut_before= Cut_before;
		cut_after= Cut_after;	
		cut_data = Cut_Data;
		cut_MC = Cut_MC;	
		basename = Basename;	
		directory = Directory;
		isrich=IsRICH;
		bins = Bins;
		file = File;
	}	
	bool ReinitializeHistos(bool refill);
	void Fill(TTree * treeMC,TTree * treeDT, Variables * vars,float (*var) (Variables * vars),float (*discr_var) (Variables * vars),bool refill=false);
	void FillEventByEventData(Variables * vars, float (*var) (Variables * vars), float (*discr_var) (Variables * vars));
	void FillEventByEventMC  (Variables * vars, float (*var) (Variables * vars), float (*discr_var) (Variables * vars));
	void SetUpBadEventSimulator(BadEventSimulator * Sim) {BadEvSim = Sim; return; };
	void LoadEventIntoBadEvSim(Variables * vars) {if(BadEvSim) BadEvSim->LoadEvent(vars);}
	void Save(FileSaver finalhistos);
	void Eval_Efficiencies();
	void SaveResults(FileSaver finalhistos);
	void Eval_Corrections();

	std::string GetName()		{return basename;}
	std::string GetDirPath()           {return directory;}
	TH1F * GetCorrectionLat_P(int lat)  {return LatCorrections_P[lat];}
	TH1F * GetGlobCorrection_P()	  {return GlobalCorrection_P;}
	TH1F * GetCorrectionLat_D(int lat)  {return LatCorrections_D[lat];}
	TH1F * GetGlobCorrection_D()	  {return GlobalCorrection_D;}

	

	TH1F * GetBeforeMC_P(int bin)	{return BeforeMC_P[bin];}
	TH1F * GetBeforeMC_D(int bin)	{return BeforeMC_D[bin];}
	TH1F * GetAfterMC_P(int bin)	{return AfterMC_P[bin];}
	TH1F * GetAfterMC_D(int bin)	{return AfterMC_D[bin];}
	TH1F * GetBeforeDataLat(int bin,int lat){return BeforeData[lat][bin];}
	TH1F * GetAfterDataLat(int bin,int lat){return AfterData[lat][bin];}
	TH1F * GetBeforeDataGlob(int bin,int lat){return BeforeDataGlob[bin];}
	TH1F * GetAfterDataGlob(int bin,int lat){return AfterDataGlob[bin];}
	
};

bool EffCorrTemplate::ReinitializeHistos(bool refill){
	TFile * input = file.GetFile();
	bool checkifsomeismissing=false;
	bool allfound=true;
	for(int bin=0;bin<bins.size();bin++){
			
			std::string name = (basename+"_MC_Before_"+to_string(bin));
			if(input) BeforeMC_P[bin] = (TH1F*)input->Get(((directory+"/"+basename+"_MC")+"/"+name).c_str());
			if((!input)||(!BeforeMC_P[bin])||refill) { BeforeMC_P[bin] = new TH1F(name.c_str(),name.c_str(),30,0,3); checkifsomeismissing=true;}	

			name = (basename+"_MC_After_"+to_string(bin));
			if(input) AfterMC_P[bin] = (TH1F*)input->Get(((directory+"/"+basename+"_MC")+"/"+name).c_str()); 
                        if((!input)||(!AfterMC_P[bin])||refill) { AfterMC_P[bin] = new TH1F(name.c_str(),name.c_str(),30,0,3); checkifsomeismissing=true;}

			name = (basename+"_MCD_Before_"+to_string(bin));
			if(input) BeforeMC_D[bin] = (TH1F*)input->Get(((directory+"/"+basename+"_MC")+"/"+name).c_str());
                        if((!input)||(!BeforeMC_D[bin])||refill){ BeforeMC_D[bin] = new TH1F(name.c_str(),name.c_str(),30,0,3); checkifsomeismissing=true;}

			name = (basename+"_MCD_After_"+to_string(bin));
			if(input) AfterMC_D[bin] = (TH1F*)input->Get(((directory+"/"+basename+"_MC")+"/"+name).c_str());
                        if((!input)||(!AfterMC_D[bin])||refill) { AfterMC_D[bin] = new TH1F(name.c_str(),name.c_str(),30,0,3); checkifsomeismissing=true;}
		
			name = (basename+"_DT_Glob_Before_"+to_string(bin));
			if(input) BeforeDataGlob[bin] = (TH1F*)input->Get((directory+"/"+basename+"_Glob"+"/"+name).c_str());
                        if((!input)||(!BeforeDataGlob[bin])||refill) { BeforeDataGlob[bin] = new TH1F(name.c_str(),name.c_str(),30,0,3); checkifsomeismissing=true;}
		
			name = (basename+"_DT_Glob_After_"+to_string(bin));
			if(input) AfterDataGlob[bin] = (TH1F*)input->Get((directory+"/"+basename+"_Glob"+"/"+name).c_str());
                        if((!input)||(!AfterDataGlob[bin])||refill) { AfterDataGlob[bin] = new TH1F(name.c_str(),name.c_str(),30,0,3); checkifsomeismissing=true;}
					
			}
		for(int lat=0;lat<10;lat++) {
			for(int bin=0;bin<bins.size();bin++){
			
				string name = (basename+"_DT_"+to_string(lat)+"_Before_"+to_string(bin));
				if(input) BeforeData[lat][bin]=(TH1F*)input->Get((directory+"/"+basename+"_lat/lat"+to_string(lat)+"/"+name).c_str());
				if((!input)||(!BeforeData[lat][bin])||refill) { BeforeData[lat][bin] = new TH1F(name.c_str(),name.c_str(),30,0,3);checkifsomeismissing=true;}
				
				name = (basename+"_DT_"+to_string(lat)+"_After_"+to_string(bin));	
				if(input) AfterData[lat][bin] = (TH1F*)input->Get((directory+"/"+basename+"_lat/lat"+to_string(lat)+"/"+name).c_str());
				if((!input)||(!AfterData[lat][bin])||refill) { AfterData[lat][bin] = new TH1F(name.c_str(),name.c_str(),30,0,3); checkifsomeismissing=true;}
			}
		}
	if(checkifsomeismissing||refill) { Refill=true; allfound=false;}
	return allfound;	
}



void EffCorrTemplate::Fill(TTree * treeMC,TTree * treeDT, Variables * vars,float (*var) (Variables * vars),float (*discr_var) (Variables * vars),bool refill){
	if(ReinitializeHistos(refill)){
		cout<<basename.c_str()<<" Filling ... (Data)"<< endl;
		vars->ReadBranches(treeDT);

		for(int i=0;i<treeDT->GetEntries()/FRAC;i++){
			UpdateProgressBar(i, treeDT->GetEntries()/FRAC);
			treeDT->GetEvent(i);
			vars->Update();
			FillEventByEventData(vars,var,discr_var);
		}
		
		cout<<basename.c_str()<<" Filling ... (MC)"<< endl;
		vars->ReadBranches(treeMC);

		for(int i=0;i<treeMC->GetEntries();i++){
			if(i%(int)FRAC!=0) continue;
			UpdateProgressBar(i, treeMC->GetEntries());
			treeMC->GetEvent(i);
			vars->Update();
			if(BadEvSim) BadEvSim->LoadEvent(vars);
			FillEventByEventMC(vars,var,discr_var);
		}
	}

	return;
}

void EffCorrTemplate::FillEventByEventData(Variables * vars, float (*var) (Variables * vars), float (*discr_var) (Variables * vars)){
	if(!Refill) return;
	int kbin;
	kbin = 	bins.GetBin(discr_var(vars));
	if(ApplyCuts((cut_before+"&"+cut_data).c_str(),vars)&&kbin>0)		{BeforeData[GetLatitude(vars)][kbin]->Fill(var(vars),vars->PrescaleFactor);  BeforeDataGlob[kbin]->Fill(var(vars),vars->PrescaleFactor);}
	if(ApplyCuts((cut_after +"&"+cut_data).c_str(),vars)&&kbin>0)		{AfterData[GetLatitude(vars)][kbin]->Fill(var(vars),vars->PrescaleFactor);   AfterDataGlob[kbin]->Fill(var(vars),vars->PrescaleFactor);}
	return;	

}

void EffCorrTemplate::FillEventByEventMC(Variables * vars, float (*var) (Variables * vars), float (*discr_var) (Variables * vars)){
	if(!Refill) return;
	float betasmear=discr_var(vars);

	if(BadEvSim)
		betasmear=BadEvSim->SimulateBadEvents(betasmear);

	int kbin =  bins.GetBin(betasmear);
	float mass = vars->R/betasmear * pow((1-pow(betasmear,2)),0.5);
	if(ApplyCuts((cut_before+"&"+cut_MC +"&IsPurePMC").c_str(),vars)&&kbin>0)  BeforeMC_P[kbin] ->Fill(mass,vars->mcweight);		
	if(ApplyCuts((cut_after +"&"+cut_MC +"&IsPurePMC").c_str(),vars)&&kbin>0)  AfterMC_P[kbin]  ->Fill(mass,vars->mcweight);
	if(ApplyCuts((cut_before+"&"+cut_MC +"&IsPureDMC").c_str(),vars)&&kbin>0)  BeforeMC_D[kbin] ->Fill(mass,vars->mcweight);		
	if(ApplyCuts((cut_after +"&"+cut_MC +"&IsPureDMC").c_str(),vars)&&kbin>0)  AfterMC_D[kbin]  ->Fill(mass,vars->mcweight);

	return;	
}


void EffCorrTemplate::Save(FileSaver finalhistos){

	for(int bin=0;bin<bins.size();bin++){
		finalhistos.Add(BeforeMC_P[bin]);
		finalhistos.Add(AfterMC_P[bin]);
		finalhistos.Add(BeforeMC_D[bin]);
		finalhistos.Add(AfterMC_D[bin]);
	}	
	finalhistos.writeObjsInFolder((directory+"/"+basename+"_MC").c_str());

	for(int lat=0;lat<10;lat++){
		for(int bin=0;bin<bins.size();bin++){
			finalhistos.Add(BeforeData[lat][bin]);
			finalhistos.Add(AfterData[lat][bin]);
		}	
		finalhistos.writeObjsInFolder((directory+"/"+basename+"_lat/lat"+to_string(lat)).c_str());
	}
	for(int bin=0;bin<bins.size();bin++){
		finalhistos.Add(BeforeDataGlob[bin]);
                finalhistos.Add(AfterDataGlob[bin]);
	}
	finalhistos.writeObjsInFolder((directory+"/"+basename+"_Glob").c_str());

	return;
}

TH1F * EvalEfficiency(TH1F * MCEfficiency,std::vector<TH1F*> Before,std::vector<TH1F*> After){

	TH1F * MCBefore=(TH1F*)MCEfficiency->Clone();
	TH1F * MCAfter =(TH1F*)MCEfficiency->Clone();
	for(int bins=0;bins<MCEfficiency->GetNbinsX();bins++){
		MCBefore->SetBinContent(bins+1,Before[bins]->Integral());
		MCBefore->SetBinError(bins+1,pow(Before[bins]->Integral(),0.2));	
		MCAfter ->SetBinContent(bins+1,After[bins]->Integral());
		MCAfter ->SetBinError(bins+1,pow(After[bins]->Integral(),0.2));	
	}	
	MCAfter->Sumw2();
	MCBefore->Sumw2();
	MCAfter->Divide(MCBefore);
	MCEfficiency=(TH1F*)MCAfter->Clone();
	return MCEfficiency;
}

float GetStatError(TFit * FitAfter,TFit * FitBefore){
	float stat_error=0;
	return pow(stat_error,0.5);	

}

void EffCorrTemplate::Eval_Efficiencies(){


	MCEfficiency_P=EvalEfficiency(MCEfficiency_P,BeforeMC_P,AfterMC_P);
	MCEfficiency_D=EvalEfficiency(MCEfficiency_D,BeforeMC_D,AfterMC_D);

	for(int lat=0;lat<10;lat++){
		FitBefore.push_back(std::vector<TFit*>());
                FitAfter.push_back(std::vector<TFit*>());	
		for(int bin=0;bin<bins.size();bin++){
			TFit * fit_before = new TFit(BeforeMC_P[bin],BeforeMC_D[bin],BeforeData[lat][bin],BeforeData[lat][bin]);		
			TFit * fit_after  = new TFit(AfterMC_P[bin],AfterMC_D[bin],AfterData[lat][bin],AfterData[lat][bin]);		
			FitBefore[lat].push_back(fit_before);	
			FitAfter[lat].push_back(fit_after);
		}
	}
	FitBefore.push_back(std::vector<TFit*>());
        FitAfter.push_back(std::vector<TFit*>());	
	for(int bin=0;bin<bins.size();bin++){
			TFit * fit_before = new TFit(BeforeMC_P[bin],BeforeMC_D[bin],BeforeDataGlob[bin],BeforeDataGlob[bin]);		
			TFit * fit_after  = new TFit(AfterMC_P[bin],AfterMC_D[bin],AfterDataGlob[bin],AfterDataGlob[bin]);		
			FitBefore[10].push_back(fit_before);	
			FitAfter[10].push_back(fit_after);
	}	

	float constrainmin[3] = {0.0001,0.0001,0.0001};
	float constrainmax[3] = {1,1,1};
	for(int bin=0;bin<bins.size();bin++){
		Do_TemplateFIT(FitBefore[10][bin],0.1,3,constrainmin,constrainmax,false,false,false);
		Do_TemplateFIT(FitAfter [10][bin],0.1,3,constrainmin,constrainmax,false,false,false);
	
		if(FitBefore[10][bin]->Tfit_outcome==0&&FitAfter[10][bin]->Tfit_outcome==0){ 
			GlobEfficiency_P->SetBinContent(bin+1,FitAfter[10][bin]->Templ_P->Integral()/FitBefore[10][bin]->Templ_P->Integral());	
			GlobEfficiency_D->SetBinContent(bin+1,FitAfter[10][bin]->Templ_D->Integral()/FitBefore[10][bin]->Templ_D->Integral());
			float stat_error=GetStatError(FitAfter[10][bin],FitBefore[10][bin]);
			GlobEfficiency_P->SetBinError(bin+1,stat_error);	
			GlobEfficiency_D->SetBinError(bin+1,stat_error);
		}	
	}
	for(int lat=0;lat<10;lat++)
		for(int bin=0;bin<bins.size();bin++){
		//Do_TemplateFIT(FitBefore[lat][bin],0.1,3);
                //Do_TemplateFIT(FitAfter [lat][bin],0.1,3);

		if(FitBefore[lat][bin]->Tfit_outcome==0&&FitAfter[lat][bin]->Tfit_outcome==0){
                        LatEfficiency_P->SetBinContent(bin+1,FitAfter[lat][bin]->Templ_P->Integral()/FitBefore[lat][bin]->Templ_P->Integral());
                        LatEfficiency_D->SetBinContent(bin+1,FitAfter[lat][bin]->Templ_D->Integral()/FitBefore[lat][bin]->Templ_D->Integral());
                        float stat_error=GetStatError(FitAfter[lat][bin],FitBefore[lat][bin]);
                        LatEfficiency_P->SetBinError(bin+1,stat_error);
                        LatEfficiency_D->SetBinError(bin+1,stat_error);
                }
		else{
                        cout<<"Integral: "<<lat<<" "<<bin<<" "<<FitAfter[lat][bin]->Data->Integral()<<" "<<FitBefore[lat][bin]->Data->Integral()<<endl;
			LatEfficiency_P->SetBinContent(bin+1,lat+1,FitAfter[lat][bin]->Data->Integral()/FitBefore[lat][bin]->Data->Integral());
                        LatEfficiency_D->SetBinContent(bin+1,lat+1,FitAfter[lat][bin]->Data->Integral()/FitBefore[lat][bin]->Data->Integral());
                        float stat_error=pow(FitAfter[lat][bin]->Data->Integral(),0.5)/FitAfter[lat][bin]->Data->Integral() + 
					 pow(FitBefore[lat][bin]->Data->Integral(),0.5)/FitBefore[lat][bin]->Data->Integral() ;
                        LatEfficiency_P->SetBinError(bin+1,stat_error);
                        LatEfficiency_D->SetBinError(bin+1,stat_error);
                }
	}	
	return;	
}


void EffCorrTemplate::SaveResults(FileSaver finalhistos){

	finalhistos.Add(MCEfficiency_P);
	finalhistos.Add(MCEfficiency_D);
	finalhistos.Add(GlobEfficiency_P);
	finalhistos.Add(GlobEfficiency_D);
	
	finalhistos.writeObjsInFolder((directory+"/"+basename).c_str());

	for(int lat=0;lat<11;lat++){
		for(int bin=0;bin<bins.size();bin++){
			finalhistos.Add(FitBefore[lat][bin]->Templ_P);
			finalhistos.Add(FitBefore[lat][bin]->Templ_D);
			finalhistos.Add(FitBefore[lat][bin]->Data);
		}
		if(lat!=10) finalhistos.writeObjsInFolder((directory+"/"+basename+"/TFit/Before/Lat"+to_string(lat)).c_str());
		else finalhistos.writeObjsInFolder((directory+"/"+basename+"/TFit/Before/Glob").c_str());	 
		for(int bin=0;bin<bins.size();bin++){
			finalhistos.Add(FitAfter[lat][bin]->Templ_P);
			finalhistos.Add(FitAfter[lat][bin]->Templ_D);
			finalhistos.Add(FitAfter[lat][bin]->Data);
		}
		if(lat!=10)finalhistos.writeObjsInFolder((directory+"/"+basename+"/TFit/After/Lat"+to_string(lat)).c_str());
		else finalhistos.writeObjsInFolder((directory+"/"+basename+"/TFit/After/Glob").c_str());
	}

	

	for(int lat=0;lat<10;lat++) 
		if(LatCorrections_P[lat]){
			finalhistos.Add(LatCorrections_P[lat]); 	
			finalhistos.writeObjsInFolder((directory+"/"+basename).c_str());
		}
	if(GlobalCorrection_P) {finalhistos.Add(GlobalCorrection_P);
		finalhistos.writeObjsInFolder((directory+"/"+basename).c_str());	
	}
	for(int lat=0;lat<10;lat++) 
		if(LatCorrections_D[lat]){
			finalhistos.Add(LatCorrections_D[lat]); 	
			finalhistos.writeObjsInFolder((directory+"/"+basename).c_str());
		}
	if(GlobalCorrection_D) {finalhistos.Add(GlobalCorrection_D);
		finalhistos.writeObjsInFolder((directory+"/"+basename).c_str());	
	}
	

	return;
}


void EffCorrTemplate::Eval_Corrections(){

	for(int lat=0;lat<10;lat++) {
		LatCorrections_P[lat] = ProjectionXtoTH1F(LatEfficiency_P,(basename + "_CorrP_lat"+to_string(lat)).c_str(),lat+1,lat+1);
		LatCorrections_P[lat]->SetName((basename + "_CorrP_lat"+to_string(lat)).c_str());
		LatCorrections_P[lat]->Sumw2();
		LatCorrections_P[lat]->Divide(MCEfficiency_P);
	}


	GlobalCorrection_P = (TH1F *) GlobEfficiency_P->Clone();
	GlobalCorrection_P-> SetName((basename + "_CorrP_glob").c_str());
	GlobalCorrection_P->Sumw2();
	GlobalCorrection_P->Divide(MCEfficiency_P);	

	for(int lat=0;lat<10;lat++) {
		LatCorrections_D[lat] = ProjectionXtoTH1F(LatEfficiency_D,(basename + "_CorrD_lat"+to_string(lat)).c_str(),lat+1,lat+1);
		LatCorrections_D[lat]->SetName((basename + "_CorrD_lat"+to_string(lat)).c_str());
		LatCorrections_D[lat]->Sumw2();
		LatCorrections_D[lat]->Divide(MCEfficiency_D);
	}


	GlobalCorrection_D = (TH1F *) GlobEfficiency_D->Clone();
	GlobalCorrection_D-> SetName((basename + "_CorrD_glob").c_str());
	GlobalCorrection_D->Sumw2();
	GlobalCorrection_D->Divide(MCEfficiency_D);	
	
	return;
}



