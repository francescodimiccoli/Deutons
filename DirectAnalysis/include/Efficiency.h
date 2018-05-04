#ifndef EFFICIENCY_H
#define EFFICIENCY_H

#include "BetaSmearing.h"
#include "Cuts.h"
#include "DBarReader.h"

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
	BadEventSimulator * BadEvSim=0x0;

	TH1F * Stat_Error=0x0;
	TH1F * Syst_Error=0x0;
	
	bool Refill = true;
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

	Efficiency(FileSaver  File, std::string Basename,std::string Directory, Binning Bins){
		file=File;
		bins=Bins;
		basename=Basename;
		directory=Directory;
	
		ReadFile();
	}

	void ReadFile(){
		before =(TH1 *) file.Get((directory+"/"+basename+"/"+basename+"_before").c_str());
                after  =(TH1 *) file.Get((directory+"/"+basename+"/"+basename+"_after").c_str());
                Eff    =(TH1 *) file.Get((directory+"/"+basename+"/"+basename+"_Eff").c_str());
        	Stat_Error = (TH1F *) file.Get((directory+"/"+basename+"/"+basename+"_Stat_Error").c_str());
		Syst_Error = (TH1F *) file.Get((directory+"/"+basename+"/"+basename+"_Syst_Error").c_str());
	}

	bool ReinitializeHistos(bool refill);
	void Fill( TTree * treeMC, Variables * vars, float (*discr_var) (Variables * vars),bool refill=false);
	void Fill(DBarReader readerMC, Variables * vars, float (*discr_var) (Variables * vars),bool refill);
	void FillEventByEventMC(Variables * vars, float (*var) (Variables * vars), float (*discr_var) (Variables * vars));
	void FillEventByEventData(Variables * vars, float (*var) (Variables * vars), float (*discr_var) (Variables * vars));

	void SetUpBadEventSimulator(BadEventSimulator * Sim) {BadEvSim = Sim; return; };
	BadEventSimulator * GetBadEventSimulator() {return BadEvSim;};
	void LoadEventIntoBadEvSim(Variables * vars) {if(BadEvSim) BadEvSim->LoadEvent(vars);}	
	
	void Save(FileSaver finalhistos);
	void SaveResults(FileSaver finalhistos);
	void Eval_Efficiency();
	void Eval_FittedEfficiency();

	void CloneEfficiency(Efficiency * Second);
	void ComposeEfficiency(Efficiency * Second);
	void Eval_StatError();
	void Eval_SystError(Efficiency * First,  Efficiency * Second);
	
	TH1 * GetEfficiency() {return Eff;}
	TGraphErrors * GetFittedEfficiency() {return FittedEff;}
	TH1 * GetBefore()     {return before;}
	TH1 * GetAfter()      {return after;}

	TH1F * GetSyst_Error() {return Syst_Error;}
	TH1F * GetStat_Error() {return Stat_Error;}

	std::string GetCut_Before()	{return cut_before;}	
	std::string GetCut_After()	{return cut_after;}
};

void Efficiency::Eval_StatError(){
	Stat_Error = (TH1F *) Eff->Clone();
	Stat_Error->SetName((basename+"_Stat_Error").c_str());
	Stat_Error->SetTitle((basename+"_Stat_Error").c_str());
	for(int i=0; i<Stat_Error->GetNbinsX();i++){
		if(Eff->GetBinContent(i+1)>0)
			Stat_Error->SetBinContent(i+1,Eff->GetBinError(i+1)/Eff->GetBinContent(i+1));
			Stat_Error->SetBinError(i+1,0);
	}
};

void Efficiency::Eval_SystError(Efficiency * First,  Efficiency * Second){
	Syst_Error = (TH1F *) Eff->Clone();
	Syst_Error->SetName((basename+"_Syst_Error").c_str());
	Syst_Error->SetTitle((basename+"_Syst_Error").c_str());
	TH1F * A = (TH1F *) First->GetEfficiency()->Clone();
	TH1F * B = (TH1F *) Second->GetEfficiency()->Clone();
	A->Smooth(2);
	B->Smooth(2);
        for(int i=0; i<Syst_Error->GetNbinsX();i++){
		if(Second->GetEfficiency()->GetBinContent(i+1)>0){
		
                	Syst_Error->SetBinContent(i+1,fabs(1 - A->GetBinContent(i+1)/B->GetBinContent(i+1))/pow(3,0.5));
                	Syst_Error->SetBinError(i+1,0);
			Eff->SetBinError(i+1, pow(pow(Stat_Error->GetBinContent(i+1),2) + pow(Syst_Error->GetBinContent(i+1),2),0.5)*Eff->GetBinContent(i+1) );
		}
	}
};


void Efficiency::ComposeEfficiency(Efficiency * Second){
	TH1F * after1= (TH1F*) after->Clone();
	TH1F * after2= (TH1F*) (TH1F*) Second->GetAfter()->Clone();

	if(after1->GetEntries()>after2->GetEntries()){
		after = (TH1F*) Second->GetAfter()->Clone();
	}
	else{
		before = (TH1F*) Second->GetBefore()->Clone();
	}

	if(Eff){
		Eff->Sumw2();
		if(Second->GetEfficiency()) Eff->Multiply(Second->GetEfficiency());
	}
	else {
		Eval_Efficiency();	
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

bool Efficiency::ReinitializeHistos(bool refill){
	bool allfound;
	if(( (TH1F*)  file.Get((directory+"/"+basename+"/"+basename+"_before").c_str()) &&
                                (TH1F*)  file.Get((directory+"/"+basename+"/"+basename+"_after").c_str()) ) && !refill) {
                before =(TH1F *) file.Get((directory+"/"+basename+"/"+basename+"_before").c_str());
                after  =(TH1F *) file.Get((directory+"/"+basename+"/"+basename+"_after").c_str());
		Refill=false;
		allfound =true;;
	}
	else allfound=false;
	return allfound;
}

void Efficiency::Fill(TTree * tree, Variables * vars, float (*discr_var) (Variables * vars),bool refill){
	if(ReinitializeHistos(refill)) return;
	cout<<basename.c_str()<<" Filling ..."<< endl;
	vars->ReadBranches(tree);
	for(int i=0;i<tree->GetEntries();i++){
		if(i%(int)FRAC!=0) continue;
		UpdateProgressBar(i, tree->GetEntries());
		tree->GetEvent(i);
		vars->Update();
		if(IsMC(vars)) { 
			if(BadEvSim) BadEvSim->LoadEvent(vars);		
			FillEventByEventMC(vars,discr_var,discr_var);
			}
		else	       FillEventByEventData(vars,discr_var,discr_var);	
	}

	return;
}

void Efficiency::Fill(DBarReader readerMC, Variables * vars, float (*discr_var) (Variables * vars),bool refill){
        if(ReinitializeHistos(refill)) return;
	if(readerMC.GetTree()->GetNbranches()>11) {Fill(readerMC.GetTree(),vars,discr_var,refill); return;}	
        cout<<basename.c_str()<<" Filling ..."<< endl;
	for(int i=0;i<readerMC.GetTreeEntries();i++){
		if(i%(int)FRAC!=0) continue; // WTF ?!
		UpdateProgressBar(i, readerMC.GetTreeEntries());
		readerMC.FillVariables(i,vars);
		vars->Update();

		if(IsMC(vars)) {
			if(BadEvSim) BadEvSim->LoadEvent(vars);
			FillEventByEventMC(vars,discr_var,discr_var);
		}
		else           FillEventByEventData(vars,discr_var,discr_var);
	}

	return;
}



void Efficiency::FillEventByEventMC(Variables * vars, float (*var) (Variables * vars), float (*discr_var) (Variables * vars)){
	if(!Refill) return;
	int kbin;
	float beta=discr_var(vars);;
		
	if(BadEvSim) beta=BadEvSim->SimulateBadEvents(beta);

	if(bins.IsUsingBetaEdges()) kbin = bins.GetBin(beta);
	else kbin =bins.GetBin(discr_var(vars));
	if(kbin>0){
			if(ApplyCuts(cut_before,vars)) { before->Fill(kbin,vars->mcweight);}
			if(ApplyCuts(cut_after ,vars)) { after ->Fill(kbin,vars->mcweight);}
	}
	return;
}


void Efficiency::FillEventByEventData(Variables * vars, float (*var) (Variables * vars), float (*discr_var) (Variables * vars)){
	if(!Refill) return;
	int kbin =bins.GetBin(discr_var(vars));
	if(kbin>0){
		if(after->GetNbinsY()==1){
			if(ApplyCuts(cut_before,vars)) before->Fill(kbin);
			if(ApplyCuts(cut_after ,vars)) after ->Fill(kbin);
		}
		else{
			if(ApplyCuts(cut_before,vars)) before->Fill(kbin,GetLatitude(vars));
			if(ApplyCuts(cut_after ,vars)) after ->Fill(kbin,GetLatitude(vars));
		}

	}
	return;
}

void Efficiency::Save(FileSaver finalhistos){
	finalhistos.Add(before);
	finalhistos.Add(after); 	
 	if(Stat_Error) finalhistos.Add(Stat_Error);
	if(Syst_Error) finalhistos.Add(Syst_Error);
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
 	if(Stat_Error) finalhistos.Add(Stat_Error);
	if(Syst_Error) finalhistos.Add(Syst_Error);
	finalhistos.writeObjsInFolder((directory+"/"+basename).c_str());
	if(fitrequested){
		cout<<FittedEff<<endl;
		finalhistos.Add(FittedEff);
		finalhistos.writeObjsInFolder((directory+"/"+basename).c_str());
	}
}

#endif

