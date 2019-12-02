#include "Efficiency.h"


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
	bool allfound = false;
	TFile* File;
	if(IsExtern){
		File = TFile::Open((outdir+"/ExternalMCEff.root").c_str());
		if(!File)  return allfound;
	} 
	
	else{
		if(!file.CheckFile()) return allfound;
		File = file.GetFile();
	}
	if(( (TH1F*)  File->Get((directory+"/"+basename+"/"+basename+"_before").c_str()) &&
			          (TH1F*)  File->Get((directory+"/"+basename+"/"+basename+"_after").c_str()) ) && !refill) {
                before =(TH1F *) File->Get((directory+"/"+basename+"/"+basename+"_before").c_str());
                after  =(TH1F *) File->Get((directory+"/"+basename+"/"+basename+"_after").c_str());
		cout<<" FOUND: "<<(directory+"/"+basename+"/"+basename+"_after").c_str()<<endl;;
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
		else	     { FillEventByEventData(vars,discr_var,discr_var);};	
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

	if(ApplyCuts("IsData",vars)) return;
	int kbin;
	kbin =bins.GetBin(var(vars));
	if(kbin>=0){
			float weight =1;
			if(!notweighted) weight = vars->mcweight;
			if(ApplyCuts(cut_before,vars)){ before->Fill(kbin,weight);}
	}

	kbin =bins.GetBin(discr_var(vars));

	if(kbin>=0){
			float weight =1;
			if(!notweighted) weight = vars->mcweight;
			if(ApplyCuts(cut_after ,vars)) after ->Fill(kbin,weight);
	}

return;
}


void Efficiency::FillEventByEventData(Variables * vars, float (*var) (Variables * vars), float (*discr_var) (Variables * vars)){
	if(!Refill) return;
	if(!ApplyCuts("IsData",vars)) return;
	int kbin =bins.GetBin(discr_var(vars));
	if(kbin>=0){
		if(after->GetNbinsY()==1){
			if(ApplyCuts(cut_before,vars)) before->Fill(kbin,vars->PrescaleFactor);
			if(ApplyCuts(cut_after ,vars)) after ->Fill(kbin,vars->PrescaleFactor);
		}
		else{
			if(ApplyCuts(cut_before,vars)) ((TH2*)before)->Fill(kbin,GetLatitude(vars),vars->PrescaleFactor);
			if(ApplyCuts(cut_after ,vars)) ((TH2*)after) ->Fill(kbin,GetLatitude(vars),vars->PrescaleFactor);
		}

	}
	return;
}

void Efficiency::Save(){
	cout<<"Saving histo with "<<before->GetEntries()<<"entries"<<endl;
	finalhistos.Add(before);
	finalhistos.Add(after); 	
//	if(Stat_Error) finalhistos.Add(Stat_Error);
//	if(Syst_Error) finalhistos.Add(Syst_Error);
        finalhistos.writeObjsInFolder((directory+"/"+basename).c_str());
}


void Efficiency::Eval_Efficiency(){
	Eff = (TH1 *)after->Clone();
	Eff -> Sumw2();
	for(int i=0;i<Eff->GetNbinsX();i++) Eff->SetBinError(i+1,0);
	Eff -> Divide(before);
	Eff ->SetName((basename+"_Eff").c_str());
	Eff ->SetTitle((basename+" Efficiency").c_str());
	return;
}

void Efficiency::Eval_TrigEfficiency(){
	TH1 * Unbias = (TH1 *) before->Clone();
	Unbias->Sumw2();
	Unbias->Add(after,-1);
	for(int i=0; i<Unbias->GetNbinsX();i++) 
		Unbias->SetBinError(i+1,2*pow(Unbias->GetBinContent(i+1),0.5));

	Unbias->Scale(100);
	
	TH1 * Denominator = (TH1F *) after->Clone();
	Denominator->Sumw2();
	Denominator->Add(Unbias);

	Eff = (TH1 *)after->Clone();	
	Eff -> Sumw2();
        Eff -> Divide(Denominator);
        Eff ->SetName((basename+"_Eff").c_str());
        Eff ->SetTitle((basename+" Efficiency").c_str());

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
	finalhistos.Add(before);
        finalhistos.Add(after);
	if(Eff)	finalhistos.Add(Eff); 	
 	//if(Stat_Error) finalhistos.Add(Stat_Error);
	//if(Syst_Error) finalhistos.Add(Syst_Error);
	finalhistos.writeObjsInFolder((directory+"/"+basename).c_str());
	if(fitrequested){
		cout<<FittedEff<<endl;
		finalhistos.Add(FittedEff);
		finalhistos.writeObjsInFolder((directory+"/"+basename).c_str());
	}
}


