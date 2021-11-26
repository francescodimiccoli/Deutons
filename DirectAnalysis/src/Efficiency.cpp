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
	//MCWeighting

	float weight =1;
	if(notweighted) weight=1;

	else if(notweighted_phtr) {
		if(((int)vars->joinCutmask&1)!=1 &&  (((int)vars->joinCutmask&4)==4)) weight = 1/100.;	
	}		
	else{ 
		weight = vars->mcweight;
		//if(bins.IsUsingBetaEdges()) weight *= vars->GetCutoffCleaningWeight(GetRFromBeta(bins.getParticle().getMass(),var(vars)),vars->Momento_gen,BetacutoffCut);
		//else weight *= vars->GetCutoffCleaningWeight(vars->RInner,vars->Momento_gen,RcutoffCut);
		weight *= vars->GetTimeDepWeight(vars->Momento_gen);				      
	}

	if(kbin>=0) if(ApplyCuts(cut_before,vars)){before->Fill(kbin,weight);}

	kbin =bins.GetBin(discr_var(vars));

	if(kbin>=0)	if(ApplyCuts(cut_after ,vars)) after ->Fill(kbin,weight);


	return;
}


void Efficiency::FillEventByEventData(Variables * vars, float (*var) (Variables * vars), float (*discr_var) (Variables * vars)){
	if(!Refill) return;
	if(!ApplyCuts("IsData",vars)) return;
	int kbin =bins.GetBin(discr_var(vars));
	if(kbin>=0){
		//se la richiesta di essere primario è attiva, escludo gli eventi con cutoff intra-bin
		/*if(cut_before.find("Primary")!=std::string::npos){
			if(bins.IsUsingBetaEdges()) { if(BetacutoffCut * vars->Rcutoff_IGRFRTI > bins.RigBin(kbin)) return; }
			else { if(RcutoffCut * vars->Rcutoff_IGRFRTI > bins.RigBin(kbin)) return;}
		}*/
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
	Eff -> Divide(before);

	TH1F * Err = (TH1F *)after->Clone();
	Err->Reset();
	
	for(int i=0;i<Err->GetNbinsX();i++) {
		float b = before->GetBinContent(i+1);
		float a = after->GetBinContent(i+1);
		float db = before->GetBinError(i+1);
		float da=after->GetBinError(i+1);
		if(a!=0&&b!=0 && (da*da/(b*b) + pow(a,2)/pow(b,4)*db*db - 2*db*da*a/pow(b,3))>0) 
			Err->SetBinContent(i+1,pow(da*da/(b*b) + pow(a,2)/pow(b,4)*db*db - 2*db*da*a/pow(b,3) ,0.5));
	}
	Err->Smooth();
	for(int i=0;i<Eff->GetNbinsX();i++) Eff->SetBinError(i+1,Err->GetBinContent(i+1));

	Eff ->SetName((basename+"_Eff").c_str());
	Eff ->SetTitle((basename+" Efficiency").c_str());
	cout<<"EFFIC CALC "<<basename<<" "<<before->GetEntries()<<" "<<after->GetEntries()<<endl;
	return;
}

void Efficiency::Eval_TrigEfficiency(){
	TH1 * Unbias = (TH1 *) before->Clone();
	Unbias->Sumw2();
	Unbias->Add(after,-1);
	Unbias->Scale(100);
	for(int i=0; i<Unbias->GetNbinsX();i++) 
		Unbias->SetBinError(i+1,0.5*pow(Unbias->GetBinContent(i+1),0.5));


	TH1 * Denominator = (TH1F *) after->Clone();
	Denominator->Sumw2();
	Denominator->Add(Unbias);

	Eff = (TH1 *)after->Clone();	
	Eff -> Sumw2();
        Eff -> Divide(Denominator);
       	for(int i=0;i<Eff->GetNbinsX();i++) {
		float a = Unbias->GetBinContent(i+1);
		float b = after->GetBinContent(i+1);
		float da = Unbias->GetBinError(i+1);
		float db=  after->GetBinError(i+1);
		if(a!=0&&b!=0) 
			Eff->SetBinError(i+1,100*a/pow((b+100*a),2)*db+(100*b/pow((b+100*a),2))*da);
	}


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
	if(before->GetNbinsY()<=1){
	if(before)	before = ConvertBinnedHisto((TH1F*)before,before->GetTitle(),bins,true);	
	if(after)	after = ConvertBinnedHisto((TH1F*)after,after->GetTitle(),bins,true);	
	if(Eff)		Eff= ConvertBinnedHisto((TH1F*)Eff,Eff->GetTitle(),bins,true);
	}

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


