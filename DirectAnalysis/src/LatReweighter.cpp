#include "LatReweighter.h"
#include "filesaver.h"

double fitfunc(double *x, double *p){
	float value;
	value = p[int((x[0]/150)*500/3)];
	return value;
}


LatReweighter::LatReweighter(FileSaver FinalHisto,std::string Basename){
	Binning bins(proton);
	Bins = bins;
	Bins.Reset();
	Bins.setBinsFromRigidity(100,0.5,50,ResponseTOF,0.00347548,5.8474);
	Bins.UseREdges();

	AllOrbit = (TH1F*) FinalHisto.Get(("Spectra/"+Basename+"_AllOrbit").c_str());
	HighLat  = (TH1F*) FinalHisto.Get(("Spectra/"+Basename+"_HighLat").c_str());
	Weights  = (TH1F*) FinalHisto.Get(("Spectra/"+Basename+"_Weights").c_str());
	ExposureTime = 	(TH1F*) FinalHisto.Get(("Spectra/"+Basename+"_ExpTime").c_str());
	if(!ExposureTime)	
		ExposureTime = 	(TH1F*) FinalHisto.Get(("Spectra/"+Basename+"_Weights").c_str());
	WeightModel = (TF1*) FinalHisto.Get("Spectra/weightmodel");

	TFile *f = TFile::Open((workdir + "/LatWeights/First2Yrs.root").c_str());	
	Reference = (TH1F*) f->Get(("Spectra/"+Basename+"_HighLat").c_str());

	cout<<"******************************************************************"<<endl;
	cout<<"Reweighter found, based on "<<Weights->GetEntries()<<" Events"<<endl;
	cout<<AllOrbit<<" "<<HighLat<<" "<<Weights<<" "<<ExposureTime<<" "<<Reference<<endl;

	cout<<"******************************************************************"<<endl;
	basename = Basename;
};


/*
void LatReweighter::LoopOnData(TTree * treeDT, Variables * vars, bool Refill){
	if(!Refill) return;
	cout<<" DATA Filling from NTuples ..."<< endl;
	vars->ReadBranches(treeDT);
	for(int i=0;i<treeDT->GetEntries()/FRAC;i++){
		UpdateProgressBar(i, treeDT->GetEntries()/FRAC);
		vars->ResetVariables();
		treeDT->GetEvent(i);
		//vars->Update();
		if(ApplyCuts(cut,vars)){
			AllOrbit->Fill(vars->R,vars->PrescaleFactor);
			if(fabs(vars->Latitude)>0.95) HighLat->Fill(vars->R,vars->PrescaleFactor);
		}
	}
}
void LatReweighter::LoopOnData(DBarReader readerDT, Variables * vars,bool Refill){
	if(!Refill) return;
	if(readerDT.GetTree()->GetNbranches()>11) {LoopOnData(readerDT.GetTree(),vars,Refill); return;}
	else 
		cout<<" DATA Filling ..."<< endl;
	for(int i=0;i<readerDT.GetTreeEntries()/FRAC;i++){
		UpdateProgressBar(i, readerDT.GetTreeEntries()/FRAC);
		readerDT.FillVariables(i,vars);
		vars->Update();
	//	vars->PrintCurrentState();
		if(ApplyCuts(cut,vars)){
			AllOrbit->Fill(vars->R,vars->PrescaleFactor);
			if(fabs(vars->Latitude)>0.95) HighLat->Fill(vars->R,vars->PrescaleFactor);
		}
	}
}


void LatReweighter::LoopOnRTI(DBarReader reader, Variables * vars,bool Refill){
	    if(!Refill) return;
		    cout<<" Exposure Time Filling ..."<< endl;
	    for(int i=0;i<reader.GetTreeEntries();i++){
		    UpdateProgressBar(i, reader.GetTreeEntries());
		    reader.FillVariables(i,vars);
		    vars->Update();
			    if(vars->good_RTI==0&&vars->isinsaa==0) UpdateZoneLivetime(vars->Livetime_RTI,vars->Rcutoff_RTI,ExposureTime,Bins);
	    }
    }
*/

void LatReweighter::CalculateWeights(){
	Weights = (TH1F *)HighLat->Clone();
	Weights->Scale(Reference->GetBinContent(Reference->FindBin(25))/(float)HighLat->GetBinContent(HighLat->FindBin(25)));
	Weights->Divide(Reference);
	Weights->SetName((basename+"_Weights").c_str());
	Weights->SetTitle((basename+"_Weights").c_str());
	ExposureTime->SetName((basename+"_ExpTime").c_str());
	ExposureTime->SetTitle((basename+"_ExpTime").c_str());


	ModelWithSpline();
};

float LatReweighter::GetCleaningWeight( float Rmeas, float Rgen, float SF){
	float min;
	if(Rmeas/SF<Rgen) min=Rmeas/SF;
	else min=Rgen;

	if(min>=40) min = 49;
	if(Rgen>=0) Rgen = 49;

	int kbin = Bins.GetBin(fabs(min));
	int kbin_gen = Bins.GetBin(fabs(Rgen));
	if(ExposureTime->GetEntries()>1)
		return ExposureTime->GetBinContent(kbin+1)/ExposureTime->GetBinContent(kbin_gen+1);
	return 1;
}

float LatReweighter::GetTimeDepWeight(float R){
	return WeightModel->Eval(R);
}


void LatReweighter::Save(FileSaver finalhisto){
	finalhisto.Add(AllOrbit);
	finalhisto.Add(HighLat);
	finalhisto.Add(Weights);
	finalhisto.Add(ExposureTime);
	finalhisto.writeObjsInFolder("Spectra");	
}

void LatReweighter::SaveResults(FileSaver finalhisto){
	finalhisto.Add(ExposureTime);
	finalhisto.Add(Weights);
	finalhisto.Add(WeightModel);
	finalhisto.writeObjsInFolder("Spectra");	

}	

void LatReweighter::ModelWithSpline(){
	
	//regularization of histo
	int nodes = Weights->GetNbinsX()/3+1;
	double p[nodes];
	double x[nodes];

	WeightModel = new TF1("weightmodel",fitfunc,0,150,nodes);
	
	int j=0;
	for(int i=0;i<Weights->GetNbinsX();i=i+5){
		p[j]=Weights->GetBinContent(i);
		cout<<p[j]<<endl;
		j++;
	}
	cout<<"**********************"<<endl;
	for(int i=0;i<nodes;i++) {
		WeightModel->SetParameter(i,p[i]);
		if(i>=28) WeightModel->FixParameter(i,(Weights->GetBinContent(82)+Weights->GetBinContent(83)+Weights->GetBinContent(84)+Weights->GetBinContent(85)+Weights->GetBinContent(86))/5);
	}
	Weights->Fit("weightmodel");

}
