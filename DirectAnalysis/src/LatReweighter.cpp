#include "LatReweighter.h"
#include "filesaver.h"

LatReweighter::LatReweighter(FileSaver FinalHisto,std::string Basename){
	AllOrbit = (TH1F*) FinalHisto.Get(("Spectra/"+Basename+"_AllOrbit").c_str());
	HighLat  = (TH1F*) FinalHisto.Get(("Spectra/"+Basename+"_HighLat").c_str());
	Weights  = (TH1F*) FinalHisto.Get(("Spectra/"+Basename+"_Weights").c_str());
	cout<<"******************************************************************"<<endl;
	cout<<"Reweighter found, based on "<<Weights->GetEntries()<<" Events"<<endl;
	cout<<"******************************************************************"<<endl;
	basename = Basename;
};

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
		//vars->PrintCurrentState();
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


void LatReweighter::CalculateWeights(){
	Weights = (TH1F *)HighLat->Clone();
	Weights->Scale(AllOrbit->GetBinContent(AllOrbit->FindBin(25))/(float)HighLat->GetBinContent(HighLat->FindBin(25)));
	Weights->Divide(AllOrbit);
	Weights->SetName((basename+"_Weights").c_str());
	Weights->SetTitle((basename+"_Weights").c_str());

	for(int i=0;i<Bins.size();i++)
		for(int j=0;j<Bins.size();j++){
			if(i<=j) ExposureMatrix->SetBinContent(i+1,j+1,(ExposureTime->GetBinContent(i+1)));
			else ExposureMatrix->SetBinContent(i+1,j+1,(ExposureTime->GetBinContent(j+1)));	
		}

};

float LatReweighter::GetWeight( float R){
	if(Weights->GetEntries()>1)
	return Weights->GetBinContent(Weights->FindBin(R)+1);
	else return 1;
}

void LatReweighter::Save(FileSaver finalhisto){
	finalhisto.Add(AllOrbit);
	finalhisto.Add(HighLat);
	finalhisto.Add(Weights);
	finalhisto.Add(ExposureTime);
	finalhisto.Add(ExposureMatrix);
	finalhisto.writeObjsInFolder("Spectra");	
}

