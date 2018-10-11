#ifndef FLUX_H
#define FLUX_H

#include "Livetime.h"
#include "Cuts.h"
#include "binning.h"
#include <fstream>
#include "rundb.h"

struct MCPar{
float Rmin,Rmax,Trigrate,gen_factor;
long int tot_ev,tot_trig;
std::string filename;
void Eval_trigrate();
};

void MCPar::Eval_trigrate(){
	std::vector<float> events;
	std::vector<float> triggers;
	std::ifstream infile;
	rundb* rdb= new rundb();
	int ret=rdb->readdb(("/cvmfs/ams.cern.ch/Offline/AMSDataDir/DataManagement/DataSetsDesc/"+filename).c_str());
	//rdb->Print();
    	rdb->Summary();	
	cout<<"********MC infos:*******"<<endl;
	cout<<"Trig. Rate: "<<rdb->GetTrigRate()<<endl;		
        cout<<"Total Gen: "<<rdb->GetTotTrigg()<<endl; 
	cout<<"Total events: "<<rdb->GetTotEvents()<<endl;
	cout<<"************************"<<endl;
        Trigrate=rdb->GetTrigRate();
	tot_ev = rdb->GetTotEvents();	
	tot_trig = rdb->GetTotTrigg();
}


	
class Flux{

	private:
	Efficiency * MCEfficiency;
	TH1F * Counts=0x0;
	TH1F * Eff_Acceptance=0x0;
	TH1F * ExposureTime=0x0;
	Binning bins;
	std::string basename;	
	std::string exposurename;

	TH1F * FluxEstim=0x0;
	MCPar param;

	TH1F * Acc_StatErr=0x0;
	TH1F * Acc_SystErr=0x0;
	TH1F * Counts_Err=0x0;

	std::vector<TH1F*> EfficiencyCorrections;

	public:


	Flux(FileSaver File, FileSaver FileRes, std::string Basename,std::string Effname, std::string EffDir,std::string CountsName,std::string ExposureName,Binning Bins){
		MCEfficiency = new Efficiency(FileRes,Effname,EffDir,Bins);
		if(FileRes.CheckFile()) Counts = (TH1F *) FileRes.Get((CountsName).c_str());	 
		ExposureTime = (TH1F *) File.Get(("Fluxes/"+Basename+"/"+ExposureName).c_str());
		cout<<("Fluxes/"+Basename+"/"+ExposureName).c_str()<<" "<<ExposureTime<<endl;

		bins = Bins;		
		basename = Basename;
		exposurename = ExposureName;
	}
	Flux(FileSaver FileRes, std::string Basename, std::string Effname, std::string EffDir,std::string CountsName,std::string ExposureName, Binning Bins){
		MCEfficiency = new Efficiency(FileRes,Effname,EffDir,Bins);
		if(FileRes.CheckFile()) {
			Counts = (TH1F *) FileRes.Get((CountsName).c_str());
			ExposureTime = (TH1F *) FileRes.Get(("Fluxes/"+Basename+"/"+ExposureName).c_str());
			Eff_Acceptance = (TH1F *) FileRes.Get(("Fluxes/"+Basename+"/Eff. Acceptance").c_str());	
			Acc_StatErr= (TH1F *) FileRes.Get(("Fluxes/"+Basename+"/Acceptance Stat. Error").c_str());
			Acc_SystErr= (TH1F *) FileRes.Get(("Fluxes/"+Basename+"/Acceptance Syst. Error").c_str());
			Counts_Err= (TH1F *) FileRes.Get(("Fluxes/"+Basename+"/Counts Error").c_str());
			FluxEstim = (TH1F *) FileRes.Get(("Fluxes/"+Basename+"/"+Basename+"_Flux").c_str());
		}
		bins = Bins;		
		basename = Basename;
	}
	bool ReinitializeHistos(bool refill) {
		if(!ExposureTime||refill) { 
			ExposureTime = new TH1F(exposurename.c_str(),exposurename.c_str(),bins.size(),0,bins.size());
			Counts=0x0;
			return false;
		}
		else return true;	
	}
	void Set_MCPar(float rmin, float rmax, float Gen_factor, std::string Filename);
	void Eval_ExposureTime(Variables * vars, TTree * treeDT,FileSaver finalhistos,bool refill);
	void ApplyEfficCorr(TH1F * Correction);
	void Eval_Flux();
	void SaveResults(FileSaver finalhistos);
	void ChangeName (std::string newname) {basename = newname; return;}
	Binning GetBins(){return bins;}
	std::string GetName(){return basename;}

	
	TH1F * GetExposureTime(){return ExposureTime;}
	TH1F * GetFlux(){return FluxEstim;}
	TH1F * GetEffAcceptance(){return Eff_Acceptance;}

	TH1F * GetAcc_StatErr() { return  Acc_StatErr;}
	TH1F * GetAcc_SystErr()	{return Acc_SystErr;}
	TH1F * GetCounts_Err() {return Counts_Err;}

	TH1F * Eval_FluxRatio(Flux * Denominator,std::string name);
};

void Flux::Set_MCPar(float rmin, float rmax, float Gen_factor, std::string Filename){
	cout<<"Setting MC parameters"<<endl;
	param.Rmin=rmin;
	param.Rmax=rmax;
	param.filename = Filename;
	param.gen_factor = Gen_factor;
	param.Eval_trigrate();
}


void Flux::ApplyEfficCorr(TH1F * Correction){
	cout<<"Correction: "<<basename<<" "<<Correction->GetEntries()<<endl;
	EfficiencyCorrections.push_back(Correction);
	return;
}

TH1F * EvalEffAcc(Efficiency* Eff, Binning bins, MCPar param){
	
	//float normalization = param.tot_ev/FRAC*(pow(param.Trigrate,-1));
	float normalization = param.tot_trig/FRAC;
	TH1F * EffAcc = (TH1F *) (Eff->GetAfter())->Clone();	
	for(int i=0;i<bins.size();i++){
		if(EffAcc -> GetBinContent(i+1)>0){
			float range;
			range = log(param.Rmax)-log(param.Rmin);
			   	
			float gen_bins= normalization*(log(bins.RigBins()[i+1])-log(bins.RigBins()[i]))/range; 
	
			EffAcc -> SetBinContent(i+1,EffAcc -> GetBinContent(i+1)/gen_bins); 
			EffAcc -> SetBinError(i+1,0);//pow(GenSpectrum -> GetBinContent(i+1),0.5)/gen_bins); 			
		}
	}
	EffAcc ->SetName("Eff. Acceptance");
	EffAcc ->SetTitle("Eff. Acceptance"); 
	EffAcc -> Sumw2();
	EffAcc -> Scale(47.78/param.gen_factor);


	cout<<"*********** Eff. Acceptance parameters ************"<<endl;
	cout<<"TOT. Ev. gen: "<<normalization<<endl;
	cout<<"Trig. rate: "<<param.Trigrate<<endl;
	return EffAcc;
}

void Flux::Eval_Flux(){
	cout<<"Counts "<<Counts<<endl;
	
	if(Counts>0) {			FluxEstim = (TH1F *) Counts->Clone();
					Counts_Err  = (TH1F *) Counts->Clone();
					for(int i=0;i<Counts_Err->GetNbinsX();i++) {
						if(Counts_Err->GetBinContent(i+1)>0&&Counts_Err->GetBinError(i+1)>0)
							Counts_Err->SetBinContent(i+1,Counts_Err->GetBinError(i+1)/Counts_Err->GetBinContent(i+1));
						else Counts_Err->SetBinContent(i+1,0);	
						Counts_Err->SetBinError(i+1,0);
					}
					Counts_Err->SetName("Counts Error");
        				Counts_Err->SetTitle("Counts Error");

	}
	else { FluxEstim = new TH1F("dummy","dummy",10,0,10);
		return;
	}	
	FluxEstim -> SetName((basename+"_Flux").c_str());
	FluxEstim -> SetTitle((basename+"_Flux").c_str());
	FluxEstim -> Sumw2();
	
	if(MCEfficiency->GetStat_Error()){	
		Acc_StatErr = (TH1F *) MCEfficiency->GetStat_Error()->Clone();
		Acc_SystErr = (TH1F *) MCEfficiency->GetSyst_Error()->Clone();
		Acc_StatErr->SetName("Acceptance Stat. Error");
		Acc_StatErr->SetTitle("Acceptance Stat. Error");	
		Acc_SystErr->SetName("Acceptance Syst. Error");
		Acc_SystErr->SetTitle("Acceptance Syst. Error");	
	}

	if(ExposureTime) FluxEstim -> Divide(ExposureTime);
	
	Eff_Acceptance = (TH1F *) MCEfficiency->GetAfter()->Clone();
	if(Eff_Acceptance){
	
	         Eff_Acceptance = EvalEffAcc(MCEfficiency,bins,param);	
		 for(int i=0;i<EfficiencyCorrections.size();i++) Eff_Acceptance -> Multiply(EfficiencyCorrections[i]); 
	        FluxEstim -> Divide(Eff_Acceptance);
		
	}

	for(int i=0;i<bins.size();i++){
		FluxEstim->SetBinError(i+1,FluxEstim->GetBinError(i+1)/(bins.EkPerMasBins()[i+1]-bins.EkPerMasBins()[i]));
		FluxEstim->SetBinContent(i+1,FluxEstim->GetBinContent(i+1)/(bins.EkPerMasBins()[i+1]-bins.EkPerMasBins()[i]));
	}
	return;		
}

void Flux::SaveResults(FileSaver finalhistos){
	std::cout<<"Results: "<<FluxEstim<<" "<<ExposureTime<<" "<<Eff_Acceptance<<std::endl;

	if(FluxEstim) finalhistos.Add(FluxEstim); 	
	if(ExposureTime) finalhistos.Add(ExposureTime);
	if(Eff_Acceptance)finalhistos.Add(Eff_Acceptance);
	if(Acc_StatErr) finalhistos.Add(Acc_StatErr);
	if(Acc_SystErr) finalhistos.Add(Acc_SystErr);
	if(Counts_Err) finalhistos.Add(Counts_Err);
	finalhistos.writeObjsInFolder(("Fluxes/"+basename).c_str());
}


void Flux::Eval_ExposureTime(Variables * vars, TTree * treeDT,FileSaver finalhistos,bool refill){
		
	if(!ExposureTime||refill){
	
		ExposureTime = new TH1F(exposurename.c_str(),exposurename.c_str(),bins.size(),0,bins.size());

		int ActualTime=0;
		vars->ReadBranches(treeDT);
		for(int i=0;i<treeDT->GetEntries()/FRAC;i++){
			UpdateProgressBar(i, treeDT->GetEntries()/FRAC);
			treeDT->GetEvent(i);
			vars->Update();
			if((int)vars->U_time!=ActualTime) {
				UpdateZoneLivetime(vars->Livetime,vars->Rcutoff,ExposureTime,bins);
				ActualTime=vars->U_time;
			}
		}
		finalhistos.Add(ExposureTime); 	
		finalhistos.writeObjsInFolder(("Fluxes/"+basename).c_str());

	}
	return;
}


TH1F * Flux::Eval_FluxRatio(Flux * Denominator,std::string name){

	TH1F * Numerator = (TH1F*) FluxEstim -> Clone();
	Numerator->SetName(name.c_str());
	Numerator->SetTitle(name.c_str());
	Numerator->Sumw2();

	Numerator->Divide(Denominator->GetFlux());
	for(int i=0;i<Numerator->GetNbinsX();i++){
		float A = FluxEstim -> GetBinContent(i+1);
		float B = Denominator->GetFlux() -> GetBinContent(i+1);	
		float sigA = FluxEstim -> GetBinError(i+1);
		float sigB = Denominator->GetFlux() -> GetBinError(i+1);	
		Numerator->SetBinError(i+1,pow(pow(sigA/A,2)+pow(sigB/B,2)+2*sigA/A*sigB/B,0.5)*A/B);
	}

	return Numerator;
}

#endif
