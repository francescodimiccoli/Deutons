#ifndef FLUX_H
#define FLUX_H

#include "Livetime.h"
#include "Cuts.h"
#include "binning.h"

struct MCPar{
float Rmin,Rmax,Trigrate;
};


class Flux{

	private:
	Efficiency * MCEfficiency;
	TH1F * Counts=0x0;
	TH1F * Geom_Acceptance=0x0;
	TH1F * Triggers=0x0;
	TH1F * ExposureTime=0x0;
	Binning bins;
	std::string basename;	
	std::string exposurename;
	std::string geomname;

	TH1F * FluxEstim=0x0;
	MCPar param;

	public:


	Flux(FileSaver File, FileSaver FileRes, std::string Basename, std::string Effname, std::string EffDir,std::string CountsName,std::string ExposureName, std::string GeomName,Binning Bins){
		MCEfficiency = new Efficiency(FileRes,Effname,EffDir,Bins);
		if(FileRes.CheckFile()) Counts = (TH1F *) FileRes.Get((CountsName).c_str());	 
		ExposureTime = (TH1F *) File.Get(("Fluxes/"+Basename+"/"+ExposureName).c_str());
		Geom_Acceptance = (TH1F *) File.Get(("Fluxes/"+Basename+"/"+GeomName).c_str());
		Triggers = (TH1F *) File.Get(("Fluxes/"+Basename+"/Triggers").c_str());

		bins = Bins;		
		basename = Basename;
		exposurename = ExposureName;
		geomname = GeomName;
	}
	Flux(FileSaver FileRes, std::string Basename, std::string Effname, std::string EffDir,std::string CountsName,std::string ExposureName, std::string GeomName, Binning Bins){
		MCEfficiency = new Efficiency(FileRes,Effname,EffDir,Bins);
		if(FileRes.CheckFile()) {
			Counts = (TH1F *) FileRes.Get((CountsName).c_str());
			ExposureTime = (TH1F *) FileRes.Get(("Fluxes/"+Basename+"/"+ExposureName).c_str());
			Geom_Acceptance = (TH1F *) FileRes.Get(("Fluxes/"+Basename+"/"+GeomName).c_str());
			FluxEstim = (TH1F *) FileRes.Get(("Fluxes/"+Basename+"/"+Basename+"_Flux").c_str());
		}
		bins = Bins;		
		basename = Basename;
	}
	bool ReinitializeHistos(bool refill) {
		if(!ExposureTime||!Triggers||!Geom_Acceptance||refill) { 
			ExposureTime = new TH1F(exposurename.c_str(),exposurename.c_str(),bins.size(),0,bins.size());
                	Triggers = new TH1F("Triggers","Triggers",ForAcceptance.size(),0,ForAcceptance.size());
			Geom_Acceptance = new TH1F(geomname.c_str(),geomname.c_str(),bins.size(),0,bins.size());
			Counts=0x0;
			return false;
		}
		else return true;	
	}
	void Set_MCPar(float rmin, float rmax, float trigrate);
	void Eval_ExposureTime(Variables * vars, TTree * treeDT,FileSaver finalhistos,bool refill);
	void Eval_GeomAcceptance(TTree * treeMC, FileSaver finalhistos,std::string cut,bool refill, bool IsRigBin=false);
	void ApplyEfficCorr(TH1F * Correction);
	void Eval_Flux();
	void SaveResults(FileSaver finalhistos);
	Binning GetBins(){return bins;}
	std::string GetName(){return basename;}

	
	TH1F * GetExposureTime(){return ExposureTime;}
	TH1F * GetFlux(){return FluxEstim;}
	TH1F * GetTriggerCounts() {return Triggers;}
	TH1F * GetGenAcceptance(){return Geom_Acceptance;}

	TH1F * Eval_FluxRatio(Flux * Denominator,std::string name);
};

void Flux::Set_MCPar(float rmin, float rmax, float trigrate){
	param.Rmin=rmin;
	param.Rmax=rmax;
	param.Trigrate=trigrate;
}

void Flux::ApplyEfficCorr(TH1F * Correction){
	cout<<"Correction: "<<basename<<" "<<Correction<<endl;
	MCEfficiency->GetEfficiency()->Multiply(Correction);
	return;
}


void Flux::Eval_Flux(){
	cout<<"Counts "<<Counts<<endl;
	
	if(Counts>0) {cout<<"1"<<endl;	FluxEstim = (TH1F *) Counts->Clone();
	}
	else { FluxEstim = new TH1F("dummy","dummy",10,0,10);
	return;
	}
//	cout<<"Flux: "<<MCEfficiency->GetEfficiency()->GetEntries()<<endl;
	FluxEstim -> SetName((basename+"_Flux").c_str());
	FluxEstim -> SetTitle((basename+"_Flux").c_str());

	FluxEstim -> Sumw2();
	FluxEstim -> Divide(MCEfficiency->GetEfficiency());
	if(ExposureTime) FluxEstim -> Divide(ExposureTime);
	
	if(Geom_Acceptance){
		float totaltriggers = Triggers->Integral();
		for(int i=0;i<bins.size();i++){
			float gen_bins= totaltriggers*(pow(param.Trigrate,-1))*(log(bins.RigBins()[i+1])-log(bins.RigBins()[i]))/(log(param.Rmax)-log(param.Rmin)); 
			if(Geom_Acceptance -> GetBinContent(i+1)>0){
				Geom_Acceptance -> SetBinContent(i+1,Geom_Acceptance -> GetBinContent(i+1)/gen_bins); 			
				Geom_Acceptance -> SetBinError(i+1,pow(Geom_Acceptance -> GetBinContent(i+1),0.5)/gen_bins); 			
				}
			}
		 Geom_Acceptance -> Sumw2();
		 Geom_Acceptance -> Scale(47.78);
		 FluxEstim -> Divide(Geom_Acceptance);
		}

	for(int i=0;i<bins.size();i++){
		FluxEstim->SetBinError(i+1,FluxEstim->GetBinError(i+1)/(bins.EkPerMasBins()[i+1]-bins.EkPerMasBins()[i]));
		FluxEstim->SetBinContent(i+1,FluxEstim->GetBinContent(i+1)/(bins.EkPerMasBins()[i+1]-bins.EkPerMasBins()[i]));
	}
	return;		
}

void Flux::SaveResults(FileSaver finalhistos){
	finalhistos.Add(FluxEstim); 	
	finalhistos.Add(ExposureTime);
	finalhistos.Add(Geom_Acceptance);
	finalhistos.writeObjsInFolder(("Fluxes/"+basename).c_str());
}
void Flux::Eval_GeomAcceptance(TTree * treeMC,FileSaver finalhistos,std::string cut,bool refill, bool IsRigBin){

	if((!Geom_Acceptance||!Triggers||refill)){
			Geom_Acceptance = new TH1F(geomname.c_str(),geomname.c_str(),bins.size(),0,bins.size());
			Triggers = new TH1F("Triggers","Triggers",ForAcceptance.size(),0,ForAcceptance.size());		
			Variables * vars = new Variables;

			vars->ReadBranches(treeMC);
			int kbin;
			for(int i=0;i<treeMC->GetEntries()/FRAC;i++){
				if(i%(int)FRAC!=0) continue;
		                UpdateProgressBar(i, treeMC->GetEntries());
                		treeMC->GetEvent(i);
                		vars->Update();
				if(!IsRigBin) kbin=bins.GetBin(GetBetaGen(vars));
				else	      kbin=bins.GetBin(GetGenMomentum(vars));
				if(ApplyCuts(cut,vars)&&kbin>0) Geom_Acceptance->Fill(kbin);
				kbin=ForAcceptance.GetBin(GetGenMomentum(vars));
				if(ApplyCuts(cut,vars)&&kbin>0) Triggers->Fill(kbin);
			}
			/*
			Efficiency * Geom = new Efficiency(finalhistos,"","",bins,cut.c_str() ,cut.c_str());
			if(!IsRigBin) Geom->Fill(RawMC, vars,GetBetaGen,true);
			else Geom->Fill(RawMC, vars,GetGenMomentum,true);

			Efficiency * Trig = new Efficiency(finalhistos,"","",ForAcceptance,cut.c_str() ,cut.c_str());
			Trig->Fill(RawMC, vars,GetGenMomentum,true);

			Geom_Acceptance = (TH1F *)Geom->GetBefore()->Clone();
			Geom_Acceptance -> SetName(geomname.c_str());
			Geom_Acceptance -> SetTitle(geomname.c_str());

			Triggers = (TH1F *)Trig->GetBefore()->Clone();
			Triggers -> SetName("Triggers");
			Triggers -> SetTitle("Triggers");
			*/
			finalhistos.Add(Geom_Acceptance); 	
			finalhistos.Add(Triggers); 	
			finalhistos.writeObjsInFolder(("Fluxes/"+basename).c_str());

	}

	return;
};

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

	return Numerator;
}

#endif
