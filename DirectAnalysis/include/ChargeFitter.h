#ifndef CHARGEFITER_H
#define CHARGEFITTER_H

#include "TH1.h"
#include "binning.h"

class ChargeFitter{

	private:
	Binning bins;
	std::vector<TH1F*> DataDistrib;
	std::vector<TH1F*> TemplatesQ1;
	std::vector<TH1F*> TemplatesQ2;
	std::vector<TH1F*> TemplatesD;
	std::vector<TH1F*> RawFragments;
	std::vector<TH1F*> Q1Mass;
	std::vector<TH1F*> PMassMC;
	std::vector<TH1F*> DMassMC;
	std::vector<TH1F*> TMassMC;
	
	std::vector<TH1F*> TemplatesQ2MC;
	std::vector<TH1F*> FragmentsMC;

	std::vector<TH1F*> FractionFits;

	BadEventSimulator * BadEvSim=0x0;

	std::string CutData,CutTempl1,CutTempl2,CutFragm,CutMass;
	std::string basename;
	FileSaver file;

	public:
	ChargeFitter(FileSaver File, std::string Basename,Binning Bins, std::string cutdata, std::string cuttempl1, std::string cuttempl2,std::string cutfragm,std::string cutmass,int nbinsx,int xmin, int xmax);
	void FillEventByEventData(Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars));
	void FillEventByEventMC(Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars));
	void SetUpBadEventSimulator(BadEventSimulator * Sim) {BadEvSim = Sim; return; };
	void LoadEventIntoBadEvSim(Variables * vars) {if(BadEvSim) BadEvSim->LoadEvent(vars);}
	bool ReinitializeHistos(bool refill);
	void Save(FileSaver finalhisto);
	void FitDistributions(float rangemin, float rangemax);
	void SaveResults(FileSaver finalResults);

	void CreateFragmentsMass(float cutmin);

};

ChargeFitter::ChargeFitter(FileSaver File, std::string Basename,Binning Bins, std::string cutdata, std::string cuttempl1, std::string cuttempl2,std::string cutfragm,std::string cutmass, int nbinsx,int xmin, int xmax){

	file = File;
	bins=Bins;
	CutData=cutdata;
	CutTempl1=cuttempl1;
	CutTempl2=cuttempl2;
	CutFragm=cutfragm;
	CutMass=cutmass;
	basename=Basename;

	for(int bin=0; bin<bins.size();bin++){
		DataDistrib.push_back(new TH1F((basename+"_Data_bin_"+to_string(bin)).c_str(),(basename+"Data_bin_"+to_string(bin)).c_str(),nbinsx,xmin,xmax ));
		TemplatesQ1.push_back(new TH1F((basename+"_Q1_bin_"+to_string(bin)).c_str(),(basename+"Q1_bin_"+to_string(bin)).c_str(),nbinsx,xmin,xmax ));
		TemplatesQ2.push_back(new TH1F((basename+"_Q2_bin_"+to_string(bin)).c_str(),(basename+"Q2_bin_"+to_string(bin)).c_str(),nbinsx,xmin,xmax ));
		TemplatesQ2MC.push_back(new TH1F((basename+"_Q2MC_bin_"+to_string(bin)).c_str(),(basename+"Q2MC_bin_"+to_string(bin)).c_str(),nbinsx,xmin,xmax ));
		TemplatesD.push_back(new TH1F((basename+"_D_bin_"+to_string(bin)).c_str(),(basename+"D_bin_"+to_string(bin)).c_str(),nbinsx,xmin,xmax ));
		RawFragments.push_back(new TH1F((basename+"_Fragm_bin_"+to_string(bin)).c_str(),(basename+"Fragm_bin_"+to_string(bin)).c_str(),nbinsx,xmin,xmax ));
		Q1Mass.push_back(new TH1F((basename+"_Mass_bin_"+to_string(bin)).c_str(),(basename+"Mass_bin_"+to_string(bin)).c_str(),nbinsx,xmin,xmax ));
		PMassMC.push_back(new TH1F((basename+"_PMassMC_bin_"+to_string(bin)).c_str(),(basename+"PMassNC_bin_"+to_string(bin)).c_str(),nbinsx,xmin,xmax ));
		DMassMC.push_back(new TH1F((basename+"_DMassMC_bin_"+to_string(bin)).c_str(),(basename+"DMassMC_bin_"+to_string(bin)).c_str(),nbinsx,xmin,xmax ));
		TMassMC.push_back(new TH1F((basename+"_TMassMC_bin_"+to_string(bin)).c_str(),(basename+"TMassMC_bin_"+to_string(bin)).c_str(),nbinsx,xmin,xmax ));
		
		
		FragmentsMC.push_back(new TH1F((basename+"_FragmMC_bin_"+to_string(bin)).c_str(),(basename+"FragmMC_bin_"+to_string(bin)).c_str(),nbinsx,xmin,xmax ));

		FractionFits.push_back(new TH1F((basename+"_fit_bin_"+to_string(bin)).c_str(),(basename+"fit_bin_"+to_string(bin)).c_str(),nbinsx,xmin,xmax ));
	}
}

bool ChargeFitter::ReinitializeHistos(bool refill){
	bool allfound=false;
	TFile * File = file.GetFile();	
	if(!File) return false;
	for(int bin=0; bin<bins.size();bin++){

		if(( (TH1F*)  File->Get((basename +"/L1Distrib./"+basename+"_Data_bin_"+to_string(bin)).c_str())) && 
		   ( (TH1F*)  File->Get((basename +"/Q1Templates/"+basename+"_Q1_bin_"+to_string(bin)).c_str())) && 
		   ( (TH1F*)  File->Get((basename +"/Q2Templates/"+basename+"_Q2_bin_"+to_string(bin)).c_str())) &&  
		   ( (TH1F*)  File->Get((basename +"/Q2TemplatesMC/"+basename+"_Q2MC_bin_"+to_string(bin)).c_str())) &&  
		   ( (TH1F*)  File->Get((basename +"/Fragments/"+basename+"_Fragm_bin_"+to_string(bin)).c_str())) &&  
		   ( (TH1F*)  File->Get((basename +"/FragmentsMC/"+basename+"_FragmMC_bin_"+to_string(bin)).c_str())) &&  
		   ( (TH1F*)  File->Get((basename +"/Mass/"+basename+"_Mass_bin_"+to_string(bin)).c_str()))

		&&!refill){ 
		
			DataDistrib[bin] = (TH1F*)  File->Get((basename +"/L1Distrib./"+basename+"_Data_bin_"+to_string(bin)).c_str());
			TemplatesQ1[bin] = (TH1F*) File->Get((basename +"/Q1Templates/"+basename+"_Q1_bin_"+to_string(bin)).c_str());
			TemplatesQ2[bin] = (TH1F*) File->Get((basename +"/Q2Templates/"+basename+"_Q2_bin_"+to_string(bin)).c_str());
			TemplatesQ2MC[bin] = (TH1F*) File->Get((basename +"/Q2TemplatesMC/"+basename+"_Q2MC_bin_"+to_string(bin)).c_str());
			RawFragments[bin] = (TH1F*) File->Get((basename +"/Fragments/"+basename+"_Fragm_bin_"+to_string(bin)).c_str());
			Q1Mass[bin] = (TH1F*) File->Get((basename +"/Mass/"+basename+"_Mass_bin_"+to_string(bin)).c_str());
			PMassMC[bin] = (TH1F*) File->Get((basename +"/PMassMC/"+basename+"_PMassMC_bin_"+to_string(bin)).c_str());
			DMassMC[bin] = (TH1F*) File->Get((basename +"/DMassMC/"+basename+"_DMassMC_bin_"+to_string(bin)).c_str());
			TMassMC[bin] = (TH1F*) File->Get((basename +"/TMassMC/"+basename+"_TMassMC_bin_"+to_string(bin)).c_str());

			FragmentsMC[bin] = (TH1F*) File->Get((basename +"/FragmentsMC/"+basename+"_FragmMC_bin_"+to_string(bin)).c_str());

			allfound =true;
		}
	}
	return allfound;
}


void ChargeFitter::FillEventByEventData(Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars)){
	int kbin;
        kbin =  bins.GetBin(var(vars)); //normal filling
	float mass = (vars->R/var(vars))*pow((1-pow(var(vars),2)),0.5);
        if(ApplyCuts(CutData,vars)&&kbin>0)
		DataDistrib[kbin]->Fill(GetL1Q(vars),vars->PrescaleFactor);	
	if(ApplyCuts(CutMass,vars)&&kbin>0)
                Q1Mass[kbin]->Fill(mass,vars->PrescaleFactor);
	if(ApplyCuts(CutFragm,vars)&&kbin>0){
		RawFragments[kbin]->Fill(mass,vars->PrescaleFactor);
		}
	
	kbin =  bins.GetBin(discr_var(vars)); //filling with beta @L1
	if(ApplyCuts(CutTempl1,vars)&&kbin>0){
  	      TemplatesQ1[kbin]->Fill(GetL2Q(vars),vars->PrescaleFactor);
	}
	if(ApplyCuts(CutTempl2,vars)&&kbin>0)
                TemplatesQ2[kbin]->Fill(GetL2Q(vars),vars->PrescaleFactor);
}


void ChargeFitter::FillEventByEventMC(Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars)){
	int kbin;
	kbin =  bins.GetBin(var(vars));
	float betasmear = var(vars);
	if(BadEvSim) betasmear=BadEvSim->SimulateBadEvents(betasmear);
	float mass = (vars->R/betasmear)*pow((1-pow(betasmear,2)),0.5);
	std::string CutMassMC = CutMass + "&IsHeliumMC";
	if(ApplyCuts(CutMassMC,vars)&&kbin>0)
                FragmentsMC[kbin]->Fill(mass,vars->mcweight);
	CutMassMC = CutMass + "&IsPurePMC";
	if(ApplyCuts(CutMassMC,vars)&&kbin>0)
                PMassMC[kbin]->Fill(mass,vars->mcweight);
	CutMassMC = CutMass + "&IsPurePMC";
	if(ApplyCuts(CutMassMC,vars)&&kbin>0)
                DMassMC[kbin]->Fill((1.875/0.938)*mass,vars->mcweight);
	CutMassMC = CutMass + "&IsPurePMC";
	if(ApplyCuts(CutMassMC,vars)&&kbin>0)
                TMassMC[kbin]->Fill((2.793/0.938)*mass,vars->mcweight);
	
	kbin =  bins.GetBin(discr_var(vars)); //filling with beta @L1
	std::string CutTempl2MC = CutData + "&IsHeliumMC";
	if(ApplyCuts(CutTempl2MC,vars)&&kbin>0)
                TemplatesQ2MC[kbin]->Fill(GetL2Q(vars),vars->mcweight);	
}


void ChargeFitter::Save(FileSaver finalhisto){

	for(int bin=0; bin<bins.size();bin++)
		finalhisto.Add(DataDistrib[bin]);
	finalhisto.writeObjsInFolder((basename+"/L1Distrib."));
	for(int bin=0; bin<bins.size();bin++)
		finalhisto.Add(TemplatesQ1[bin]);
	finalhisto.writeObjsInFolder((basename+"/Q1Templates"));

	for(int bin=0; bin<bins.size();bin++)
		finalhisto.Add(TemplatesQ2[bin]);
	finalhisto.writeObjsInFolder((basename+"/Q2Templates"));

	for(int bin=0; bin<bins.size();bin++)
                finalhisto.Add(TemplatesQ2MC[bin]);
        finalhisto.writeObjsInFolder((basename+"/Q2TemplatesMC"));
	for(int bin=0; bin<bins.size();bin++)
		finalhisto.Add(RawFragments[bin]);
	finalhisto.writeObjsInFolder((basename+"/Fragments"));
	for(int bin=0; bin<bins.size();bin++)
		finalhisto.Add(Q1Mass[bin]);
	finalhisto.writeObjsInFolder((basename+"/Mass"));
	
	for(int bin=0; bin<bins.size();bin++)
		finalhisto.Add(FragmentsMC[bin]);
	finalhisto.writeObjsInFolder((basename+"/FragmentsMC"));

	for(int bin=0; bin<bins.size();bin++)
		finalhisto.Add(PMassMC[bin]);
	finalhisto.writeObjsInFolder((basename+"/PMassMC"));
	for(int bin=0; bin<bins.size();bin++)
		finalhisto.Add(DMassMC[bin]);
	finalhisto.writeObjsInFolder((basename+"/DMassMC"));
	for(int bin=0; bin<bins.size();bin++)
		finalhisto.Add(TMassMC[bin]);
	finalhisto.writeObjsInFolder((basename+"/TMassMC"));

}


void ChargeFitter::FitDistributions(float rangemin, float rangemax){

	for(int bin=2; bin<bins.size();bin++){
		TObjArray *Tpl;
        	Tpl = new TObjArray(2);	
		Tpl->Add(TemplatesQ1[bin]);
		Tpl->Add(TemplatesQ2[bin]);		

		TFractionFitter * Fit = new TFractionFitter(DataDistrib[bin], Tpl ,"q");
		Fit -> SetRangeX(DataDistrib[bin] -> FindBin(rangemin), DataDistrib[bin] -> FindBin(rangemax));
		
		float outcome = 1;
		bool fitcondition = (DataDistrib[bin] -> GetEntries()>0 && TemplatesQ1[bin]-> GetEntries()>0 && TemplatesQ2[bin]-> GetEntries()>0);
		if(fitcondition) outcome=Fit->Fit();
		if(outcome == 0){
			TH1F * Result = (TH1F *) Fit-> GetPlot();
			FractionFits[bin] = Result;
			float itot= Result->Integral();
			double w1,e1 = 0;
			double w2,e2 = 0;
		
			////////////////
			Fit ->GetResult(0,w1,e1);
			Fit ->GetResult(1,w2,e2);
				
			float i1 = TemplatesQ1[bin] ->Integral(TemplatesQ1[bin] -> FindBin(rangemin), TemplatesQ1[bin] -> FindBin(rangemax));
			float i2 = TemplatesQ2[bin] ->Integral(TemplatesQ2[bin] -> FindBin(rangemin), TemplatesQ2[bin] -> FindBin(rangemax));

			cout<<w1<<" "<<w2<<": chisquare= "<<Fit-> GetChisquare()/(float) (Fit -> GetNDF())<<endl;

			TemplatesQ1[bin] -> Scale(w1*itot/i1);
			TemplatesQ2[bin] -> Scale(w2*itot/i2);
	
		} 	
			
	}
};

void ChargeFitter::CreateFragmentsMass(float cutmin){
	for(int bin=2; bin<bins.size();bin++){
		float falsefragments = TemplatesQ1[bin]->Integral(TemplatesQ1[bin]->FindBin(cutmin)+1,TemplatesQ1[bin]->FindBin(2.3)+1);
		Q1Mass[bin]->Scale((1/Q1Mass[bin]->Integral())*falsefragments);
		RawFragments[bin]->Add(Q1Mass[bin],-1);
	
	}	

}

void ChargeFitter::SaveResults(FileSaver finalResults){
	Save(finalResults);
	for(int bin=0; bin<bins.size();bin++)
                finalResults.Add(FractionFits[bin]);
        finalResults.writeObjsInFolder((basename+"/Fits"));
}
#endif
