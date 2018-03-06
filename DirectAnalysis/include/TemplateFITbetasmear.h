#ifndef TEMPLATEFITBETASMEAR_H
#define TEMPLATEFITBETASMEAR_H

#include "BadEventSimulator.h"

using namespace std;

// not used
TSpline3 * ExtractCutoffWeight(TH1F * ExposureTime){

        ExposureTime->Scale(1/ExposureTime->GetBinContent(ExposureTime->GetMaximumBin()));

        double x[ExposureTime->GetNbinsX()];
        double y[ExposureTime->GetNbinsX()];

        for(int i=0;i<ExposureTime->GetNbinsX();i++){
                x[i]=ExposureTime->GetBinCenter(i+1);
                y[i]=ExposureTime->GetBinContent(i+1);
        }

        TSpline3 * CutoffWeight = new TSpline3("CutoffWeight",x,y,ExposureTime->GetNbinsX());
        CutoffWeight->SetName("CutoffWeight");
        return CutoffWeight;

}


struct TFit {
   TH1F * Templ_P =0x0;
   TH1F * Templ_D =0x0;
   TH1F * Templ_He=0x0;

   TH1F * Data =0x0;
   TH1F * DataPrim=0x0;	


   float wheightP,wheightD,wheightHe;
   float ContribP,ContribD,ContribHe;
   float errP,errD,errHe;
   float fitrangemin = 0.6;
   float fitrangemax = 4;
		
   TFractionFitter *Tfit;
   int Tfit_outcome=-1;
   
   float DCounts=0;
   float PCounts=0;
   float ChiSquare=500;			
   float StatErr=0;	

   TH1F * Templ_DPrim;
   TH1F * Templ_PPrim;
   
   float DCountsPrim=0;		
   float PCountsPrim=0;	



   TFit(){}	
   TFit(TH1F * templ_P, TH1F * templ_D, TH1F * data, TH1F * dataPrim) { Templ_P= templ_P; Templ_D=templ_D; Data=data; DataPrim=dataPrim; }

};


struct Systpar{
	int steps;
	float sigma;
	float shift;
};

struct BestChi {
	int i=0;
	int j=0;
	float chimin=0;
	void FindMinimum(TH2F * Histo) {
		float Best = 9999999;
		for(int x=0;x<Histo->GetNbinsX();x++)
			for(int y=0;y<Histo->GetNbinsY();y++)
				if(Histo->GetBinContent(x+1,y+1)<Best){
					Best=Histo->GetBinContent(x+1,y+1);
				        i=x;
					j=y;
					chimin=Best;	
				}
	}
};

class TemplateFIT {

	private:
	std::vector<std::vector<std::vector<TFit *>>> fits;
	std::vector<BestChi *> BestChiSquare;	
	std::vector<TH1F *> TransferFunction;


	std::vector<TH2F *> DCountsSpread;
	std::vector<TH1F *> WeightedDCounts;
	std::vector<TH2F *> TFitChisquare;
	TF1 * HeContModel;
	TH1F * MCHeContRatio;
	TH1F * MeasuredHeContRatio;

	TH1F * StatError;
	TH1F * SystError;
	TH1F * HeContError;

	TH1F * ProtonCounts;
	TH1F * DeuteronCounts;

	TH1F * ProtonCountsPrim;
	TH1F * DeuteronCountsPrim;


	TH1F * BestChiSquares;
	TH1F * OriginalChiSquares;

	TH1F * BestFitSigma;
	TH1F * BestFitShift;

	Binning bins;
        std::string cut;
	std::string cutprimary;
	std::string discr_var;

	std::string basename;

	TH1F * Exposure_Time;


	Systpar systpar;
	bool fitDisabled=false;
	bool isrich=false;
	BadEventSimulator * BadEvSim=0x0;
	bool IsFitNoise = false;
	float constrainmin[3];
	float constrainmax[3];


	public:	
	//standard constructor
	TemplateFIT(std::string Basename,Binning Bins, std::string Cut, int Nbins, float Xmin, float Xmax, bool IsRich=false ,int steps=11,float sigma=80,float shift=40,TH1F * ExposureTime=0x0){
		
		for(int bin=0;bin<Bins.size();bin++){
			fits.push_back(std::vector<std::vector<TFit *>>());
			for(int i=0;i<steps;i++){
				fits[bin].push_back(std::vector<TFit *>());
				for(int j=0;j<steps;j++){

					TFit * fit = new TFit;
					string named    =Basename + "_Data_" +to_string(bin)+" "+to_string(i)+" "+to_string(j);
					string namedprim=Basename + "_DataPrim_" +to_string(bin)+" "+to_string(i)+" "+to_string(j);
					string nameP    =Basename + "_MCP_"      +to_string(bin)+" "+to_string(i)+" "+to_string(j);
					string nameD    =Basename + "_MCD_"      +to_string(bin)+" "+to_string(i)+" "+to_string(j);
					string nameHe   =Basename + "_MCHe_"     +to_string(bin)+" "+to_string(i)+" "+to_string(j);


					fit->Templ_P =  new TH1F(nameP.c_str(),nameP.c_str(),Nbins,Xmin,Xmax);
					fit->Templ_D =  new TH1F(nameD.c_str(),nameD.c_str(),Nbins,Xmin,Xmax);
					fit->Templ_He=  new TH1F(nameHe.c_str(),nameHe.c_str(),Nbins,Xmin,Xmax);
					fit->Data    =  new TH1F(named.c_str(),named.c_str(),Nbins,Xmin,Xmax);
					fit->DataPrim=  new TH1F(namedprim.c_str(),namedprim.c_str(),Nbins,Xmin,Xmax);
					fits[bin][i].push_back(fit);
				}
			}
		}
	
		basename=Basename;
		cut = Cut;
		cutprimary=Cut+"&IsPrimary";
		bins=Bins;
		
		StatError  = new TH1F("StatError","StatError",bins.size(),0,bins.size()) ;
        	SystError  = new TH1F("SystError","SystError",bins.size(),0,bins.size()) ;

		ProtonCounts    = new TH1F("Proton Counts","Proton Counts",bins.size(),0,bins.size()) ;
        	DeuteronCounts  = new TH1F("Deuteron Counts","Deuteron Counts",bins.size(),0,bins.size()) ;

		ProtonCountsPrim    = new TH1F("Primary Proton Counts","Primary Proton Counts",bins.size(),0,bins.size()) ;
        	DeuteronCountsPrim  = new TH1F("Primary Deuteron Counts","Primary Deuteron Counts",bins.size(),0,bins.size()) ;
	
	
		BestFitSigma  = new TH1F("Best FIt Sigma","Best FIt Sigma",bins.size(),0,bins.size()) ;
        	BestFitShift  = new TH1F("Best Fit Shift","Best FIt Shift",bins.size(),0,bins.size()) ;

		BestChiSquares     = new TH1F("Best ChiSquare","Best ChiSquare",bins.size(),0,bins.size()) ;
        	OriginalChiSquares = new TH1F("Original ChiSquare","Original CHiSquare",bins.size(),0,bins.size()) ;

		Exposure_Time=(TH1F*)ExposureTime;

		isrich = IsRich;

		systpar.sigma=sigma;
		systpar.shift=shift;
		systpar.steps=steps;

		SetFitConstraints(0.0001,1,0.0001,1,0.000085,0.0012);
	}

	//reading constructor

	TemplateFIT(FileSaver  File, std::string Basename,Binning Bins, bool IsRich=false, int steps=11,float sigma=60,float shift=40,TH1F * ExposureTime=0x0){

		TFile * file = File.GetFile();

		for(int bin=0;bin<Bins.size();bin++){
			fits.push_back(std::vector<std::vector<TFit *>>());
			for(int i=0;i<steps;i++){
				fits[bin].push_back(std::vector<TFit *>());
				for(int j=0;j<steps;j++){

					TFit * fit = new TFit;
					string named    =Basename + "/Bin "+ to_string(bin)+"/Data/" + Basename + "_Data_" +to_string(bin)+" "+to_string(0)+" "+to_string(5);
					string namedprim=Basename + "/Bin "+ to_string(bin)+"/Data/" + Basename + "_DataPrim_" +to_string(bin)+" "+to_string(0)+" "+to_string(5);
					string nameP    =Basename + "/Bin "+ to_string(bin)+"/TemplateP/" + Basename + "_MCP_"      +to_string(bin)+" "+to_string(i)+" "+to_string(j);
					string nameD    =Basename + "/Bin "+ to_string(bin)+"/TemplateD/" + Basename + "_MCD_"      +to_string(bin)+" "+to_string(i)+" "+to_string(j);
					string nameHe    =Basename + "/Bin "+ to_string(bin)+"/TemplateHe/" + Basename + "_MCHe_"      +to_string(bin)+" "+to_string(i)+" "+to_string(j);


					fit->Templ_P =  (TH1F *)file->Get(nameP.c_str());
					fit->Templ_D =  (TH1F *)file->Get(nameD.c_str());
					fit->Templ_He=  (TH1F *)file->Get(nameHe.c_str());
					fit->Data    =  (TH1F *)file->Get(named.c_str());
					fit->DataPrim=  (TH1F *)file->Get(namedprim.c_str());

					fits[bin][i].push_back(fit);
				}
			}
		}	



		basename=Basename;

		bins=Bins; 

		StatError  = new TH1F("StatError","StatError",bins.size(),0,bins.size()) ;
		SystError  = new TH1F("SystError","SystError",bins.size(),0,bins.size()) ;

		ProtonCounts    = new TH1F("Proton Counts","Proton Counts",bins.size(),0,bins.size()) ;
		DeuteronCounts  = new TH1F("Deuteron Counts","Deuteron Counts",bins.size(),0,bins.size()) ;
	
		ProtonCountsPrim    = new TH1F("Primary Proton Counts","Primary Proton Counts",bins.size(),0,bins.size()) ;
        	DeuteronCountsPrim  = new TH1F("Primary Deuteron Counts","Primary Deuteron Counts",bins.size(),0,bins.size()) ;
	

		BestFitSigma  = new TH1F("Best Fit Sigma","Best Fit Sigma",bins.size(),0,bins.size()) ;
		BestFitShift  = new TH1F("Best Fit Shift","Best Fit Shift",bins.size(),0,bins.size()) ;

		BestChiSquares     = new TH1F("Best ChiSquare","Best ChiSquare",bins.size(),0,bins.size()) ;
        	OriginalChiSquares = new TH1F("Original ChiSquare","Original CHiSquare",bins.size(),0,bins.size()) ;

		isrich=IsRich;

		Exposure_Time=(TH1F*)ExposureTime;

		systpar.sigma=sigma;
		systpar.shift=shift;
		systpar.steps=steps;

		MCHeContRatio = (TH1F*) File.Get((basename+"/Fit Results/HeliumContamination/MC Expected He over P").c_str());
		if(!MCHeContRatio) MCHeContRatio = Eval_MCHeContRatio("MC Expected He over P");

		SetFitConstraints(0.0001,1,0.0001,1,0.000085,0.0012);

	}

	void SetFitConstraints(float minP, float maxP,float minD,float maxD,float minHe,float maxHe) {
		constrainmin[0]=minP; constrainmin[1]=minD; constrainmin[2]=minHe; 
		constrainmax[0]=maxP; constrainmax[1]=maxD; constrainmax[2]=maxHe;}

	void ReinitializeHistos(){}; //dummy
	void Eval_TransferFunction();
	TH1F * Eval_MCHeContRatio(std::string name);
	void Fill(TTree * treeMC,TTree * treeDT, Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars) );
	float SmearBeta(float Beta, float stepsigma, float stepshift,float R);
	float SmearBetaRICH(float Beta, float stepsigma, float stepshift);

	void FillEventByEventData(Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars));
	void FillEventByEventMC(Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars));


	void ExtractCounts(FileSaver finalhisto);
	void EvalFinalParameters();
	void EvalFinalErrors();
	void CalculateFinalPDCounts();
	void Save(FileSaver finalhisto,bool recreate=false);	
	void SaveFitResults(FileSaver finalhisto);
	void SumUpMassDistrib(FileSaver finalhisto);

	void SetUpBadEventSimulator(BadEventSimulator * Sim) {BadEvSim = Sim; return; };
	void SetFitWithNoiseMode(){IsFitNoise = true; if(BadEvSim)  BadEvSim->SetFrequency(1); return;}

	void SetSystematicParameters(int steps,float sigma,float shift){ systpar.steps=steps; systpar.shift=shift; systpar.sigma=sigma; return;};
	void SetFitRange(float min, float max){ 
		for(int bin=0;bin<bins.size();bin++) for(int i=0;i<systpar.steps;i++) for(int j=0;j<systpar.steps;j++) {
			fits[bin][i][j]->fitrangemin=min; 
			fits[bin][i][j]->fitrangemax=max; }
			return;
		}
	void SetFitRangeByQuantiles(float quant_min,float quant_max);
	void DisableFit(){fitDisabled=true;}
	void SetHeliumContamination(TF1 * HelimCont) {HeContModel=HelimCont; return;};
	BadEventSimulator * GetBadEventSimulator() {return BadEvSim;}
	void LoadEventIntoBadEvSim(Variables * vars) {if(BadEvSim) BadEvSim->LoadEvent(vars);}
	void Eval_ContError();
	std::string GetName(){return basename;}
	TH1F * GetStatError(){ return StatError;}
	TH1F * GetSystError() { return SystError;}	
	Binning  GetBinning() {return bins;}
	void  RebinAll(int f=2);	
	float GetHeContaminationWeight(int bin) { return fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->ContribHe; }
	float GetHeContaminationErr   (int bin) { return fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->errHe; }

	TH2F * GetDCountsSpread(int bin)	{ return DCountsSpread[bin];}
	TH2F * GetChiSquareSpread(int bin)      { return TFitChisquare[bin];}	
	TH1F * GetWeightedDCounts(int bin)     { return WeightedDCounts[bin];}

};

void TemplateFIT::Eval_TransferFunction(){
	for(int bin=0;bin<fits.size();bin++){
		TH1F * transferfunction = (TH1F *) fits[bin][0][5]->DataPrim->Clone();
		transferfunction->Sumw2();
		transferfunction->Divide(fits[bin][0][5]->Data);
		transferfunction->Smooth();
		TransferFunction.push_back(transferfunction);
	}
	return;
}

void TemplateFIT::Fill(TTree * treeMC,TTree * treeDT, Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars) ){
	
	cout<<basename.c_str()<<" Filling ... (Data)"<< endl;
	vars->ReadBranches(treeDT);

	for(int i=0;i<treeDT->GetEntries()/FRAC;i++){
		UpdateProgressBar(i, treeDT->GetEntries()/FRAC);
		treeDT->GetEvent(i);
		vars->Update();
		FillEventByEventData(vars,var,discr_var);
	}
	
	cout<<basename.c_str()<<" Filling ... (MC Protons)"<< endl;
	vars->ReadBranches(treeMC);

	for(int i=0;i<treeMC->GetEntries();i++){
		if(i%(int)FRAC!=0) continue;
		UpdateProgressBar(i, treeMC->GetEntries());
		treeMC->GetEvent(i);
		vars->Update();

		if(BadEvSim) BadEvSim->LoadEvent(vars);
		FillEventByEventMC(vars,var,discr_var);
	}

	return;
}


void TemplateFIT::FillEventByEventData(Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars)){

	int kbin;
	kbin = 	bins.GetBin(discr_var(vars));
	if(ApplyCuts(cut,vars)&&kbin>0){
		for(int i=0;i<systpar.steps;i++)
			for(int j=0;j<systpar.steps;j++){
				fits[kbin][i][j]->Data->Fill(var(vars),vars->PrescaleFactor);		
				//if(vars->R>0) cout<<"Prescaling: "<<vars->PrescaleFactor<<"; Beta: "<<vars->Beta<<"; R: "<<vars->R<<"; Mass "<<var(vars)<<endl;	
				if(ApplyCuts(cutprimary,vars)) fits[kbin][i][j]->DataPrim->Fill(var(vars),vars->PrescaleFactor);
				}
	}
	return;	
}

float TemplateFIT::SmearBetaRICH(float Beta, float stepsigma, float stepshift){
	float angle;
	if(Beta<0.87) angle= acos(1/(1.25*Beta))*10e4;
	else 	      angle= acos(1/(1.15*Beta))*10e4;
	float shiftstart=-systpar.shift;
	angle = angle + (shiftstart+(2*systpar.shift/(float)systpar.steps)*stepshift) + Rand->Gaus(0,(float)((2*systpar.sigma/systpar.steps)*stepsigma));
	if(Beta<0.87) return 1/(1.25*cos(angle/10e4));
	else	      return 1/(1.15*cos(angle/10e4));
}


float TemplateFIT::SmearBeta(float Beta, float stepsigma, float stepshift,float R){

	float time = 1.2/(Beta*3e-4);
	float shiftstart=-systpar.shift;

	float tailcontrolfactor=110./90.;
	if(R<2.7) tailcontrolfactor=1;//migration tail fixing

	float smeartime = (shiftstart+(2*systpar.shift/(float)systpar.steps)*stepshift) + Rand->Gaus(0,(float) tailcontrolfactor*((2*systpar.sigma/systpar.steps)*stepsigma));
	time = time + smeartime;
	return 1.2/(time*3e-4);

}



void TemplateFIT::FillEventByEventMC(Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars)){

	std::string cutP=cut+"&IsProtonMC";
	std::string cutD=cut+"&IsPureDMC";
	std::string cutHe=cut+"&IsFragmentedTMC";
	//std::string cutD=cut+"&IsDeutonMC";
	//std::string cutHe=cut+"&IsHeliumMC";


	cutHe.erase(cutHe.find("IsPreselected&"),14);
	cutHe = "IsPreselectedInner&" + cutHe;  //releasing cut for more stat. in Tritium templates

	if((ApplyCuts(cutP,vars)||ApplyCuts(cutD,vars)||ApplyCuts(cutHe,vars))){
		for(int i=0;i<systpar.steps;i++)
			for(int j=0;j<systpar.steps;j++){
				float betasmear;

				if(!ApplyCuts("IsOnlyFromToF",vars)) betasmear = SmearBetaRICH(vars->BetaRICH_new,(float)i,(float)j); 
				else 	betasmear = SmearBeta(vars->Beta,(float)i,(float)j,vars->R);    

				int kbin;
				if(!IsFitNoise){
				   if(bins.IsUsingBetaEdges()){	 
					kbin = bins.GetBin(betasmear);

					float mass = vars->R/betasmear * pow((1-pow(betasmear,2)),0.5);

					if(ApplyCuts(cutP,vars)&&kbin>0)  fits[kbin][i][j]->Templ_P->Fill(mass,vars->mcweight);		
					if(ApplyCuts(cutD,vars)&&kbin>0)  fits[kbin][i][j]->Templ_D->Fill(mass,vars->mcweight);
					if(ApplyCuts(cutHe,vars)&&kbin>0) fits[kbin][i][j]->Templ_He->Fill(mass,vars->mcweight);
					}
				    else {
				   	kbin = bins.GetBin(discr_var(vars));
					if(ApplyCuts(cutP,vars)&&kbin>0)  fits[kbin][i][j]->Templ_P->Fill(betasmear,vars->mcweight);
                                        if(ApplyCuts(cutD,vars)&&kbin>0)  fits[kbin][i][j]->Templ_D->Fill(betasmear,vars->mcweight);
                                        if(ApplyCuts(cutHe,vars)&&kbin>0) fits[kbin][i][j]->Templ_He->Fill(betasmear,vars->mcweight);

				   }	

				}
				else{
					if(bins.IsUsingBetaEdges()) {
						kbin = bins.GetBin(betasmear);
						float mass = vars->R/betasmear * pow((1-pow(betasmear,2)),0.5);		

						if(ApplyCuts(cutP,vars)&&kbin>0)  fits[kbin][i][j]->Templ_P->Fill(mass,vars->mcweight);		
						if(ApplyCuts(cutD,vars)&&kbin>0)  fits[kbin][i][j]->Templ_D->Fill(mass,vars->mcweight);

						float betabad = betasmear;
						if(BadEvSim) {betabad=BadEvSim->SimulateBadEvents(betasmear); 
							kbin = bins.GetBin(betabad);
						}
						float mass_bad = vars->R/betabad * pow((1-pow(betabad,2)),0.5);
						if(ApplyCuts(cutP,vars)&&kbin>0)  fits[kbin][i][j]->Templ_He->Fill(mass_bad,vars->mcweight);
					}					
					else {
						if(BadEvSim) betasmear=BadEvSim->SimulateBadEvents(betasmear);
						kbin = bins.GetBin(discr_var(vars));
						if(ApplyCuts(cutP,vars)&&kbin>0)  fits[kbin][i][j]->Templ_P->Fill(betasmear,vars->mcweight);
	                                        if(ApplyCuts(cutD,vars)&&kbin>0)  fits[kbin][i][j]->Templ_D->Fill(betasmear,vars->mcweight);
        	                                if(ApplyCuts(cutHe,vars)&&kbin>0) fits[kbin][i][j]->Templ_He->Fill(betasmear,vars->mcweight);
					}
				}
			}
	}
	return;	
}



void TemplateFIT::Save(FileSaver finalhisto,bool recreate){

	for(int bin=0;bin<bins.size();bin++){ 
		finalhisto.Add(fits[bin][0][5]->Data);
		finalhisto.Add(fits[bin][0][5]->DataPrim);

		finalhisto.writeObjsInFolder((basename + "/Bin "+ to_string(bin)+"/Data").c_str(),recreate);
	}

	for(int bin=0;bin<bins.size();bin++){ 
		for(int i=0;i<systpar.steps;i++)
                        for(int j=0;j<systpar.steps;j++){

				finalhisto.Add(fits[bin][i][j]->Templ_P);
		}
		finalhisto.writeObjsInFolder((basename + "/Bin "+ to_string(bin)+"/TemplateP").c_str(),recreate);
	}

	for(int bin=0;bin<bins.size();bin++){ 
		for(int i=0;i<systpar.steps;i++)
                        for(int j=0;j<systpar.steps;j++){

				finalhisto.Add(fits[bin][i][j]->Templ_D);
		}
		finalhisto.writeObjsInFolder((basename + "/Bin "+ to_string(bin)+"/TemplateD").c_str(),recreate);
	}

	for(int bin=0;bin<bins.size();bin++){ 
		for(int i=0;i<systpar.steps;i++)
                        for(int j=0;j<systpar.steps;j++){

				finalhisto.Add(fits[bin][i][j]->Templ_He);
		}
		finalhisto.writeObjsInFolder((basename + "/Bin "+ to_string(bin)+"/TemplateHe").c_str(),recreate);
	}

	return;
}

void TemplateFIT::SaveFitResults(FileSaver finalhisto){


	for(int bin=0;bin<bins.size();bin++){
	TH1F * OriginalP=(TH1F*)fits[bin][0][5]->Templ_P->Clone();
	TH1F * BestP=(TH1F*)fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->Templ_P->Clone();
	OriginalP->SetName("Original Proton MC ");
	BestP->SetName("Best #chi^{2} Mod. Proton MC ");
	finalhisto.Add(OriginalP);
	finalhisto.Add(BestP);
	finalhisto.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesP/Bin"+to_string(bin)).c_str());
	
	TH1F * OriginalD=(TH1F*)fits[bin][0][5]->Templ_D->Clone();
	TH1F * BestD=(TH1F*)fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->Templ_D->Clone();
	OriginalD->SetName("Original Deuton MC ");
	BestD->SetName("Best #chi^{2} Mod. Deuton MC ");
	finalhisto.Add(OriginalD);
	finalhisto.Add(BestD);
	finalhisto.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesD/Bin"+to_string(bin)).c_str());

	TH1F * OriginalHe=(TH1F*)fits[bin][0][5]->Templ_He->Clone();
	TH1F * BestHe=(TH1F*)fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->Templ_He->Clone();
	OriginalHe->SetName("Original Helium Fragm. MC ");
	BestHe->SetName("Best #chi^{2} Helium Fragm. MC ");
	finalhisto.Add(OriginalHe);
	finalhisto.Add(BestHe);
	finalhisto.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesHe/Bin"+to_string(bin)).c_str());
	
	}

	for(int bin=0;bin<bins.size();bin++){ 
		finalhisto.Add(fits[bin][0][5]->Data);
		finalhisto.Add(fits[bin][0][5]->DataPrim);
	finalhisto.writeObjsInFolder((basename+"/Fit Results/Data/Bin"+to_string(bin)).c_str());	
	}
	
	for(int bin=0;bin<bins.size();bin++){
                for(int i=0;i<systpar.steps;i++)
                        for(int j=0;j<systpar.steps;j++){
				if(!(i==0&&i==5))finalhisto.Add(fits[bin][i][j]->Templ_P);
			}	
	finalhisto.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesP/Bin"+to_string(bin)).c_str());
	}

	for(int bin=0;bin<bins.size();bin++){
                for(int i=0;i<systpar.steps;i++)
                        for(int j=0;j<systpar.steps;j++){
				if(!(i==0&&i==5))finalhisto.Add(fits[bin][i][j]->Templ_D);
				}
	finalhisto.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesD/Bin"+to_string(bin)).c_str());
	}
	
	for(int bin=0;bin<bins.size();bin++){
                for(int i=0;i<systpar.steps;i++)
                        for(int j=0;j<systpar.steps;j++){
				if(!(i==0&&i==5))finalhisto.Add(fits[bin][i][j]->Templ_He);
				}
	finalhisto.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesHe/Bin"+to_string(bin)).c_str());
	}

	for(int bin=0;bin<bins.size();bin++){
		if(BestChiSquare[bin]&&fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->Tfit){
			TH1F * FIT;
			if(fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->Tfit_outcome!=-1){ 
				FIT=(TH1F*)fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->Tfit-> GetPlot();
				FIT->SetName(("Fraction Fit bin" + to_string(bin)).c_str());
				finalhisto.Add(FIT);
				finalhisto.writeObjsInFolder((basename+"/Fit Results/FractionFits/Bin"+to_string(bin)).c_str());
			}
		}
	}


	for(int i=0;i<bins.size();i++) 
                finalhisto.Add(TransferFunction[i]);
        finalhisto.writeObjsInFolder((basename + "/Fit Results/TrasnferFunctions/").c_str());

	for(int i=0;i<bins.size();i++){
		finalhisto.Add(DCountsSpread[i]);
		finalhisto.writeObjsInFolder((basename+"/Fit Results/Spreads/DCounts").c_str());
	}
	for(int i=0;i<bins.size();i++){
		finalhisto.Add(TFitChisquare[i]);
		finalhisto.writeObjsInFolder((basename+"/Fit Results/Spreads/ChiSquare").c_str());	
	}
	for(int i=0;i<bins.size();i++){
		finalhisto.Add(WeightedDCounts[i]);
		finalhisto.writeObjsInFolder((basename+"/Fit Results/Spreads/Weighted D counts").c_str());	
	}

	if(HeContModel) finalhisto.Add(HeContModel);
	if(HeContError) finalhisto.Add(HeContError);
	if(MeasuredHeContRatio) finalhisto.Add(MeasuredHeContRatio);
	if(MCHeContRatio) finalhisto.Add(MCHeContRatio);
	finalhisto.writeObjsInFolder((basename+"/Fit Results/HeliumContamination").c_str());	

	finalhisto.Add(BestFitSigma);
	finalhisto.Add(BestFitShift);
	finalhisto.Add(StatError);
	finalhisto.Add(SystError);
	finalhisto.Add(ProtonCounts);
	finalhisto.Add(DeuteronCounts);
	finalhisto.Add(ProtonCountsPrim);
	finalhisto.Add(DeuteronCountsPrim);

	finalhisto.Add(BestChiSquares   ); 
	finalhisto.Add(OriginalChiSquares);

	finalhisto.writeObjsInFolder((basename+"/Fit Results/").c_str());

}


void TemplateFIT::SumUpMassDistrib(FileSaver finalhisto){

	TH1F * SummedMass = (TH1F*)fits[0][0][5]->Data ->Clone(); 
	for(int bin=0; bin<bins.size()-2;bin++){
		SummedMass->Add(fits[bin][0][5]->Data);
	}
	finalhisto.Add(SummedMass);
	finalhisto.writeObjsInFolder((basename + "/SummedData").c_str());
	return;
}

void TemplateFIT::SetFitRangeByQuantiles(float quant_min,float quant_max){
	for(int bin=0;bin<bins.size();bin++) 
		for(int i=0;i<systpar.steps;i++) 
			for(int j=0;j<systpar.steps;j++){
				Double_t q[2]= {quant_min,quant_max};
				Double_t x[2];
				fits[bin][i][j]->Data->GetQuantiles(2,x,q);
				fits[bin][i][j]->fitrangemin = x[0];
				fits[bin][i][j]->fitrangemax = x[1];
				cout<<"bin: "<<bin<<": "<<fits[bin][i][j]->fitrangemin<<" "<<fits[bin][i][j]->fitrangemax<<endl;
			}
	return;
}


void Pre_Scale(TH1F * PHisto,TH1F * DHisto, TH1F * HeHisto){
	if(DHisto>0&&PHisto>0) DHisto->Scale(0.02*PHisto->Integral());
	if(HeHisto>0&&PHisto>0) HeHisto->Scale(0.1*PHisto->Integral());
	return;
}


void Do_TemplateFIT(TFit * Fit,float fitrangemin,float fitrangemax,float constrain_min[], float constrain_max[], TF1 * HeCont=0x0, float bincenter=1){

	Pre_Scale(Fit ->  Templ_P,Fit ->  Templ_D,Fit ->  Templ_He);

	TObjArray *Tpl;
	Tpl = new TObjArray(2);
	if(Fit ->  Templ_P)  if(Fit ->  Templ_P ->GetEntries()>0) Tpl -> Add( Fit ->  Templ_P );
	if(Fit ->  Templ_D)  if(Fit ->  Templ_D ->GetEntries()>0) Tpl -> Add( Fit ->  Templ_D );
	if(Fit ->  Templ_He) if(Fit ->  Templ_He->GetEntries()>0) Tpl -> Add( Fit ->  Templ_He);

	float min=fitrangemin;
	float max=fitrangemax;

	bool fitcondition = (Fit -> Data->Integral()>500)&&(Fit -> Templ_P->Integral()>50) &&(Fit -> Templ_D->Integral()>50);
		
	
	if(fitcondition) { 
		Fit -> Tfit = new TFractionFitter(Fit -> Data, Tpl ,"q");
		Fit -> Tfit -> SetRangeX(Fit -> Data -> FindBin(min), Fit -> Data -> FindBin(max));
		
		Fit -> Tfit -> Constrain(1, constrain_min[0] ,constrain_max[0]);
                Fit -> Tfit -> Constrain(2, constrain_min[1] ,constrain_max[1]);
		if(Fit ->  Templ_He) {
			if(Fit ->  Templ_He -> Integral()>50)
				if(HeCont) Fit -> Tfit -> Constrain(3,0.4*HeCont->Eval(bincenter),1.6*HeCont->Eval(bincenter));
				else Fit -> Tfit -> Constrain(3, constrain_min[2] ,constrain_max[2]);	 
			else Fit -> Tfit -> Constrain(3, 0.000001 ,0.0000011);
		}
		if(Fit -> Tfit ) Fit -> Tfit_outcome = Fit -> Tfit -> Fit();

		for(int fit_attempt=0; fit_attempt<20; fit_attempt++) {
			cout<<fit_attempt<<endl;
			if(Fit -> Tfit_outcome == 0) break;
			else {
				cout<<fit_attempt<<endl;
				Fit -> Tfit_outcome = Fit -> Tfit -> Fit();
			}
		}

		if(Fit -> Tfit_outcome==0){
			TH1F * Result = (TH1F *) Fit-> Tfit -> GetPlot();
			float itot= Result->Integral();
			double w1,e1 = 0;
			double w2,e2 = 0;
			double w3,e3 = 0;
		
			Fit -> Tfit ->GetResult(0,w1,e1);
			Fit -> Tfit ->GetResult(1,w2,e2);
			if(Fit-> Templ_He) Fit -> Tfit ->GetResult(2,w3,e3);
		
			float i1 = Fit-> Templ_P  ->Integral(Fit->Templ_P -> FindBin(min), Fit->Templ_P -> FindBin(max));
			float i2 = Fit-> Templ_D  ->Integral(Fit->Templ_D -> FindBin(min), Fit->Templ_D -> FindBin(max));
			float i3=1;
			if(Fit-> Templ_He) i3 = Fit-> Templ_He ->Integral(Fit->Templ_He -> FindBin(min), Fit->Templ_He -> FindBin(max));

			Fit ->ContribP= w1;
			Fit ->ContribD= w2;
			Fit ->ContribHe=w3;
			
			Fit ->errP= e1;
			Fit ->errD= e2;
			Fit ->errHe=e3;

			Fit ->wheightP= w1*itot/i1;
			Fit ->wheightD= w2*itot/i2;
			Fit->wheightHe= w3*itot/i3;

			if(HeCont) cout<<"He expected: "<<HeCont->Eval(bincenter)<<": ";
			cout<<w1<<" "<<w2<<" "<<w3<<endl;

			Fit ->  Templ_P  -> Scale(Fit ->wheightP);
			Fit ->  Templ_D  -> Scale(Fit ->wheightD);
			if(Fit-> Templ_He) Fit ->  Templ_He -> Scale(Fit ->wheightHe);

			float Cov01=0;

			Cov01= Fit -> Tfit->GetFitter()->GetCovarianceMatrixElement(0,1);
			Fit -> StatErr = pow(pow(w2*e2,2)+pow(w1*e1,2)-2*Cov01*w1*w2,0.5);

			Fit -> ChiSquare = Fit -> Tfit -> GetChisquare()/(float) (Fit ->  Tfit -> GetNDF());
			Fit -> DCounts = Fit ->  Templ_D -> Integral();
			Fit -> PCounts = Fit ->  Templ_P -> Integral();
		}
		else{
			Fit ->wheightP= 0;
			Fit ->wheightD= 0;
		}
	}
	return;
}

float EvalFitProbability(float chi){
	TF1 * CumulativeChiSquare=new TF1("Chi","exp(-x)*x^-0.5",0.001,10);
	return CumulativeChiSquare->Integral(chi,100);
}

void  TemplateFIT::RebinAll(int f){
	for(int bin=0;bin<bins.size();bin++){
	 for(int sigma=0;sigma<systpar.steps;sigma++)
		 for(int shift=0;shift<systpar.steps;shift++){
			cout<<bin<<" "<<sigma<<" "<<shift<<" "<< fits[bin][sigma][shift] -> DataPrim->GetNbinsX()<<endl;
			 fits[bin][sigma][shift] ->  Templ_P->Rebin(f);
			 fits[bin][sigma][shift] ->  Templ_D->Rebin(f);
			 fits[bin][sigma][shift] ->  Templ_He->Rebin(f);
			 fits[bin][sigma][shift] -> Data->Rebin(f);
			 fits[bin][sigma][shift] -> DataPrim->Rebin(f);			
		 }
	}
};

void TemplateFIT::ExtractCounts(FileSaver finalhisto){

	Eval_TransferFunction();

	for(int bin=0;bin<bins.size();bin++){

		// Template Fits
		for(int sigma=0;sigma<systpar.steps;sigma++){
			for(int shift=0;shift<systpar.steps;shift++){
				cout<<fits[bin][sigma][shift]<<endl;

				//TEMPORARY
			/*		fits[bin][sigma][shift]->  Templ_P->Rebin(2);
					fits[bin][sigma][shift] ->  Templ_D->Rebin(2);
					fits[bin][sigma][shift] ->  Templ_He->Rebin(2);
					fits[bin][sigma][shift] -> Data->Rebin(2);
					fits[bin][sigma][shift] -> DataPrim->Rebin(2);
			*/		
				if(!fitDisabled) {
					Do_TemplateFIT(fits[bin][sigma][shift],fits[bin][sigma][shift]->fitrangemin,fits[bin][sigma][shift]->fitrangemax,
					constrainmin,constrainmax,HeContModel,bins.BetaBinCent(bin));
					cout<<bin<<endl;
				}
				fits[bin][sigma][shift]->Templ_DPrim=(TH1F*) fits[bin][sigma][shift]->Templ_D->Clone();
				fits[bin][sigma][shift]->Templ_DPrim->Multiply(TransferFunction[bin]);
				if(!fitDisabled) fits[bin][sigma][shift]->DCountsPrim = fits[bin][sigma][shift]->Templ_DPrim->Integral();
				else { 
					fits[bin][sigma][shift]->DCountsPrim = fits[bin][sigma][shift]->DataPrim->Integral();
					fits[bin][sigma][shift]->DCounts     = fits[bin][sigma][shift]->DataPrim->Integral();
				}
				fits[bin][sigma][shift]->Templ_PPrim=(TH1F*) fits[bin][sigma][shift]->Templ_P->Clone();
				fits[bin][sigma][shift]->Templ_PPrim->Multiply(TransferFunction[bin]);
				if(!fitDisabled) fits[bin][sigma][shift]->PCountsPrim = fits[bin][sigma][shift]->Templ_PPrim->Integral();
				else { 
					fits[bin][sigma][shift]->PCountsPrim = fits[bin][sigma][shift]->DataPrim->Integral();
				        fits[bin][sigma][shift]->PCounts     = fits[bin][sigma][shift]->DataPrim->Integral();
				     }
			}
		}

		// Histograms for systematic error evaluation
		TH2F * dcountsspread = new TH2F(("DCountsSpread Bin " +to_string(bin)).c_str(),("DCountsSpread Bin " +to_string(bin)).c_str(),systpar.steps,0,2*systpar.sigma,systpar.steps,-systpar.shift,systpar.shift);
		TH2F * tfitchisquare = new TH2F(("ChiSquare Bin " +to_string(bin)).c_str(),("ChiSquare Bin " +to_string(bin)).c_str(),systpar.steps,0,2*systpar.sigma,systpar.steps,-systpar.shift,systpar.shift);
		TH1F * weighteddcounts  = new TH1F(("Weighted Counts Bin " +to_string(bin)).c_str(),("Weighted Counts Bin " +to_string(bin)).c_str(),35,0.4*fits[bin][0][4]->DCounts,1.7*fits[bin][0][4]->DCounts);
		BestChi * MinimumChi = new BestChi();


		for(int sigma=0;sigma<systpar.steps;sigma++){
			for(int shift=0;shift<systpar.steps;shift++){
				dcountsspread->SetBinContent(sigma+1,shift+1,fits[bin][sigma][shift]->DCounts);
				if(fits[bin][sigma][shift]->ChiSquare>0)
					tfitchisquare->SetBinContent(sigma+1,shift+1,fits[bin][sigma][shift]->ChiSquare);
				else  tfitchisquare->SetBinContent(sigma+1,shift+1,500);

				if(	fits[bin][sigma][shift]->DCounts>0.1*fits[bin][0][4]->DCounts &&
						fits[bin][sigma][shift]->DCounts<2*fits[bin][0][4]->DCounts  )

					weighteddcounts -> Fill(fits[bin][sigma][shift]->DCounts,EvalFitProbability(fits[bin][sigma][shift]->ChiSquare));
			}
		}

		MinimumChi->FindMinimum(tfitchisquare);	


		BestChiSquare.push_back(MinimumChi);			
		DCountsSpread.push_back(dcountsspread);
		TFitChisquare.push_back(tfitchisquare);
		WeightedDCounts.push_back(weighteddcounts);


	}

	EvalFinalParameters();
	EvalFinalErrors();
	CalculateFinalPDCounts();	
	return;
}

void TemplateFIT::EvalFinalParameters(){
	for(int bin=0;bin<bins.size();bin++){
		float BinSTD=0;
		for(int i=0;i<systpar.steps;i++){
			for(int j=0;j<systpar.steps;j++){
				BinSTD+=pow(TFitChisquare[bin]->GetBinContent(i+1,j+1)-BestChiSquare[bin]->chimin,2);
			}
		}
		BinSTD=pow(BinSTD,0.5)/(systpar.steps*systpar.steps);

		BestFitSigma->SetBinContent(bin+1,(float)((2*systpar.sigma/systpar.steps)*BestChiSquare[bin]->i));
		BestFitShift->SetBinContent(bin+1,-systpar.shift+(2*systpar.shift/(float)systpar.steps)*BestChiSquare[bin]->j);

		
	
		if(BinSTD>0){
		BestFitSigma->SetBinError(bin+1,pow(pow(4/BinSTD,2)+pow(1.5*(2*systpar.sigma/(float)systpar.steps),2),0.5));
		BestFitShift->SetBinError(bin+1,pow(pow(4/BinSTD,2)+pow(1.5*(2*systpar.shift/(float)systpar.steps),2),0.5));
		}

		if(BestChiSquare[bin]->chimin<1000) 	       BestChiSquares     ->SetBinContent(bin+1,BestChiSquare[bin]->chimin);
                if(TFitChisquare[bin]->GetBinContent(1,6)<1000)OriginalChiSquares ->SetBinContent(bin+1,TFitChisquare[bin]->GetBinContent(1,6));
		BestChiSquares     ->SetBinError(bin+1,0.25);	
                OriginalChiSquares ->SetBinError(bin+1,0.25);
	}
	return;
}

void TemplateFIT::EvalFinalErrors(){

	for(int bin=0;bin<bins.size();bin++){

		if(fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->DCounts>0){
			StatError -> SetBinContent(bin+1,fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->StatErr);
			SystError -> SetBinContent(bin+1,WeightedDCounts[bin]->GetStdDev()/fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->DCounts);
			StatError -> SetBinError(bin+1,0);
			SystError -> SetBinError(bin+1,0);
		}

	}
	return;

}


void TemplateFIT::CalculateFinalPDCounts(){

		HeContError = new TH1F ("HeContError","HeContError",bins.size(),0,bins.size());	
		MeasuredHeContRatio = new TH1F("Measured He over P","Measured He over P",bins.size(),0,bins.size());
		Eval_ContError();
		StatError -> Smooth(1);
		SystError -> Smooth(1);
		
	for(int bin=0;bin<bins.size();bin++){

		float toterr= pow(pow(StatError -> GetBinContent(bin+1),2) + pow(SystError -> GetBinContent(bin+1),2) + pow(HeContError -> GetBinContent(bin+1),2),0.5);
		
		float CountsP      = fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->PCounts;
		float CountsD      = fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->DCounts;
		float CountsP_prim = fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->PCountsPrim;
		float CountsD_prim = fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->DCountsPrim;
		
		float ErrCountsP      = toterr*CountsP;
		float ErrCountsD      = toterr*CountsD;
		float ErrCountsP_prim = toterr*CountsP_prim;
		float ErrCountsD_prim = toterr*CountsP_prim;
		
		if(!IsFitNoise){
			CountsP += fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->Templ_He->Integral();
			CountsD += (CountsD/CountsP)*fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->Templ_He->Integral();
			ErrCountsP = toterr*CountsP;
			ErrCountsP_prim = toterr*CountsP_prim;
			toterr = pow(toterr,2)	- pow(HeContError -> GetBinContent(bin+1),2);
			toterr += pow(0.1*fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->Templ_He->Integral(37,62)/CountsD,2);
			toterr = pow(toterr,0.5);
			ErrCountsD = toterr*CountsD;
			ErrCountsD_prim = toterr*CountsP_prim;		
		}

		ProtonCounts->SetBinContent(bin+1,fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->PCounts);
		ProtonCounts->SetBinError(bin+1,toterr*fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->PCounts);
		
		DeuteronCounts->SetBinContent(bin+1,fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->DCounts);
                DeuteronCounts->SetBinError(bin+1,toterr*fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->DCounts);		
		
		ProtonCountsPrim->SetBinContent(bin+1,fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->PCountsPrim);
		ProtonCountsPrim->SetBinError(bin+1,toterr*fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->PCountsPrim);

		DeuteronCountsPrim->SetBinContent(bin+1,fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->DCountsPrim);
                DeuteronCountsPrim->SetBinError(bin+1,toterr*fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->DCountsPrim);		
	}
		

}

TH1F * TemplateFIT::Eval_MCHeContRatio(std::string name){
	TH1F * HeCountRatio = new TH1F(name.c_str(),name.c_str(),bins.size(),0,bins.size());
	for(int bin=0;bin<bins.size();bin++){
                TH1F * HeliumCounts = (TH1F*) fits[bin][0][5]->Templ_He->Clone();
                TH1F * ProtonCounts = (TH1F*) fits[bin][0][5]->Templ_P->Clone();	
		if(ProtonCounts->Integral()>0){
			HeCountRatio->SetBinContent(bin+1,HeliumCounts->Integral()/ProtonCounts->Integral());
			HeCountRatio->SetBinError(bin+1,pow(HeliumCounts->Integral(),0.5)/ProtonCounts->Integral());
		}
	}
	return HeCountRatio;
}

void TemplateFIT::Eval_ContError(){

	for(int bin=0;bin<bins.size();bin++){
		TH1F * HeliumCounts = (TH1F*) fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->Templ_He->Clone();
		TH1F * DeuteronCounts = (TH1F*) fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->Templ_D->Clone();
		TH1F * ProtonCounts = (TH1F*) fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->Templ_P->Clone();
		if(DeuteronCounts->Integral()>0)
		HeContError->SetBinContent(bin+1,0.2*HeliumCounts->Integral(HeliumCounts->FindBin(1.5),HeliumCounts->GetNbinsX())/DeuteronCounts->Integral());
		if(ProtonCounts->Integral()>0)
			MeasuredHeContRatio->SetBinContent(bin+1,HeliumCounts->Integral()/ProtonCounts->Integral());
			float toterr= pow(pow(StatError -> GetBinContent(bin+1),2) + pow(SystError -> GetBinContent(bin+1),2) ,0.5);
			MeasuredHeContRatio->SetBinError(bin+1,toterr*HeliumCounts->Integral()/ProtonCounts->Integral());
	}
	return;
}
#endif
