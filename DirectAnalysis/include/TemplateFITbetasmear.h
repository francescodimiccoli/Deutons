#ifndef TEMPLATEFITBETASMEAR_H
#define TEMPLATEFITBETASMEAR_H

#include "BadEventSimulator.h"
#include "LatReweighter.h" 
#include "Livetime.h"
#include "TFractionFitter.h"
#include <TGraphErrors.h>
#include "Tool.h"

class Tool;

using namespace std;

// not used
/*TSpline3 * ExtractCutoffWeight(TH1F * ExposureTime){

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

}*/

struct TFit {
   TH1F * Templ_P =0x0;
   TH1F * Templ_D =0x0;
   TH1F * Templ_DPrim=0x0;
   TH1F * Templ_PPrim=0x0;
   TH1F * Templ_He=0x0;
   TH1F * Templ_Noise=0x0; 	

   TH1F * Data =0x0;
   TH1F * DataPrim=0x0;	


   float wheightP,wheightD,wheightHe,wheightNoise;
   float ContribP,ContribD,ContribHe,ContribNoise;
   float errP,errD,errHe,errNoise;
   float fitrangemin = 0.6;
   float fitrangemax = 3.5;
		
   TFractionFitter *Tfit;
   int Tfit_outcome=-1;
   
   float DCounts=0;
   float PCounts=0;
   float TCounts=0;


   float ChiSquare=500;			
   
   float StatErrP=0;
   float StatErrD=0;
   float StatErrT=0;
	

  
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
	void FindMinimum(TH2F * Histo, TH2F * relerr) {
		float meanErr = relerr->Integral()/(Histo->GetNbinsY()*Histo->GetNbinsX());
		float Best = 9999999;
		for(int x=0;x<Histo->GetNbinsX();x++)
			for(int y=0;y<Histo->GetNbinsY();y++)
				if(Histo->GetBinContent(x+1,y+1)<Best){
					if(relerr->GetBinContent(x+1,y+1)<2*meanErr){
						Best=Histo->GetBinContent(x+1,y+1);
						i=x;
						j=y;
						chimin=Best;	
					}
				}
	}
};

void Do_TemplateFIT(TFit * Fit,float fitrangemin,float fitrangemax,float constrain_min[], float constrain_max[], bool isfitnoise, bool highmasstailconstrain, bool IsFitPrim=false );

TSpline3 * GetSplineFromHisto(TH1F * Graph, Binning bins);

class TemplateFIT : public Tool{

	private:
	std::vector<std::vector<std::vector<TFit *>>> fits;
	std::vector<BestChi *> BestChiSquare;	
	std::vector<TH1F *> TransferFunction;


	std::vector<TH2F *> DCountsSpread;
	std::vector<TH2F *> DErrSpread;
	std::vector<TH1F *> WeightedDCounts;
	std::vector<TH2F *> TFitChisquare;
	TF1 * HeContModel=0x0;
	TH1F * MCHeContRatio;
	TH1F * MeasuredHeContRatio;

	TH1F * StatErrorP;
	TH1F * StatErrorD;
	TH1F * StatErrorT;
	
	TH1F * SystError;
	TH1F * HeContError;

	TH1F * ProtonCounts;
	TH1F * DeuteronCounts;
	TH1F * TritiumCounts;


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

	TSpline3 * ExposureTime;
	LatReweighter * Latweighter = new LatReweighter("LatWeights","",500,0,150);

	Systpar systpar;
	bool fitDisabled=false;
	bool isrich=false;
	BadEventSimulator * BadEvSim=0x0;
	bool IsFitNoise = false;
	float constrainmin[3];
	float constrainmax[3];
	bool highmassconstrain=false;
	int ActualTime=0;


	FileSaver ExternalTemplates;
	bool checkfiletemplates=false;

	public:	
	//standard constructor
	TemplateFIT(std::string Basename,Binning Bins, std::string Cut, int Nbins, float Xmin, float Xmax, bool IsRich=false ,int steps=11,float sigma=60,float shift=60){
		
			ExternalTemplates.setName("AnalysisFiles/ExternalTemplates.root");		
			checkfiletemplates = ExternalTemplates.CheckFile();

			if(checkfiletemplates){ 
				cout<<"External Template file found: "<<endl; 
			}
	
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
					string nameNo   =Basename + "_MCNoise_"  +to_string(bin)+" "+to_string(i)+" "+to_string(j);



					fit->Templ_P =  new TH1F(nameP.c_str(),nameP.c_str(),Nbins,Xmin,Xmax);
					fit->Templ_D =  new TH1F(nameD.c_str(),nameD.c_str(),Nbins,Xmin,Xmax);
					fit->Templ_He=  new TH1F(nameHe.c_str(),nameHe.c_str(),Nbins,Xmin,Xmax);
					fit->Templ_Noise=  new TH1F(nameNo.c_str(),nameNo.c_str(),Nbins,Xmin,Xmax);
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
		
		StatErrorP  = new TH1F("StatErrorP","StatErrorP",bins.size(),0,bins.size()) ;
		StatErrorD  = new TH1F("StatErrorD","StatErrorD",bins.size(),0,bins.size()) ;
		StatErrorT  = new TH1F("StatErrorT","StatErrorT",bins.size(),0,bins.size()) ;

        	SystError  = new TH1F("SystError","SystError",bins.size(),0,bins.size()) ;

		ProtonCounts    = new TH1F("Proton Counts","Proton Counts",bins.size(),0,bins.size()) ;
        	DeuteronCounts  = new TH1F("Deuteron Counts","Deuteron Counts",bins.size(),0,bins.size()) ;
		TritiumCounts  = new TH1F("Tritium Counts","Tritium Counts",bins.size(),0,bins.size()) ;

		TFile * f = TFile::Open("LatWeights/ExposureModel.root");		
		TH1F * Exp = (TH1F*) f->Get("HEExposure");	
	        Exp->Scale(1/Exp->GetBinContent(Exp->GetMaximumBin()));
		ExposureTime = GetSplineFromHisto(Exp,PRB);

		ProtonCountsPrim    = new TH1F("Primary Proton Counts","Primary Proton Counts",bins.size(),0,bins.size()) ;
        	DeuteronCountsPrim  = new TH1F("Primary Deuteron Counts","Primary Deuteron Counts",bins.size(),0,bins.size()) ;
	
	
		BestFitSigma  = new TH1F("Best FIt Sigma","Best FIt Sigma",bins.size(),0,bins.size()) ;
        	BestFitShift  = new TH1F("Best Fit Shift","Best FIt Shift",bins.size(),0,bins.size()) ;

		BestChiSquares     = new TH1F("Best ChiSquare","Best ChiSquare",bins.size(),0,bins.size()) ;
        	OriginalChiSquares = new TH1F("Original ChiSquare","Original CHiSquare",bins.size(),0,bins.size()) ;


		isrich = IsRich;

		systpar.sigma=sigma;
		systpar.shift=shift;
		systpar.steps=steps;

		 SetFitConstraints(0.9,1,0.015,0.06,0.0001,0.005);
	}

	//reading constructor

	TemplateFIT(FileSaver  File, std::string Basename,Binning Bins, bool IsRich=false, int steps=11,float sigma=60,float shift=60){

		TFile * file = File.GetFile();
		TFile * externalfile;

		ExternalTemplates.setName("/afs/cern.ch/work/f/fdimicco/private/Deutons/DirectAnalysis/AnalysisFiles/ExternalTemplates.root");		

		checkfiletemplates = ExternalTemplates.CheckFile();
				
		if(checkfiletemplates){ externalfile=ExternalTemplates.GetFile();
			cout<<"External Template file found: "<<externalfile<<endl; 
		}
			
		for(int bin=0;bin<Bins.size();bin++){
			fits.push_back(std::vector<std::vector<TFit *>>());
			for(int i=0;i<steps;i++){
				fits[bin].push_back(std::vector<TFit *>());
				for(int j=0;j<steps;j++){

					TFit * fit = new TFit;
					string named    =Basename + "/Bin "+ to_string(bin)+"/Data/" + Basename + "_Data_" +to_string(bin)+" "+to_string(0)+" "+to_string(3);
					string namedprim=Basename + "/Bin "+ to_string(bin)+"/Data/" + Basename + "_DataPrim_" +to_string(bin)+" "+to_string(0)+" "+to_string(3);
					string nameP    =Basename + "/Bin "+ to_string(bin)+"/TemplateP/" + Basename + "_MCP_"      +to_string(bin)+" "+to_string(i)+" "+to_string(j);
					string nameD    =Basename + "/Bin "+ to_string(bin)+"/TemplateD/" + Basename + "_MCD_"      +to_string(bin)+" "+to_string(i)+" "+to_string(j);
					string nameHe    =Basename + "/Bin "+ to_string(bin)+"/TemplateHe/" + Basename + "_MCHe_"      +to_string(bin)+" "+to_string(i)+" "+to_string(j);
					string nameNo    =Basename + "/Bin "+ to_string(bin)+"/TemplateNoise/" + Basename + "_MCNoise_"      +to_string(bin)+" "+to_string(i)+" "+to_string(j);

					if(checkfiletemplates) {
						fit->Templ_P    =  (TH1F *)externalfile->Get(nameP.c_str());
						fit->Templ_D 	=  (TH1F *)externalfile->Get(nameD.c_str());
						fit->Templ_He	=  (TH1F *)externalfile->Get(nameHe.c_str());
						fit->Templ_Noise=  (TH1F *)externalfile->Get(nameNo.c_str());
					}

					else {
						fit->Templ_P    =  (TH1F *)file->Get(nameP.c_str());
						fit->Templ_D 	=  (TH1F *)file->Get(nameD.c_str());
						fit->Templ_He	=  (TH1F *)file->Get(nameHe.c_str());
						fit->Templ_Noise=  (TH1F *)file->Get(nameNo.c_str());
					}

					fit->Data    	=  (TH1F *)file->Get(named.c_str());
					fit->DataPrim	=  (TH1F *)file->Get(namedprim.c_str());


					fits[bin][i].push_back(fit);
				}
			}
		}	



		basename=Basename;

		bins=Bins; 

		StatErrorP  = new TH1F("StatErrorP","StatErrorP",bins.size(),0,bins.size()) ;
		StatErrorD  = new TH1F("StatErrorD","StatErrorD",bins.size(),0,bins.size()) ;
		StatErrorT  = new TH1F("StatErrorT","StatErrorT",bins.size(),0,bins.size()) ;

		SystError  = new TH1F("SystError","SystError",bins.size(),0,bins.size()) ;

		ProtonCounts    = new TH1F("Proton Counts","Proton Counts",bins.size(),0,bins.size()) ;
		DeuteronCounts  = new TH1F("Deuteron Counts","Deuteron Counts",bins.size(),0,bins.size()) ;
		TritiumCounts    = new TH1F("Tritium Counts","Deuteron Counts",bins.size(),0,bins.size()) ;
		
		TFile * f = TFile::Open("/afs/cern.ch/work/f/fdimicco/private/Deutons/DirectAnalysis/LatWeights/ExposureModel.root");		
		TH1F * Exp = (TH1F*) f->Get("HEExposure");
		Exp->Scale(1/Exp->GetBinContent(Exp->GetMaximumBin()));
		ExposureTime = GetSplineFromHisto(Exp,PRB);

		ProtonCountsPrim    = new TH1F("Primary Proton Counts","Primary Proton Counts",bins.size(),0,bins.size()) ;
        	DeuteronCountsPrim  = new TH1F("Primary Deuteron Counts","Primary Deuteron Counts",bins.size(),0,bins.size()) ;
	

		BestFitSigma  = new TH1F("Best Fit Sigma","Best Fit Sigma",bins.size(),0,bins.size()) ;
		BestFitShift  = new TH1F("Best Fit Shift","Best Fit Shift",bins.size(),0,bins.size()) ;

		BestChiSquares     = new TH1F("Best ChiSquare","Best ChiSquare",bins.size(),0,bins.size()) ;
        	OriginalChiSquares = new TH1F("Original ChiSquare","Original CHiSquare",bins.size(),0,bins.size()) ;

		isrich=IsRich;

		systpar.sigma=sigma;
		systpar.shift=shift;
		systpar.steps=steps;

		MCHeContRatio = (TH1F*) File.Get((basename+"/Fit Results/HeliumContamination/MC Expected He over P").c_str());
		//if(!MCHeContRatio) MCHeContRatio = Eval_MCHeContRatio("MC Expected He over P");

		SetFitConstraints(0.9,1,0.015,0.06,0.0001,0.005);

	}

	void SetFitConstraints(float minP, float maxP,float minD,float maxD,float minHe,float maxHe, bool highmasstailconstrain=false) {
		constrainmin[0]=minP; constrainmin[1]=minD; constrainmin[2]=minHe; 
		constrainmax[0]=maxP; constrainmax[1]=maxD; constrainmax[2]=maxHe;
		highmassconstrain = highmasstailconstrain;
	}

	bool ReinitializeHistos(bool refill){return true;}; //dummy
	void Eval_TransferFunction();
	TH1F * Eval_MCHeContRatio(std::string name);
	void Fill(TTree * treeMC,TTree * treeDT, Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars) );
	float SmearBeta(float Beta, float stepsigma, float stepshift,float R);
	float SmearBetaRICH(float Beta, float stepsigma, float stepshift);

	void FillEventByEventData(Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars));
	void FillEventByEventMC(Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars));

	float GetCutoffWeight(float particle_m, float beta, float m);
	void ConvoluteTempletesWithCutoff();
	void ExtractCounts(FileSaver finalhistos);
	void EvalFinalParameters();
	void EvalFinalErrors();
	void CalculateFinalPDCounts();
	void Save();	
	void SaveFitResults(FileSaver finalhistos);
	void SumUpMassDistrib(FileSaver finalhistos);

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
	void SetLatitudeReweighter(LatReweighter * weighter) { Latweighter = weighter;}
	void Eval_ContError();
	std::string GetName(){return basename;}
	TH1F * GetStatErrorP(){ return StatErrorP;}
	TH1F * GetStatErrorD(){ return StatErrorD;}
	TH1F * GetStatErrorT(){ return StatErrorT;}

	TH1F * GetSystError() { return SystError;}	
	Binning  GetBinning() {return bins;}
	void  RebinAll(int f=2);	
	float GetHeContaminationWeight(int bin) { return fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->ContribHe; }
	float GetHeContaminationErr   (int bin) { return fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->errHe; }

	TH2F * GetDCountsSpread(int bin)	{ return DCountsSpread[bin];}
	TH2F * GetChiSquareSpread(int bin)      { return TFitChisquare[bin];}	
	TH1F * GetWeightedDCounts(int bin)     { return WeightedDCounts[bin];}

};

#endif
