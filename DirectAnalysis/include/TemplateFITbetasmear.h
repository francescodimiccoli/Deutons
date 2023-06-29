#ifndef TEMPLATEFITBETASMEAR_H
#define TEMPLATEFITBETASMEAR_H

#include "LatReweighter.h" 
#include "Livetime.h"
#include "TFractionFitter.h"
#include <TGraphErrors.h>
#include "Tool.h"
#include "GlobalPaths.h"
#include "PlottingFunctions.h"
#include "ShiftPeak.h"

class Tool;

using namespace std;

float FindConstraintD(TH1F * dt);
	
struct TFit {
   TH1F * Templ_P =0x0;
   TH1F * Templ_D =0x0;
   TH1F * Templ_DPrim=0x0;
   TH1F * Templ_PPrim=0x0;
   TH1F * Templ_He=0x0;
   TH1F * Templ_HePrim=0x0;
  
   TH1F * Templ_Noise=0x0; 	

   TH1F * Data =0x0;
   TH1F * DataAmbient =0x0;
   TH1F * DataPrim=0x0;	

   TH1F * Residual=0x0;	


   float wheightP,wheightD,wheightHe,wheightNoise;
   float ContribP,ContribD,ContribHe,ContribNoise;
   float errP,errD,errHe,errNoise;
   float fitrangemin = 0.5;
   float fitrangemax = 7;
		
   TFractionFitter *Tfit;
   int Tfit_outcome=-1;
   
   float DCounts=0;
   float PCounts=0;
   float TCounts=0;


   float ChiSquare=500;			
   float ndf=0;
 	  
   float StatErrP=0;
   float StatErrD=0;
   float StatErrT=0;
	

  
   float DCountsPrim=0;		
   float PCountsPrim=0;	



   TFit(){}	
   TFit(TH1F * templ_P, TH1F * templ_D, TH1F * data, TH1F * dataPrim, TH1F* dataAmbient) { Templ_P= templ_P; Templ_D=templ_D; Data=data; DataPrim=dataPrim; DataAmbient = dataAmbient; }
  TFit(TH1F * templ_P, TH1F * templ_D, TH1F * data, TH1F * dataPrim) { Templ_P= templ_P; Templ_D=templ_D; Data=data; DataPrim=dataPrim; }
  void RegularizeTemplateError();
  void AddSysttoTempl();
  TFit * Clone();	
};


struct Systpar{
	int mode;
	int steps;
	float sigma;
	float shift;

	float SigmaStart;	
	float SigmaEnd;	
	float ShiftStart;	
	float ShiftEnd;	


	float GetSigmaStart();
	float GetSigmaEnd();
	float GetShiftStart();
	float GetShiftEnd();
};

struct BestChi {
	int i=0;
	int j=0;
	int z=0;
	float chimin=0;
	float centroid_x=0;
	float centroid_y=0;
	float centroid_parx=0;
	float centroid_pary=0;
	void FindMinimum(TH2F * Histo) {
			Histo->GetMinimumBin(i,j,z);
			i--;
			j--;
			chimin = Histo->GetBinContent(i+1,j+1);
		}
	void FindMinimum(TH2F * Histo, int fixed_i) {
			TH1D * s = (TH1D*) (Histo->ProjectionY("s",fixed_i+1,fixed_i+1));
			j=s->GetMinimumBin();
			j--;
			i=fixed_i;
			chimin = Histo->GetBinContent(fixed_i+1,j+1);
			}
	void FindCentroid(TH2F * histo) {
		TH2F* Histo = (TH2F*) histo->Clone("histo");
		float chimin= Histo->GetMinimum();
		for(int i=0;i<Histo->GetNbinsX();i++){
			for(int j=0;j<Histo->GetNbinsX();j++){
			if(Histo->GetBinContent(i+1,j+1)>15*chimin) Histo->SetBinContent(i+1,j+1,15*chimin); 
			}  
		}
		TH1F * Centroid_X = (TH1F*) Histo->ProjectionX();
		TH1F * Centroid_Y = (TH1F*) Histo->ProjectionY();
//		FindMinimum(Histo);
	
		float binx1 = Centroid_X->GetMinimumBin();
		float biny1 = Centroid_Y->GetMinimumBin();
		float parx1 = Centroid_X->GetBinCenter(binx1);
		float pary1 = Centroid_Y->GetBinCenter(biny1);
		float wx1 = Centroid_X->GetBinContent(binx1);
		float wy1 = Centroid_Y->GetBinContent(biny1);
/*		Centroid_X->SetBinContent(Centroid_X->GetMinimumBin(),1e7);
		Centroid_Y->SetBinContent(Centroid_Y->GetMinimumBin(),1e7);
		float binx2 = Centroid_X->GetMinimumBin();
		float biny2 = Centroid_Y->GetMinimumBin();
		float parx2 = Centroid_X->GetBinCenter(binx2);
		float pary2 = Centroid_Y->GetBinCenter(biny2);
		float wx2 = Centroid_X->GetBinContent(binx2);
		float wy2 = Centroid_Y->GetBinContent(biny2);
		Centroid_X->SetBinContent(Centroid_X->GetMinimumBin(),1e7);
		Centroid_Y->SetBinContent(Centroid_Y->GetMinimumBin(),1e7);
		float binx3 = Centroid_X->GetMinimumBin();
		float biny3 = Centroid_Y->GetMinimumBin();
		float parx3 = Centroid_X->GetBinCenter(binx3);
		float pary3 = Centroid_Y->GetBinCenter(biny3);
		float wx3 = Centroid_X->GetBinContent(binx3);
		float wy3 = Centroid_Y->GetBinContent(biny3);
		Centroid_X->SetBinContent(Centroid_X->GetMinimumBin(),1e7);
		Centroid_Y->SetBinContent(Centroid_Y->GetMinimumBin(),1e7);
		float binx4 = Centroid_X->GetMinimumBin();
		float biny4 = Centroid_Y->GetMinimumBin();
		float parx4 = Centroid_X->GetBinCenter(binx4);
		float pary4 = Centroid_Y->GetBinCenter(biny4);
		float wx4 = Centroid_X->GetBinContent(binx4);
		float wy4 = Centroid_Y->GetBinContent(biny4);
*/

		float binx2 = binx1-1;
		float biny2 = biny1-1;
		float parx2 = Centroid_X->GetBinCenter(binx2);
		float pary2 = Centroid_Y->GetBinCenter(biny2);
		float wx2 = Centroid_X->GetBinContent(binx2);
		float wy2 = Centroid_Y->GetBinContent(biny2);
		float binx3 = binx1+1;
		float biny3 = biny1+1;
		float parx3 = Centroid_X->GetBinCenter(binx3);
		float pary3 = Centroid_Y->GetBinCenter(biny3);
		float wx3 = Centroid_X->GetBinContent(binx3);
		float wy3 = Centroid_Y->GetBinContent(biny3);
		
		if(wx2==0) wx2=10000;
		if(wx3==0) wx3=10000;


		centroid_x = ((binx1/wx1)+(binx2/wx2)+(binx3/wx3))/(1/wx1+1/wx2+1/wx3);
		centroid_y = ((biny1/wy1)+(biny2/wy2)+(biny3/wy3))/(1/wy1+1/wy2+1/wy3);
		centroid_parx = ((parx1/wx1)+(parx2/wx2)+(parx3/wx3))/(1/wx1+1/wx2+1/wx3);
		centroid_pary = ((pary1/wy1)+(pary2/wy2)+(pary3/wy3))/(1/wy1+1/wy2+1/wy3);
	
	}

};
	
void Do_TemplateFIT(TFit * Fit,float fitrangemin,float fitrangemax,float constrain_min[], float constrain_max[], bool isfitnoise, bool highmasstailconstrain, bool IsFitPrim=false );

TSpline3 * GetSplineFromHisto(TH1F * Graph, Binning bins);

void AdjustErrorsforFullStat(TH1F* histo,int ntimebins);

class TemplateFIT : public Tool{

	private:
	std::vector<std::vector<std::vector<TFit *>>> fits;
	std::vector<std::vector<std::vector<TFit *>>> fits_copy;
	std::vector<TFit *>fits_best;
	BestChi * BestChiSquare;	
	std::vector<TH1F *> TransferFunction;
	double massbins[46]={0.4,0.476667,0.553333,0.63,0.706667,0.783333,0.86,0.936667,1.01333,1.09,1.16667,1.24333,1.32,1.39667,1.47333,1.55,1.62667,1.70333,1.78,1.85667,1.93333,2.01,2.08667,2.16333,2.24,2.31667,2.39333,2.47,2.54667,2.62333,2.7,2.77667,2.85333,2.93,3.00667,3.08333,3.16,3.23667,3.39,3.54333,3.77333,4.00333,4.31,4.61667,4.92333,6};

	TF1 * MCmodel;
 
	std::vector<TH2F *> DCountsSpread;
	std::vector<TH2F *> DErrSpread;
	std::vector<TH1F *> WeightedDCounts;
	std::vector<TH2F *> TFitChisquare;

	TH1F * BetaResolutionDT;
   	TH1F * BetaResolutionMC;
 	TH1F * BetaResolutionMC_tuned;


	TH2F* Global_ChiSquare;	

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
	TH1F * ProtonTail;
	TH1F * TritiumCounts;


	TH1F * ProtonCountsPrim;
	TH1F * DeuteronCountsPrim;

	TH1F * BestChiSquares;
	TH1F * OriginalChiSquares;

	TH1F * BestSigma;
	TH1F * BestShift;


	Binning bins;
        std::string cut;
	std::string cutoff;
	std::string cutprimary;
	std::string discr_var;
	std::string cutP; 
        std::string cutD ;
        std::string cutHe;


	TFile * regularizationfile=0x0;
	TF1 * Reg_sigma=0x0;
	TF1 * Reg_shift=0x0;

	std::string basename;
	int time=1429054000;

	Systpar systpar;
	bool fitDisabled=false;
	bool isrich=false;
	BadEventSimulator * BadEvSim=0x0;
	bool IsFitNoise = false;
	float constrainmin[3];
	float constrainmax[3];
	bool highmassconstrain=false;
	int ActualTime=0;
	bool IsLocalFit=false;
	bool IsLocalConstrainedFit=false;
	bool usebestonly=false;
	bool usecentroid=false;

	float template1scalefactor1=1;
	float template1scalefactor2=1;
	float template1scalefactor3=1;

	bool lowstatDmode=false;
	bool adjousttail=false;
	bool adjoustfixedtail=false;
	float Mag=0;
	float Mid=0;
	float Midbin=0;
	float Fast=0;
	int Forcesigma=-1;
	std::vector<TFit*> tailFT;
	std::vector<TH1F *> chitail;

	bool adjoustpeak=false;
	int iso=0;	
	float magpeak=0;

	float sharpenfactor=0;
	bool sharpen=false;

	bool IsExtern = false;

	bool useMCreweighting=true;
	bool useMCtuning=true;

	public:	
	//standard constructor
	TemplateFIT(std::string Basename,Binning Bins, std::string Cut, std::string CutP, std::string CutD, std::string CutHe,std::string Cutoff, int Nbins, float Xmin, float Xmax, bool IsRich,int steps,float sigma,float shift,int smearingmode){
		
			for(int bin=0;bin<Bins.size();bin++){
			fits.push_back(std::vector<std::vector<TFit *>>());
			for(int i=0;i<steps;i++){
				fits[bin].push_back(std::vector<TFit *>());
				for(int j=0;j<steps;j++){

					TFit * fit = new TFit;
					string named    =Basename + "_Data_" +to_string(bin)+" "+to_string(i)+" "+to_string(j);
					string namedamb =Basename + "_DataAmbient_" +to_string(bin)+" "+to_string(i)+" "+to_string(j);
					string namedprim=Basename + "_DataPrim_" +to_string(bin)+" "+to_string(i)+" "+to_string(j);
					string nameP    =Basename + "_MCP_"      +to_string(bin)+" "+to_string(i)+" "+to_string(j);
					string nameD    =Basename + "_MCD_"      +to_string(bin)+" "+to_string(i)+" "+to_string(j);
					string nameHe   =Basename + "_MCHe_"     +to_string(bin)+" "+to_string(i)+" "+to_string(j);
					string nameNo   =Basename + "_MCNoise_"  +to_string(bin)+" "+to_string(i)+" "+to_string(j);

				if(Nbins==0){
					fit->Templ_P =  new TH1F(nameP.c_str(),nameP.c_str(),45,massbins);
					fit->Templ_D =  new TH1F(nameD.c_str(),nameD.c_str(),45,massbins);
					fit->Templ_He=  new TH1F(nameHe.c_str(),nameHe.c_str(),45,massbins);
					fit->Templ_Noise=  new TH1F(nameNo.c_str(),nameNo.c_str(),45,massbins);
					fit->Data    =  new TH1F(named.c_str(),named.c_str(),45,massbins);
					fit->DataPrim=  new TH1F(namedprim.c_str(),namedprim.c_str(),45,massbins);
					fit->DataAmbient=  new TH1F(namedamb.c_str(),namedprim.c_str(),45,massbins);
				}
				else{
					fit->Templ_P =  new TH1F(nameP.c_str(),nameP.c_str(),Nbins,Xmin,Xmax);
					fit->Templ_D =  new TH1F(nameD.c_str(),nameD.c_str(),Nbins,Xmin,Xmax);
					fit->Templ_He=  new TH1F(nameHe.c_str(),nameHe.c_str(),Nbins,Xmin,Xmax);
					fit->Templ_Noise=  new TH1F(nameNo.c_str(),nameNo.c_str(),Nbins,Xmin,Xmax);
					fit->Data    =  new TH1F(named.c_str(),named.c_str(),Nbins,Xmin,Xmax);
					fit->DataPrim=  new TH1F(namedprim.c_str(),namedprim.c_str(),Nbins,Xmin,Xmax);
					fit->DataAmbient=  new TH1F(namedamb.c_str(),namedprim.c_str(),Nbins,Xmin,Xmax);
					}		
				fits[bin][i].push_back(fit);
				}
			}
		
			if(IsRich){
			BetaResolutionDT=new TH1F("BetaResolutionDT","BetaResolutionDT",1000,0.9,1.25);
			BetaResolutionMC=new TH1F("BetaResolutionMC","BetaResolutionMC",1000,0.9,1.25);
			BetaResolutionMC_tuned=new TH1F("BetaResolutionMC_tuned","BetaResolutionMC_tuned",1000,0.9,1.25);
			}

			else{
			BetaResolutionDT=new TH1F("BetaResolutionDT","BetaResolutionDT",1000,0.7,1.5);
			BetaResolutionMC=new TH1F("BetaResolutionMC","BetaResolutionMC",1000,0.7,1.5);
			BetaResolutionMC_tuned=new TH1F("BetaResolutionMC_tuned","BetaResolutionMC_tuned",1000,0.7,1.5);
			}


		}
	
		MCmodel = new TF1("MCmodel","pol3",0,5);
		MCmodel->SetParameter(0,-0.257527);
		MCmodel->SetParameter(1,0.737386);
		MCmodel->SetParameter(2,-0.607996);
		MCmodel->SetParameter(3,0.161477);

		IsLocalFit = false;
		IsLocalConstrainedFit=false;
		basename=Basename;
		cut = Cut;
		cutoff = Cutoff;
		cutprimary=Cut+"&IsPrimary";
		bins=Bins;
	
		cutP=CutP;
		cutD=CutD;
		cutHe=CutHe;
		cutP=cut+cutP;
		cutD=cut+cutD;
		cutHe=cut+cutHe;

		StatErrorP  = new TH1F("StatErrorP","StatErrorP",bins.size(),0,bins.size()) ;
		StatErrorD  = new TH1F("StatErrorD","StatErrorD",bins.size(),0,bins.size()) ;
		StatErrorT  = new TH1F("StatErrorT","StatErrorT",bins.size(),0,bins.size()) ;

        	SystError  = new TH1F("SystError","SystError",bins.size(),0,bins.size()) ;

		ProtonCounts    = new TH1F("Proton Counts","Proton Counts",bins.size(),0,bins.size()) ;
        	DeuteronCounts  = new TH1F("Deuteron Counts","Deuteron Counts",bins.size(),0,bins.size()) ;
		ProtonTail  = new TH1F("Proton Tail","Proton Tail",bins.size(),0,bins.size()) ;
		TritiumCounts  = new TH1F("Tritium Counts","Tritium Counts",bins.size(),0,bins.size()) ;


		ProtonCountsPrim    = new TH1F("Primary Proton Counts","Primary Proton Counts",bins.size(),0,bins.size()) ;
        	DeuteronCountsPrim  = new TH1F("Primary Deuteron Counts","Primary Deuteron Counts",bins.size(),0,bins.size()) ;
	
	

		BestChiSquares     = new TH1F("Best ChiSquare","Best ChiSquare",bins.size(),0,bins.size()) ;
        	OriginalChiSquares = new TH1F("Original ChiSquare","Original CHiSquare",bins.size(),0,bins.size()) ;


		isrich = IsRich;

		if(smearingmode==0 && !IsRich) systpar.mode=0;
		if(smearingmode==1 && !IsRich) systpar.mode=1;
		if(smearingmode==0 &&  IsRich) systpar.mode=2;
		if(smearingmode==1 &&  IsRich) systpar.mode=3;


		systpar.sigma=sigma;
		systpar.shift=shift;
		systpar.steps=steps;

		SetFitConstraints(0.9,1,0.015,0.06,0.0001,0.005);
	}

	//reading constructor

	TemplateFIT(FileSaver  File, std::string Basename,Binning Bins, bool IsRich, int steps,float sigma,float shift,int smearingmode,int ntimebins=1){

		cout<<"Templates: "<<Basename<<endl;
		TFile * filecurrent = File.GetFile();
	
		TFile * file;
		file = filecurrent;

		BetaResolutionDT = (TH1F*) filecurrent->Get((Basename+"/BetaRes/BetaResolutionDT").c_str());
		BetaResolutionMC = (TH1F*) filecurrent->Get((Basename+"/BetaRes/BetaResolutionMC").c_str());
		BetaResolutionMC_tuned = (TH1F*) filecurrent->Get((Basename+"/BetaRes/BetaResolutionMC_tuned").c_str());
		cout<<BetaResolutionMC_tuned->Integral()<<endl;

		for(int bin=0;bin<Bins.size();bin++){
			fits.push_back(std::vector<std::vector<TFit *>>());
			cout<<"Reading Bin: "<<bin<<endl;
			for(int i=0;i<steps;i++){
				fits[bin].push_back(std::vector<TFit *>());
				for(int j=0;j<steps;j++){

					TFit * fit = new TFit;
					string named    =Basename + "/Bin "+ to_string(bin)+"/Data/" + Basename + "_Data_" +to_string(bin)+" "+to_string(0)+" "+to_string(5);
					string namedamb =Basename + "/Bin "+ to_string(bin)+"/Data/" + Basename + "_DataPrim_" +to_string(bin)+" "+to_string(0)+" "+to_string(5);
					string namedprim=Basename + "/Bin "+ to_string(bin)+"/Data/" + Basename + "_DataAmbient_" +to_string(bin)+" "+to_string(0)+" "+to_string(5);
					string nameP    =Basename + "/Bin "+ to_string(bin)+"/TemplateP/" + Basename + "_MCP_"      +to_string(bin)+" "+to_string(i)+" "+to_string(j);
					string nameD    =Basename + "/Bin "+ to_string(bin)+"/TemplateD/" + Basename + "_MCD_"      +to_string(bin)+" "+to_string(i)+" "+to_string(j);
					string nameHe    =Basename + "/Bin "+ to_string(bin)+"/TemplateHe/" + Basename + "_MCHe_"      +to_string(bin)+" "+to_string(i)+" "+to_string(j);
					string nameNo    =Basename + "/Bin "+ to_string(bin)+"/TemplateNoise/" + Basename + "_MCNoise_"      +to_string(bin)+" "+to_string(i)+" "+to_string(j);

					if(i==0&&j==5){
						fit->Templ_P    =  (TH1F *)filecurrent->Get(nameP.c_str());
						fit->Templ_D 	=  (TH1F *)filecurrent->Get(nameD.c_str());
						fit->Templ_He	=  (TH1F *)filecurrent->Get(nameHe.c_str());
						fit->Templ_Noise=  (TH1F *)filecurrent->Get(nameNo.c_str());
					}

					else{
						fit->Templ_P    =  (TH1F *)file->Get(nameP.c_str());
						fit->Templ_D 	=  (TH1F *)file->Get(nameD.c_str());
						fit->Templ_He	=  (TH1F *)file->Get(nameHe.c_str());
						fit->Templ_Noise=  (TH1F *)file->Get(nameNo.c_str());
					}

					/*AdjustErrorsforFullStat(fit->Templ_P ,ntimebins);	
					  AdjustErrorsforFullStat(fit->Templ_D ,ntimebins);	
					  AdjustErrorsforFullStat(fit->Templ_He,ntimebins);	
					 */
					
					fit->Data    	=  (TH1F *)filecurrent->Get(named.c_str());
					fit->DataPrim	=  (TH1F *)filecurrent->Get(namedprim.c_str());
					fit->DataAmbient	=  (TH1F *)filecurrent->Get(namedamb.c_str());



					fits[bin][i].push_back(fit);
				}
			}
		}

		if(IsExtern){
		for(int bin=0;bin<Bins.size();bin++){
			for(int i=0;i<steps;i++){
				for(int j=0;j<steps;j++){

					string nameP    =Basename + "/Bin "+ to_string(bin)+"/TemplateP/" + Basename + "_MCP_"      +to_string(bin)+" "+to_string(0)+" "+to_string(5);
					string nameD    =Basename + "/Bin "+ to_string(bin)+"/TemplateD/" + Basename + "_MCD_"      +to_string(bin)+" "+to_string(0)+" "+to_string(5);
					string nameHe    =Basename + "/Bin "+ to_string(bin)+"/TemplateHe/" + Basename + "_MCHe_"      +to_string(bin)+" "+to_string(0)+" "+to_string(5);
		
						fits[bin][i][j]->Templ_P    	=  BuildTemplateFromExternal(fits[bin][0][5]->Templ_P,(TH1F *)file->Get(nameP.c_str()),  fits[bin][i][j]->Templ_P);
						fits[bin][i][j]->Templ_D 	=  BuildTemplateFromExternal(fits[bin][0][5]->Templ_D,(TH1F *)file->Get(nameD.c_str()),  fits[bin][i][j]->Templ_D);
						fits[bin][i][j]->Templ_He	=  BuildTemplateFromExternal(fits[bin][0][5]->Templ_He,(TH1F *)file->Get(nameHe.c_str()),fits[bin][i][j]->Templ_He);
						cout<<i<<" "<<j<<fits[bin][i][j]->Templ_P->Integral()<<endl;
				}
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
		ProtonTail  = new TH1F("Proton Tail","Proton Tail",bins.size(),0,bins.size()) ;
		TritiumCounts    = new TH1F("Tritium Counts","Deuteron Counts",bins.size(),0,bins.size()) ;
		

		ProtonCountsPrim    = new TH1F("Primary Proton Counts","Primary Proton Counts",bins.size(),0,bins.size()) ;
        	DeuteronCountsPrim  = new TH1F("Primary Deuteron Counts","Primary Deuteron Counts",bins.size(),0,bins.size()) ;

		BestChiSquares     = new TH1F("Best ChiSquare","Best ChiSquare",bins.size(),0,bins.size()) ;
        	OriginalChiSquares = new TH1F("Original ChiSquare","Original CHiSquare",bins.size(),0,bins.size()) ;

		isrich=IsRich;
		IsLocalFit = false;
		IsLocalConstrainedFit = false;
	

		if(smearingmode==0 && !IsRich) systpar.mode=0;
		if(smearingmode==1 && !IsRich) systpar.mode=1;
		if(smearingmode==0 &&  IsRich) systpar.mode=2;
		if(smearingmode==1 &&  IsRich) systpar.mode=3;

		systpar.sigma=sigma;
		systpar.shift=shift;
		systpar.steps=steps;

		cout<<"SMEARING PARAMETERS "<<Basename<<": "<<endl;
		cout<<"Shift: "<<systpar.GetShiftStart()<<" - "<<systpar.GetShiftEnd()<<endl;;
		cout<<"Sigma: "<<systpar.GetSigmaStart()<<" - "<<systpar.GetSigmaEnd()<<endl;;


		MCHeContRatio = (TH1F*) File.Get((basename+"/Fit Results/HeliumContamination/MC Expected He over P").c_str());
		//if(!MCHeContRatio) MCHeContRatio = Eval_MCHeContRatio("MC Expected He over P");

		SetFitConstraints(0.9,1,0.015,0.06,0.0001,0.005);

	}

	TH1F * BuildTemplateFromExternal(TH1F* core, TH1F *reference, TH1F* external);

	void SetFitConstraints(float minP, float maxP,float minD,float maxD,float minHe,float maxHe, bool highmasstailconstrain=false) {
		constrainmin[0]=minP; constrainmin[1]=minD; constrainmin[2]=minHe; 
		constrainmax[0]=maxP; constrainmax[1]=maxD; constrainmax[2]=maxHe;
		highmassconstrain = highmasstailconstrain;
	}

	void BuildRegularizedTemplates(int bin,float centroid_x,float centroid_y, TH2F* Chi);

	void SetRegularizationFile(TFile* file,std::string timename){
		regularizationfile=file;
		time=std::atoi(timename.substr(timename.find("-")+1,10).c_str())-4665600;
                std::cout<<"Set_Regularization: "<<time<<" "<<regularizationfile<<std::endl;
	
	}

	bool ReinitializeHistos(bool refill){return true;}; //dummy
	void Eval_TransferFunction();
	TH1F * Eval_MCHeContRatio(std::string name);
	void Fill(TTree * treeMC,TTree * treeDT, Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars) );
	float SmearBeta(float Beta,float M_gen, float stepsigma, float stepshift,float R);
	float SmearBeta_v2(float Beta, float Beta_gen,float  HighRsigmaMC, float HighRsigmaDT, float relativeshift , TF1* MCmodel,  float stepsigma, float stepshift);
	float SmearBetaRICH(float Beta, float stepsigma, float stepshift);
	float SmearBetaRICH_v2(int kbin, float Beta, float Betagen,float HighRsigmaMC,float HighRsigmaDT, float relativeshift, float stepsigma, float stepshift);
	float SmearRRICH_R(float R, float R_gen, float stepsigma, float stepshift);

	void FillEventByEventData(Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars));
	void FillEventByEventMC(Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars));

	float GetCutoffWeight(float particle_m, float beta, float m);
	void ConvoluteTempletesWithCutoff();
	void ExtractCounts(FileSaver finalhistos,int force_shift=-1);
	void EvalFinalParameters();
	void EvalFinalErrors();
	void CalculateFinalPDCounts();
	void Save();	
	void SaveFitResults(FileSaver finalhistos);
	void SumUpMassDistrib(FileSaver finalhistos);
	void SimpleExtractPrimaries();
	void SetUpBadEventSimulator(BadEventSimulator * Sim) {BadEvSim = Sim; return; };
	void SetFitWithNoiseMode(){IsFitNoise = true; if(BadEvSim)  BadEvSim->SetFrequency(1); return;}

	void SetSystematicParameters(int steps,float sigma,float shift){ systpar.steps=steps; systpar.shift=shift; systpar.sigma=sigma; return;};
	void SetFitRange(float min, float max){ 
		for(int bin=0;bin<bins.size();bin++) for(int i=0;i<systpar.steps;i++) for(int j=0;j<systpar.steps;j++) {
			fits[bin][i][j]->fitrangemin=min; 
			fits[bin][i][j]->fitrangemax=max; 
			}
			return;
		}
	void SetFitRangeByQuantiles(float quant_min,float quant_max);
	void DisableFit(){fitDisabled=true;}
	void SetHeliumContamination(TF1 * HelimCont) {HeContModel=HelimCont; return;};
	BadEventSimulator * GetBadEventSimulator() {return BadEvSim;}
	void LoadEventIntoBadEvSim(Variables * vars) {if(BadEvSim) BadEvSim->LoadEvent(vars);}
	void Eval_ContError();
	std::string GetName(){return basename;}
	TH1F * GetStatErrorP(){ return StatErrorP;}
	TH1F * GetStatErrorD(){ return StatErrorD;}
	TH1F * GetStatErrorT(){ return StatErrorT;}

	void SetTemplateScaleFactor(float scale1,float scale2, float scale3) {template1scalefactor1=scale1; template1scalefactor2=scale2; template1scalefactor3=scale3;}
	TH1F * GetSystError() { return SystError;}	
	Binning  GetBinning() {return bins;}
	void  RebinAll(int f=2);	
	float GetHeContaminationWeight(int bin) { return fits[bin][BestChiSquare->i][BestChiSquare->j]->ContribHe; }
	float GetHeContaminationErr   (int bin) { return fits[bin][BestChiSquare->i][BestChiSquare->j]->errHe; }
	void SetLocalFit(){IsLocalFit = true; }
	void SetLocalConstrainedFit(){IsLocalConstrainedFit = true; }
	void SetLowStatDMode(){lowstatDmode=true;}
	void SetAdjoustTailMode(float magn,float midpoint, float enmidpoint, float fast, int forcesigma=-1){Mag=magn; Mid=midpoint; Midbin=enmidpoint; Fast=fast; Forcesigma=forcesigma; adjousttail=true;}
	void SetFixModTailMode(float magn,float midpoint, float enmidpoint, float fast, int forcesigma=-1){Mag=magn; Mid=midpoint; Midbin=enmidpoint; Fast=fast; Forcesigma=forcesigma; adjoustfixedtail=true;}
	void SetAdjoustPeakMode(int isotope, float mag){iso=isotope; magpeak=mag; adjoustpeak=true;}
	void SetSharpenMode(float SharpenFactor) {sharpen=true; sharpenfactor=SharpenFactor;}

	TH2F * GetDCountsSpread(int bin)	{ return DCountsSpread[bin];}
	TH2F * GetChiSquareSpread(int bin)      { return TFitChisquare[bin];}	
	TH1F * GetWeightedDCounts(int bin)     { return WeightedDCounts[bin];}
	
	void SetAsExtern() {IsExtern=true;}
	void SetNotWeightedMC() {useMCreweighting=false;}
	void SetNotTunedMC() {useMCtuning=false;}
	void SetUseBestOnly() {usebestonly=true;}
	void UseCentroid() {usecentroid=true;}


};

#endif
