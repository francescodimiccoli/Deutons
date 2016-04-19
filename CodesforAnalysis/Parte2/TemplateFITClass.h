#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <cstring>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <cstdlib>
#include <string>
#include <array>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF2.h"
#include "TVector3.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSpline.h"
#include "TFractionFitter.h"
#include "THStack.h"
#include "TNtuple.h"
#include "TObject.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"



using namespace std;

struct TFit{
	TH1F * Templ_P ;
	TH1F * Templ_D ;
	TH1F * Templ_He;
	TH1F * Data;
	TFractionFitter *Tfit;	
	int Tfit_outcome;
};


class TemplateFIT
{

private:
	std::vector<std::vector<TFit *>> fits;

public:
	// Templates
	TH1 * TemplateP	;	
	TH1 * TemplateD	;
	TH1 * TemplateHe;
	
	//Data
	TH1 * DATA	 ;
	
	//Counts
	TH1 * PCounts	;
        TH1 * DCounts	;	
	
	int nbins;
	bool Geomag;

	//creation constructors
	TemplateFIT(std::string basename ,int Nbins ,float val_min , float val_max){

		TemplateP   = 	new TH2F((basename + "_P"   ).c_str(),(basename + "_P"   ).c_str(),100,val_min,val_max,Nbins,0,Nbins);
                TemplateD   =	new TH2F((basename + "_D"   ).c_str(),(basename + "_D"   ).c_str(),100,val_min,val_max,Nbins,0,Nbins);
                TemplateHe  =	new TH2F((basename + "_He"  ).c_str(),(basename + "_He"  ).c_str(),100,val_min,val_max,Nbins,0,Nbins);

		DATA        =	new TH2F((basename + "_Data"  ).c_str(),(basename + "_Data"  ).c_str(),100,val_min,val_max,Nbins,0,Nbins);
		
		PCounts	   =	new TH1F((basename + "_PCounts"  	 ).c_str(),(basename + "_PCounts"   	).c_str(),Nbins,0,Nbins);
		DCounts	   =	new TH1F((basename + "_DCounts"  	 ).c_str(),(basename + "_DCounts"   	).c_str(),Nbins,0,Nbins);
	
		nbins = Nbins;
		
		Geomag = false;	
	}

	TemplateFIT(std::string basename ,int Nbins ,float val_min , float val_max, int n){

                TemplateP   =   new TH2F((basename + "_P"   ).c_str(),(basename + "_P"   ).c_str(),100,val_min,val_max,Nbins,0,Nbins);
                TemplateD   =   new TH2F((basename + "_D"   ).c_str(),(basename + "_D"   ).c_str(),100,val_min,val_max,Nbins,0,Nbins);
                TemplateHe  =   new TH2F((basename + "_He"  ).c_str(),(basename + "_He"  ).c_str(),100,val_min,val_max,Nbins,0,Nbins);

                DATA 	    =   new TH3F((basename + "_Data").c_str(),(basename + "_Data").c_str(),100,val_min,val_max,Nbins,0,Nbins,n,0,n);

                PCounts     =    new TH2F((basename + "_PCounts"      ).c_str(),(basename + "_PCounts"   ).c_str(),Nbins,0,Nbins,n,0,n);
                DCounts     =    new TH2F((basename + "_DCounts"      ).c_str(),(basename + "_DCounts"   ).c_str(),Nbins,0,Nbins,n,0,n);

                nbins = Nbins;

		Geomag = true;
        }



	//reading constructor
	TemplateFIT(TFile * file , std::string basename_MC , std::string basename_data, float val_min , float val_max){

		TemplateP   =	(TH1 *)file->Get((basename_MC   + "_P"     ).c_str());	
                TemplateD   =	(TH1 *)file->Get((basename_MC   + "_D"     ).c_str());
                TemplateHe  =	(TH1 *)file->Get((basename_MC   + "_He"    ).c_str());
	                     
                DATA        =	(TH1 *)file->Get((basename_data + "_Data"     ).c_str());
		
		nbins =  TemplateP -> GetNbinsY();
			
		PCounts	   =	new TH1F((basename_data + "_PCounts"  	 ).c_str(),(basename_data + "_PCounts"   	).c_str(),nbins,0,nbins);
		DCounts	   =	new TH1F((basename_data + "_DCounts"  	 ).c_str(),(basename_data + "_DCounts"   	).c_str(),nbins,0,nbins);
		
		Geomag = false;
		
		fits.push_back(std::vector<TFit *>());
	}

	TemplateFIT(TFile * file , std::string basename_MC, std::string basename_data, float val_min , float val_max,int n){

		TemplateP   =	(TH1 *)file->Get((basename_MC   + "_P"     ).c_str());	
                TemplateD   =	(TH1 *)file->Get((basename_MC   + "_D"     ).c_str());
                TemplateHe  =	(TH1 *)file->Get((basename_MC   + "_He"    ).c_str());
	                     
	        DATA        =	(TH1 *)file->Get((basename_data + "_Data"   ).c_str());
		
		nbins =  TemplateP -> GetNbinsY();
			
        	PCounts = 	new TH2F((basename_data + "_PCounts"      ).c_str(),(basename_data + "_PCounts"   ).c_str(),nbins,0,nbins,n,0,n);
        	DCounts = 	new TH2F((basename_data + "_DCounts"      ).c_str(),(basename_data + "_DCounts"   ).c_str(),nbins,0,nbins,n,0,n);
		
		Geomag = true;	
		
		for(int n=0; n<nbins;n++) fits.push_back(std::vector<TFit *>());
	}



	//Methods

	void Write();
	
	TH1F * Extract_Bin(TH1 * Histo, int bin,int lat=0);
	
	void Do_TemplateFIT(TFit * Fit, int lat=0);
	
	int GetFitOutcome(int bin,int lat=0){if(fits[bin][lat]) return fits[bin][lat]->Tfit_outcome; else {cout<<"Fit not yet performed: bin nr. "<<bin<<endl; return -1;}};

	double GetFitWheights(int par, int bin,int lat=0);

	double GetFitErrors(int par,int bin,int lat=0);

	TH1F * GetResult_P (int bin, int lat=0){ TH1F *res =(TH1F*)fits[bin][lat] -> Templ_P -> Clone() ; res -> Scale(GetFitWheights(0,bin,lat)); return res;};	
	TH1F * GetResult_D (int bin, int lat=0){ TH1F *res =(TH1F*)fits[bin][lat] -> Templ_D -> Clone() ; res -> Scale(GetFitWheights(1,bin,lat)); return res;};
	TH1F * GetResult_He(int bin, int lat=0){ TH1F *res =(TH1F*)fits[bin][lat] -> Templ_He-> Clone() ; res -> Scale(GetFitWheights(2,bin,lat)); return res;};

	TH1F * GetResult_Data(int bin,int lat=0)	      { return fits[bin][lat] -> Data; };
	
	void TemplateFits();
	
	void TemplateFitPlot(TCanvas * c, std::string var_name,int bin,int lat=0);
};
