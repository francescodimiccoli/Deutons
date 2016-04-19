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
	std::vector<TFit *> fits;

public:
	// Templates
	TH1 * TemplateP	;	
	TH1 * TemplateD	;
	TH1 * TemplateHe;
	
	//Data
	TH1 * Data_Prim	 ;
	TH1 * Data_Geomag;
	
	//Counts
	TH1 * PCounts	;
        TH1 * DCounts	;	
	TH1 * PCountsgeo;	
        TH1 * DCountsgeo;	
	
	int nbins;
	//creation constructors
	TemplateFIT(std::string basename ,int Nbins ,float val_min , float val_max,int n){

		TemplateP   = 	new TH2F((basename + "_P"   ).c_str(),(basename + "_P"   ).c_str(),100,val_min,val_max,Nbins,0,Nbins);
                TemplateD   =	new TH2F((basename + "_D"   ).c_str(),(basename + "_D"   ).c_str(),100,val_min,val_max,Nbins,0,Nbins);
                TemplateHe  =	new TH2F((basename + "_He"  ).c_str(),(basename + "_He"  ).c_str(),100,val_min,val_max,Nbins,0,Nbins);

	Data_Prim   =	new TH2F((basename + "_Data_Prim"  ).c_str(),(basename + "_Data_Prim"  ).c_str(),100,val_min,val_max,Nbins,0,Nbins);
		Data_Geomag =	new TH3F((basename + "_Data_Geomag").c_str(),(basename + "_Data_Geomag").c_str(),100,val_min,val_max,Nbins,0,Nbins,n,0,n);		
		
		
		PCounts	   =	new TH1F((basename + "_PCounts"  	 ).c_str(),(basename + "_PCounts"   	).c_str(),Nbins,0,Nbins);
		DCounts	   =	new TH1F((basename + "_DCounts"  	 ).c_str(),(basename + "_DCounts"   	).c_str(),Nbins,0,Nbins);
        	PCountsgeo = 	new TH2F((basename + "_PCounts_geo"      ).c_str(),(basename + "_PCounts_geo"   ).c_str(),Nbins,0,Nbins,n,0,n);
        	DCountsgeo = 	new TH2F((basename + "_DCounts_geo"      ).c_str(),(basename + "_DCounts_geo"   ).c_str(),Nbins,0,Nbins,n,0,n);
	
		nbins = Nbins;

		
	}

	//reading constructor
	TemplateFIT(TFile * file , std::string basename , float val_min , float val_max,int n){

		TemplateP   =	(TH1 *)file->Get((basename + "_P"     ).c_str());	
                TemplateD   =	(TH1 *)file->Get((basename + "_D"     ).c_str());
                TemplateHe  =	(TH1 *)file->Get((basename + "_He"    ).c_str());
	                     
                Data_Prim   =	(TH1 *)file->Get((basename + "_Data_Prim"     ).c_str());
	        Data_Geomag =	(TH1 *)file->Get((basename + "_Data_Geomag"   ).c_str());
		
		nbins =  TemplateP -> GetNbinsY();
			
		PCounts	   =	new TH1F((basename + "_PCounts"  	 ).c_str(),(basename + "_PCounts"   	).c_str(),nbins,0,nbins);
		DCounts	   =	new TH1F((basename + "_DCounts"  	 ).c_str(),(basename + "_DCounts"   	).c_str(),nbins,0,nbins);
        	PCountsgeo = 	new TH2F((basename + "_PCounts_geo"      ).c_str(),(basename + "_PCounts_geo"   ).c_str(),nbins,0,nbins,n,0,n);
        	DCountsgeo = 	new TH2F((basename + "_DCounts_geo"      ).c_str(),(basename + "_DCounts_geo"   ).c_str(),nbins,0,nbins,n,0,n);
		
	}
	void Write();
	
	TH1F * Extract_Bin_histos(TH1 * Histo, int bin);

	TH1F * Extract_Bin_histos_geo(TH1 * Histo, int bin, int lat);
	
	void Do_TemplateFIT(TFit * Fit);
	
	int GetFitOutcome(int bin){if(fits[bin]) return fits[bin]->Tfit_outcome; else {cout<<"Fit not yet performed: bin nr. "<<bin<<endl; return -1;}}

	double GetFitWheights(int par, int bin);

	double GetFitErrors(int par,int bin);

	TH1F * GetResult_P (int bin){ TH1F *res =(TH1F*)fits[bin] -> Templ_P -> Clone() ; res -> Scale(GetFitWheights(0,bin)); return res;}	
	TH1F * GetResult_D (int bin){ TH1F *res =(TH1F*)fits[bin] -> Templ_D -> Clone() ; res -> Scale(GetFitWheights(1,bin)); return res;}
	TH1F * GetResult_He(int bin){ TH1F *res =(TH1F*)fits[bin] -> Templ_He-> Clone() ; res -> Scale(GetFitWheights(2,bin)); return res;}

	TH1F * GetResult_Data(int bin)	      { return fits[bin] -> Data; }
	
	void TemplateFits();
	
	void TemplateFitPlot(TCanvas * c, std::string var_name,int bin);
};
