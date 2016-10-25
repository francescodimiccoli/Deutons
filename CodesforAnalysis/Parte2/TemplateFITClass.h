



using namespace std;

struct TFit {
   TH1F * Templ_P ;
   TH1F * Templ_D ;
   TH1F * Templ_He;
   TH1F * Data;
   float wheightP,wheightD,wheightHe;
   TFractionFitter *Tfit;
   int Tfit_outcome;
};


class TemplateFIT {

   private:
      std::vector<std::vector<TFit *>> fits;
	int binstemplate = 50;

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

      bool TemplateFITenabled = true;

      //Fit constraints
      std::vector<float>  lowP ,lowD ,lowHe    ;
      std::vector<float>  highP,highD,highHe   ;
	
      float tolerance=0;
      float minrange,maxrange=0;
      //creation constructors
      //standard
      TemplateFIT(std::string basename ,int Nbins ,float val_min , float val_max,int mc_types = 0, int n=0)
      {

         TemplateP   = 	new TH2F((basename + "_P"   ).c_str(),(basename + "_P"   ).c_str(),binstemplate,val_min,val_max,Nbins,0,Nbins);
         TemplateHe  =   new TH2F((basename + "_He"  ).c_str(),(basename + "_He"  ).c_str(),binstemplate,val_min,val_max,Nbins,0,Nbins);

         if(mc_types == 0)
            TemplateD   =	new TH2F((basename + "_D"   ).c_str(),(basename + "_D"     ).c_str(),binstemplate,val_min,val_max,Nbins,0,Nbins);
         if(mc_types >  0 )
            TemplateD   =   new TH3F((basename + "_D"   ).c_str(),(basename + "_D"     ).c_str(),binstemplate,val_min,val_max,Nbins,0,Nbins,mc_types,0,mc_types);
         if(n == 0) {
            DATA        =	new TH2F((basename + "_Data").c_str(),(basename + "_Data"  ).c_str(),binstemplate,val_min,val_max,Nbins,0,Nbins);
            Geomag = false;
         }
         if(n > 0) {
            DATA        =   new TH3F((basename + "_Data").c_str(),(basename + "_Data"  ).c_str(),binstemplate,val_min,val_max,Nbins,0,Nbins,n,0,n);
            Geomag = true;
         }
         PCounts	   =	new TH1F((basename + "_PCounts"  	 ).c_str(),(basename + "_PCounts"   	).c_str(),Nbins,0,Nbins);
         DCounts	   =	new TH1F((basename + "_DCounts"  	 ).c_str(),(basename + "_DCounts"   	).c_str(),Nbins,0,Nbins);

         nbins = Nbins;
      }



      //reading constructors
      //standard
      TemplateFIT(TFile * file , std::string basename_MC , std::string basename_data)
      {

         TemplateP   =	(TH1 *)file->Get((basename_MC   + "_P"     ).c_str());
         TemplateD   =	(TH1 *)file->Get((basename_MC   + "_D"     ).c_str());
         TemplateHe  =	(TH1 *)file->Get((basename_MC   + "_He"    ).c_str());

         DATA        =	(TH1 *)file->Get((basename_data + "_Data"     ).c_str());

         nbins =  TemplateP -> GetNbinsY();

         PCounts	   =	new TH1F((basename_data + "_PCounts"  	 ).c_str(),(basename_data + "_PCounts"   	).c_str(),nbins,0,nbins);
         DCounts	   =	new TH1F((basename_data + "_DCounts"  	 ).c_str(),(basename_data + "_DCounts"   	).c_str(),nbins,0,nbins);

         Geomag = false;

         fits.push_back(std::vector<TFit *>());
	 minrange = 0;
	 maxrange = 100;
      }
      //geom. zones
      TemplateFIT(TFile * file , std::string basename_MC, std::string basename_data, int n)
      {

         TemplateP   =	(TH1 *)file->Get((basename_MC   + "_P"     ).c_str());
         TemplateD   =	(TH1 *)file->Get((basename_MC   + "_D"     ).c_str());
         TemplateHe  =	(TH1 *)file->Get((basename_MC   + "_He"    ).c_str());

         DATA        =	(TH1 *)file->Get((basename_data + "_Data"   ).c_str());

         nbins =  TemplateP -> GetNbinsY();

         PCounts = 	new TH2F((basename_data + "_PCounts"      ).c_str(),(basename_data + "_PCounts"   ).c_str(),nbins,0,nbins,n,0,n);
         DCounts = 	new TH2F((basename_data + "_DCounts"      ).c_str(),(basename_data + "_DCounts"   ).c_str(),nbins,0,nbins,n,0,n);

         Geomag = true;

         for(int lat=0; lat<n; lat ++) fits.push_back(std::vector<TFit *>());
     	 minrange = 0;
	 maxrange = 100;

	 }



      //Methods

      void Write();

      TH1F * Extract_Bin(TH1 * Histo, int bin,int third_dim=0,bool reverse=false);
      void SetFitConstraints(float LowP=0,float HighP=1, float LowD=0,float HighD=1, float LowHe=0,float HighHe=1);
      void SetFitConstraints(TH1F * ContHe, float LowP=0,float HighP=1, float LowD=0,float HighD=1);
      void SetTolerance(float tol);
      void SetFitRange(float min,float max); 
      void Do_TemplateFIT(TFit * Fit, int bin, int lat=0);
      int GetFitOutcome(uint bin,uint lat=0)
      {
         if(lat < fits.size() && bin < fits[lat].size())
				return fits[lat][bin]->Tfit_outcome;
         cout<<"Fit not yet performed: bin nr. "<<bin<<endl;
         return -2;
      };
      void PrintResults(int bin,int lat=0);
      double GetFitWheights(int par, int bin,int lat=0);
      double GetFitFraction(int par, int bin,int lat=0);
      double GetFitErrors(int par,int bin,int lat=0);
      TH1F * GetResult_P (int bin, int lat=0) {TH1F * res = (TH1F*)fits[lat][bin]->Templ_P ->Clone(); res -> Scale(fits[lat][bin]->wheightP); return res;};
      TH1F * GetResult_D (int bin, int lat=0) {TH1F * res = (TH1F*)fits[lat][bin]->Templ_D ->Clone(); res -> Scale(fits[lat][bin]->wheightD); return res;};
      TH1F * GetResult_He(int bin, int lat=0) {TH1F * res = (TH1F*)fits[lat][bin]->Templ_He ->Clone(); res -> Scale(fits[lat][bin]->wheightHe); return res;};
      TH1F * GetResult_Data(int bin,int lat=0) {return (TH1F*)fits[lat][bin] -> Data;};
      void DisableFit();
      void TemplateFits(int mc_type=1);
      void TemplateFitPlot(TVirtualPad * c, std::string var_name,int bin,int lat=0);
};
