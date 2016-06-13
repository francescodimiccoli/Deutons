using namespace std;



class Flux {
   private:
      //Delta E;
      TH1 * DeltaE;
      std::string suffixname;
      Binning bins;

   public:

      TH1 * Counts;      ///<counts
      TH1 * Acceptance;  ///<Acceptance
      TH1 * Exposure;    ///<Exposure Time
      TH1 * Fluxes;      ///<Fluxes

      std::string name;

      ///creation constructor
      Flux (std::string basename, std::string suffix, Binning B);
      Flux (std::string basename, std::string suffix, Binning B, int n);
      ///reading constructor
      Flux (TFile * file, std::string basename, std::string suffix, std::string dirname, std::string acceptname, Binning B);
      
      void Write();
      void Set_Exposure_Time (TH2 * Exposure );
      void Set_Exposure_Time (TH1 * Tempi );
      TH1 * ExtractParticularMC_cs (TH1 * Histo, int lat_zones=1 ,int mc_type=0);
      void Set_DeltaE (int n);
      void Add_SystFitError (TH1* syst_err);
      void Eval_Flux (int n,bool deutons=1,int mc_type=0);
};


void Flux::Write()
{
   if (Counts->GetEntries() >0) Counts->Write();
   return;
}


void Flux::Set_Exposure_Time (TH2 * hExposure )
{
   Exposure  = (TH1 *) hExposure -> ProjectionX ("",0,10) -> Clone();
   return;
}


void Flux::Set_Exposure_Time (TH1 * Tempi)
{
   string hname=name + "TimeZone_" + suffixname;
   Exposure  = new TH2F ( hname.c_str(), hname.c_str(), nbinsr,  0,nbinsr,   Tempi->GetNbinsX(),0,Tempi->GetNbinsX() );
      for (int R=0; R<Exposure ->GetNbinsX() ; R++ )
      for (int lat =1; lat< Tempi->GetNbinsX(); lat++) {
         ((TH2 *) Exposure)  -> SetBinContent (R+1,lat+1,Tempi -> GetBinContent (lat) );
         ((TH2 *) Exposure)  -> SetBinError (R+1,lat+1,10);
      }
   return;
}



TH1 * Flux::ExtractParticularMC_cs (TH1 * Histo, int lat_zones, int mc_type)
{
   if (lat_zones == 1) {
      TH1F * Slice = new TH1F ("","",Histo->GetNbinsX(),0,Histo->GetNbinsX() );
      for (int i = 0; i< Histo->GetNbinsX(); i++) {
         Slice->SetBinContent (i+1,Histo->GetBinContent (i+1,mc_type+1) );
         Slice->SetBinError (i+1,Histo->GetBinError (i+1,mc_type+1) );
      }
      return Slice;
   } else {
      TH2F * Slice = new TH2F ("","",Histo->GetNbinsX(),0,Histo->GetNbinsX(),lat_zones,0,lat_zones);
      for (int lat=0; lat<lat_zones; lat++) {
         for (int i = 0; i< Histo->GetNbinsX(); i++) {
            Slice->SetBinContent (i+1,lat+1,Histo->GetBinContent (i+1,lat+1,mc_type+1) );
            Slice->SetBinError (i+1,lat+1,Histo->GetBinError (i+1,lat+1,mc_type+1) );
         }
      }
      return Slice;
   }
}



void Flux::Set_DeltaE (int n)
{
   string hname=name + "DeltaE_" + suffixname;
   
   if (n>1) {
      DeltaE   = new TH2F (  hname.c_str(), hname.c_str(), bins.size(),  0, bins.size(),  n,0,n);
      DeltaE   ->Sumw2();
      for (int iR=0; iR<DeltaE->GetNbinsX(); iR++)
         for (int lat =0; lat<DeltaE->GetNbinsX(); lat++)
            DeltaE->SetBinContent (iR+1,lat+1, bins.EkBin (iR) - bins.EkBin (iR-1) );
   } else { // 1D
      DeltaE   = new TH1F ( hname.c_str(), hname.c_str() ,bins.size(),  0, bins.size()  );
      DeltaE   ->Sumw2();
      for (int iR=1; iR<DeltaE->GetNbinsX(); iR++)
         DeltaE->SetBinContent (iR+0,bins.EkBin (iR) - bins.EkBin (iR-1) );
   }
}


void Flux::Add_SystFitError (TH1* syst_err)
{
   if (Counts)
      for (int iR = 0; iR < Counts    ->GetNbinsX(); iR++)
         Counts  ->   SetBinError (iR+1,Counts-> GetBinError (iR+1)+fabs (syst_err->GetBinContent (iR+1) ) ) ;
   return;

}


void Flux::Eval_Flux (int n,bool deutons,int mc_type)
{
   Set_DeltaE (deutons);

   if (Counts) {
      Fluxes = (TH1 *) Counts -> Clone();
      if (deutons)  Fluxes->Divide (ExtractParticularMC_cs ( Acceptance , n,mc_type) );
      else          Fluxes->Divide (Acceptance);
      Fluxes->Divide (Exposure);
      Fluxes->Divide (DeltaE);
   }

   return;
}


Flux::Flux (std::string basename, std::string suffix, Binning inbin)
{
   string hname=basename + "Counts_" + suffix;
   Counts	 = new TH1F ( hname.c_str(), hname.c_str(),inbin.size(),0,inbin.size() );
   name = basename;
   suffixname = suffix;
   bins=inbin;
}

Flux::Flux (std::string basename, std::string suffix, Binning inbin, int n)
{
   string hname=basename + "Counts_" + suffix;
   Counts   = new TH2F ( hname.c_str(), hname.c_str(), inbin.size(),  0, inbin.size(),  n,0,n);
   name = basename;
   suffixname = suffix;
   bins=inbin;
}

Flux::Flux (TFile * file, std::string basename, std::string suffix, std::string dirname, std::string acceptname, Binning inbin)
{
   Counts	 = (TH1 *) file->Get ( (basename + "Counts_" + suffix       ).c_str() );
   Acceptance  = (TH1 *) file->Get ( ("/" + dirname +"/" +acceptname + "_" + suffix  ).c_str() );
   name = basename;
   suffixname = suffix;
   bins=inbin;
}
