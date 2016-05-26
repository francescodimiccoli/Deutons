using namespace std;


void TimeGeo (TH2 * Exposure, TH1 * Tempi)
{
   for (int R=0; R<Exposure ->GetNbinsX() ; R++ )
      for (int lat =1; lat< Tempi->GetNbinsX(); lat++) {
         Exposure  -> SetBinContent (R+1,lat+1,Tempi -> GetBinContent (lat) );
         Exposure  -> SetBinError (R+1,lat+1,10);
      }
   return;
}



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
      Flux (std::string basename, std::string suffix, Binning B)
      {
         string hname=basename + "Counts_" + suffix;
         Counts	 = new TH1F ( hname.c_str(), hname.c_str(),B.size(),0,B.size());
         name = basename;
         suffixname = suffix;
         bins=Bin;
      }

      Flux (std::string basename, std::string suffix, Binning B, int n)
      {
         string hname=basename + "Counts_" + suffix;
         Counts   = new TH2F ( hname.c_str(), hname.c_str(), B.size(),  0, B.size(),  n,0,n);
         name = basename;
         suffixname = suffix;
         bins=Bin;
      }

      //reading constructor
      Flux (TFile * file, std::string basename, std::string suffix, std::string dirname, std::string acceptname, Binning B)
      {
         Counts	 = (TH1 *) file->Get ( (basename + "Counts_" + suffix       ).c_str() );
         Acceptance  = (TH1 *) file->Get ( ("/" + dirname +"/" +acceptname + "_" + suffix  ).c_str() );
         name = basename;
         suffixname = suffix;
         bins=Bin;
      }

      void Write();

      void Set_Exposure_Time (TH2 * ExposureR ,TH2 * ExposureTOF,TH2 * ExposureNaF,TH2 * ExposureAgl );
      void Set_Exposure_Time (TH1 * Tempi );
      TH1 * ExtractParticularMC_cs (TH1 * Histo, int lat_zones=1 ,int mc_type=0);
      void Set_DeltaE (int n,bool deutons=1);
      void Add_SystFitError (int n, TH1* syst_errR,TH1* syst_errTOF,TH1* syst_errNaF,TH1* syst_errAgl);
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
   TimeGeo ( (TH2 *) Exposure  ,Tempi);

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



void Flux::Set_DeltaE (int n,bool deutons)
{
   if (n>1) {
      string hname=name + "DeltaE_" + suffixname;
      DeltaE   = new TH2F (  hname.c_str(), hname.c_str(), bins.size(),  0, bins.size(),  n,0,n);
      DeltaE   ->Sumw2();


      if (deutons) {
         for (int iR=0; iR<DeltaE_R->GetNbinsX(); iR++)
            for (int lat =1; lat<DeltaE_R->GetNbinsX(); lat++) DeltaE_R->SetBinContent (iR+1,lat+1,deltaencindeut[iR]);
      } else { //only in R bins p/d have different DeltaE
         for (int iR=0; iR<DeltaE_R->GetNbinsX(); iR++)
            for (int lat =1; lat<DeltaE_R->GetNbinsX(); lat++) DeltaE_R->SetBinContent (iR+1,lat+1,deltaencinprot[iR]);
      }
      for (int lat =0; lat<DeltaE->GetNbinsX(); lat++)
         for (int iR=0; iR<DeltaE->GetNbinsX(); iR++)
            DeltaE->SetBinContent (iR+1,lat+1, bins.EkBin(i) - bins.EkBin(i-1));


   } else { // 1D
      DeltaE_R   = new TH1F ( (name + "DeltaE_R"   ).c_str(), (name + "DeltaE_R"   ).c_str(),nbinsr,  0,nbinsr  );


      DeltaE_R   ->Sumw2();



      for (int iR=1; iR<DeltaE_R->GetNbinsX(); iR++)  DeltaE_R->SetBinContent (iR+0,bins.EkBin(i) - bins.EkBin(i-1));

   }
}
void Flux::Add_SystFitError (int n, TH1* syst_errR,TH1* syst_errTOF,TH1* syst_errNaF,TH1* syst_errAgl)
{
   if (Counts_R  ) for (int iR = 0; iR < Counts_R  ->GetNbinsX(); iR++) Counts_R  -> SetBinError (iR+1,Counts_R  -> GetBinError (iR+1)+fabs (syst_errR->GetBinContent (iR+1) ) ) ;
   if (Counts_TOF) for (int iR = 0; iR < Counts_TOF  ->GetNbinsX(); iR++) Counts_TOF  -> SetBinError (iR+1,Counts_TOF  -> GetBinError (iR+1)+fabs (syst_errTOF->GetBinContent (iR+1) ) ) ;
   if (Counts_NaF) for (int iR = 0; iR < Counts_NaF  ->GetNbinsX(); iR++) Counts_NaF  -> SetBinError (iR+1,Counts_NaF  -> GetBinError (iR+1)+fabs (syst_errNaF->GetBinContent (iR+1) ) ) ;
   if (Counts_Agl) for (int iR = 0; iR < Counts_Agl  ->GetNbinsX(); iR++) Counts_Agl  -> SetBinError (iR+1,Counts_Agl  -> GetBinError (iR+1)+fabs (syst_errAgl->GetBinContent (iR+1) ) ) ;
}


void Flux::Eval_Flux (int n,bool deutons,int mc_type)
{
   Flux::Set_DeltaE (n,deutons);

   if (Counts_R  ) 	Flux_R	= (TH1 *) Counts_R   -> Clone();
   if (Counts_TOF)	Flux_TOF= (TH1 *) Counts_TOF -> Clone();
   if (Counts_NaF)	Flux_NaF= (TH1 *) Counts_NaF -> Clone();
   if (Counts_Agl)	Flux_Agl= (TH1 *) Counts_Agl -> Clone();

   if (deutons) {
      if (Counts_R  ) 	Flux_R	->Divide (ExtractParticularMC_cs ( Acceptance_R      ,n,mc_type) );
      if (Counts_TOF)	Flux_TOF->Divide (ExtractParticularMC_cs ( Acceptance_TOF    ,n,mc_type) );
      if (Counts_NaF)	Flux_NaF->Divide (ExtractParticularMC_cs ( Acceptance_NaF    ,n,mc_type) );
      if (Counts_Agl)	Flux_Agl->Divide (ExtractParticularMC_cs ( Acceptance_Agl    ,n,mc_type) );
   }

   else {
      if (Counts_R  ) 	Flux_R	->Divide (Acceptance_R  	);
      if (Counts_TOF)	Flux_TOF->Divide (Acceptance_TOF	);
      if (Counts_NaF)	Flux_NaF->Divide (Acceptance_NaF	);
      if (Counts_Agl)	Flux_Agl->Divide (Acceptance_Agl	);
   }

   if (Counts_R  )	Flux_R	->Divide (Exposure_R 		);
   if (Counts_TOF)  Flux_TOF->Divide (Exposure_TOF          );
   if (Counts_NaF)  Flux_NaF->Divide (Exposure_NaF          );
   if (Counts_Agl)  Flux_Agl->Divide (Exposure_Agl          );

   if (Counts_R  )	Flux_R	->Divide (DeltaE_R 		);
   if (Counts_TOF)  Flux_TOF->Divide (DeltaE_TOF            );
   if (Counts_NaF)  Flux_NaF->Divide (DeltaE_NaF            );
   if (Counts_Agl)  Flux_Agl->Divide (DeltaE_Agl            );

   return;
}






