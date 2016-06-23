using namespace std;

TH1 * ExposureTime(TH2 * esposizionegeo);
void TimeGeo(TH2 * Exposure, TH1 * Tempi);

class Flux
{
   private:
      //Delta E;
      TH1 * DeltaE_R  ;
      TH1 * DeltaE_TOF;
      TH1 * DeltaE_NaF;
      TH1 * DeltaE_Agl;

   public:
      //counts
      TH1 * Counts_R	;
      TH1 * Counts_TOF;
      TH1 * Counts_NaF;
      TH1 * Counts_Agl;

      //Acceptance
      TH1 * Acceptance_R   ;
      TH1 * Acceptance_TOF ;
      TH1 * Acceptance_NaF ;
      TH1 * Acceptance_Agl ;

      //Exposure Time
      TH1 * Exposure_R  ;
      TH1 * Exposure_TOF;
      TH1 * Exposure_NaF;
      TH1 * Exposure_Agl;

      //Fluxes
      TH1 * Flux_R  ;
      TH1 * Flux_TOF;
      TH1 * Flux_NaF;
      TH1 * Flux_Agl;

      std::string name;


      //creation constructor
      Flux(std::string basename) {
         Counts_R	 = new TH1F((basename + "Counts_R"      ).c_str(),(basename + "Counts_R"      ).c_str(),nbinsr,0,nbinsr);
         Counts_TOF      = new TH1F((basename + "Counts_TOF"    ).c_str(),(basename + "Counts_TOF"    ).c_str(),nbinsToF,0,nbinsToF);
         Counts_NaF      = new TH1F((basename + "Counts_NaF"    ).c_str(),(basename + "Counts_NaF"    ).c_str(),nbinsNaF,0,nbinsNaF);
         Counts_Agl      = new TH1F((basename + "Counts_Agl"    ).c_str(),(basename + "Counts_Agl"    ).c_str(),nbinsAgl,0,nbinsAgl);
         name = basename;
      }

      Flux(std::string basename,int n) {
         Counts_R    = new TH2F((basename + "Counts_R"   ).c_str(),(basename + "Counts_R"   ).c_str(),nbinsr,  0,nbinsr,  n,0,n);
         Counts_TOF  = new TH2F((basename + "Counts_TOF" ).c_str(),(basename + "Counts_TOF" ).c_str(),nbinsToF,0,nbinsToF,n,0,n);
         Counts_NaF  = new TH2F((basename + "Counts_NaF" ).c_str(),(basename + "Counts_NaF" ).c_str(),nbinsNaF,0,nbinsNaF,n,0,n);
         Counts_Agl  = new TH2F((basename + "Counts_Agl" ).c_str(),(basename + "Counts_Agl" ).c_str(),nbinsAgl,0,nbinsAgl,n,0,n);
         name = basename;
      }

      //reading constructor
      Flux(TFile * file, std::string basename, std::string dirname, std::string acceptname,int n) {
         Counts_R	 =(TH1 *)file->Get((basename + "Counts_R"        ).c_str());
         Counts_TOF      =(TH1 *)file->Get((basename + "Counts_TOF"      ).c_str());
         Counts_NaF      =(TH1 *)file->Get((basename + "Counts_NaF"      ).c_str());
         Counts_Agl      =(TH1 *)file->Get((basename + "Counts_Agl"      ).c_str());

         Acceptance_R  =(TH1 *)file->Get(("/" + dirname +"/" +acceptname + "_R"   ).c_str());
         Acceptance_TOF=(TH1 *)file->Get(("/" + dirname +"/" +acceptname + "_TOF" ).c_str());
         Acceptance_NaF=(TH1 *)file->Get(("/" + dirname +"/" +acceptname + "_NaF" ).c_str());
         Acceptance_Agl=(TH1 *)file->Get(("/" + dirname +"/" +acceptname + "_Agl" ).c_str());
         name = basename;
      }

      void Write();

      void Set_Exposure_Time(TH2 * ExposureR ,TH2 * ExposureTOF,TH2 * ExposureNaF,TH2 * ExposureAgl );
      void Set_Exposure_Time(TH1 * Tempi );
      TH1 * ExtractParticularMC_cs(TH1 * Histo, int lat_zones=1 ,int mc_type=0);
      void Set_DeltaE(int n,bool deutons=1);
      void Add_SystFitError(int n, TH1* syst_errR,TH1* syst_errTOF,TH1* syst_errNaF,TH1* syst_errAgl);
      void Eval_Flux(int n,bool deutons=1,int mc_type=0);
};


void Flux::Write()
{
   if(Counts_R       ->GetEntries()>0) Counts_R       ->Write();
   if(Counts_TOF     ->GetEntries()>0) Counts_TOF     ->Write();
   if(Counts_NaF     ->GetEntries()>0) Counts_NaF     ->Write();
   if(Counts_Agl     ->GetEntries()>0) Counts_Agl     ->Write();

   return;
}

void Flux::Set_Exposure_Time(TH2 * ExposureR ,TH2 * ExposureTOF,TH2 * ExposureNaF,TH2 * ExposureAgl ) {

   Exposure_R  = ExposureTime(ExposureR  );
   Exposure_TOF= ExposureTime(ExposureTOF);
   Exposure_NaF= ExposureTime(ExposureNaF);
   Exposure_Agl= ExposureTime(ExposureAgl);

   return;
}

void Flux::Set_Exposure_Time(TH1 * Tempi) {

   Exposure_R  = new TH2F ((name + "TimeZone_R"   ).c_str(),(name + "TimeZone_R"   ).c_str(), nbinsr,  0,nbinsr,   Tempi->GetNbinsX(),0,Tempi->GetNbinsX());
   Exposure_TOF= new TH2F ((name + "TimeZone_TOF" ).c_str(),(name + "TimeZone_TOF" ).c_str(), nbinsToF,0,nbinsToF, Tempi->GetNbinsX(),0,Tempi->GetNbinsX());
   Exposure_NaF= new TH2F ((name + "TimeZone_NaF" ).c_str(),(name + "TimeZone_NaF" ).c_str(), nbinsNaF,0,nbinsNaF, Tempi->GetNbinsX(),0,Tempi->GetNbinsX());
   Exposure_Agl= new TH2F ((name + "TimeZone_Agl" ).c_str(),(name + "TimeZone_Agl" ).c_str(), nbinsAgl,0,nbinsAgl, Tempi->GetNbinsX(),0,Tempi->GetNbinsX());

   TimeGeo((TH2 *)Exposure_R  ,Tempi);
   TimeGeo((TH2 *)Exposure_TOF,Tempi);
   TimeGeo((TH2 *)Exposure_NaF,Tempi);
   TimeGeo((TH2 *)Exposure_Agl,Tempi);

   return;
}

TH1 * Flux::ExtractParticularMC_cs(TH1 * Histo, int lat_zones, int mc_type){
	if(lat_zones == 1){
		TH1F * Slice = new TH1F("","",Histo->GetNbinsX(),0,Histo->GetNbinsX());
		for(int i = 0; i< Histo->GetNbinsX();i++){
			Slice->SetBinContent(i+1,Histo->GetBinContent(i+1,mc_type+1));
			Slice->SetBinError(i+1,Histo->GetBinError(i+1,mc_type+1));
			}
		return Slice;
	}
	else{
		TH2F * Slice = new TH2F("","",Histo->GetNbinsX(),0,Histo->GetNbinsX(),lat_zones,0,lat_zones);
		for(int lat=0;lat<lat_zones;lat++){
			for(int i = 0; i< Histo->GetNbinsX();i++){
				Slice->SetBinContent(i+1,lat+1,Histo->GetBinContent(i+1,lat+1,mc_type+1));
				Slice->SetBinError(i+1,lat+1,Histo->GetBinError(i+1,lat+1,mc_type+1));
				}
		}
		return Slice;	
	}
}



void Flux::Set_DeltaE(int n,bool deutons) {
   if(n>1) {
      DeltaE_R   = new TH2F ((name + "DeltaE_R"   ).c_str(),(name + "DeltaE_R"   ).c_str(),nbinsr,  0,nbinsr,  n,0,n);
      DeltaE_TOF = new TH2F ((name + "DeltaE_TOF" ).c_str(),(name + "DeltaE_TOF" ).c_str(),nbinsToF,0,nbinsToF,n,0,n);
      DeltaE_NaF = new TH2F ((name + "DeltaE_NaF" ).c_str(),(name + "DeltaE_NaF" ).c_str(),nbinsNaF,0,nbinsNaF,n,0,n);
      DeltaE_Agl = new TH2F ((name + "DeltaE_Agl" ).c_str(),(name + "DeltaE_Agl" ).c_str(),nbinsAgl,0,nbinsAgl,n,0,n);

	DeltaE_R   ->Sumw2();
	DeltaE_TOF ->Sumw2();
	DeltaE_NaF ->Sumw2();
	DeltaE_Agl ->Sumw2();	
      
    if(deutons) {
         for(int iR=0; iR<DeltaE_R->GetNbinsX(); iR++)
            for(int lat =1; lat<DeltaE_R->GetNbinsX(); lat++) DeltaE_R->SetBinContent(iR+1,lat+1, DRB.EkBinCent(iR+1)-DRB.EkBinCent(iR));
      }
      else { //only in R bins p/d have different DeltaE
         for(int iR=0; iR<DeltaE_R->GetNbinsX(); iR++)
            for(int lat =1; lat<DeltaE_R->GetNbinsX(); lat++) DeltaE_R->SetBinContent(iR+1,lat+1, PRB.EkBinCent(iR+1)-PRB.EkBinCent(iR));
      }
      for(int iR=0; iR<DeltaE_TOF->GetNbinsX(); iR++)
         for(int lat =0; lat<DeltaE_TOF->GetNbinsX(); lat++) DeltaE_TOF->SetBinContent(iR+1,lat+1, ToFDB.EkBinCent(iR+1)-ToFPB.EkBinCent(iR));
      for(int iR=0; iR<DeltaE_NaF->GetNbinsX(); iR++)
         for(int lat =0; lat<DeltaE_NaF->GetNbinsX(); lat++) DeltaE_NaF->SetBinContent(iR+1,lat+1, NaFDB.EkBinCent(iR+1)-NaFPB.EkBinCent(iR));
      for(int iR=0; iR<DeltaE_Agl->GetNbinsX(); iR++)
         for(int lat =0; lat<DeltaE_TOF->GetNbinsX(); lat++) DeltaE_Agl->SetBinContent(iR+1,lat+1, AglDB.EkBinCent(iR+1)-AglDB.EkBinCent(iR));

   }
   else {
      DeltaE_R   = new TH1F ((name + "DeltaE_R"   ).c_str(),(name + "DeltaE_R"   ).c_str(),nbinsr,  0,nbinsr  );
      DeltaE_TOF = new TH1F ((name + "DeltaE_TOF" ).c_str(),(name + "DeltaE_TOF" ).c_str(),nbinsToF,0,nbinsToF);
      DeltaE_NaF = new TH1F ((name + "DeltaE_NaF" ).c_str(),(name + "DeltaE_NaF" ).c_str(),nbinsNaF,0,nbinsNaF);
      DeltaE_Agl = new TH1F ((name + "DeltaE_Agl" ).c_str(),(name + "DeltaE_Agl" ).c_str(),nbinsAgl,0,nbinsAgl);
      
 	DeltaE_R   ->Sumw2();     
	 DeltaE_TOF ->Sumw2();
	 DeltaE_NaF ->Sumw2();
	 DeltaE_Agl ->Sumw2();

      if(deutons) {
         for(int iR=0; iR<DeltaE_R->GetNbinsX(); iR++)  DeltaE_R->SetBinContent(iR+1,DRB.EkBin(iR+1)-DRB.EkBin(iR));
      }
      else { //only in R bins p/d have different DeltaE
         for(int iR=0; iR<DeltaE_R->GetNbinsX(); iR++)  DeltaE_R->SetBinContent(iR+1,PRB.EkBin(iR+1)-PRB.EkBin(iR));
      }
      for(int iR=0; iR<DeltaE_TOF->GetNbinsX(); iR++) DeltaE_TOF->SetBinContent(iR+1,ToFPB.EkBin(iR+1)-ToFPB.EkBin(iR));
      for(int iR=0; iR<DeltaE_NaF->GetNbinsX(); iR++) DeltaE_NaF->SetBinContent(iR+1,NaFPB.EkBin(iR+1)-NaFPB.EkBin(iR));
      for(int iR=0; iR<DeltaE_Agl->GetNbinsX(); iR++) DeltaE_Agl->SetBinContent(iR+1, AglPB.EkBin(iR+1)-AglPB.EkBin(iR));

   }
}
void Flux::Add_SystFitError(int n, TH1* syst_errR,TH1* syst_errTOF,TH1* syst_errNaF,TH1* syst_errAgl){
		if(Counts_R  ) for(int iR = 0; iR < Counts_R  ->GetNbinsX();iR++) Counts_R  -> SetBinError(iR+1,Counts_R  -> GetBinError(iR+1)+fabs(syst_errR->GetBinContent(iR+1)) ) ;	
		if(Counts_TOF) for(int iR = 0; iR < Counts_TOF  ->GetNbinsX();iR++) Counts_TOF  -> SetBinError(iR+1,Counts_TOF  -> GetBinError(iR+1)+fabs(syst_errTOF->GetBinContent(iR+1)) ) ;		
		if(Counts_NaF) for(int iR = 0; iR < Counts_NaF  ->GetNbinsX();iR++) Counts_NaF  -> SetBinError(iR+1,Counts_NaF  -> GetBinError(iR+1)+fabs(syst_errNaF->GetBinContent(iR+1)) ) ;
		if(Counts_Agl) for(int iR = 0; iR < Counts_Agl  ->GetNbinsX();iR++) Counts_Agl  -> SetBinError(iR+1,Counts_Agl  -> GetBinError(iR+1)+fabs(syst_errAgl->GetBinContent(iR+1)) ) ;
}


void Flux::Eval_Flux(int n,bool deutons,int mc_type) {
	Flux::Set_DeltaE(n,deutons);

	if(Counts_R  ) 	Flux_R	= (TH1 *) Counts_R   -> Clone();
	if(Counts_TOF)	Flux_TOF= (TH1 *) Counts_TOF -> Clone();
	if(Counts_NaF)	Flux_NaF= (TH1 *) Counts_NaF -> Clone();
	if(Counts_Agl)	Flux_Agl= (TH1 *) Counts_Agl -> Clone();

	if(deutons){	
		if(Counts_R  ) 	Flux_R	->Divide (ExtractParticularMC_cs( Acceptance_R      ,n,mc_type));
		if(Counts_TOF)	Flux_TOF->Divide (ExtractParticularMC_cs( Acceptance_TOF    ,n,mc_type));	
		if(Counts_NaF)	Flux_NaF->Divide (ExtractParticularMC_cs( Acceptance_NaF    ,n,mc_type));
		if(Counts_Agl)	Flux_Agl->Divide (ExtractParticularMC_cs( Acceptance_Agl    ,n,mc_type));
	}
	
	else{
		if(Counts_R  ) 	Flux_R	->Divide (Acceptance_R  	);
		if(Counts_TOF)	Flux_TOF->Divide (Acceptance_TOF	);
		if(Counts_NaF)	Flux_NaF->Divide (Acceptance_NaF	);
		if(Counts_Agl)	Flux_Agl->Divide (Acceptance_Agl	);
	}

	if(Counts_R  )	Flux_R	->Divide (Exposure_R 		);
	if(Counts_TOF)  Flux_TOF->Divide (Exposure_TOF          );
	if(Counts_NaF)  Flux_NaF->Divide (Exposure_NaF          );
	if(Counts_Agl)  Flux_Agl->Divide (Exposure_Agl          );

	if(Counts_R  )	Flux_R	->Divide (DeltaE_R 		);
	if(Counts_TOF)  Flux_TOF->Divide (DeltaE_TOF            );
	if(Counts_NaF)  Flux_NaF->Divide (DeltaE_NaF            );
	if(Counts_Agl)  Flux_Agl->Divide (DeltaE_Agl            );

	return;
}


TH1 * ExposureTime(TH2 * esposizionegeo) {
   return (TH1 *) esposizionegeo -> ProjectionX("",0,10) -> Clone();
}

void TimeGeo(TH2 * Exposure, TH1 * Tempi) {
	//HistInfo(Tempi);
	for(int R=0; R<Exposure ->GetNbinsX() ; R++ )
		for(int lat =1; lat< Tempi->GetNbinsX(); lat++){
			Exposure  -> SetBinContent(R+1,lat+1,Tempi -> GetBinContent(lat));
			Exposure  -> SetBinError(R+1,lat+1,10);
		}	
	return;
}
