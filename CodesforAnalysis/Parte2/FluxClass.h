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

	std::string name;	
	//creation constructor
	Flux(std::string basename){
		Counts_R	= new TH1F((basename + "Counts_R"      ).c_str(),(basename + "Counts_R"      ).c_str(),43,0,43);
		Counts_TOF      = new TH1F((basename + "Counts_TOF"    ).c_str(),(basename + "Counts_TOF"    ).c_str(),18,0,18);
		Counts_NaF      = new TH1F((basename + "Counts_NaF"    ).c_str(),(basename + "Counts_NaF"    ).c_str(),18,0,18);
		Counts_Agl      = new TH1F((basename + "Counts_Agl"    ).c_str(),(basename + "Counts_Agl"    ).c_str(),18,0,18);
		name = basename;
	}
	
	Flux(std::string basename,int n){
		Counts_R    = new TH2F((basename + "Counts_R"   ).c_str(),(basename + "Counts_R"   ).c_str(),43,0,43,n,0,n);
		Counts_TOF  = new TH2F((basename + "Counts_TOF" ).c_str(),(basename + "Counts_TOF" ).c_str(),18,0,18,n,0,n);
		Counts_NaF  = new TH2F((basename + "Counts_NaF" ).c_str(),(basename + "Counts_NaF" ).c_str(),18,0,18,n,0,n);
		Counts_Agl  = new TH2F((basename + "Counts_Agl" ).c_str(),(basename + "Counts_Agl" ).c_str(),18,0,18,n,0,n);	
		name = basename;
	}

	//reading constructor
	Flux(TFile * file, std::string basename, std::string dirname, std::string acceptname,int n){
		Counts_R	=(TH1 *)file->Get((basename + "Counts_R"        ).c_str());
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

};


void Flux::Write()
{
	if(Counts_R       ->GetEntries()>0) Counts_R       ->Write();       
	if(Counts_TOF     ->GetEntries()>0) Counts_TOF     ->Write(); 
	if(Counts_NaF     ->GetEntries()>0) Counts_NaF     ->Write(); 
	if(Counts_Agl     ->GetEntries()>0) Counts_Agl     ->Write(); 

	return;
}

void Flux::Set_Exposure_Time(TH2 * ExposureR ,TH2 * ExposureTOF,TH2 * ExposureNaF,TH2 * ExposureAgl ){

	Exposure_R  = ExposureTime(ExposureR  ); 
	Exposure_TOF= ExposureTime(ExposureTOF);
	Exposure_NaF= ExposureTime(ExposureNaF);
	Exposure_Agl= ExposureTime(ExposureAgl);

	return;
}

void Flux::Set_Exposure_Time(TH1 * Tempi){
		
	Exposure_R  = new TH2F ((name + "TimeZone_R"   ).c_str(),(name + "TimeZone_R"   ).c_str(),43,0,43,Tempi->GetNbinsX(),0,Tempi->GetNbinsX());
	Exposure_TOF= new TH2F ((name + "TimeZone_TOF" ).c_str(),(name + "TimeZone_TOF" ).c_str(),18,0,18,Tempi->GetNbinsX(),0,Tempi->GetNbinsX());
	Exposure_NaF= new TH2F ((name + "TimeZone_NaF" ).c_str(),(name + "TimeZone_NaF" ).c_str(),18,0,18,Tempi->GetNbinsX(),0,Tempi->GetNbinsX());
	Exposure_Agl= new TH2F ((name + "TimeZone_Agl" ).c_str(),(name + "TimeZone_Agl" ).c_str(),18,0,18,Tempi->GetNbinsX(),0,Tempi->GetNbinsX());	

	TimeGeo((TH2 *)Exposure_R  ,Tempi);
	TimeGeo((TH2 *)Exposure_TOF,Tempi);
	TimeGeo((TH2 *)Exposure_NaF,Tempi);
	TimeGeo((TH2 *)Exposure_Agl,Tempi);
	
	return;
}










TH1 * ExposureTime(TH2 * esposizionegeo) {
	return (TH1 *) esposizionegeo -> ProjectionX("",0,10) -> Clone();
}

void TimeGeo(TH2 * Exposure, TH1 * Tempi){
	for(int R=0;R<Exposure ->GetNbinsX() ; R++ ) 
			for(int lat =0; lat< Tempi->GetNbinsX();lat++) 
				Exposure  -> SetBinContent(R+1,lat+1,Tempi -> GetBinContent(lat));			

	return;
}
