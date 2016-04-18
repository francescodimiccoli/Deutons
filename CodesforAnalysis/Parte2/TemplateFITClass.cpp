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
	
	double GetFitWheights(int par, int bin);
	
	int GetFitOutcome(int bin){if(fits[bin]) return fits[bin]->Tfit_outcome; else {cout<<"Fit not yet performed: bin nr. "<<bin<<endl; return -1;}}

	TH1F * GetResult_P (int bin){ TH1F *res =(TH1F*)fits[bin] -> Templ_P -> Clone() ; res -> Scale(GetFitWheights(0,bin)); return res;}	
	TH1F * GetResult_D (int bin){ TH1F *res =(TH1F*)fits[bin] -> Templ_D -> Clone() ; res -> Scale(GetFitWheights(1,bin)); return res;}
	TH1F * GetResult_He(int bin){ TH1F *res =(TH1F*)fits[bin] -> Templ_He-> Clone() ; res -> Scale(GetFitWheights(2,bin)); return res;}

	TH1F * GetResult_Data(int bin)	      { return fits[bin] -> Data; }
	
	void TemplateFits();
	
};

void TemplateFIT::Write(){

	if(TemplateP  -> GetEntries() > 0)  TemplateP  -> Write();
        if(TemplateD  -> GetEntries() > 0)  TemplateD  -> Write();
        TemplateHe -> Write();
        
	if(Data_Prim  -> GetEntries() > 0)  Data_Prim  -> Write();
        if(Data_Geomag-> GetEntries() > 0)  Data_Geomag-> Write();
	return;
}


TH1F * TemplateFIT::Extract_Bin_histos(TH1 * Histo, int bin){
	TH1F * Slice = new TH1F("","",Histo->GetNbinsX(),0,Histo->GetXaxis()->GetBinLowEdge(101));
	for(int i = 0; i< Histo->GetNbinsX();i++)
		Slice->SetBinContent(i+1,Histo->GetBinContent(i+1,bin+1));
	return Slice;
}

TH1F * TemplateFIT::Extract_Bin_histos_geo(TH1 * Histo, int bin, int lat){
        TH1F * Slice = (TH1F *)((TH3F*)Histo) -> ProjectionX ("",bin+1,bin+1,lat+1,lat+1) -> Clone();
        return Slice;
}

void TemplateFIT::Do_TemplateFIT(TFit * Fit){
	TObjArray *Tpl;
	Tpl = new TObjArray(3);
	Tpl -> Add( Fit ->  Templ_P );
	Tpl -> Add( Fit ->  Templ_D );
	Tpl -> Add( Fit ->  Templ_He);
	
	Fit -> Tfit = new TFractionFitter(Fit -> Data, Tpl ,"q");
	Fit -> Tfit_outcome = 1;//fit -> Fit();
	fits.push_back(Fit);
	
	return;
}

double TemplateFIT::GetFitWheights(int par, int bin){
	if(GetFitOutcome(bin)==-1) return  1;
	if(GetFitOutcome(bin)>0)   return  1;
	if(GetFitOutcome(bin)==0){
		double w1,e1=0;
		fits[bin]-> Tfit ->GetResult(par,w1,e1);
		TH1F * Result = (TH1F*)fits[bin] -> Tfit -> GetPlot();
		float itot= Result->Integral();
                float i1;
		if(par == 0) i1 = fits[bin]-> Templ_P ->Integral();
		if(par == 1) i1 = fits[bin]-> Templ_D ->Integral();
		if(par == 2) i1 = fits[bin]-> Templ_He ->Integral();
		return w1/(i1*itot); 
		}	
}

void TemplateFIT::TemplateFits(){
	
	
	for(int bin=0; bin<nbins ; bin++){
		TFit * Fit = new TFit;
		Fit->Templ_P =  (TH1F *)TemplateFIT::Extract_Bin_histos(TemplateP, bin);	
		Fit->Templ_D =  (TH1F *)TemplateFIT::Extract_Bin_histos(TemplateD, bin);	
		Fit->Templ_He=  (TH1F *)TemplateFIT::Extract_Bin_histos(TemplateHe,bin);
		Fit->Data    =  (TH1F *)TemplateFIT::Extract_Bin_histos(Data_Prim ,bin);
		TemplateFIT::Do_TemplateFIT(Fit);

		TH1F * ResultPlot_P  = GetResult_P (bin);		
		TH1F * ResultPlot_D  = GetResult_D (bin);
		TH1F * ResultPlot_He = GetResult_He(bin);
		
		PCounts -> SetBinContent(bin+1,ResultPlot_P->Integral());
		DCounts -> SetBinContent(bin+1,ResultPlot_D->Integral());
	}
	return;
}



