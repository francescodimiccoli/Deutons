using namespace std;


class TemplateFIT
{
private:

	TH1 * ResultP;
	TH1 * ResultD;
	TH1 * ResultHe;
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
		
		ResultP  =	TemplateP ; 
	        ResultD  =      TemplateD ; 
        	ResultHe =      TemplateHe; 
		
		PCounts	   =	new TH1F((basename + "_PCounts"  	 ).c_str(),(basename + "_PCounts"   	).c_str(),Nbins,0,Nbins);
		DCounts	   =	new TH1F((basename + "_PCounts"  	 ).c_str(),(basename + "_PCounts"   	).c_str(),Nbins,0,Nbins);
        	PCountsgeo = 	new TH2F((basename + "_PCounts_geo"      ).c_str(),(basename + "_PCounts_geo"   ).c_str(),Nbins,0,Nbins,n,0,n);
        	DCountsgeo = 	new TH2F((basename + "_PCounts_geo"      ).c_str(),(basename + "_PCounts_geo"   ).c_str(),Nbins,0,Nbins,n,0,n);
	
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
		
		ResultP  =	TemplateP ; 
	        ResultD  =      TemplateD ; 
        	ResultHe =      TemplateHe; 	
	
		nbins = TemplateP -> GetNbinsY();
	
		PCounts	   =	new TH1F((basename + "_PCounts"  	 ).c_str(),(basename + "_PCounts"   	).c_str(),nbins,0,nbins);
		DCounts	   =	new TH1F((basename + "_PCounts"  	 ).c_str(),(basename + "_PCounts"   	).c_str(),nbins,0,nbins);
        	PCountsgeo = 	new TH2F((basename + "_PCounts_geo"      ).c_str(),(basename + "_PCounts_geo"   ).c_str(),nbins,0,nbins,n,0,n);
        	DCountsgeo = 	new TH2F((basename + "_PCounts_geo"      ).c_str(),(basename + "_PCounts_geo"   ).c_str(),nbins,0,nbins,n,0,n);

	}
	void Write();
	
	TH1F * Extract_Bin_histos(TH1 * Histo, int bin);

	TH1F * Extract_Bin_histos_geo(TH1 * Histo, int bin, int lat);
	
	TFractionFitter * Do_TemplateFIT(int bin, TH1F * PMC, TH1F *DMC, TH1F *HeMC, TH1F* Data);
	
	void TemplateFits();
	
	TH1F * GetResult_P (int bin){ return Extract_Bin_histos(ResultP ,bin);	}	
	TH1F * GetResult_D (int bin){ return Extract_Bin_histos(ResultD ,bin);	}
	TH1F * GetResult_He(int bin){ return Extract_Bin_histos(ResultHe,bin);	}

	TH1F * GetResult_Data(int bin){ return Extract_Bin_histos(Data_Prim,bin);  }
	TH1F * GetResult_Data(int bin,int lat){ return Extract_Bin_histos_geo(Data_Prim,bin,lat);  }
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
	return (TH1F *)((TH2*)Histo) -> ProjectionX ("",bin+1,bin+1);	
	    
}

TH1F * TemplateFIT::Extract_Bin_histos_geo(TH1 * Histo, int bin, int lat){
        return (TH1F *)((TH3*)Histo) -> ProjectionX ("",bin+1,bin+1,lat+1,lat+1) -> Clone();

}

TFractionFitter * TemplateFIT::Do_TemplateFIT(int bin, TH1F * PMC, TH1F *DMC, TH1F *HeMC, TH1F* Data){
	TObjArray *Tpl;
	Tpl -> Add(PMC);
	Tpl -> Add(DMC);
	Tpl -> Add(HeMC);
	TFractionFitter * fit = new TFractionFitter(Data,Tpl,"q");
	return fit;
}



void TemplateFIT::TemplateFits(){
	for(int bin=0; bin<nbins ; bin++){
		TH1F * Templ_P =  TemplateFIT::Extract_Bin_histos(TemplateP, bin);	
		TH1F * Templ_D =  TemplateFIT::Extract_Bin_histos(TemplateD, bin);	
		TH1F * Templ_He=  TemplateFIT::Extract_Bin_histos(TemplateHe,bin);
		TH1F * Data    =  TemplateFIT::Extract_Bin_histos(TemplateHe,bin);

		TFractionFitter * fit = TemplateFIT::Do_TemplateFIT(bin, Templ_P, Templ_D, Templ_He, Data);

		Templ_P ->Scale(1);
		Templ_D ->Scale(1);
		Templ_He->Scale(1);	

		for(int R=0;R<Templ_P->GetNbinsX();R++){
			ResultP ->SetBinContent(R+1,bin +1,Templ_P -> GetBinContent(R+1));	
			ResultD ->SetBinContent(R+1,bin +1,Templ_P -> GetBinContent(R+1));	
			ResultHe->SetBinContent(R+1,bin +1,Templ_P -> GetBinContent(R+1));
		}
		PCounts -> SetBinContent(bin+1,Templ_P->Integral());
		DCounts -> SetBinContent(bin+1,Templ_D->Integral());
	}
	return;
}



