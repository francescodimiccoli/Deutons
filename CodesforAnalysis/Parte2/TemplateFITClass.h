using namespace std;


class TemplateFIT
{
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

	//creation constructors
	TemplateFIT(std::string basename , float val_min , float val_max,int n){

		TemplateP   = 	new TH2F((basename + "_P"   ).c_str(),(basename + "_P"   ).c_str(),100,val_min,val_max,18,0,18);
                TemplateD   =	new TH2F((basename + "_D"   ).c_str(),(basename + "_D"   ).c_str(),100,val_min,val_max,18,0,18);
                TemplateHe  =	new TH2F((basename + "_He"  ).c_str(),(basename + "_He"  ).c_str(),100,val_min,val_max,18,0,18);
	
		Data_Prim   =	new TH2F((basename + "_Data_Prim"  ).c_str(),(basename + "_Data_Prim"  ).c_str(),100,val_min,val_max,18,0,18);
		Data_Geomag =	new TH3F((basename + "_Data_Geomag").c_str(),(basename + "_Data_Geomag").c_str(),100,val_min,val_max,18,0,18,n,0,n);		
	}

	//reading constructor
	TemplateFIT(TFile * file , std::string basename , float val_min , float val_max,int n){

		TemplateP   =	(TH1 *)file->Get((basename + "_P"     ).c_str());	
                TemplateD   =	(TH1 *)file->Get((basename + "_D"     ).c_str());
                TemplateHe  =	(TH1 *)file->Get((basename + "_He"    ).c_str());
	                     
                Data_Prim   =	(TH1 *)file->Get((basename + "_Data_Prim"     ).c_str());
	        Data_Geomag =	(TH1 *)file->Get((basename + "_Data_Geomag"   ).c_str());
	}

	void Write();


};

void TemplateFIT::Write(){

	if(TemplateP  -> GetEntries() > 0)  TemplateP  -> Write();
        if(TemplateD  -> GetEntries() > 0)  TemplateD  -> Write();
        if(TemplateHe -> GetEntries() > 0)  TemplateHe -> Write();
                                                     
        if(Data_Prim  -> GetEntries() > 0)  Data_Prim  -> Write();
        if(Data_Geomag-> GetEntries() > 0)  Data_Geomag-> Write();
	return;
}
