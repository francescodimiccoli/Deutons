	TH2F * esposizionegeo_R    = new TH2F("esposizionegeo_R"  ,"esposizionegeo_R"  ,43,0,43,11,0,11);
	TH2F * esposizionedgeoTOF  = new TH2F("esposizionedgeoTOF","esposizionedgeoTOF",18,0,18,11,0,11);
	TH2F * esposizionedgeoNaF  = new TH2F("esposizionedgeoNaF","esposizionedgeoNaF",18,0,18,11,0,11);
	TH2F * esposizionedgeoAgl  = new TH2F("esposizionedgeoAgl","esposizionedgeoAgl",18,0,18,11,0,11);
	TH2F * esposizionepgeoTOF  = new TH2F("esposizionepgeoTOF","esposizionepgeoTOF",18,0,18,11,0,11);
	TH2F * esposizionepgeoNaF  = new TH2F("esposizionepgeoNaF","esposizionepgeoNaF",18,0,18,11,0,11);
	TH2F * esposizionepgeoAgl  = new TH2F("esposizionepgeoAgl","esposizionepgeoAgl",18,0,18,11,0,11);

	TH1F * Tempi = new TH1F("Tempi","Tempi",11,0,11);

void UpdateZoneLivetime (Binning bins, int zonageo, TH2F * esposizionegeo);

 void ExposureTime_Fill(float Zona){
	
	if((int)Tup.U_time!=ActualTime) {
		Tempi -> SetBinContent((int)Zona, Tempi -> GetBinContent((int)Zona)+Tup.Livetime );
	//protons exposure time	
	UpdateZoneLivetime( PRB,  (int)Zona,esposizionegeo_R );
	UpdateZoneLivetime( ToFPB,(int)Zona,esposizionepgeoTOF );
	UpdateZoneLivetime( NaFPB,(int)Zona,esposizionepgeoNaF );
	UpdateZoneLivetime( AglPB,(int)Zona,esposizionepgeoAgl );
	
	//deutone exposure time
	UpdateZoneLivetime( ToFDB,(int)Zona,esposizionedgeoTOF );
	UpdateZoneLivetime( NaFDB,(int)Zona,esposizionedgeoNaF );
	UpdateZoneLivetime( AglDB,(int)Zona,esposizionedgeoAgl );
	
	ActualTime = (int)Tup.U_time;
	}
	return;
}	
	

void ExposureTime_Write(){

	esposizionegeo_R    -> Write();
	esposizionedgeoTOF  -> Write();
	esposizionedgeoNaF  -> Write();
	esposizionedgeoAgl  -> Write();
	esposizionepgeoTOF  -> Write();
	esposizionepgeoNaF  -> Write();
	esposizionepgeoAgl  -> Write();
	Tempi -> Write();
	return;

}



void UpdateZoneLivetime (Binning bins, int zonageo, TH2F * esposizionegeo){

	for(int i=0;i<bins.RigBins().size()-1;i++)
                        if(bins.RigBins()[i+1]>=SF*Tup.Rcutoff){
                                esposizionegeo -> SetBinContent(i+1, zonageo, esposizionegeo -> GetBinContent(i+1, zonageo) + Tup.Livetime) ;
	}
	return;

}

