#include "PlottingFunctions/DeutonFlux_Plot.h"

using namespace std;


void DeutonFlux(string filename) {
	
	cout<<"*************** DEUTONS FLUXES CALCULATION *******************"<<endl;

        cout<<"*** Reading  P1 file ****"<<endl;
        TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");


	TH2F * esposizionegeo_R    = (TH2F*)inputHistoFile->Get(      "esposizionegeo"        );
	TH2F * esposizionedgeoTOF  = (TH2F*)inputHistoFile->Get(      "esposizionedgeo"       );
	TH2F * esposizionedgeoNaF  = (TH2F*)inputHistoFile->Get(      "esposizionedgeoNaF"    );
	TH2F * esposizionedgeoAgl  = (TH2F*)inputHistoFile->Get(      "esposizionedgeoAgl"    );
	TH2F * esposizionepgeoTOF  = (TH2F*)inputHistoFile->Get(      "esposizionepgeo"       );
	TH2F * esposizionepgeoNaF  = (TH2F*)inputHistoFile->Get(      "esposizionepgeoNaF"    );
	TH2F * esposizionepgeoAgl  = (TH2F*)inputHistoFile->Get(      "esposizionepgeoAgl"    );



	Flux * D_Flux         = new Flux(inputHistoFile, "D_Flux"         ,"Results","Corr_AcceptanceD",1);
	Flux * D_Flux_geo     = new Flux(inputHistoFile, "D_Flux_geo"     ,"Results","Geomag_AcceptanceD",11);

	Flux * D_Flux_Dist         = new Flux(inputHistoFile, "D_Flux_Dist"         ,"Results","Corr_AcceptanceD",1);
	Flux * D_Flux_geo_Dist     = new Flux(inputHistoFile, "D_Flux_geo_Dist"     ,"Results","Geomag_AcceptanceD",11);

	Flux * P_Flux         = new Flux(inputHistoFile, "P_Flux"         ,"Results","Corr_AcceptanceP",1);
	Flux * P_Flux_Dist    = new Flux(inputHistoFile, "P_Flux_Dist"    ,"Results","Corr_AcceptanceP",1);

	cout<<"*************** DEUTONS FLUXES CALCULATION *******************"<<endl;
	
	//Evaluation additional Fit Error
	
	TH1F * SystR  ;
        TH1F * SystTOF;
        TH1F * SystNaF;
        TH1F * SystAgl;
	
	if(D_Flux -> Counts_R  ) SystR   = (TH1F*) D_Flux -> Counts_R   ->Clone();
	if(D_Flux -> Counts_TOF) SystTOF = (TH1F*) D_Flux -> Counts_TOF ->Clone();	
	if(D_Flux -> Counts_NaF) SystNaF = (TH1F*) D_Flux -> Counts_NaF ->Clone();
	if(D_Flux -> Counts_Agl) SystAgl = (TH1F*) D_Flux -> Counts_Agl ->Clone();
	
	if(D_Flux -> Counts_R  )SystR   -> Add(	(TH1F*) D_Flux_Dist -> Counts_R  ->Clone()	,-1); 
        if(D_Flux -> Counts_TOF)SystTOF -> Add(	(TH1F*) D_Flux_Dist -> Counts_TOF->Clone()	,-1);
        if(D_Flux -> Counts_NaF)SystNaF -> Add(	(TH1F*) D_Flux_Dist -> Counts_NaF->Clone()	,-1);
        if(D_Flux -> Counts_Agl)SystAgl -> Add(	(TH1F*) D_Flux_Dist -> Counts_Agl->Clone()	,-1);

	D_Flux      -> Add_SystFitError(1,SystR,SystTOF,SystNaF,SystAgl);
	D_Flux_Dist -> Add_SystFitError(1,SystR,SystTOF,SystNaF,SystAgl);


	 cout<<"*************** DEUTONS FLUXES CALCULATION *******************"<<endl;

	//

	D_Flux         -> Set_Exposure_Time (esposizionegeo_R,esposizionedgeoTOF,esposizionedgeoNaF,esposizionedgeoAgl);
	D_Flux_geo     -> Set_Exposure_Time (Tempi);

	D_Flux_Dist    -> Set_Exposure_Time (esposizionegeo_R,esposizionedgeoTOF,esposizionedgeoNaF,esposizionedgeoAgl);
	D_Flux_geo_Dist-> Set_Exposure_Time (Tempi);

	P_Flux         -> Set_Exposure_Time (esposizionegeo_R,esposizionepgeoTOF,esposizionepgeoNaF,esposizionepgeoAgl);
	P_Flux_Dist    -> Set_Exposure_Time (esposizionegeo_R,esposizionepgeoTOF,esposizionepgeoNaF,esposizionepgeoAgl);

	enum {protons, deutons};

	D_Flux         -> Eval_Flux(1 , deutons, 2 );
	D_Flux_geo     -> Eval_Flux(11, deutons, 2 );

	D_Flux_Dist    -> Eval_Flux(1 , deutons, 2 );
	D_Flux_geo_Dist-> Eval_Flux(11, deutons, 2 );

	P_Flux         -> Eval_Flux(1 , protons );
	P_Flux_Dist    -> Eval_Flux(1 , protons );

	 cout<<"*************** DEUTONS FLUXES CALCULATION *******************"<<endl;



	//Fit on Mass
	TH1F * DeutonsPrimaryFlux_TOF 		= (TH1F *)D_Flux     -> Flux_TOF ;
	TH1F * DeutonsPrimaryFlux_NaF 		= (TH1F *)D_Flux     -> Flux_NaF ;
	TH1F * DeutonsPrimaryFlux_Agl 		= (TH1F *)D_Flux     -> Flux_Agl ;

	TH2F * DeutonsGeomagFlux_TOF  		= (TH2F *)D_Flux_geo -> Flux_TOF ;
	TH2F * DeutonsGeomagFlux_NaF 		= (TH2F *)D_Flux_geo -> Flux_NaF ;
	TH2F * DeutonsGeomagFlux_Agl  		= (TH2F *)D_Flux_geo -> Flux_Agl ;

	  cout<<"*** Plotting ...  ****"<<endl;
	

	//Fit on Distance
	TH1F * DeutonsPrimaryFlux_Dist_TOF 	= (TH1F *)D_Flux_Dist-> Flux_TOF ;
	TH1F * DeutonsPrimaryFlux_Dist_NaF 	= (TH1F *)D_Flux_Dist-> Flux_NaF ;
	TH1F * DeutonsPrimaryFlux_Dist_Agl 	= (TH1F *)D_Flux_Dist-> Flux_Agl ;

	TH2F * DeutonsGeomagFlux_Dist_TOF  	= (TH2F *)D_Flux_geo_Dist-> Flux_TOF ;
	TH2F * DeutonsGeomagFlux_Dist_NaF 	= (TH2F *)D_Flux_geo_Dist-> Flux_NaF ;
	TH2F * DeutonsGeomagFlux_Dist_Agl  	= (TH2F *)D_Flux_geo_Dist-> Flux_Agl ;

 cout<<"*** Plotting ...  ****"<<endl;

	//for D/P ratio: P Flux
	//Fit on Mass
	TH1F * ProtonsPrimaryFlux_TOF 		= (TH1F *)P_Flux     -> Flux_TOF ;
	TH1F * ProtonsPrimaryFlux_NaF 		= (TH1F *)P_Flux     -> Flux_NaF ;
	TH1F * ProtonsPrimaryFlux_Agl 		= (TH1F *)P_Flux     -> Flux_Agl ;

	//Fit on Distance
	TH1F * ProtonsPrimaryFlux_Dist_TOF 	= (TH1F *)P_Flux_Dist-> Flux_TOF ;
	TH1F * ProtonsPrimaryFlux_Dist_NaF 	= (TH1F *)P_Flux_Dist-> Flux_NaF ;
	TH1F * ProtonsPrimaryFlux_Dist_Agl 	= (TH1F *)P_Flux_Dist-> Flux_Agl ;

 cout<<"*** Plotting ...  ****"<<endl;

		cout<< D_Flux -> Counts_R   <<endl;
	        cout<< D_Flux -> Counts_TOF <<endl;
                cout<< D_Flux -> Counts_NaF <<endl;
	

	//D/P ratio
	//Fit on Mass	
	TH1F * DP_ratioTOF = (TH1F *)D_Flux     -> Flux_TOF -> Clone();
	TH1F * DP_ratioNaF = (TH1F *)D_Flux     -> Flux_NaF -> Clone();
	TH1F * DP_ratioAgl = (TH1F *)D_Flux     -> Flux_Agl -> Clone();


	cout<< D_Flux -> Counts_R   <<endl;
        cout<< D_Flux -> Counts_TOF <<endl;
        cout<< D_Flux -> Counts_NaF <<endl;
	
	DP_ratioTOF ->  Divide((TH1F *)P_Flux     -> Flux_TOF			);
	DP_ratioNaF ->  Divide((TH1F *)P_Flux     -> Flux_NaF			);
	DP_ratioAgl ->  Divide((TH1F *)P_Flux     -> Flux_Agl			);





 cout<<"*** Plotting ...  ****"<<endl;

	//Fit on Distance
	TH1F * DP_ratioTOF_Dist = (TH1F *)D_Flux_Dist-> Flux_TOF -> Clone();
	TH1F * DP_ratioNaF_Dist = (TH1F *)D_Flux_Dist-> Flux_NaF -> Clone();
	TH1F * DP_ratioAgl_Dist = (TH1F *)D_Flux_Dist-> Flux_Agl -> Clone();

	DP_ratioTOF_Dist ->  Divide((TH1F *)P_Flux_Dist-> Flux_TOF                   );
	DP_ratioNaF_Dist ->  Divide((TH1F *)P_Flux_Dist-> Flux_NaF                   );
	DP_ratioAgl_Dist ->  Divide((TH1F *)P_Flux_Dist-> Flux_Agl                   );	
 cout<<"*************** DEUTONS FLUXES CALCULATION *******************"<<endl;

 cout<<"*** Plotting ...  ****"<<endl;

	DeutonsPrimaryFlux_TOF 	    	->SetName("DeutonsPrimaryFlux_TOF"); 
	DeutonsPrimaryFlux_NaF 	  	->SetName("DeutonsPrimaryFlux_NaF");
	DeutonsPrimaryFlux_Agl 	  	->SetName("DeutonsPrimaryFlux_Agl");
 cout<<"*** Plotting ...  ****"<<endl;
	DeutonsGeomagFlux_TOF  	  	->SetName("DeutonsGeomagFlux_TOF");
	DeutonsGeomagFlux_NaF 	  	->SetName("DeutonsGeomagFlux_NaF");
	DeutonsGeomagFlux_Agl  	  	->SetName("DeutonsGeomagFlux_Agl");
 cout<<"*** Plotting ...  ****"<<endl;
	DeutonsPrimaryFlux_Dist_TOF	->SetName("DeutonsPrimaryFlux_Dist_TOF");
	DeutonsPrimaryFlux_Dist_NaF	->SetName("DeutonsPrimaryFlux_Dist_NaF");
	DeutonsPrimaryFlux_Dist_Agl	->SetName("DeutonsPrimaryFlux_Dist_Agl");
 cout<<"*** Plotting ...  ****"<<endl;
	DeutonsGeomagFlux_Dist_TOF 	->SetName("DeutonsGeomagFlux_Dist_TOF");
	DeutonsGeomagFlux_Dist_NaF 	->SetName("DeutonsGeomagFlux_Dist_NaF");
	DeutonsGeomagFlux_Dist_Agl 	->SetName("DeutonsGeomagFlux_Dist_Agl");
 cout<<"*** Plotting ...  ****"<<endl;
	ProtonsPrimaryFlux_TOF 	  	->SetName("ProtonsPrimaryFlux_TOF");
	ProtonsPrimaryFlux_NaF 	  	->SetName("ProtonsPrimaryFlux_NaF");
	ProtonsPrimaryFlux_Agl 	  	->SetName("ProtonsPrimaryFlux_Agl");
 cout<<"*** Plotting ...  ****"<<endl;
	ProtonsPrimaryFlux_Dist_TOF	->SetName("ProtonsPrimaryFlux_Dist_TOF");
	ProtonsPrimaryFlux_Dist_NaF	->SetName("ProtonsPrimaryFlux_Dist_NaF");
	ProtonsPrimaryFlux_Dist_Agl	->SetName("ProtonsPrimaryFlux_Dist_Agl");
 cout<<"*** Plotting ...  ****"<<endl;
	DP_ratioTOF 	 -> SetName("DP_ratioTOF"		);
	DP_ratioNaF 	 -> SetName("DP_ratioNaF"		);
	DP_ratioAgl 	 -> SetName("DP_ratioAgl"		);	
 cout<<"*** Plotting ...  ****"<<endl;
	DP_ratioTOF_Dist -> SetName("DP_ratioTOF_Dist"	);
	DP_ratioNaF_Dist -> SetName("DP_ratioNaF_Dist"	);
	DP_ratioAgl_Dist -> SetName("DP_ratioAgl_Dist"	);


 cout<<"*************** DEUTONS FLUXES CALCULATION *******************"<<endl;


	finalHistos.Add(DeutonsPrimaryFlux_TOF 	   );
        finalHistos.Add(DeutonsPrimaryFlux_NaF 	    );
        finalHistos.Add(DeutonsPrimaryFlux_Agl 	    );
        finalHistos.Add(DeutonsGeomagFlux_TOF  	    );
        finalHistos.Add(DeutonsGeomagFlux_NaF 	    );
        finalHistos.Add(DeutonsGeomagFlux_Agl  	    );
        finalHistos.Add(DeutonsPrimaryFlux_Dist_TOF );
        finalHistos.Add(DeutonsPrimaryFlux_Dist_NaF );
        finalHistos.Add(DeutonsPrimaryFlux_Dist_Agl );
        finalHistos.Add(DeutonsGeomagFlux_Dist_TOF  );
        finalHistos.Add(DeutonsGeomagFlux_Dist_NaF  );
        finalHistos.Add(DeutonsGeomagFlux_Dist_Agl  );
        finalHistos.Add(ProtonsPrimaryFlux_TOF 	    );
        finalHistos.Add(ProtonsPrimaryFlux_NaF 	    );
        finalHistos.Add(ProtonsPrimaryFlux_Agl 	    );
        finalHistos.Add(ProtonsPrimaryFlux_Dist_TOF );
        finalHistos.Add(ProtonsPrimaryFlux_Dist_NaF );
        finalHistos.Add(ProtonsPrimaryFlux_Dist_Agl );
	finalHistos.Add(DP_ratioTOF 	);
        finalHistos.Add(DP_ratioNaF 	);
        finalHistos.Add(DP_ratioAgl 	);
        finalHistos.Add(DP_ratioTOF_Dist);
        finalHistos.Add(DP_ratioNaF_Dist);
        finalHistos.Add(DP_ratioAgl_Dist);

        finalHistos.writeObjsInFolder("Results/Fluxes");

        cout<<"*** Plotting ...  ****"<<endl;


	DeutonFlux_Plot(DeutonsPrimaryFlux_TOF 	   ,
                        DeutonsPrimaryFlux_NaF 	   ,
                        DeutonsPrimaryFlux_Agl 	   ,
                        DeutonsGeomagFlux_TOF  	   ,
                        DeutonsGeomagFlux_NaF 	   ,
                        DeutonsGeomagFlux_Agl  	   ,
                        DeutonsPrimaryFlux_Dist_TOF,
                        DeutonsPrimaryFlux_Dist_NaF,
                        DeutonsPrimaryFlux_Dist_Agl,
                        DeutonsGeomagFlux_Dist_TOF ,
                        DeutonsGeomagFlux_Dist_NaF ,
                        DeutonsGeomagFlux_Dist_Agl ,
                        ProtonsPrimaryFlux_TOF 	   ,
                        ProtonsPrimaryFlux_NaF 	   ,
                        ProtonsPrimaryFlux_Agl 	   ,
                        ProtonsPrimaryFlux_Dist_TOF,
                        ProtonsPrimaryFlux_Dist_NaF,
                        ProtonsPrimaryFlux_Dist_Agl,
			DP_ratioTOF 		,
			DP_ratioNaF 	,
                        DP_ratioAgl 	,
                        DP_ratioTOF_Dist,
                        DP_ratioNaF_Dist,
                        DP_ratioAgl_Dist,
			D_Flux,
			P_Flux
	
	);

	return;

}

