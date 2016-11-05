#include "PlottingFunctions/ProtonFlux_Plot.h"

using namespace std;

Flux * P_Flux         = new Flux("P_Flux"    	  );
Flux * P_Flux_geo     = new Flux("P_Flux_geo"     ,11);
Flux * P_Flux_geo_prim= new Flux("P_Flux_geo_prim",11);

//DvsMC check
Flux * P_Flux_pre = new Flux("P_Flux_pre" );
Flux * P_Flux_sel = new Flux("P_Flux_sel" );

void ProtonFlux_Fill(int zona) {

	if(!trgpatt.IsPhysical()) return;
	if(Tup.Beta<=0||Tup.R<=0) return;

	int Kbin=PRB.GetRBin(Tup.R);
	

	if(Distcut && Likcut) {
		P_Flux_geo-> Counts_R -> Fill(Kbin,zona);
		if(Tup.R>1.2*Tup.Rcutoff) {
			P_Flux -> Counts_R-> Fill(Kbin);
			P_Flux_geo_prim -> Counts_R -> Fill(Kbin,zona);
		}
	}
	
	if(Tup.Beta<=0||Tup.R<=0||Tup.R<1.2*Tup.Rcutoff) return;

	if(Herejcut && ProtonsMassWindow ) P_Flux_pre -> Counts_R -> Fill(Kbin);
	if(Distcut  && Likcut )  P_Flux_sel -> Counts_R -> Fill(Kbin);
		

	return;
}


void ProtonFlux_Write() {
	P_Flux      	->Write();
	P_Flux_geo  	->Write();
	P_Flux_geo_prim ->Write();
	P_Flux_pre  	->Write();
	P_Flux_sel  	->Write();

	return;
}


void ProtonFlux(string filename) {


	cout<<"*************** PROTONS FLUXES CALCULATION *******************"<<endl;

	cout<<"*** Reading  P1 file ****"<<endl;
        TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");


	TH2F * esposizionegeo_R    = (TH2F*)inputHistoFile->Get(        "esposizionegeo_R"        );
	TH2F * esposizionepgeoTOF  = (TH2F*)inputHistoFile->Get(	"esposizionepgeoTOF"	);
	TH2F * esposizionepgeoNaF  = (TH2F*)inputHistoFile->Get(	"esposizionepgeoNaF"	);
	TH2F * esposizionepgeoAgl  = (TH2F*)inputHistoFile->Get(	"esposizionepgeoAgl"	);



	Flux * P_Flux         = new Flux(inputHistoFile, "P_Flux"    	 ,"Results","Corr_AcceptanceP",1);
	Flux * P_Flux_geo     = new Flux(inputHistoFile, "P_Flux_geo"     ,"Results","Geomag_AcceptanceP",11);
	Flux * P_Flux_geo_prim= new Flux(inputHistoFile, "P_Flux_geo_prim","Results","Geomag_AcceptanceP",11);

	Flux * P_Flux_pre     = new Flux(inputHistoFile, "P_Flux_pre" 	 ,"Results","Corr_AcceptancePreP",1);
	Flux * P_Flux_sel     = new Flux(inputHistoFile, "P_Flux_sel"     ,"Results","Corr_AcceptanceP",1);

	cout<<"*************** PROTONS FLUXES CALCULATION *******************"<<endl;

	Tempi = (TH1F *)inputHistoFile->Get("Tempi");

	P_Flux         -> Set_Exposure_Time (esposizionegeo_R,esposizionepgeoTOF,esposizionepgeoNaF,esposizionepgeoAgl);
	P_Flux_geo     -> Set_Exposure_Time (Tempi);
	P_Flux_geo_prim-> Set_Exposure_Time (Tempi);

	P_Flux_pre     -> Set_Exposure_Time (esposizionegeo_R,esposizionepgeoTOF,esposizionepgeoNaF,esposizionepgeoAgl);
	P_Flux_sel     -> Set_Exposure_Time (esposizionegeo_R,esposizionepgeoTOF,esposizionepgeoNaF,esposizionepgeoAgl);

	enum {protons, deutons};

	P_Flux         -> Eval_Flux(1 , protons);
	P_Flux_geo     -> Eval_Flux(11, protons);
	P_Flux_geo_prim-> Eval_Flux(11, protons);

	P_Flux_pre     -> Eval_Flux(1 , protons);
	P_Flux_sel     -> Eval_Flux(1 , protons);

	TH1F * ProtonsPrimaryFlux = (TH1F *)P_Flux     -> Flux_R ;
	TH2F * ProtonsGeomagFlux  = (TH2F *)P_Flux_geo -> Flux_R ;
	TH1F * P_pre_PrimaryFlux  = (TH1F *)P_Flux_pre -> Flux_R ; //P flux with pre-selections alone
	TH1F * P_sel_PrimaryFlux  = (TH1F *)P_Flux_sel -> Flux_R ; //P flux with full-set selections


	ProtonsPrimaryFlux  ->SetName("ProtonsPrimaryFlux" 	);
	ProtonsGeomagFlux   ->SetName("ProtonsGeomagFlux"  	);
	P_pre_PrimaryFlux   ->SetName("P_pre_PrimaryFlux"  	);
	P_sel_PrimaryFlux   ->SetName("P_sel_PrimaryFlux"  	);

	finalHistos.Add( ProtonsPrimaryFlux);
        finalHistos.Add( ProtonsGeomagFlux );
        finalHistos.Add( P_pre_PrimaryFlux );
        finalHistos.Add( P_sel_PrimaryFlux );

	finalHistos.writeObjsInFolder("Results");

        cout<<"*** Plotting ...  ****"<<endl;

	ProtonFlux_Plot( ProtonsPrimaryFlux,
                         ProtonsGeomagFlux ,
                         P_pre_PrimaryFlux ,
                         P_sel_PrimaryFlux ,
			 Tempi		
	);


	return;

}



