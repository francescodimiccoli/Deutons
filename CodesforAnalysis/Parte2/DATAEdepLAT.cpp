#include "PlottingFunctions/DATAEdepLAT_Plot.h"


TH2F * EdepUTOF_lat  =new TH2F("EdepUTOF_lat" ,"EdepUTOF_lat" ,1000,0,40,11,0,11);
TH2F * EdepLTOF_lat  =new TH2F("EdepLTOF_lat" ,"EdepLTOF_lat" ,1000,0,40,11,0,11);
TH2F * EdepTrack_lat =new TH2F("EdepTrack_lat","EdepTrack_lat",1000,0,40,11,0,11);


void DATAEdepLAT_Fill(int zona){

	if(!trgpatt.IsPhysical()) return;
        if(Tup.Beta<=0||Tup.R<=0) return;	
	if(!Herejcut) return;
	if(!ProtonsMassWindow) return;
	

	EdepUTOF_lat	-> Fill(fabs(EdepTOFbeta->Eval(Tup.Beta)-Tup.EdepTOFU)/(pow(EdepTOFbeta->Eval(Tup.Beta),2)*etofu ->Eval(Tup.Beta)),zona); 
        EdepLTOF_lat 	-> Fill(fabs(EdepTOFbeta->Eval(Tup.Beta)-Tup.EdepTOFD)/(pow(EdepTOFbeta->Eval(Tup.Beta),2)*etofd ->Eval(Tup.Beta)),zona); 
        EdepTrack_lat   -> Fill(fabs(EdepTrackbeta->Eval(Tup.Beta)-Tup.EdepTrack)/(pow(EdepTrackbeta->Eval(Tup.Beta),2)*etrack->Eval(Tup.Beta)),zona); 


	return;
}


void DATAEdepLAT_Write(){

	EdepUTOF_lat  ->Write();
	EdepLTOF_lat  ->Write();
	EdepTrack_lat ->Write();

}



void DATAEdepLAT(string filename){

	cout<<"****************************** DATA LATITUDE E. dep.  **************************************"<<endl;

        cout<<"*** Reading  P1 file ****"<<endl;
        TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");
	

	TH2F * EdepUTOF_lat  =(TH2F *)inputHistoFile->Get("EdepUTOF_lat"); 
        TH2F * EdepLTOF_lat  =(TH2F *)inputHistoFile->Get("EdepLTOF_lat"); 
        TH2F * EdepTrack_lat =(TH2F *)inputHistoFile->Get("EdepTrack_lat"); 
	

	cout<<"****************************** DATA LATITUDE E. dep.  **************************************"<<endl;

	TH1F *EdepUTOF_coll[11];
	TH1F *EdepLTOF_coll[11];
	TH1F *EdepTrack_coll[11];


	for(int i=1;i<11;i++){
		EdepUTOF_coll[i]=  ProjectionXtoTH1F(EdepUTOF_lat,  ("Geo. Zone" + to_string(i)).c_str(), i, i+1);	
		EdepLTOF_coll[i]=  ProjectionXtoTH1F(EdepLTOF_lat,  ("Geo. Zone" + to_string(i)).c_str(), i, i+1);	
		EdepTrack_coll[i]= ProjectionXtoTH1F(EdepTrack_lat, ("Geo. Zone" + to_string(i)).c_str(), i, i+1);	
		
	}


	//normalization
	
	for(int i=1;i<11;i++){
	
		EdepUTOF_coll[i] -> Scale(1/EdepUTOF_coll[i]->GetEntries());
                EdepLTOF_coll[i] -> Scale(1/EdepLTOF_coll[i]->GetEntries());
                EdepTrack_coll[i]-> Scale(1/EdepTrack_coll[i]->GetEntries());
	}


	
	cout<<"*** Plotting ...  ****"<<endl;
	
	DATAEdepLAT_Plot(
	
		EdepUTOF_coll,
                EdepLTOF_coll,
                EdepTrack_coll


	);	


}

