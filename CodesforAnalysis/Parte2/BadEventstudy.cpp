#include "PlottingFunctions/BadEventStudy_Plot.h"


using namespace std;

TH2F * EdepUtofvsLtof_TOF = new TH2F("EdepUtofvsLtof_TOF","EdepUtofvsLtof_TOF",1000,-50,100,1000,-50,100);
TH2F * EdepUtofvsLtof_NaF = new TH2F("EdepUtofvsLtof_NaF","EdepUtofvsLtof_NaF",1000,-50,100,1000,-50,100);
TH2F * EdepUtofvsLtof_Agl = new TH2F("EdepUtofvsLtof_Agl","EdepUtofvsLtof_Agl",1000,-50,100,1000,-50,100);

TH2F * EdepUtofvsTrack_TOF = new TH2F("EdepUtofvsTrack_TOF","EdepUtofvsTrack_TOF",1000,-50,100,1000,-50,100);
TH2F * EdepUtofvsTrack_NaF = new TH2F("EdepUtofvsTrack_NaF","EdepUtofvsTrack_NaF",1000,-50,100,1000,-50,100);
TH2F * EdepUtofvsTrack_Agl = new TH2F("EdepUtofvsTrack_Agl","EdepUtofvsTrack_Agl",1000,-50,100,1000,-50,100);

TH2F * EdepLtofvsTrack_TOF = new TH2F("EdepLtofvsTrack_TOF","EdepLtofvsTrack_TOF",1000,-50,100,1000,-50,100);
TH2F * EdepLtofvsTrack_NaF = new TH2F("EdepLtofvsTrack_NaF","EdepLtofvsTrack_NaF",1000,-50,100,1000,-50,100);
TH2F * EdepLtofvsTrack_Agl = new TH2F("EdepLtofvsTrack_Agl","EdepLtofvsTrack_Agl",1000,-50,100,1000,-50,100);


TH3F * EdepUtofvsLtofD_TOF = new TH3F("EdepUtofvsLtofD_TOF","EdepUtofvsLtofD_TOF",1000,-50,100,1000,-50,100,6,0,6);
TH3F * EdepUtofvsLtofD_NaF = new TH3F("EdepUtofvsLtofD_NaF","EdepUtofvsLtofD_NaF",1000,-50,100,1000,-50,100,6,0,6);
TH3F * EdepUtofvsLtofD_Agl = new TH3F("EdepUtofvsLtofD_Agl","EdepUtofvsLtofD_Agl",1000,-50,100,1000,-50,100,6,0,6);

TH3F * EdepUtofvsTrackD_TOF = new TH3F("EdepUtofvsTrackD_TOF","EdepUtofvsTrackD_TOF",1000,-50,100,1000,-50,100,6,0,6);
TH3F * EdepUtofvsTrackD_NaF = new TH3F("EdepUtofvsTrackD_NaF","EdepUtofvsTrackD_NaF",1000,-50,100,1000,-50,100,6,0,6);
TH3F * EdepUtofvsTrackD_Agl = new TH3F("EdepUtofvsTrackD_Agl","EdepUtofvsTrackD_Agl",1000,-50,100,1000,-50,100,6,0,6);

TH3F * EdepLtofvsTrackD_TOF = new TH3F("EdepLtofvsTrackD_TOF","EdepLtofvsTrackD_TOF",1000,-50,100,1000,-50,100,6,0,6);
TH3F * EdepLtofvsTrackD_NaF = new TH3F("EdepLtofvsTrackD_NaF","EdepLtofvsTrackD_NaF",1000,-50,100,1000,-50,100,6,0,6);
TH3F * EdepLtofvsTrackD_Agl = new TH3F("EdepLtofvsTrackD_Agl","EdepLtofvsTrackD_Agl",1000,-50,100,1000,-50,100,6,0,6);

TH2F * EdepUtofvsLtofHe_TOF = new TH2F("EdepUtofvsLtofHe_TOF","EdepUtofvsLtofHe_TOF",1000,-50,100,1000,-50,100);
TH2F * EdepUtofvsLtofHe_NaF = new TH2F("EdepUtofvsLtofHe_NaF","EdepUtofvsLtofHe_NaF",1000,-50,100,1000,-50,100);
TH2F * EdepUtofvsLtofHe_Agl = new TH2F("EdepUtofvsLtofHe_Agl","EdepUtofvsLtofHe_Agl",1000,-50,100,1000,-50,100);

TH2F * EdepUtofvsTrackHe_TOF = new TH2F("EdepUtofvsTrackHe_TOF","EdepUtofvsTrackHe_TOF",1000,-50,100,1000,-50,100);
TH2F * EdepUtofvsTrackHe_NaF = new TH2F("EdepUtofvsTrackHe_NaF","EdepUtofvsTrackHe_NaF",1000,-50,100,1000,-50,100);
TH2F * EdepUtofvsTrackHe_Agl = new TH2F("EdepUtofvsTrackHe_Agl","EdepUtofvsTrackHe_Agl",1000,-50,100,1000,-50,100);

TH2F * EdepLtofvsTrackHe_TOF = new TH2F("EdepLtofvsTrackHe_TOF","EdepLtofvsTrackHe_TOF",1000,-50,100,1000,-50,100);
TH2F * EdepLtofvsTrackHe_NaF = new TH2F("EdepLtofvsTrackHe_NaF","EdepLtofvsTrackHe_NaF",1000,-50,100,1000,-50,100);
TH2F * EdepLtofvsTrackHe_Agl = new TH2F("EdepLtofvsTrackHe_Agl","EdepLtofvsTrackHe_Agl",1000,-50,100,1000,-50,100);






void BadEventStudy_Fill(){

	float mass = 0;
	//cuts
	if(Tup.Beta<=0||Tup.R<=0) return;
	if(!(Likcut)) return;
	if(!(Betastrongcut)) return;

	float eutof=-(EdepTOFbeta->Eval(Tup.Beta)-Tup.EdepTOFU)/(pow(EdepTOFbeta->Eval(Tup.Beta),2)*etofu->Eval(Tup.Beta));
	float eltof=-(EdepTOFbeta->Eval(Tup.Beta)-Tup.EdepTOFD)/(pow(EdepTOFbeta->Eval(Tup.Beta),2)*etofd->Eval(Tup.Beta));	
	float etrac=-(EdepTrackbeta->Eval(Tup.Beta)-Tup.EdepTrack)/(pow(EdepTrackbeta->Eval(Tup.Beta),2)*etrack->Eval(Tup.Beta));	

	if(cmask.isOnlyFromToF()){
		mass = ((Tup.R/Tup.Beta)*pow((1-pow(Tup.Beta,2)),0.5));
		if(mass>1.875&&mass<2.4){
			if(Massa_gen<1){	
				EdepUtofvsLtof_TOF -> Fill(eutof,eltof);
				EdepUtofvsTrack_TOF-> Fill(eutof,etrac);
				EdepLtofvsTrack_TOF-> Fill(eltof,etrac);
			}
			if(Massa_gen>1&&Massa_gen<2){	
				EdepUtofvsLtofD_TOF -> Fill(eutof,eltof,ReturnMCGenType());
				EdepUtofvsTrackD_TOF-> Fill(eutof,etrac,ReturnMCGenType());
				EdepLtofvsTrackD_TOF-> Fill(eltof,etrac,ReturnMCGenType());
			}
		}
		if(Massa_gen>3){	
				EdepUtofvsLtofHe_TOF -> Fill(eutof,eltof);
				EdepUtofvsTrackHe_TOF-> Fill(eutof,etrac);
				EdepLtofvsTrackHe_TOF-> Fill(eltof,etrac);
			}
	}

	if(cmask.isFromNaF()){
		mass = ((Tup.R/Tup.BetaRICH)*pow((1-pow(Tup.BetaRICH,2)),0.5));
		if(mass>1.875&&mass<2.4){
			if(Massa_gen<1){	
				EdepUtofvsLtof_NaF -> Fill(eutof,eltof);
				EdepUtofvsTrack_NaF-> Fill(eutof,etrac);
				EdepLtofvsTrack_NaF-> Fill(eltof,etrac);
			}
			if(Massa_gen>1&&Massa_gen<2){	
				EdepUtofvsLtofD_NaF -> Fill(eutof,eltof,ReturnMCGenType());
				EdepUtofvsTrackD_NaF-> Fill(eutof,etrac,ReturnMCGenType());
				EdepLtofvsTrackD_NaF-> Fill(eltof,etrac,ReturnMCGenType());
			}
		}
		if(Massa_gen>3){	
				EdepUtofvsLtofHe_NaF -> Fill(eutof,eltof);
				EdepUtofvsTrackHe_NaF-> Fill(eutof,etrac);
				EdepLtofvsTrackHe_NaF-> Fill(eltof,etrac);
			}
	}

	if(cmask.isFromAgl()){
		mass = ((Tup.R/Tup.BetaRICH)*pow((1-pow(Tup.BetaRICH,2)),0.5));
		if(mass>1.875&&mass<2.4){
			if(Massa_gen<1){	
				EdepUtofvsLtof_Agl -> Fill(eutof,eltof);
				EdepUtofvsTrack_Agl-> Fill(eutof,etrac);
				EdepLtofvsTrack_Agl-> Fill(eltof,etrac);
			}
			if(Massa_gen>1&&Massa_gen<2){	
				EdepUtofvsLtofD_Agl -> Fill(eutof,eltof,ReturnMCGenType());
				EdepUtofvsTrackD_Agl-> Fill(eutof,etrac,ReturnMCGenType());
				EdepLtofvsTrackD_Agl-> Fill(eltof,etrac,ReturnMCGenType());
			}
		}
			if(Massa_gen>3){	
				EdepUtofvsLtofHe_Agl -> Fill(eutof,eltof);
				EdepUtofvsTrackHe_Agl-> Fill(eutof,etrac);
				EdepLtofvsTrackHe_Agl-> Fill(eltof,etrac);
			}

		
	}



}



void BadEventStudy_Write(){

	EdepUtofvsLtof_TOF  ->Write();
        EdepUtofvsLtof_NaF  ->Write();
        EdepUtofvsLtof_Agl  ->Write();
                           
        EdepUtofvsTrack_TOF ->Write();
        EdepUtofvsTrack_NaF ->Write();
        EdepUtofvsTrack_Agl ->Write();
                           
        EdepLtofvsTrack_TOF ->Write();
        EdepLtofvsTrack_NaF ->Write();
        EdepLtofvsTrack_Agl ->Write();

	EdepUtofvsLtofD_TOF  ->Write();
        EdepUtofvsLtofD_NaF  ->Write();
        EdepUtofvsLtofD_Agl  ->Write();
                           
        EdepUtofvsTrackD_TOF ->Write();
        EdepUtofvsTrackD_NaF ->Write();
        EdepUtofvsTrackD_Agl ->Write();
                           
        EdepLtofvsTrackD_TOF ->Write();
        EdepLtofvsTrackD_NaF ->Write();
        EdepLtofvsTrackD_Agl ->Write();

	EdepUtofvsLtofHe_TOF  ->Write();
        EdepUtofvsLtofHe_NaF  ->Write();
        EdepUtofvsLtofHe_Agl  ->Write();
                           
        EdepUtofvsTrackHe_TOF ->Write();
        EdepUtofvsTrackHe_NaF ->Write();
        EdepUtofvsTrackHe_Agl ->Write();
                           
        EdepLtofvsTrackHe_TOF ->Write();
        EdepLtofvsTrackHe_NaF ->Write();
        EdepLtofvsTrackHe_Agl ->Write();

}


void BadEventStudy(std::string filename){

	cout<<"*** Bad Event study ***"<<endl;

	cout<<"*** Reading  P1 file ****"<<endl;
	TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");	

		TH2F * EdepUtofvsLtof_TOF   = 	(TH2F *) inputHistoFile->Get("EdepUtofvsLtof_TOF");
		TH2F * EdepUtofvsLtof_NaF   =   (TH2F *) inputHistoFile->Get("EdepUtofvsLtof_NaF");
		TH2F * EdepUtofvsLtof_Agl   =   (TH2F *) inputHistoFile->Get("EdepUtofvsLtof_Agl");
		TH2F * EdepUtofvsTrack_TOF  =   (TH2F *) inputHistoFile->Get("EdepUtofvsTrack_TOF");
		TH2F * EdepUtofvsTrack_NaF  =   (TH2F *) inputHistoFile->Get("EdepUtofvsTrack_NaF");
		TH2F * EdepUtofvsTrack_Agl  =   (TH2F *) inputHistoFile->Get("EdepUtofvsTrack_Agl");
		TH2F * EdepLtofvsTrack_TOF  =   (TH2F *) inputHistoFile->Get("EdepLtofvsTrack_TOF");
		TH2F * EdepLtofvsTrack_NaF  =   (TH2F *) inputHistoFile->Get("EdepLtofvsTrack_NaF");
		TH2F * EdepLtofvsTrack_Agl  =   (TH2F *) inputHistoFile->Get("EdepLtofvsTrack_Agl");
		TH3F * EdepUtofvsLtofD_TOF  =   (TH3F *) inputHistoFile->Get("EdepUtofvsLtofD_TOF");
		TH3F * EdepUtofvsLtofD_NaF  =   (TH3F *) inputHistoFile->Get("EdepUtofvsLtofD_NaF");
		TH3F * EdepUtofvsLtofD_Agl  =   (TH3F *) inputHistoFile->Get("EdepUtofvsLtofD_Agl");
		TH3F * EdepUtofvsTrackD_TOF =   (TH3F *) inputHistoFile->Get("EdepUtofvsTrackD_TOF");
		TH3F * EdepUtofvsTrackD_NaF =   (TH3F *) inputHistoFile->Get("EdepUtofvsTrackD_NaF");
		TH3F * EdepUtofvsTrackD_Agl =   (TH3F *) inputHistoFile->Get("EdepUtofvsTrackD_Agl");
		TH3F * EdepLtofvsTrackD_TOF =   (TH3F *) inputHistoFile->Get("EdepLtofvsTrackD_TOF");
		TH3F * EdepLtofvsTrackD_NaF =   (TH3F *) inputHistoFile->Get("EdepLtofvsTrackD_NaF");
		TH3F * EdepLtofvsTrackD_Agl =   (TH3F *) inputHistoFile->Get("EdepLtofvsTrackD_Agl");
		TH2F * EdepUtofvsLtofHe_TOF =   (TH2F *) inputHistoFile->Get("EdepUtofvsLtofHe_TOF");
		TH2F * EdepUtofvsLtofHe_NaF =   (TH2F *) inputHistoFile->Get("EdepUtofvsLtofHe_NaF");
		TH2F * EdepUtofvsLtofHe_Agl =   (TH2F *) inputHistoFile->Get("EdepUtofvsLtofHe_Agl");
		TH2F * EdepUtofvsTrackHe_TOF=   (TH2F *) inputHistoFile->Get("EdepUtofvsTrackHe_TOF");
		TH2F * EdepUtofvsTrackHe_NaF=   (TH2F *) inputHistoFile->Get("EdepUtofvsTrackHe_NaF");
		TH2F * EdepUtofvsTrackHe_Agl=   (TH2F *) inputHistoFile->Get("EdepUtofvsTrackHe_Agl");
		TH2F * EdepLtofvsTrackHe_TOF=   (TH2F *) inputHistoFile->Get("EdepLtofvsTrackHe_TOF");
		TH2F * EdepLtofvsTrackHe_NaF=   (TH2F *) inputHistoFile->Get("EdepLtofvsTrackHe_NaF");
		TH2F * EdepLtofvsTrackHe_Agl=   (TH2F *) inputHistoFile->Get("EdepLtofvsTrackHe_Agl");

	cout<<"*** Bad Event study ***"<<endl;
	


	cout<<"*** Plotting ****"<<endl;
	
	BadEventStudy_plot(
		EdepUtofvsLtof_TOF  , 	
                EdepUtofvsLtof_NaF   ,
                EdepUtofvsLtof_Agl   ,
                EdepUtofvsTrack_TOF  ,
                EdepUtofvsTrack_NaF  ,
                EdepUtofvsTrack_Agl  ,
                EdepLtofvsTrack_TOF  ,
                EdepLtofvsTrack_NaF  ,
		EdepLtofvsTrack_Agl  ,
                EdepUtofvsLtofD_TOF  ,
                EdepUtofvsLtofD_NaF  ,
                EdepUtofvsLtofD_Agl  ,
                EdepUtofvsTrackD_TOF ,
                EdepUtofvsTrackD_NaF ,
                EdepUtofvsTrackD_Agl ,
                EdepLtofvsTrackD_TOF ,
                EdepLtofvsTrackD_NaF ,
                EdepLtofvsTrackD_Agl ,
                EdepUtofvsLtofHe_TOF ,
                EdepUtofvsLtofHe_NaF ,
                EdepUtofvsLtofHe_Agl ,
                EdepUtofvsTrackHe_TOF,
                EdepUtofvsTrackHe_NaF,
                EdepUtofvsTrackHe_Agl,
                EdepLtofvsTrackHe_TOF,
                EdepLtofvsTrackHe_NaF,
                EdepLtofvsTrackHe_Agl
	);	




}



