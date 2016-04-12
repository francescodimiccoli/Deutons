#include "TH2.h"
#include "TH3.h"
#include <TVector3.h>
#include <fstream>
#include <sstream>
#include <vector>
#include "TCanvas.h"
#include "TLegend.h"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <cstdlib>
#include <stdio.h>
#include <stdarg.h>
#include <TSpline.h>
#include "TFractionFitter.h"
#include "THStack.h"
#include "TNtuple.h"
#include "TObject.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include <cstring>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include <math.h>


extern const int nbinsr=43;


#include "Functions_auto.h"
#include "Parte2/Definitions.h"
#include "Parte2/FitError.h"
#include "Parte2/EfficiencyClass.h"
#include "Parte2/LATcorrClass.h"
#include "Parte2/ACCEPTANCEClass.h"
#include "Parte2/FluxClass.h"
#include "Parte2/MCpreeff.h"
#include "Parte2/MCUnbiaseff.h"
#include "Parte2/Hecut.h"
#include "Parte2/SlidesforPlot.h"
#include "Parte2/MCQualeff.h"
#include "Parte2/Cuts.h"
#include "Parte2/MCTrackeff.h"
#include "Parte2/MC_do_preSeleff.h"
#include "Parte2/MigrationMatrix.h"
#include "Parte2/MCFullSeteff.h"
#include "Parte2/DATAUnbiaseff.h"
#include "Parte2/CorrelazionePreselezioni.h"
#include "Parte2/DATApreSeleff.h"
#include "Parte2/DATAQualeff.h"
#include "Parte2/CorrLAT.h"
#include "Parte2/Acceptance.h"
#include "Parte2/ProtonFlux.h"
/*#include "Parte2/DVSMCpreSeleff.h"
#include "Parte2/DVSMCQualeff2.h"
#include "Parte2/DVSMCTrackeff.h"
#include "Parte2/Deutons.h"
#include "Parte2/DeutonsDist.h"
#include "Parte2/MCMC.h"
#include "Parte2/DeutonsFlux.h"
#include "Parte2/DeutonsFlux_Dist.h"*/
#include "FillIstogram.h"
using namespace std;

int main(int argc, char * argv[])
{
	cout<<"Month _ Indx _ Frac "<<endl;
	frac=argv[3];
	FRAC=atoi(argv[3]);
	INDX=atoi(argv[2]);
	mese=argv[1];
	cout<<"****************************** R BINS ***************************************"<<endl;
	for(int i=0;i<nbinsr+1;i++)
	{
		float temp=i+14;
		bin[i]=0.1*pow(10,temp/(9.5*2));
		if(i<nbinsr) {R_cent[i]=0.1*pow(10,(temp+0.5)/(9.5*2));
			encindeut[i]=pow(((1+pow((R_cent[i]/1.875),2))),0.5)-1;
			encinprot[i]=pow(((1+pow((R_cent[i]/0.938),2))),0.5)-1;
		}
		cout<<bin[i]<<endl;
	}
	for(int i=0;i<nbinsr;i++) {
		deltaencinprot[i]=(pow(((1+pow((bin[i+1]/0.938),2))),0.5)-1)-(pow(((1+pow((bin[i]/0.938),2))),0.5)-1);
		deltaencindeut[i]=(pow(((1+pow((bin[i+1]/1.875),2))),0.5)-1)-(pow(((1+pow((bin[i]/1.875),2))),0.5)-1);
	}

	cout<<"**************************** BETA BINS TOF***********************************"<<endl;
	float B=0.4;
	float B1=0;
	float B2=0;
	float E=0.1;
	int binnum=0;
	float a=(log(1)-log(0.1))/18;
	float E2=exp(log(0.1)+1.5*a);
	while(B1<0.85){
		E=exp(log(0.1)+binnum*a);
		E2=exp(log(0.1)+(binnum+0.5)*a);
		B1=sqrt(1-1/(pow(E+1,2)));
		B2=sqrt(1-1/(pow(E2+1,2)));
		Betabins[binnum]=B1;
		BetabinsR_P[binnum]=0.938*pow(pow(B1,2)/(1-pow(B1,2)),0.5);
		BetabinsR_D[binnum]=1.875*pow(pow(B1,2)/(1-pow(B1,2)),0.5);
		Betacent[binnum]=B2;
		Ekincent[binnum]=1/pow(1-pow(B2,2),0.5)-1;
		binnum++;
	}
	for(int i=0;i<18;i++) deltaencinTOF[i]=(1/pow(1-pow(Betabins[i+1],2),0.5)-1)-(1/pow(1-pow(Betabins[i],2),0.5)-1);
	for(int i=0;i<18;i++) cout<<Betabins[i]<<" "<<BetabinsR_D[i]<<" "<<BetabinsR_P[i]<<endl;
	string TitoliTOF[18];
	for(int i=0;i<18;i++){
		ostringstream ss;
		ss<<Betabins[i];
		TitoliTOF[i]= ss.str();
		cout<<TitoliTOF[i]<<", ";
	}
	cout<<endl;

	cout<<"**************************** BETA BINS NaF***********************************"<<endl;
	a=(log(4.025)-log(0.666))/18;
	E2=exp(log(0.666)+1.5*a);
	binnum=0;
	while(B1<0.98){
		E=exp(log(0.666)+binnum*a);
		E2=exp(log(0.666)+(binnum+0.5)*a);
		B1=sqrt(1-1/(pow(E+1,2)));
		B2=sqrt(1-1/(pow(E2+1,2)));
		BetabinsNaF[binnum]=B1;
		BetabinsNaFR_P[binnum]=0.938*pow(pow(B1,2)/(1-pow(B1,2)),0.5);
		BetabinsNaFR_D[binnum]=1.875*pow(pow(B1,2)/(1-pow(B1,2)),0.5);
		BetacentNaF[binnum]=B2;
		EkincentNaF[binnum]=1/pow(1-pow(B1,2),0.5)-1;
		binnum++;
	}
	for(int i=0;i<18;i++) deltaencinNaF[i]=(1/pow(1-pow(BetabinsNaF[i+1],2),0.5)-1)-(1/pow(1-pow(BetabinsNaF[i],2),0.5)-1);
	for(int i=0;i<18;i++) cout<<BetabinsNaF[i]<<" "<<BetabinsNaFR_D[i]<<" "<<EkincentNaF[i]<<" "<<deltaencinNaF[i]<<endl;
	string TitoliNaF[18];
	for(int i=0;i<18;i++){
		ostringstream ss;
		ss<<BetabinsNaF[i];
		TitoliNaF[i]= ss.str();
		cout<<TitoliNaF[i]<<", ";
	}
	cout<<endl;
	cout<<"**************************** BETA BINS Agl***********************************"<<endl;
	a=(log(9.01)-log(2.57))/18;
	E2=exp(log(2.57)+1.5*a);
	binnum=0;
	while(B1<0.995){
		E=exp(log(2.57)+binnum*a);
		E2=exp(log(2.57)+(binnum+0.5)*a);
		B1=sqrt(1-1/(pow(E+1,2)));
		B2=sqrt(1-1/(pow(E2+1,2)));
		BetabinsAgl[binnum]=B1;
		BetabinsAglR_P[binnum]=0.938*pow(pow(B1,2)/(1-pow(B1,2)),0.5);
		BetabinsAglR_D[binnum]=1.875*pow(pow(B1,2)/(1-pow(B1,2)),0.5);
		BetacentAgl[binnum]=B2;
		EkincentAgl[binnum]=1/pow(1-pow(B2,2),0.5)-1;
		binnum++;
	}
	for(int i=0;i<18;i++) deltaencinAgl[i]=(1/pow(1-pow(BetabinsAgl[i+1],2),0.5)-1)-(1/pow(1-pow(BetabinsAgl[i],2),0.5)-1);
	for(int i=0;i<18;i++) cout<<BetabinsAgl[i]<<" "<<BetabinsAglR_D[i]<<" "<<EkincentAgl[i]<<" "<<deltaencinAgl[i]<<endl;
	string TitoliAgl[18];
	for(int i=0;i<18;i++){
		ostringstream ss;
		ss<<BetabinsAgl[i];
		TitoliAgl[i]= ss.str();
		cout<<TitoliAgl[i]<<", ";
	}
	cout<<endl;

	////////// BINNAGGIO IN BETA
	/*for(int i=0;i<19;i++){
	  BetaD[i]=Betabins[i];
	  BetaP[i]=Betabins[i];
	  BetaNaFD[i]=BetabinsNaF[i];
	  BetaNaFP[i]=BetabinsNaF[i];
	  BetaAglD[i]=BetabinsAgl[i];
	  BetaAglP[i]=BetabinsAgl[i];
	  }*/	
	////////////////////////////
	/////////// BINNAGGIO IN RIGIDITA'
	for(int i=0;i<20;i++){
		BetaD[i]=BetabinsR_D[i];
		BetaP[i]=BetabinsR_P[i];
		BetaNaFD[i]=BetabinsNaFR_D[i];
		BetaNaFP[i]=BetabinsNaFR_P[i];
		BetaAglD[i]=BetabinsAglR_D[i];
		BetaAglP[i]=BetabinsAglR_P[i];
		cout<<"TOF: "<<BetaD[i]<<" "<<BetaP[i]<<" NaF "<<BetaNaFD[i]<<" "<<BetaNaFP[i]<<" Agl "<<BetaAglD[i]<<" "<<BetaAglP[i]<<endl;
	}
	////////////////////////////

	cout<<"************************ ISTOGRAM FILLING **************************************************************"<<endl;
	FillIstogram(INDX,frac,mese);
	string	nomefile="../Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
	TFile *file1 =TFile::Open(nomefile.c_str());

	cout<<"************************* ANALYSIS **********************************************************************"<<endl;
	if(INDX==2){
		MCpreeff(file1);
		MCUnbiaseff(file1);
		Hecut(file1);
		SlidesforPlot(file1);
		MCQualeff(file1);
		MCTrackeff(file1);
		MCFullseteff(file1);
		MC_do_preSeleff(file1);
		Correlazione_Preselezioni(file1);
		MigrationMatrix(file1);
		DATAUnbiaseff(file1);
		DATApreSeleff(file1);
		//DVSMCTrackeff(file1);
		DATAQualeff(file1);
		/*DVSMCpreSeleff(file1);
		  DVSMCQualeff2(file1);
		  DeutonsTemplFits(file1);
		  DeutonsTemplFits_Dist(file1);
		  DeutonFlux(file1);
		  DeutonFlux_Dist(file1);*/
	}
	cout<<"************************* RESULTS  **************************************************************"<<endl;
	
	if(INDX==2){	
		CorrLAT();
		Acceptance();
		ProtonFlux();
	}
	cout<<"************************** OUTPUT **************************************************************"<<endl;
	return 1;
}
