#include "Functions_auto.cpp"


#include "Parte2/Definitions.cpp"
#include "Parte2/FitError.cpp"
#include "Parte2/EfficiencyClass.cpp"
#include "Parte2/LATcorrClass.cpp"
#include "Parte2/ACCEPTANCEClass.cpp"
#include "Parte2/FluxClass.cpp"
#include "Parte2/TemplateFITClass.cpp"
#include "Parte2/DatavsMCClass.h"
#include "Parte2/MCpreeff.cpp"
#include "Parte2/MCUnbiaseff.cpp"
#include "Parte2/Hecont.cpp"
#include "Parte2/SlidesforPlot.cpp"
#include "Parte2/Qualcutoptimization.cpp"
#include "Parte2/MCQualeff.cpp"
#include "Parte2/Cuts.cpp"
#include "Parte2/MCTrackeff.cpp"
#include "Parte2/Eff_Factorizationtest.cpp"
#include "Parte2/MigrationMatrix.cpp"
#include "Parte2/MCFullSeteff.cpp"
#include "Parte2/DATAUnbiaseff.cpp"
#include "Parte2/CorrelazionePreselezioni.cpp"
#include "Parte2/DATApreSeleff.cpp"
#include "Parte2/DATAQualeff.cpp"
#include "Parte2/DATARICHeff.cpp"
#include "Parte2/CorrLAT.cpp"
#include "Parte2/Acceptance.cpp"
#include "Parte2/ProtonFlux.cpp"
#include "Parte2/DeutonsCountsExtraction.cpp"
#include "Parte2/DeutonsCountsExtraction_Dist.cpp"
#include "Parte2/MCMC.cpp"
#include "Parte2/DVSMCQualeff.cpp"
#include "Parte2/DVSMCRICHeff.cpp"
#include "Parte2/DVSMCQualeff_D.cpp"
#include "Parte2/DeutonsFlux.cpp"
#include "Parte2/DVSMCPreSeleff.cpp"
#include "Parte2/DVSMCPreSeleff_D.cpp"
//#include "Parte2/DVSMCTrackeff.cpp"*/


#include "FillIstogram.cpp"


using namespace std;

int main(int argc, char * argv[])
{
	cout<<"Month _ Indx _ Frac _ output"<<endl;
	cout<<argc<<endl;
	if(argc == 1 ) {
		cout<<"No Month specified: running 2012_05"<<endl;
		mese = "2012_05";
		cout<<"No Mode specified: running Mode 2"<<endl;
		INDX = 0;
		cout<<"No fraction specified: running 35"<<endl;
		frac = "35";
		cout<<"No output path specified: writing locally"<<endl;
                outputpath = "../";
	}
	if(argc == 2 ) {
                mese=argv[1];
		cout<<"No Mode specified: running Mode 2"<<endl;
                INDX = 2;
                cout<<"No fraction specified: running 35"<<endl;
                frac = "35"; 
                cout<<"No output path specified: writing locally"<<endl;
                outputpath = "../";
        }
	if(argc == 3 ) {
                mese=argv[1];
		INDX=atoi(argv[2]);
		cout<<"No fraction specified: running 35"<<endl;
                frac = "35"; 
                cout<<"No output path specified: writing locally"<<endl;
                outputpath = "../";
        }
	if(argc == 4) {
		mese=argv[1];
                INDX=atoi(argv[2]);
		frac= argv[3];
		cout<<"No output path specified: writing locally"<<endl;
		outputpath = "../";
	}
	if(argc == 5) {
		mese=argv[1];
                INDX=atoi(argv[2]);
                frac=argv[3];
		outputpath=argv[4];
	}
	cout<<"****************************** INPUT PAR. ***********************************"<<endl;
	cout<<endl;
	cout<<"Month: "<<mese<<endl; 
	cout<<endl;
        cout<<"Mode: "<<INDX<<endl;
	cout<<endl;
        cout<<"Fraction: "<<frac<<endl;
	cout<<endl;
        cout<<"Output dir: "<<outputpath<<"Histos/"<<mese<<endl;
	cout<<endl;
	cout<<"****************************** R BINS ***************************************"<<endl;
	for(int i=0;i<nbinsr+1;i++)
	{
		float temp=i+14;
		bin[i]=0.1*pow(10,temp/(9.5*2));
		Rbins[i]=bin[i];
		if(i<nbinsr) {R_cent[i]=0.1*pow(10,(temp+0.5)/(9.5*2));
			encindeut[i]=pow(((1+pow((R_cent[i]/1.875),2))),0.5)-1;
			encinprot[i]=pow(((1+pow((R_cent[i]/0.938),2))),0.5)-1;
		}
	}
	for(int i=0;i<Rbins.size();i++) {
		deltaencinprot[i]=(pow(((1+pow((Rbins[i+1]/0.938),2))),0.5)-1)-(pow(((1+pow((Rbins[i]/0.938),2))),0.5)-1);
		deltaencindeut[i]=(pow(((1+pow((Rbins[i+1]/1.875),2))),0.5)-1)-(pow(((1+pow((Rbins[i]/1.875),2))),0.5)-1);
		cout<<R_cent[i]<<endl;
	}

	cout<<"**************************** BETA BINS TOF***********************************"<<endl;
	float B=0.4;
	float B1=0;
	float B2=0;
	float E=0.1;
	int binnum=0;
	float a=(log(1)-log(0.1))/nbinsToF;
	float E2=exp(log(0.1)+1.5*a);
	while(B1<0.85){
		E=exp(log(0.1)+binnum*a);
		E2=exp(log(0.1)+(binnum+0.5)*a);
		B1=sqrt(1-1/(pow(E+1,2)));
		B2=sqrt(1-1/(pow(E2+1,2)));
		Betabins[binnum]=B1;
		if (binnum<BetabinsR_P.size()) BetabinsR_P[binnum]=0.938*pow(pow(B1,2)/(1-pow(B1,2)),0.5);
		if (binnum<BetabinsR_D.size()) BetabinsR_D[binnum]=1.875*pow(pow(B1,2)/(1-pow(B1,2)),0.5);
		Betacent[binnum]=B2;
		Ekincent[binnum]=1/pow(1-pow(B2,2),0.5)-1;
		binnum++;
	}

	string TitoliTOF[nbinsToF];
	for(int i=0;i<nbinsToF ;i++) {
		deltaencinTOF[i]=(1/pow(1-pow(Betabins[i+1],2),0.5)-1)-(1/pow(1-pow(Betabins[i],2),0.5)-1);
		
		ostringstream ss;
		ss<<Betabins[i];
		TitoliTOF[i]= ss.str();
		cout<<BetabinsR_P[i]<<" "<<BetabinsR_D[i]<<endl;
	}
	cout<<endl;

	cout<<"**************************** BETA BINS NaF***********************************"<<endl;
	a=(log(4.025)-log(0.666))/nbinsNaF;
	E2=exp(log(0.666)+1.5*a);
	binnum=0;
	while(B1<0.98){
		E=exp(log(0.666)+binnum*a);
		E2=exp(log(0.666)+(binnum+0.5)*a);
		B1=sqrt(1-1/(pow(E+1,2)));
		B2=sqrt(1-1/(pow(E2+1,2)));
		if (binnum<BetabinsNaF   .size()) BetabinsNaF   [binnum]=B1;
		if (binnum<BetabinsNaFR_P.size()) BetabinsNaFR_P[binnum]=0.938*pow(pow(B1,2)/(1-pow(B1,2)),0.5);
		if (binnum<BetabinsNaFR_D.size()) BetabinsNaFR_D[binnum]=1.875*pow(pow(B1,2)/(1-pow(B1,2)),0.5);
		BetacentNaF[binnum]=B2;
		EkincentNaF[binnum]=1/pow(1-pow(B1,2),0.5)-1;
		binnum++;
	}

	string TitoliNaF[nbinsNaF];
	
	for(int i=0;i<nbinsNaF;i++) {	
		deltaencinNaF[i]=(1/pow(1-pow(BetabinsNaF[i+1],2),0.5)-1)-(1/pow(1-pow(BetabinsNaF[i],2),0.5)-1);
		
		ostringstream ss;
		ss<<BetabinsNaF[i];
		TitoliNaF[i]= ss.str();
	}
	cout<<endl;
	cout<<"**************************** BETA BINS Agl***********************************"<<endl;
	
	a=(log(9.01)-log(2.57))/nbinsAgl ;
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

	string TitoliAgl[nbinsAgl];
	for(int i=0;i<nbinsAgl; i++) {
		deltaencinAgl[i]=(1/pow(1-pow(BetabinsAgl[i+1],2),0.5)-1)-(1/pow(1-pow(BetabinsAgl[i],2),0.5)-1);
		
		ostringstream ss;
		ss<<BetabinsAgl[i];
		TitoliAgl[i]= ss.str();
	}
	cout<<endl;

	/*//////// BINNAGGIO IN BETA
     BetaD   =Betabins   ;
	  BetaP   =Betabins   ;
	  BetaNaFD=BetabinsNaF;
	  BetaNaFP=BetabinsNaF;
	  for(int i=0;i<19;i++){
	    BetaAglD[i]=BetabinsAgl[i];
	    BetaAglP[i]=BetabinsAgl[i];
	  }*/	
	////////////////////////////
	/////////// BINNAGGIO IN RIGIDITA'
	BetaD=BetabinsR_D;
	BetaP=BetabinsR_P;
	BetaNaFD=BetabinsNaFR_D;
	BetaNaFP=BetabinsNaFR_P;
	for(int i=0;i<20;i++){	
		BetaAglD[i]=BetabinsAglR_D[i];
		BetaAglP[i]=BetabinsAglR_P[i];
	}
	////////////////////////////

	cout<<"************************ ISTOGRAM FILLING **************************************************************"<<endl;
	FillIstogram(INDX,frac,mese);
	string	nomefile="../Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
	TFile *file1 =TFile::Open(nomefile.c_str());

	cout<<"************************* ANALYSIS **********************************************************************"<<endl;
	if(INDX==2){
		Hecut(file1);
                SlidesforPlot(file1);
		DistanceCut(file1);
		Correlazione_Preselezioni(file1);
		
		MCpreeff(file1);
                MCUnbiaseff(file1);		
		MCQualeff(file1);
		FluxFactorizationtest(file1);
		MCTrackeff(file1);
		MCFullseteff(file1);
		MigrationMatrix(file1);
		DATAUnbiaseff(file1);
		DATApreSeleff(file1);
		//DVSMCTrackeff(file1);
		DATAQualeff(file1);
		DATARICHeff(file1);
		if(frac=="tot") DeutonsTemplFits();
		if(frac=="tot") DeutonsTemplFits_Dist();
	}
	cout<<"************************* RESULTS  **************************************************************"<<endl;
	
	if(INDX==2){	
		CorrLAT();
		DVSMCPreSeleff();
		DVSMCPreSeleffD();
		DVSMCRICHeff();
		DVSMCQualeff2();
		DVSMCQualeffD();
		Acceptance();
		ProtonFlux();
		if(frac=="tot") DeutonFlux();
	}
	cout<<"************************** OUTPUT **************************************************************"<<endl;
	return 1;
}


