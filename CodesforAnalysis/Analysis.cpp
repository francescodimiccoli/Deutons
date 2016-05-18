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

	float ekmin=0.1, ekmax=1;
	ToFDB.Setbins(nbinsToF, ekmin, ekmax);
	ToFPB.Setbins(nbinsToF, ekmin, ekmax);
	string TitoliTOF[nbinsToF];

	cout<<"**************************** BETA BINS NaF***********************************"<<endl;

	ekmin=0.666, ekmax=4.025;
	NaFDB.Setbins(nbinsNaF, ekmin, ekmax);
	NaFPB.Setbins(nbinsNaF, ekmin, ekmax);

	string TitoliNaF[nbinsNaF];
	cout<<endl;
	cout<<"**************************** BETA BINS Agl***********************************"<<endl;

	ekmin=2.57, ekmax=9.01;
	AglDB.Setbins(nbinsAgl, ekmin, ekmax);	
	AglPB.Setbins(nbinsAgl, ekmin, ekmax);

	string TitoliAgl[nbinsAgl];

	cout<<endl;


	
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


