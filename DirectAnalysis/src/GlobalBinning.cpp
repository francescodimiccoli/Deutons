#include "Globals.h"


int Ev_Num;
int Timebeg;
float FRAC =1;
int OFFSET =0;

int nbinsr=43;
int nbinsToF=18;
int nbinsNaF=18;
int nbinsAgl=18;

float ToFsmearSigma= 79.8;
float ToFsmearShift= -1.5;

TRandom3 * Rand= new TRandom3(time(0));

Particle proton(0.9382720813, 1, 1);  // proton mass 938 MeV
Particle deuton(1.8756129   , 1, 2);  // deuterium mass 1876 MeV, Z=1, A=2


Binning ToFDB(deuton);
Binning ToFPB(proton);
Binning NaFDB(deuton);
Binning NaFPB(proton);
Binning AglDB(deuton);
Binning AglPB(proton);
Binning ForEffCorr(proton);
Binning ForEffCorr_D(deuton);
	
Binning DRB(deuton);
Binning PRB(proton);

Binning ForAcceptance(proton);

//resolution binning
Binning ToFRigB(proton);
Binning NaFRigB(proton);
Binning AglRigB(proton);

Binning PResB(proton);

TMVA::Reader *readerTOF;
TMVA::Reader *readerNaF;
TMVA::Reader *readerAgl;

TF1 * ResponseTOF = new TF1("ResponseTOF","x*(1-[0]/x^[1]) - [2]",0,1);
TF1 * ResponseNaF = new TF1("ResponseNaF","x*(1-[0]/x^[1]) - [2]",0,1);
TF1 * ResponseAgl = new TF1("ResponseAgl","x*(1-[0]/x^[1]) - [2]",0,1);

void SetBins(){	

	DRB.Reset();
	PRB.Reset();
	ToFDB.Reset();
	ToFPB.Reset();
	NaFDB.Reset();
	NaFPB.Reset();
	AglDB.Reset();
	AglPB.Reset();
	ForAcceptance.Reset();
	ToFRigB.Reset();
	NaFRigB.Reset();
	AglRigB.Reset();
	PResB  .Reset();


	cout<<"H.E. bins"<<endl;
	DRB.setBinsFromRigidity(nbinsr, 0.5, 100,ResponseTOF,0.00347548,5.8474); 
	PRB.setBinsFromRigidity(nbinsr, 0.5, 100,ResponseTOF,0.00347548,5.8474);
	ForAcceptance.setBinsFromRigidity(2*nbinsr,0.5,250,ResponseTOF,0.00347548,5.8474);

	cout<<"TOF bins"<<endl;
	float ekmin=0.1, ekmax=0.82;
	float betamin=0.55; float betamax=0.853;
	ToFDB.setBinsFromRDatacard ("/data1/home/fdimicco/Deutons/DirectAnalysis/bindatacard.data", betamin, betamax ,ResponseTOF,0.00347548,5.8474);
	ToFPB.setBinsFromRDatacard ("/data1/home/fdimicco/Deutons/DirectAnalysis/bindatacard.data", betamin, betamax ,ResponseTOF,0.00347548,5.8474);
	ToFRigB.setBinsFromBeta (nbinsToF, betamin, betamax ,ResponseTOF,0.00347548,5.8474);

	cout<<"NaF bins"<<endl;
	betamin=0.85, betamax=0.977;
	NaFDB.setBinsFromRDatacard ("/data1/home/fdimicco/Deutons/DirectAnalysis/bindatacard.data", betamin, betamax ,ResponseNaF,-0.000859132,-30.5065);
	NaFPB.setBinsFromRDatacard ("/data1/home/fdimicco/Deutons/DirectAnalysis/bindatacard.data", betamin, betamax ,ResponseNaF,-0.000859132,-30.5065);
	NaFRigB.setBinsFromEkPerMass(nbinsNaF, ekmin, ekmax,ResponseNaF,-0.000859132,-30.5065);

	cout<<"Agl bins"<<endl;
	betamin=0.97, betamax=0.995;
	AglDB.setBinsFromRDatacard ("/data1/home/fdimicco/Deutons/DirectAnalysis/bindatacard.data", betamin, betamax ,ResponseAgl,0,67.8521);
	AglPB.setBinsFromRDatacard ("/data1/home/fdimicco/Deutons/DirectAnalysis/bindatacard.data", betamin, betamax ,ResponseAgl,0,67.8521);
	AglRigB.setBinsFromEkPerMass(nbinsAgl, ekmin, ekmax,ResponseAgl,4.28781e-05,67.8521);


	//PResB.setBinsFromRigidity(60, 0.5, 100,ResponseTOF,0.00347548,5.8474);	
	//ToFResB.setBinsFromRigidity(25, 1,8,ResponseTOF,0.00347548,5.8474);
	//NaFResB.setBinsFromRigidity(25, 3,13,ResponseNaF,-0.000859132,-30.5065);
	//AglResB.setBinsFromRigidity(25, 6,25,ResponseAgl,4.28781e-05,67.8521);

	PRB.Print();

	cout<<"**TOF**"<<endl;
	ToFPB.Print();
	ToFDB.Print();

	cout<<"**NaF**"<<endl;
	NaFPB.Print();
	NaFDB.Print();

	cout<<"**Agl**"<<endl;
	AglPB.Print();
	AglDB.Print();



	return;
}

void SetUpUsualBinning(){

	SetBins();

	ToFPB.UseBetaEdges();
	NaFPB.UseBetaEdges();
	AglPB.UseBetaEdges();

	ToFDB.UseBetaEdges();
	NaFDB.UseBetaEdges();
	AglDB.UseBetaEdges();

	ToFRigB.UseREdges();
	NaFRigB.UseREdges();
	AglRigB.UseREdges();



	PRB.UseREdges();
	ForAcceptance.UseREdges();
	cout<<endl;
	return;
}

void SetUpEffCorrBinning(){

	SetBins();

	ToFPB.UseREdges();
	NaFPB.UseREdges();
	AglPB.UseREdges();

	ToFDB.UseBetaEdges();
	NaFDB.UseBetaEdges();
	AglDB.UseBetaEdges();

	ToFRigB.UseREdges();
	NaFRigB.UseREdges();
	AglRigB.UseREdges();


	PRB.UseREdges();
	ForAcceptance.UseREdges();
	cout<<endl;
	return;
}

void SetUpTOIBinning(){

	SetBins();

	ToFPB.UseBetaTOIEdges();
	NaFPB.UseBetaTOIEdges();
	AglPB.UseBetaTOIEdges();

	ToFDB.UseBetaTOIEdges();
	NaFDB.UseBetaTOIEdges();
	AglDB.UseBetaTOIEdges();

	ToFRigB.UseRTOIEdges();
	NaFRigB.UseRTOIEdges();
	AglRigB.UseRTOIEdges();


	PRB.UseRTOIEdges();
	ForAcceptance.UseRTOIEdges();
	cout<<endl;
	return;
}


void SetUpRigTOIBinning(){
	SetBins();

	ToFPB.UseRTOIEdges();
	NaFPB.UseRTOIEdges();
	AglPB.UseRTOIEdges();

	ToFDB.UseRTOIEdges();
	NaFDB.UseRTOIEdges();
	AglDB.UseRTOIEdges();

	ToFRigB.UseRTOIEdges();
	NaFRigB.UseRTOIEdges();
	AglRigB.UseRTOIEdges();

	PRB.UseRTOIEdges();
	ForAcceptance.UseRTOIEdges();
	cout<<endl;
	return;

}

void UpdateProgressBar(int currentevent, int totalentries)
{
        int newratio =(int)1000*(currentevent/    (float)totalentries);
        int oldratio =(int)1000*((currentevent-1)/(float)totalentries);
        if(newratio>oldratio)
                std::cout<<'\r' << "Progress : "<< (float)(newratio+1)/10 << " %"<< std::flush; //+1 pour finir a 100%
}

