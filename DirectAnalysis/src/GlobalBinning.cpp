#include "Globals.h"
#include "GlobalPaths.h"

std::string workdir ="/afs/cern.ch/user/f/fdimicco/Work/Deutons/DirectAnalysis";
std::string outdir ="/eos/ams/user/f/fdimicco/AnalysisFiles/";

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


Binning ForEffCorr(proton);
Binning ForEffCorr_D(deuton);
Binning HefragmToF(deuton);
Binning HefragmNaF(deuton);
Binning HefragmAgl(deuton);

	
Binning DRB(deuton);
Binning PRB(proton);

Binning ForAcceptance(proton);


Binning PResB(proton);

TMVA::Reader *readerTOF;
TMVA::Reader *readerNaF;
TMVA::Reader *readerAgl;

TF1 * ResponseTOF = new TF1("ResponseTOF","x*(1-[0]/x^[1]) - [2]",0,1);
TF1 * ResponseNaF = new TF1("ResponseNaF","x*(1-[0]/x^[1]) - [2]",0,1);
TF1 * ResponseAgl = new TF1("ResponseAgl","x*(1-[0]/x^[1]) - [2]",0,1);


RangeMerger Global;
RangeMerger GlobalRig;


void SetBins(){	

	Global.Reset();
	DRB.Reset();
	PRB.Reset();
	ForAcceptance.Reset();
	GlobalRig.Reset();
	PResB  .Reset();
	HefragmToF.Reset();
	HefragmNaF.Reset();
	HefragmAgl.Reset();

	cout<<"H.E. bins"<<endl;
	DRB.setBinsFromRigidity(nbinsr, 0.5, 100,ResponseTOF,0.00347548,5.8474); 
	PRB.setBinsFromRigidity(nbinsr, 0.5, 100,ResponseTOF,0.00347548,5.8474);
	ForAcceptance.setBinsFromRigidity(3*nbinsr,0.5,1250,ResponseTOF,0.00347548,5.8474);

	cout<<"Global Bins"<<endl;
	Global.setBinsFromRDatacard((workdir+"/bindatacard_mod.data").c_str(),ResponseTOF,ResponseNaF,ResponseAgl);
	GlobalRig.setBinsFromRDatacard((workdir+"/bindatacard_mod.data").c_str(),ResponseTOF,ResponseNaF,ResponseAgl);
	
	cout<<"TOF bins"<<endl;
	float ekmin=0.1, ekmax=0.82;
	float betamin=0.55; float betamax=0.853;
	HefragmToF.setBinsFromEkPerMass (4,0.15,0.504,ResponseTOF,0.00347548,5.8474);

	cout<<"NaF bins"<<endl;
	betamin=0.85, betamax=0.977;
	HefragmNaF.setBinsFromEkPerMass (1,1.5,3,ResponseNaF,0.00347548,5.8474);


	cout<<"Agl bins"<<endl;
	betamin=0.97, betamax=0.995;
	HefragmAgl.setBinsFromEkPerMass (2,2.6,11.9,ResponseAgl,0.00347548,5.8474);

	//PResB.setBinsFromRigidity(60, 0.5, 100,ResponseTOF,0.00347548,5.8474);	
	//ToFResB.setBinsFromRigidity(25, 1,8,ResponseTOF,0.00347548,5.8474);
	//NaFResB.setBinsFromRigidity(25, 3,13,ResponseNaF,-0.000859132,-30.5065);
	//AglResB.setBinsFromRigidity(25, 6,25,ResponseAgl,4.28781e-05,67.8521);

	cout<<"**TOF**"<<endl;
	Global.GetToFPBins().Print();
	Global.GetToFDBins().Print();


	cout<<"**NaF**"<<endl;
	Global.GetNaFPBins().Print();
	Global.GetNaFDBins().Print();


	cout<<"**Agl**"<<endl;
	Global.GetAglPBins().Print();
	Global.GetAglDBins().Print();


	cout<<"**Global**"<<endl;
	Global.GetGlobalPBins().Print();
	Global.GetGlobalDBins().Print();


	return;
}

void SetUpUsualBinning(){

	SetBins();

	Global.UseBetaEdges();
	GlobalRig.UseREdges();

	HefragmToF.UseBetaEdges();
	HefragmNaF.UseBetaEdges();
	HefragmAgl.UseBetaEdges();

	PRB.UseREdges();
	ForAcceptance.UseREdges();
	cout<<endl;
	return;
}

void SetUpTOIBinning(){

	SetBins();

	Global.UseBetaTOIEdges();
	GlobalRig.UseRTOIEdges();


	PRB.UseRTOIEdges();
	ForAcceptance.UseRTOIEdges();
	cout<<endl;
	return;
}


void SetUpRigTOIBinning(){
	SetBins();

	Global.UseRTOIEdges();
	GlobalRig.UseRTOIEdges();

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

