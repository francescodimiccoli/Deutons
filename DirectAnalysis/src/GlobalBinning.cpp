#include "Globals.h"
#include "GlobalPaths.h"

std::string workdir ="/afs/cern.ch/user/f/fdimicco/Work/Deutons/DirectAnalysis";
std::string outdir ="/eos/ams/user/f/fdimicco/AnalysisFiles/";

int Ev_Num;
int Timebeg;
float FRAC =1;
int OFFSET =0;

int nbinsr=44;
int nbinsToF=18;
int nbinsNaF=18;
int nbinsAgl=18;

float ToFsmearSigma= 79.8;
float ToFsmearShift= -1.5;

float RcutoffCut = 1.2;
float BetacutoffCut = 1.2;

float con_min=0.4;
float con_max=1.4;

TRandom3 * Rand= new TRandom3(time(0));

Particle proton(0.9382720813, 1, 1);  // proton mass 938 MeV
Particle deuton(1.8756129   , 1, 2);  // deuterium mass 1876 MeV, Z=1, A=2


Binning ForEffCorr(proton);
Binning ForEffCorr_D(deuton);
Binning HefragmToF(deuton);
Binning HefragmNaF(deuton);
Binning HefragmAgl(deuton);


Binning UnfoldingToF(proton);
Binning UnfoldingNaF(proton);
Binning UnfoldingAgl(proton);

Binning UnfoldingToF_D(deuton);
Binning UnfoldingNaF_D(deuton);
Binning UnfoldingAgl_D(deuton);

	
Binning DRB(deuton);
Binning PRB(proton);

Binning ForAcceptance(proton);
Binning ForCutoff(proton);  
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
	ForCutoff.Reset();
	UnfoldingToF.Reset();
	UnfoldingNaF.Reset();
	UnfoldingAgl.Reset();
	UnfoldingToF_D.Reset();
	UnfoldingNaF_D.Reset();
	UnfoldingAgl_D.Reset();
	
	cout<<"H.E. bins"<<endl;
	DRB.setBinsFromRDatacard ((workdir+"/bindatacard_PMIT.data").c_str(), 0.1, 0.9999999 ,ResponseTOF,0.00347548,5.8474); 
	PRB.setBinsFromRDatacard ((workdir+"/bindatacard_PMIT.data").c_str(), 0.1, 0.9999999 ,ResponseTOF,0.00347548,5.8474);
	ForAcceptance.setBinsFromRigidity(3*nbinsr,0.5,1250,ResponseTOF,0.00347548,5.8474);

	cout<<"Global Bins"<<endl;
	Global.setBinsFromRDatacard((workdir+"/bindatacard_Hepaper.data").c_str(),ResponseTOF,ResponseNaF,ResponseAgl);
	GlobalRig.setBinsFromRDatacard((workdir+"/bindatacard_Hepaper.data").c_str(),ResponseTOF,ResponseNaF,ResponseAgl);
	
	cout<<"TOF bins"<<endl;
	float ekmin=0.1, ekmax=0.82;
	float betamin=0.55; float betamax=0.853;
	HefragmToF.setBinsFromEkPerMass (4,0.15,0.504,ResponseTOF,0.00347548,5.8474);
	UnfoldingToF.setBinsFromRDatacard((workdir+"/bindatacard_Hepaper.data").c_str(), 0.555, 0.91 ,ResponseTOF,0.00347548,5.8474);
	UnfoldingToF_D.setBinsFromRDatacard((workdir+"/bindatacard_Hepaper.data").c_str(), 0.555, 0.91 ,ResponseTOF,0.00347548,5.8474);


	cout<<"NaF bins"<<endl;
	betamin=0.85, betamax=0.977;
	HefragmNaF.setBinsFromEkPerMass (1,1.5,3,ResponseNaF,0.00347548,5.8474);
//	UnfoldingNaF.setBinsFromRDatacard((workdir+"/bindatacard_Hepaper.data").c_str(), 0.813, 0.9878 ,ResponseNaF,0.00347548,5.8474);
//	UnfoldingNaF_D.setBinsFromRDatacard((workdir+"/bindatacard_Hepaper.data").c_str(), 0.813, 0.9878 ,ResponseNaF,0.00347548,5.8474);
	UnfoldingNaF.setBinsFromRDatacard((workdir+"/bindatacard_Hepaper.data").c_str(), 0.75, 0.9878 ,ResponseNaF,0.00347548,5.8474);
	UnfoldingNaF_D.setBinsFromRDatacard((workdir+"/bindatacard_Hepaper.data").c_str(), 0.75, 0.9878 ,ResponseNaF,0.00347548,5.8474);



	cout<<"Agl bins"<<endl;
	betamin=0.97, betamax=0.995;
	HefragmAgl.setBinsFromEkPerMass (2,2.6,11.9,ResponseAgl,0.00347548,5.8474);
//	UnfoldingAgl.setBinsFromRDatacard((workdir+"/bindatacard_Hepaper.data").c_str(), 0.953,0.9977 ,ResponseNaF,0.00347548,5.8474);
//	UnfoldingAgl_D.setBinsFromRDatacard((workdir+"/bindatacard_Hepaper.data").c_str(), 0.953,0.9977 ,ResponseNaF,0.00347548,5.8474);
	UnfoldingAgl.setBinsFromRDatacard((workdir+"/bindatacard_Hepaper.data").c_str(), 0.945,0.9977 ,ResponseNaF,0.00347548,5.8474);
	UnfoldingAgl_D.setBinsFromRDatacard((workdir+"/bindatacard_Hepaper.data").c_str(), 0.945,0.9977 ,ResponseNaF,0.00347548,5.8474);




	ForCutoff.setBinsFromRigidity(100,0.5,50,ResponseTOF,0.00347548,5.8474);
	ForCutoff.UseREdges();

	cout<<"**H.E.**"<<endl;
	PRB.Print();
	cout<<"**TOF**"<<endl;
	GlobalRig.GetToFPBins().Print();
	GlobalRig.GetToFDBins().Print();


	cout<<"**NaF**"<<endl;
	GlobalRig.GetNaFPBins().Print();
	GlobalRig.GetNaFDBins().Print();


	cout<<"**Agl**"<<endl;
	GlobalRig.GetAglPBins().Print();
	GlobalRig.GetAglDBins().Print();


	cout<<"**Global**"<<endl;
	GlobalRig.GetGlobalPBins().Print();
	GlobalRig.GetGlobalDBins().Print();

	cout<<"*** For Unfolding ***"<<endl;

	cout<<"**TOF**"<<endl;
	UnfoldingToF.Print();
	
	cout<<"**NaF**"<<endl;
	UnfoldingNaF.Print();
	
	cout<<"**Agl**"<<endl;
	UnfoldingAgl.Print();
	
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

