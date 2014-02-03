#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include <TVector3.h>
#include <fstream>
#include <sstream>
#include <math.h>
#include <cstring>
#include <vector>
#include "TMath.h"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

double bin[24];
int U_time;
float Latitude;
int zonageo;
float Rcutoff;
float CaricaTOF=0;
float CaricaTRD=0;
float CaricaTrack=0;
float ProbQ=0;
float Qbest=0;
std::vector<float> * Endep = 0;
//float Endep[4];
int layernonusati=0;
int NAnticluster=0;
int NTRDSegments=0;
int NTofClusters=0;
int NTofClustersusati=0;
double Rup=0;
double Rdown=0;
double R=0;
std::vector<float> * chiq = 0;
std::vector<float> * R_ = 0;
//float chiq[6];
//float R_[6];
float residuoY=0;
float residuoX=0;
int fuoriY=0;
int fuoriX=0;
int clusterTOFfuori=0;
int clusterTrackfuori=0;
float Beta=0;
float Betacorr=0;
float BetaRICH=-1;
float Massa=0;
float EdepTRD=0;
int NTRDclusters=0;
float Livetime;
int a,b,c,d,e,f,g,h,l,m,n,o,p,q,r,s,t,u,v=0;
int a1,b1,c1,d1,e1=0;
int alpha,Gamma,delta,alpha1,gamma1,delta1=0;
int a2=0;
std::vector<float> * TRDclusters = 0;
std::vector<double> * ResiduiX;
std::vector<double> * ResiduiY;
float Chisquare=0; 
//float TRDclusters[1000];
float endepostatrack=0;
int NTrackHits=0;
std::vector<float> * clusterTrack = 0;
int protoni[23][11];
int deutoni[23][11];
int background[23];
int preselezionate[23][11];
int quality[23][11];
double geomag[12]={0,0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.3};
int tbeg,tend;
double tempozona[11]={0,0,0,0,0,0,0,0,0,0,0};
int contasecondi[11]={0,0,0,0,0,0,0,0,0,0,0};
int zona;
float Time[11]={0,0,0,0,0,0,0,0,0,0,0};
float Massa_gen,Momento_gen,Beta_gen;
int efficienzagen[23];


bool Quality(TTree *albero,int i)
{
	int k = albero->GetEvent(i);
	 bool selection = true;
//QUALITY
if(NAnticluster>=1)  selection = false; else f++;
if(NTRDSegments<1)  selection = false; else g++;
if((NTofClusters-NTofClustersusati)>1)  selection = false; else h++;
if(fabs(Rup-Rdown)/R>0.2)  selection = false; else l++;
if(ProbQ<0.43) selection = false; else m++;
//CONTROLLOFIT
a2=0;
for(int j=0;j<6;j++)
        if((*chiq)[j]<0||(*chiq)[j]>15||fabs((*R_)[j]-R)/R>0.2) a2++;
if(a2>0) selection = false; else n++;
fuoriX=0;
fuoriY=0;
for(int layer=0;layer<7;layer++) {
	if((*ResiduiX)[layer]<0) continue;
	else {
		residuoX=residuoX+(*ResiduiX)[layer];
		if((*ResiduiX)[layer]>0.0025) fuoriX++;
	     }
		if((*ResiduiY)[layer]<0) continue;
	else {
		residuoY=residuoY+(*ResiduiY)[layer];
		if((*ResiduiY)[layer]>0.0025) fuoriY++;
	     }	
}
if((7-layernonusati)!=0){
	residuoX=residuoX/(7-layernonusati);
	residuoY=residuoY/(7-layernonusati);
}
if(Chisquare>2.5) selection = false; else r++;
if(layernonusati>2) selection = false; else s++;
//CARICA
if(fabs(CaricaTOF-Qbest)>0.1) selection = false; else t++;
if(fabs(CaricaTRD-Qbest)>0.17) selection = false; else v++;
if(fabs(CaricaTrack-Qbest)>0.1) selection = false; else u++;
//=true: Disattiva selezione/////
selection=selection;
////////////////////////////////
return selection;

}

bool Protoni (TTree *albero,int i)
{
        bool isprot=true;
        clusterTOFfuori=0;
        clusterTrackfuori=0;
        int k = albero->GetEvent(i);
        float EndepTOF=((*Endep)[0]+(*Endep)[1]+(*Endep)[2]+(*Endep)[3])/4;
        //TOF
        for (int j=0;j<4;j++)
        if (fabs((*Endep)[j]-fverap->Eval(R))>1) clusterTOFfuori++;
        if(fabs(EndepTOF-fverap->Eval(R))>1) isprot=false;
        else a++;
        if(clusterTOFfuori>0) isprot=false;
        else b++;
        if(fabs(EndepTOF-EdepTOFbeta->Eval(Betacorr))>1) isprot=false;
	else alpha++;
	//TRACKER
        for (int j=0;j<NTrackHits;j++)
        if(fabs((*clusterTrack)[j]-fveratrp->Eval(R))>40) clusterTrackfuori++;
        if(fabs(endepostatrack/NTrackHits-fveratrp->Eval(R))>20) isprot=false;
        else c++;
        if(clusterTrackfuori>3) isprot=false;
        else d++;
	if(fabs(endepostatrack/NTrackHits-EdepTrackbeta->Eval(Betacorr))>20) isprot=false;
	else Gamma++;
        //TRD
        if(fabs(EdepTRD/NTRDclusters-fveraTRDp->Eval(R))>2) isprot=false;
        else e++;
	if(fabs(EdepTRD/NTRDclusters-EdepTRDbeta->Eval(Betacorr))>1.2) isprot=false;
	else delta++;
	//=true: Disattiva selezione/////
	//isprot=true;
	////////////////////////////////
        return isprot;
}

bool Deutoni (TTree *albero,int i)
{
        bool isprot=true;
        clusterTOFfuori=0;
        clusterTrackfuori=0;
        int k = albero->GetEvent(i);
        float EndepTOF=((*Endep)[0]+(*Endep)[1]+(*Endep)[2]+(*Endep)[3])/4;
        //TOF
        for (int j=0;j<4;j++)
        if (fabs((*Endep)[j]-fvera->Eval(R))>1) clusterTOFfuori++;
        if(fabs(EndepTOF-fvera->Eval(R))>1) isprot=false;
        else a1++;
        if(clusterTOFfuori>0) isprot=false;
        else b1++;
	if(fabs(EndepTOF-EdepTOFbeta->Eval(Betacorr))>1) isprot=false;
        else alpha1++;
        //TRACKER
        for (int j=0;j<NTrackHits;j++)
        if(fabs((*clusterTrack)[j]-fveratr->Eval(R))>40) clusterTrackfuori++;
        if(fabs(endepostatrack/NTrackHits-fveratr->Eval(R))>20) isprot=false;
        else c1++;
        if(clusterTrackfuori>3) isprot=false;
        else d1++;
        if(fabs(endepostatrack/NTrackHits-EdepTrackbeta->Eval(Betacorr))>20) isprot=false;
        else gamma1++;
	//TRD
        if(fabs(EdepTRD/NTRDclusters-fveraTRD->Eval(R))>2) isprot=false;
        else e1++;
	if(fabs(EdepTRD/NTRDclusters-EdepTRDbeta->Eval(Betacorr))>1.2) isprot=false;
        else delta1++;
	///SELEZIONE SU RvsBETA
	if(Betacorr>f1->Eval(R)&&BetaRICH<0) isprot=false;
        if (Betacorr<0.5&&Massa<1.7699) isprot= false;
        if((Betacorr>=0.5&&Betacorr<0.6)&&Massa<1.84) isprot= false;
        if((Betacorr>=0.6&&Betacorr<0.7)&&Massa<1.754) isprot=false;
        if((Betacorr>=0.7&&Betacorr<0.85)&&Massa<1.754) isprot=false;
        if(Massa>=2.5||Massa<1.67) isprot=false;

	//=true: Disattiva selezione/////
	//isprot=true;
	////////////////////////////////
        return isprot;
}

