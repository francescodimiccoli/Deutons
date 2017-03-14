#include "TFile.h"
#include "TH2.h"
#include "TF2.h"
#include <TVector3.h>
#include "TMath.h"
#include <TFile.h>
#include "TFile.h"
#include "TH2.h"
#include "TF2.h"
#include <TVector3.h>
#include "TMath.h"
#include <TMVA/Reader.h>
#include <TMVA/Tools.h>

double bin[44];
float tempobin[44]= {0};
float tempobingeo[44][11];
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
std::vector<float> * trtot_edep = 0;
std::vector<float> * trtrack_edep = 0;
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
float residuoY=0;
float residuoX=0;
int fuoriY=0;
int fuoriX=0;
int clusterTOFfuori=0;
int clusterTrackfuori=0;
int clusterTRDfuori=0;
float Beta=0;
float BetaR=0;
float Beta_corr=0;
float BetaRICH=-1;
float BetaRICH_new=-1;
int RICHmask=0;
int RICHmask_new=0;
float Richtotused=0;
float RichPhEl=0;
float Massa=0;
float EdepTRD=0;
int NTRDclusters=0;
float Livetime;
int a,b,c,d,e,f,g,h,l,m,n,o,p,q,r,s,t,u,v=0;
int a1,b1,c1,d1,e1=0;
int alpha,Gamma,delta,alpha1,gamma1,delta1=0;
int a2=0;
float EndepTOF;
float Chisquare=0;
float endepostatrack=0;
int NTrackHits=0;
int clusterTrack = 0;
int clustertotTrack = 0;
int protoni[43][11];
int prescaled_protoni[43][11];
int prescaled_selezioni[10][43][11];
int prescaled_protoniMC[43];
int prescaled_selezioniMC[10][43];
int primari[43][11];
int primari_RICH_MC[43];
int primari_RICH[43];
int PROTONI[43];
int DEUTONI[43];
int primari_D[43][11];
int primari_bkgnd[43][11];
int deutoni[43][11];
int background[43];
int Template_D[43];
int background_He[43];
int background_el[43];
int background_pos[43];
int preselezionate[43][11]= {{0}};
int preselezionate_D[43][6]= {0};
int pres_unbias[43][11];
int preselected[43][11];
int quality[43][11];
int quality_D[43];
int Unbiascounts[43];
int Unbias_geo[43][11];
double geomag[12]= {0,0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.3};
int tbeg,tend;
double tempozona[11]= {0};
int contasecondi[11]= {0};
int zona;
float Time[11]= {0};
float Massa_gen,Momento_gen,Beta_gen;
int efficienzagen[43];
int efficienzagen_D[43][6]= {{0}};
int efficienzagen_He[43];
int efficienzagen_el[43];
int efficienzagen_pos[43];
float Etot=0;
float R_corr=0;
int rig;
int bet;
int Particle_ID=0;
double probP,probD;
double matricervsbeta[100][100];
double matricervsbetad[100][100];
float ThetaTOF=0;
float PhiTOF=0;
float ThetaTrack=0;
float PhiTrack=0;
float ThetaTRD=0;
float PhiTRD=0;
float PhiProjTOF=0;
float PhiProjTrack=0;
float PhiProjTRD=0;
float CooTOF[3]= {0,0,0};
float CooTrack[3]= {0,0,0};
float CooTRD[3]= {0,0,0};
float RGDT=0;
float BT=0;
float distR=0;
float distB=0;
float distETOFU=0;
float distETOFD=0;
float distETrack=0;
float distETRD=0;
float DistTOF=0;
float DistTOFU=0;
float DistTOFD=0;
float DistTrack=0;
float DistTRD=0;
float Dist=0;
float Dist_P=0;
float Dist1=1000000;
float Dist2=1000000;
float Dist3=1000000;
float CooTOF_p[3]= {0,0,0};
float CooTrack_p[3]= {0,0,0};
float CooTRD_p[3]= {0,0,0};
float Rmin=0;
float Rmin_P=0;
float RminTOF=0;
float RminTrack=0;
float RminTRD=0;
float DistTOF_P=DistTOF;
float DistTrack_P=DistTrack;
float DistTRD_P=DistTRD;
float PhiProjTOF_P=0;
float PhiProjTrack_P=0;
float PhiProjTRD_P=0;
float pi=3.14159265358979312e+00;
float EdepTOFU=0;
float EdepTOFD=0;
float EdepTrack=0;
float E_depTRD=0;
int DR1,DR2,DR3=0;
float PSCAL=0;
float PSCALTOF=0;
float PSCALTrack=0;
float PSCALTRD=0;
float PSCALTOF2=0;
float IsCharge1=0;
float PSCALTrack2=0;
float PSCALTRD2=0;
float PSCALTOF3=0;
float PSCALTrack3=0;
float PSCALTRD3=0;
float PSCALTOF4=0;
float PSCALTrack4=0;
float PSCALTRD4=0;
float passo=0.01;
float Dist5D=0;
float Dist5D_P=0;
long double L_fakeD=0;
long double L_trueD=0;
long double LDiscriminant=0;
double LSignal=0;
double LBkgnd=0;
float Unbias=0;
long double trueD=0;
long double fakeD=0;
float Pres_Unbias=0;
float Preselected=0;
int IsPrescaled=0;
int CUTMASK=0;
float Beta_pre=0;
float R_pre=0;
float TOF_Up_Down=0;
float DiffTrackEdep=0;
float Track_Up_Down=0;
float NTofUsed=0;
float diffR=0;
float ResX[7],ResY[7],trtoted[9],trtred[9],EdepLayer[4];
float nAnticluster;
float Layernonusati;
float XTOF,XTrack,XTRD,YTOF,YTrack,YTRD,ZTOF,ZTrack,ZTRD=0;
float PDiscriminant=0;
int Cutmask=0;
int UnbiasPre=9;
int noR=33;
float Velocity=0;
int seconds=0;
float qL1,qInner,qLtof,qUtof=0;
double  R_L1=0;

TMVA::Reader *reader;
Float_t BDT_response;
float EdepECAL=0;




using namespace std;

void BDTreader()
{
   TMVA::Tools::Instance();
   reader = new TMVA::Reader( "V:Color:!Silent" );
   reader->AddVariable("NAnticluster", &nAnticluster);
   reader->AddVariable("Chisquare",&Chisquare);
   reader->AddVariable("layernonusati",&Layernonusati);
   reader->AddVariable("NTofUsed := NTofClusters - NTofClustersusati",&NTofUsed);
   reader->AddVariable("diffR := TMath::Abs(Rup-Rdown)/R",&diffR);
   reader->AddVariable("TOF_Up_Down := TMath::Abs(Endep[2]+Endep[3]-Endep[0]-Endep[1])", &TOF_Up_Down);
   reader->BookMVA("BDTmethod", "/storage/gpfs_ams/ams/users/fdimicco/Deutons/TMVA/QualityBDT_BDT.weights.xml");
}

void Eval_BDT(TTree *albero,int i){
   for(int layer=0; layer<9; layer++) trtoted[layer]=(*trtot_edep)[layer];
   for(int layer=0; layer<9; layer++) trtred[layer]=(*trtrack_edep)[layer];
   for(int layer=0; layer<4; layer++)  EdepLayer[layer]=(*Endep)[layer];
   diffR=fabs(Rup-Rdown)/R;
   NTofUsed=NTofClusters-NTofClustersusati;
   Layernonusati=layernonusati;
   nAnticluster=NAnticluster;
   BDT_response=reader->EvaluateMVA("BDTmethod");
   
}

bool Quality(TTree *albero,int i)
{
   albero->GetEvent(i);
   bool selection = true;
   TOF_Up_Down=fabs(((*Endep)[2]+(*Endep)[3])-((*Endep)[0]+(*Endep)[1]));
   DiffTrackEdep=0;
   for(int layer=0; layer<1; layer++) DiffTrackEdep+=fabs((*trtot_edep)[layer]-(*trtrack_edep)[layer]);

   ///////Matter or Antimatter
   if (R<0) return false;
   ////////////////////////////
   
   R=fabs(R);
   /////////////////////////// Calcolo LIKELIHOOD ////////////////////////////

   {
	   float var[9]= {0,0,0,0,0,0,0,0,0};
	   var[0]=NAnticluster;
	   var[1]=(NTofClusters-NTofClustersusati);
	   if(R!=0) var[2]=(fabs(Rup-Rdown)/R);
	   else var[2]=1;
	   var[3]=layernonusati;
	   var[4]=fuoriX;
	   var[5]=Chisquare;
	   var[6]=TOF_Up_Down;
	   var[7]=Track_Up_Down;
	   var[8]=DiffTrackEdep;
	   double Ltrue=1;
	   double Lfalse=1;
	   for(int m=0; m<6; m++) {
		   if(m!=4){
			   Lfalse=Lfalse*Bkgnd[m]->Eval(var[m]);
			   Ltrue=Ltrue*Signal[m]->Eval(var[m]);
		   }
	   }
	   LDiscriminant=Ltrue/(Ltrue+Lfalse);
   }
   if((((int)Cutmask)>>11)==512||(((int)Cutmask)>>11)==0) {
	   float var[9]= {0,0,0,0,0,0,0,0,0};
	   var[0]=NAnticluster;
	   var[1]=(NTofClusters-NTofClustersusati);
	   if(R!=0) var[2]=(fabs(Rup-Rdown)/R);
	   else var[2]=1;
	   var[3]=layernonusati;
	   var[4]=fuoriX;
	   var[5]=Chisquare;
	   var[6]=Richtotused;
	   var[7]=RichPhEl;
	   var[8]=DiffTrackEdep;
	   double Ltrue=1;
	   double Lfalse=1;
	   for(int m=0; m<8; m++) {
		   if(m!=4){
			   if((((int)Cutmask)>>11)==512) {
				   Lfalse=Lfalse*BkgndNaF[m]->Eval(var[m]);
				   Ltrue=Ltrue*SignalNaF[m]->Eval(var[m]);
			   }
			   if((((int)Cutmask)>>11)==0) {
				   Lfalse=Lfalse*BkgndAgl[m]->Eval(var[m]);
				   Ltrue=Ltrue*SignalAgl[m]->Eval(var[m]);
			   }
		   }	
	   }
	   LDiscriminant=Ltrue/(Ltrue+Lfalse);
   }

   ///////////////////////////////////////////////////////////////////////////
   
   /////////////////////////// Variables ///////////////////////////
   EdepTrack=0;
   EdepTOFU=((*Endep)[0]+(*Endep)[1])/2;
   EdepTOFD=((*Endep)[2]+(*Endep)[3])/2;

   for(int layer=1; layer<8; layer++) EdepTrack+=(*trtrack_edep)[layer];
   EdepTrack=EdepTrack/7;
   EndepTOF=((*Endep)[0]+(*Endep)[1]+(*Endep)[2]+(*Endep)[3])/4;
   
   R_corr=R;
   Massa=pow(fabs(pow(fabs(R_corr)*pow((1-pow(Beta,2)),0.5)/Beta,2)),0.5);
   
   if(qInner<1.5&&qUtof<1.5&&qLtof<1.5) IsCharge1=1;
   
   Velocity=0;
   Velocity=Beta;
   if((((int)Cutmask)>>11)==0||(((int)Cutmask)>>11)==512) Velocity=BetaRICH_new;
   
   //////////////////////////////////////////////////////////////////////

   return true;

}

void Eval_Distance(float RG,float M, TF1 *RBETA)
{
	RGDT=0;
	BT=0;
	distR=0;
	distB=0;
	distETOFU=0;
	distETOFD=0;
	distETrack=0;
	distETRD=0;
	DistTOFU=0;
	DistTOFD=0;
	DistTrack=0;
	DistTRD=0;
	DR1=0;
	DR2=0;
	DR3=0;
	for(int P=0; P<3; P++) CooTOF[P]=0;
	for(int P=0; P<3; P++) CooTrack[P]=0;
	for(int P=0; P<3; P++) CooTRD[P]=0;
	//////////////////// CALCOLO DISTANZA /////////////////////////
	Dist1=1000000;
	Dist2=1000000;
	Dist3=1000000;
	for(int z=0; z<1e6; z++) {
		passo=0.075;
		BT=RBETA->Eval(RGDT);
		distR=(RGDT-RG)/(pow(RGDT,2)*Rig->Eval(RGDT));
		distB=(BT-Beta)/(pow(BT,2)*beta->Eval(BT));
		if((((int)Cutmask)>>11)==512) distB=(BT-BetaRICH_new)/(pow(BT,2)*betaNaF->Eval(BT));
		if((((int)Cutmask)>>11)==0) distB=(BT-BetaRICH_new)/(pow(BT,2)*betaAgl->Eval(BT));
		distETOFU=(EdepTOFbeta->Eval(BT)-EdepTOFU)/(pow(EdepTOFbeta->Eval(BT),2)*etofu->Eval(BT));
		distETrack=(EdepTrackbeta->Eval(BT)-EdepTrack)/(pow(EdepTrackbeta->Eval(BT),2)*etrack->Eval(BT));
		distETOFD=(EdepTOFbeta->Eval(BT)-EdepTOFD)/(pow(EdepTOFbeta->Eval(BT),2)*etofd->Eval(BT));

		//std::cout<<Dist<<" "<<z<<" : "<<R<<" "<<Rmin<<" "<<RGDT<<std::endl;

		Dist=pow(pow(distR,2)+pow(distB,2)+pow(distETrack,2)+pow(distETOFU,2)+pow(distETOFD,2),0.5);
		if(Dist<Dist1) {
			DR1=0;
			Dist1=Dist;
			CooTOF[0]=distETOFU;
			CooTOF[1]=distB;
			CooTOF[2]=distR;
			CooTrack[0]=distETrack;
			CooTrack[1]=distB;
			CooTrack[2]=distR;
			CooTRD[0]=distETOFD;
			CooTRD[1]=distB;
			CooTRD[2]=distR;
			Rmin=RGDT;
		} else DR1++;
		if(!((((int)Cutmask)>>11)==512||(((int)Cutmask)>>11)==0)) {
			if(DR1>25) break;
			RGDT=RGDT+passo;
		}
		if(((((int)Cutmask)>>11)==512||(((int)Cutmask)>>11)==0)) {
			if(fabs(Dist-Dist2)<0.0001) DR2++;
			if(DR2>4||DR1>3000) {
				//	if(R<20&&(((int)Cutmask)>>11)==0) cout<<" "<<Dist1<<" R "<<R<<" Rmin "<<Rmin<<" dR "<<CooTOF[2]<<" dB "<<CooTOF[1]<<" "<<distETOFU<<" "<<distETrack<<" "<<distETOFD<<endl;
				break;
			}
			RGDT=RGDT+passo;
		}
      Dist2=Dist;
      if(z>1e5) std::cout<<"cazzo"<<std::endl;
   }
   Dist5D=Dist1;
     //{cout<<Rmin<<" : "<<R<<" "<<Dist5D<<endl;
    // cout<<endl;}
   /////////////////////////////////////////////////////////////////

   PhiTOF=acos(CooTOF[2]/pow(pow(CooTOF[0],2)+pow(CooTOF[2],2),0.5));
   if(CooTOF[0]<0) PhiTOF=2*pi-PhiTOF;
   ThetaTOF=acos(CooTOF[1]/Dist1);
   if(CooTOF[2]<0) ThetaTOF=2*pi-ThetaTOF;
   PhiTrack=acos(CooTrack[2]/pow(pow(CooTrack[0],2)+pow(CooTrack[2],2),0.5));
   if(CooTrack[0]<0) PhiTrack=2*pi-PhiTrack;
   ThetaTrack=acos(CooTrack[1]/Dist2);
   if(CooTrack[2]<0) ThetaTrack=2*pi-ThetaTrack;
   PhiTRD=acos(CooTRD[2]/pow(pow(CooTRD[0],2)+pow(CooTRD[2],2),0.5));
   if(CooTRD[0]<0) PhiTRD=2*pi-PhiTRD;
   ThetaTRD=acos(CooTRD[1]/Dist3);
   if(CooTRD[2]<0) ThetaTRD=2*pi-ThetaTRD;
   DistTOF=pow(pow(CooTOF[0],2)+pow(CooTOF[1],2)+pow(CooTOF[2],2),0.5);
   DistTrack=pow(pow(CooTrack[0],2)+pow(CooTrack[1],2)+pow(CooTrack[2],2),0.5);
   DistTRD=pow(pow(CooTRD[0],2)+pow(CooTRD[1],2)+pow(CooTRD[2],2),0.5);
   /////////////////////////// VETTORE TANGENTE ///////////////////////////////////
   float TanTOF[3]= {0,0,0};
   float TanTrack[3]= {0,0,0};
   float TanTRD[3]= {0,0,0};
   //passo=passo*5;
   //Rmin=RminTOF;
   TanTOF[0]=(EdepTOFbeta->Eval(RBETA->Eval(Rmin+passo))-EdepTOFbeta->Eval(RBETA->Eval(Rmin)))/(pow(EdepTOFbeta->Eval(RBETA->Eval(Rmin)),2)*etofu->Eval(RBETA->Eval(Rmin)));
   //Rmin=RminTrack;
   TanTrack[0]=(EdepTrackbeta->Eval(RBETA->Eval(Rmin+passo))-EdepTrackbeta->Eval(RBETA->Eval(Rmin)))/(pow(EdepTrackbeta->Eval(RBETA->Eval(Rmin)),2)*etrack->Eval(RBETA->Eval(Rmin)));
   //Rmin=RminTRD;
   //TanTRD[0]=(EdepTRDbeta->Eval(RBETA->Eval(Rmin+passo))-EdepTRDbeta->Eval(RBETA->Eval(Rmin)))/(pow(EdepTRDbeta->Eval(RBETA->Eval(Rmin)),2)*etrd->Eval(RBETA->Eval(Rmin)));
   TanTRD[0]=(EdepTOFbeta->Eval(RBETA->Eval(Rmin+passo))-EdepTOFbeta->Eval(RBETA->Eval(Rmin)))/(pow(EdepTOFbeta->Eval(RBETA->Eval(Rmin)),2)*etofd->Eval(RBETA->Eval(Rmin)));
   //Rmin=RminTOF;

   TanTOF[1]=(RBETA->Eval(Rmin+passo)-RBETA->Eval(Rmin))/(pow(RBETA->Eval(Rmin),2)*beta->Eval(RBETA->Eval(Rmin)));
   //Rmin=RminTrack;
   TanTrack[1]=(RBETA->Eval(Rmin+passo)-RBETA->Eval(Rmin))/(pow(RBETA->Eval(Rmin),2)*beta->Eval(RBETA->Eval(Rmin)));
   //Rmin=RminTRD;
   //TanTRD[1]=(RBETA->Eval(Rmin+passo)-RBETA->Eval(Rmin))/(pow(RBETA->Eval(Rmin),2)*beta->Eval(RBETA->Eval(Rmin)));
   TanTRD[1]=(RBETA->Eval(Rmin+passo)-RBETA->Eval(Rmin))/(pow(RBETA->Eval(Rmin),2)*beta->Eval(RBETA->Eval(Rmin)));

   TanTOF[2]=passo/(pow(Rmin,2)*Rig->Eval(Rmin));
   TanTrack[2]=passo/(pow(Rmin,2)*Rig->Eval(Rmin));
   //TanTRD[2]=passo/(pow(RminTRD,2)*Rig->Eval(RminTRD));
   TanTRD[2]=passo/(pow(Rmin,2)*Rig->Eval(Rmin));

   for(int i=0; i<3; i++) {TanTOF[i]=3*TanTOF[i]; TanTrack[i]=3*TanTrack[i]; TanTRD[i]=3*TanTRD[i];}
   PSCALTOF=CooTOF[0]*TanTOF[0]+CooTOF[1]*TanTOF[1]+CooTOF[2]*TanTOF[2];
   PSCALTrack=CooTrack[0]*TanTrack[0]+CooTrack[1]*TanTrack[1]+CooTrack[2]*TanTrack[2];
   PSCALTRD=CooTRD[0]*TanTRD[0]+CooTRD[1]*TanTRD[1]+CooTRD[2]*TanTRD[2];

   float normaCooTOF=0;
   float normaCooTrack=0;
   float normaCooTRD=0;
   float normaTanTOF=0;
   float normaTanTrack=0;
   float normaTanTRD=0;
   for(int i=0; i<3; i++) {
      normaCooTOF=normaCooTOF+CooTOF[i]*CooTOF[i];
      normaCooTrack=normaCooTrack+CooTrack[i]*CooTrack[i];
      normaCooTRD=normaCooTRD+CooTRD[i]*CooTRD[i];
      normaTanTOF=normaTanTOF+TanTOF[i]*TanTOF[i];
      normaTanTrack=normaTanTrack+TanTrack[i]*TanTrack[i];
      normaTanTRD=normaTanTRD+TanTRD[i]*TanTRD[i];
   }

   normaCooTOF=pow(normaCooTOF,0.5);
   normaCooTrack=pow(normaCooTrack,0.5);
   normaCooTRD=pow(normaCooTRD,0.5);
   normaTanTOF=pow(normaTanTOF,0.5);
   normaTanTrack=pow(normaTanTrack,0.5);
   normaTanTRD=pow(normaTanTRD,0.5);

   PSCAL=CooTOF[0]*TanTOF[0]+CooTrack[0]*TanTrack[0]+CooTRD[0]*TanTRD[0]+ CooTOF[1]*TanTOF[1]+CooTOF[2]*TanTOF[2] ;
   PSCALTOF=acos(PSCALTOF/(DistTOF*normaTanTOF));
   PSCALTrack=acos(PSCALTrack/(DistTrack*normaTanTrack));
   PSCALTRD=acos(PSCALTRD/(DistTRD*normaTanTRD));

   float normaCoo=0;
   float normaTan=0;
   normaCoo=CooTOF[0]*CooTOF[0]+CooTrack[0]*CooTrack[0]+CooTRD[0]*CooTRD[0]+CooTOF[1]*CooTOF[1]+CooTOF[2]*CooTOF[2];
   normaTan=TanTOF[0]*TanTOF[0]+TanTrack[0]*TanTrack[0]+TanTRD[0]*TanTRD[0]+TanTOF[1]*TanTOF[1]+TanTOF[2]*TanTOF[2];
   normaCoo=pow(normaCoo,0.5);
   normaTan=pow(normaTan,0.5);
   PSCAL=acos(PSCAL/(Dist*normaTan));

   //cout<<PSCAL<<" "<<PSCALTOF<<" "<<PSCALTrack<<" "<<PSCALTRD<<endl;
   //////////////////////////////////////////////////////////////////////////
   //
   //
   ////////////////////////////// PROIEZIONE P. PERP. (TRot.)////////////////
   float AsseBetaTOF[3]= {0,1,0};
   float AsseBetaTrack[3]= {0,1,0};
   float AsseBetaTRD[3]= {0,1,0};
   float PVettTOF[3]= {0,0,0};
   float PVettTrack[3]= {0,0,0};
   float PVettTRD[3]= {0,0,0};
   PVettTOF[0]=AsseBetaTOF[1]*TanTOF[2]-AsseBetaTOF[2]*TanTOF[1];
   PVettTOF[1]=AsseBetaTOF[0]*TanTOF[2]-AsseBetaTOF[2]*TanTOF[0];
   PVettTOF[2]=AsseBetaTOF[0]*TanTOF[1]-AsseBetaTOF[1]*TanTOF[0];
   PVettTrack[0]=AsseBetaTrack[1]*TanTrack[2]-AsseBetaTrack[2]*TanTrack[1];
   PVettTrack[1]=AsseBetaTrack[0]*TanTrack[2]-AsseBetaTrack[2]*TanTrack[0];
   PVettTrack[2]=AsseBetaTrack[0]*TanTrack[1]-AsseBetaTrack[1]*TanTrack[0];
   PVettTRD[0]=AsseBetaTRD[1]*TanTRD[2]-AsseBetaTRD[2]*TanTRD[1];
   PVettTRD[1]=AsseBetaTRD[0]*TanTRD[2]-AsseBetaTRD[2]*TanTRD[0];
   PVettTRD[2]=AsseBetaTRD[0]*TanTRD[1]-AsseBetaTRD[1]*TanTRD[0];
   PSCALTOF2=TanTOF[0]*AsseBetaTOF[0]+TanTOF[1]*AsseBetaTOF[1]+TanTOF[2]*AsseBetaTOF[2];
   PSCALTrack2=TanTrack[0]*AsseBetaTrack[0]+TanTrack[1]*AsseBetaTrack[1]+TanTrack[2]*AsseBetaTrack[2];
   PSCALTRD2=TanTRD[0]*AsseBetaTRD[0]+TanTRD[1]*AsseBetaTRD[1]+TanTRD[2]*AsseBetaTRD[2];
   TVector3 R_TOF(PVettTOF[0],PVettTOF[1],PVettTOF[2]);
   TVector3 R_Track(PVettTrack[0],PVettTrack[1],PVettTrack[2]);
   TVector3 R_TRD(PVettTRD[0],PVettTRD[1],PVettTRD[2]);
   TVector3 X_TOF(0,0,1);
   TVector3 X_Track(0,0,1);
   TVector3 X_TRD(0,0,1);
   X_TOF.Rotate((acos(PSCALTOF2/normaTanTOF)),R_TOF);
   X_Track.Rotate((acos(PSCALTrack2/normaTanTrack)),R_Track);
   X_TRD.Rotate((acos(PSCALTRD2/normaTanTRD)),R_TRD);
   TVector3 Y_TOF(1,0,0);
   TVector3 Y_Track(1,0,0);
   TVector3 Y_TRD(1,0,0);
   Y_TOF.Rotate((acos(PSCALTOF2/normaTanTOF)),R_TOF);
   Y_Track.Rotate((acos(PSCALTrack2/normaTanTrack)),R_Track);
   Y_TRD.Rotate((acos(PSCALTRD2/normaTanTRD)),R_TRD);
   TVector3 Z_TOF(0,1,0);
   TVector3 Z_Track(0,1,0);
   TVector3 Z_TRD(0,1,0);
   Z_TOF.Rotate((acos(PSCALTOF2/normaTanTOF)),R_TOF);
   Z_Track.Rotate((acos(PSCALTrack2/normaTanTrack)),R_Track);
   Z_TRD.Rotate((acos(PSCALTRD2/normaTanTRD)),R_TRD);
   PSCALTOF2=CooTOF[0]*X_TOF(0)+CooTOF[1]*X_TOF(1)+CooTOF[2]*X_TOF(2);
   PSCALTrack2=CooTrack[0]*X_Track(0)+CooTrack[1]*X_Track(1)+CooTrack[2]*X_Track(2);
   PSCALTRD2=CooTRD[0]*X_TRD(0)+CooTRD[1]*X_TRD(1)+CooTRD[2]*X_TRD(2);
   PSCALTOF3=CooTOF[0]*Y_TOF(0)+CooTOF[1]*Y_TOF(1)+CooTOF[2]*Y_TOF(2);
   PSCALTrack3=CooTrack[0]*Y_Track(0)+CooTrack[1]*Y_Track(1)+CooTrack[2]*Y_Track(2);
   PSCALTRD3=CooTRD[0]*Y_TRD(0)+CooTRD[1]*Y_TRD(1)+CooTRD[2]*Y_TRD(2);
   PSCALTOF4=CooTOF[0]*Z_TOF(0)+CooTOF[1]*Z_TOF(1)+CooTOF[2]*Z_TOF(2);
   PSCALTrack4=CooTrack[0]*Z_Track(0)+CooTrack[1]*Z_Track(1)+CooTrack[2]*Z_Track(2);
   PSCALTRD4=CooTRD[0]*Z_TRD(0)+CooTRD[1]*Z_TRD(1)+CooTRD[2]*Z_TRD(2);
   ////////////////////////////////////////////////////////////////
   //cout<<"Diff "<<(fabs(PSCALTOF2-PSCALTrack2)+fabs(PSCALTOF2-PSCALTRD2)+fabs(PSCALTRD2-PSCALTrack2))/fabs(PSCALTrack2)<<" "<<Dist<<endl;
   RminTOF=Rmin;
   RminTrack=Rmin;
   RminTRD=Rmin;
   if(M<1) {DistTOF_P=DistTOF; DistTrack_P=DistTrack; DistTRD_P=DistTRD; Dist5D_P=Dist5D; }
   return;
}



bool Protoni (TTree *albero,int i)
{

   clusterTOFfuori=0;
   clusterTrackfuori=0;
   clusterTRDfuori=0;
   albero->GetEvent(i);
   Eval_Distance(R,0.938,protons);

   return true;
}

bool Deutoni (TTree *albero,int i)
{
   clusterTOFfuori=0;
   clusterTrackfuori=0;
   clusterTRDfuori=0;
   albero->GetEvent(i);
   Eval_Distance(R,1.875,deutons);
   
   return true;
}


