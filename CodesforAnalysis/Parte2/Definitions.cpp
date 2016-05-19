using namespace std;

#include "Binning.h"

//// Input-Output Variables

string inputpath="/storage/gpfs_ams/ams/users/fdimicco/Deutons";
string outputpath;
string mese;
string frac;


//// Global Variables
float R=0;
float Beta=0;
float BetaRICH=0;
float RminTOF=0;
float RminTrack=0;
float RminTRD=0;
float XTOF=0;
float XTrack=0;
float XTRD=0;
float YTOF=0;
float YTrack=0;
float YTRD=0;
float Rcutoff=0;
float LDiscriminant=0;
float Massa_gen=0;
float MC_type=0;
float Massa=0;
float BDT_response=0;
float D_TOF,D_Track,D_TRD,Discr=0;
float Zona=0;
float CUTMASK=0;
float IsPrescaled=0;
float Latitude=0;
float EdepL1=0;
float Rmin=0;
float X=0;
float YTOFU=0;
float YTOFD=0;
int Cutmask=0;
float Dist5D=0;
float Dist5D_P=0;
float Momento_gen=0;
float Ev_Num=0;
float Trig_Num=0;
float R_pre=0;
float Beta_gen=0;
float Beta_pre=0;
float EdepECAL=0;
float EdepTOFU=0;
float EdepTOFD=0;
float EdepTrack=0;
float BetaRICH_new=0;
double geomag[12]= {0,0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.3};
float Rcut[11]= {18,18,16,14,12,10,8,6,4,2,1};
int INDX=0;
int FRAC=0;
double geomagC[11]= {0,0.05,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.15};
int notpassed[3]= {155,139,11};
int passed[3]= {187,155,139};
float B=0.4;
float B1=0;
float B2=0;
float E=0.1;
int binnum=1;
float distcut=4;
float ddiscrcut=0.15;
bool qualcut=true;
float a= (log (0.9)-log (0.1) ) /nbinsbeta;
float E2=exp (log (0.1)+1.5*a);


float Betacent      [nbinsbeta] = {0};
float Ekincent      [nbinsbeta] = {0};

float BetacentNaF   [nbinsNaF]  = {0};
float EkincentNaF   [nbinsNaF]  = {0};
float BetacentAgl   [nbinsAgl]  = {0};
float EkincentAgl   [nbinsAgl]  = {0};

float Unbias=0;


std::array <float, nbinsr+1> Rbins;


double R_cent[nbinsr];
float encinprot     [nbinsr];
float encindeut     [nbinsr];
float deltaencinprot[nbinsr];
float deltaencindeut[nbinsr];
float deltaencinTOF [nbinsToF];
float deltaencinNaF [nbinsNaF];
float deltaencinAgl [nbinsAgl];
float Var=0;
float Var2=0;
float Var3=0;
int Norm[11]= {0};

//cuts
bool Likcut=false;
bool Distcut=false;
bool Herejcut=false;
bool Betastrongcut=false;

TH1F * Esposizione[10];
TH1F * Tempi;
TH2F * esposizionegeo;
TH2F * esposizionepgeo;
TH2F * esposizionepgeoNaF;
TH2F * esposizionepgeoAgl;
TH2F * esposizionedgeo;
TH2F * esposizionedgeoNaF;
TH2F * esposizionedgeoAgl;

//retrieve MC particle species
float ReturnMass_Gen()
{
   float Mass_gen=0;
   if ( ( ( (int) MC_type) &0xFF    ) >0)      Mass_gen = 0.938;
   if ( ( ( (int) MC_type) &0xFF00  ) >0)      Mass_gen = 1.875;
   if ( ( ( (int) MC_type) &0xFF0000) >0)      Mass_gen = 3.725;
   return Mass_gen;
}

//retrieve MC cross section type
int ReturnMCGenType()
{
   int mc_type=-1;
   int cursor=0;
   if (ReturnMass_Gen() <1&&ReturnMass_Gen() >0) cursor=0 ;
   if (ReturnMass_Gen() <2&&ReturnMass_Gen() >1) cursor=8 ;
   if (ReturnMass_Gen() <4&&ReturnMass_Gen() >3) cursor=16;
   for (int i=2; i<8; i++) {
      if ( ( ( ( (int) MC_type) >> (cursor+i) ) & 1 ) ==1) mc_type=i-2;
   }
   if (mc_type == -1) std::cout<<"ERROR: MC cross section type not found"<<std::endl;
   return mc_type;
}


void FillBinMGen (TH1* h, int bin)
{
   int mass = ReturnMCGenType();
   h->Fill (bin, mass);
   return;
}

void FillBinMGen (TH3* h, int bin, int S)
{
   int mass = ReturnMCGenType();
   h->Fill (bin, mass, S);
   return;
}



DBinning ToFDB;
PBinning ToFPB;
DBinning NaFDB;
PBinning NaFPB;
DBinning AglDB;
PBinning AglPB;

Binning RB;



