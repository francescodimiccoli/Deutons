using namespace std;

#include "Binning.h"

//// Input-Output Variables

string inputpath="/storage/gpfs_ams/ams/users/fdimicco/Deutons";
string outputpath;
string mese;
string frac;



extern const int nbinsr=43;
extern const int nbinsbeta=18;
extern const int nbinsToF=18;
extern const int nbinsNaF=18;
extern const int nbinsAgl=18;



TF1 *protons = new TF1("f1","pow((pow(x,2)/pow(0.938,2)/(1 + pow(x,2)/pow(0.938,2))),0.5)",0.1,100);
TF1 *deutons = new TF1("f1","pow((pow(x,2)/pow(1.875,2)/(1 + pow(x,2)/pow(1.875,2))),0.5)",0.1,100);


////////////// DEFINIZIONE SPLINES //////////////////
TSpline3 *Rig;
TSpline3 *beta;
TF1 *betaNaF;
TF1 *betaAgl;
TSpline3 *eL1;
TSpline3 *etofu;
TSpline3 *etrack;
TSpline3 *etofd;
TSpline3 *EdepL1beta;
TSpline3 *EdepTOFbeta;
TSpline3 *EdepTrackbeta;
TSpline3 *EdepTOFDbeta;
TSpline3 *Corr_L1;
TSpline3 *Corr_TOFU ;
TSpline3 *Corr_Track;
TSpline3 *Corr_TOFD; 


///  Variables retrieved from Root ntuples
struct Tuplevar {
float BDT_response;
float Beta;
float Beta_pre;
float BetaRICH;
float Dist5D;
float Dist5D_P;
float EdepECAL;
float EdepL1;
float EdepTOFD;
float EdepTOFU;
float EdepTrack;
float Ev_Num;
int   Cutmask;
float Latitude;
float LDiscriminant;
float MC_type;
float Momento_gen;
float R;
float Rcutoff;
float Rmin;
float R_pre;
float Trig_Num;
float Unbias;
};

//// Global Variables





float Massa_gen=0;
float Massa=0;
float Zona=0;
float IsPrescaled=0;
float X=0;
float Beta_gen=0;

float BetaRICH_new=0;
float Rcut[11]= {18,18,16,14,12,10,8,6,4,2,1};
int INDX=0;
int FRAC=0;
double geomagC[11]= {0,0.05,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.15};
int notpassed[3]= {155,139,11};
int passed[3]= {187,155,139};





TH1F* Tempi;



float encinprot     [nbinsr];
float encindeut     [nbinsr];
float RUsed=0;
int Norm[11]= {0};

//cuts
bool Likcut=false;
bool Distcut=false;
bool Herejcut=false;
bool Betastrongcut=false;



Tuplevar Tup;

//retrieve MC particle species
float ReturnMass_Gen()
{
   float Mass_gen=0;
   if ( ( ( (int) Tup.MC_type) &0xFF    ) >0)      Mass_gen = 0.938;
   if ( ( ( (int) Tup.MC_type) &0xFF00  ) >0)      Mass_gen = 1.875;
   if ( ( ( (int) Tup.MC_type) &0xFF0000) >0)      Mass_gen = 3.725;
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
      if ( ( ( ( (int) Tup.MC_type) >> (cursor+i) ) & 1 ) ==1) mc_type=i-2;
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


