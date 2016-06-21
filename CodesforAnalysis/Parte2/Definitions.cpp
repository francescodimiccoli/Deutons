using namespace std;

#include "binning.h"

//// Input-Output Variables

string mese;
TFile* inputHistoFile;
TFile* fileFinalPlots;


extern const int nbinsr=43;
extern const int nbinsToF=18;
extern const int nbinsNaF=18;
extern const int nbinsAgl=18;



TF1 *protons = new TF1("f1","pow((pow(x,2)/pow(0.938,2)/(1 + pow(x,2)/pow(0.938,2))),0.5)",0.1,100);
TF1 *deutons = new TF1("f1","pow((pow(x,2)/pow(1.875,2)/(1 + pow(x,2)/pow(1.875,2))),0.5)",0.1,100);


enum mode {BUILDALL, BUILDSEPD, READ};

////////////// DEFINIZIONE SPLINES //////////////////
TSpline3 *Rig;
TSpline3 *beta;
TF1 *betaNaF;
TF1 *betaAgl;
TSpline3 *eL1;
TSpline3 *etofu;
TSpline3 *etrack;
TSpline3 *EdepL1beta;
TSpline3 *EdepTOFbeta;
TSpline3 *EdepTrackbeta;



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
   float Cutmask;
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

float Rcut[11]= {18,18,16,14,12,10,8,6,4,2,1};
double geomagC[11]= {0,0.05,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.15}; // Only for drawing, but in various files


TH1F* Tempi;
float RUsed=0;

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



class FileSaver {
   public:
      FileSaver() : filename("") {fArr=new TObjArray();}
      FileSaver (string fname) : filename (fname) {fArr=new TObjArray();}
      void writeObjsInFolder(string folder);
      void setName(string fname) {filename=fname;}
      void Add(TObject* obj) {fArr->Add(obj);}
      string setName() {return filename;}
   private:
      string filename;
      TObjArray* fArr;
      TFile* fileFinalPlots;
};

void FileSaver::writeObjsInFolder(string folder)
{
   cout<<"*** Updating Results file in "<< folder << endl;
   fileFinalPlots-> Open(filename.data(), "UPDATE");
   if (!fileFinalPlots->GetDirectory(folder.data()))
      fileFinalPlots->mkdir(folder.data());
   fileFinalPlots->cd   (folder.data());
   for (int i = 0; i <= fArr->GetLast(); i++)
      fArr->At(i)->Write();
   fileFinalPlots->Write();
   fileFinalPlots->Flush();
   fileFinalPlots->Close();
   fArr->Clear();
   return;
}

FileSaver finalPlots;

Particle proton(0.9382720813, 1, 1);  // proton mass 938 MeV
Particle deuton(1.8756129   , 1, 2);  // deuterium mass 1876 MeV, Z=1, A=2

Binning ToFDB(deuton);
Binning ToFPB(proton);
Binning NaFDB(deuton);
Binning NaFPB(proton);
Binning AglDB(deuton);
Binning AglPB(proton);

Binning DRB(deuton);
Binning PRB(proton);

Cutmask cmask;
