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

// Rigidity cutoff safety factor

float SF = 1.2;

TF1 *protons = new TF1("f1","pow((pow(x,2)/pow(0.938,2)/(1 + pow(x,2)/pow(0.938,2))),0.5)",0.1,100);
TF1 *deutons = new TF1("f1","pow((pow(x,2)/pow(1.875,2)/(1 + pow(x,2)/pow(1.875,2))),0.5)",0.1,100);

TF1 *RBeta = new TF1("R(Beta)_for_protons","pow((pow(0.938,2)*(pow(x,2)/(1-pow(x,2)))),0.5)",0.1,0.999999999999999999999999999999);

enum mode {BUILDALL, BUILDSEPD, READ};

enum {Betaedges,Redges};
////////////// DEFINIZIONE SPLINES //////////////////
TSpline3 *Rig;
TSpline3 *beta;
TF1 *betaNaF;
TF1 *betaAgl;
TSpline3 *eL1;
TSpline3 *etofu;
TSpline3 *etofd;
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
   float R_pre;
   float Trig_Num;
   float PhysBPatt;
   float mcweight;
   float U_time;
   float Livetime;	
   float qL1;
   float qInner;
   float qUtof;
   float qLtof;
   float Rmin; 	
   float R_L1;
};

//// Global Variables




int    ActualTime=0;
float Massa_gen=0;

float Rcut[11]= {18,18,16,14,12,10,8,6,4,2,1};
double geomagC[11]= {0,0.05,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.15}; // Only for drawing, but in various files


float RUsed=0;

//cuts
bool Likcut=false;
bool Distcut=false;
bool Herejcut=false;
bool Betastrongcut=false;
bool IsHeL1=false;
bool IsPfromHeL1=false;
bool ProtonsMassWindow=false;
bool ProtonsMassThres=false;

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


void FillBinMGen (TH1* h, int bin,float weight)
{
   int mass = ReturnMCGenType();
   ((TH2*)h)->Fill (bin, mass,weight);
   return;
}

void FillBinMGen (TH3* h, int bin, int S)
{
   int mass = ReturnMCGenType();
   ((TH3*)h)->Fill (bin, mass, S);
   return;
}



class FileSaver {
   public:
      FileSaver(bool Isfinal=false) : filename("") {fArr=new TObjArray(); IsFinal = Isfinal;}
      FileSaver (string fname,bool Isfinal=false) : filename (fname) {fArr=new TObjArray(); IsFinal = Isfinal;}
      void writeObjsInFolder(string folder, bool recreate = false);
      void writeObjs();
      void setName(string fname) {filename=fname;}
      void Add(TObject* obj) {fArr->Add(obj);}
      string setName() {return filename;}
   private:
      string filename;
      TObjArray* fArr;
      bool IsFinal;	
};

void FileSaver::writeObjsInFolder(string folder, bool recreate)
{
	cout<<"*** Updating "<<filename.c_str()<<" file in "<< folder << endl;
	TFile * fileFinalPlots;

	if(recreate) fileFinalPlots=TFile::Open(filename.c_str(), "RECREATE");
	else fileFinalPlots=TFile::Open(filename.c_str(), "UPDATE");

	if (!fileFinalPlots->GetDirectory(folder.c_str()))
		fileFinalPlots->mkdir(folder.c_str());
	fileFinalPlots->cd   (folder.c_str());
	for (int i = 0; i <= fArr->GetLast(); i++){
		fArr->At(i)->Write();
	}
	fileFinalPlots->Flush();
	fileFinalPlots->Write();
	if(IsFinal) fileFinalPlots->Close(); 
	fArr->Clear();
	return;
}


bool Isfinal=true;
FileSaver finalPlots(Isfinal);
FileSaver finalHistos;

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
TriggPatt trgpatt;


void HistInfo(TH1* histo) {
   cout << ">>> Class " << histo->ClassName () << " >> Info :  ";
  histo->Print();
   return;
   }

TH1F* TH1DtoTH1F(TH1D* hd) {
   TH1F*  hf=(TH1F*)hd->Clone();
   hf->SetName(hd->GetTitle());
    return hf;
}



TH1F* ProjectionXtoTH1F(TH2F* h2, string title, int binmin, int binmax) {
   TH1D* hd=h2->ProjectionX(title.data(),binmin, binmax);
   TH1F* hf=TH1DtoTH1F(hd);
   return hf;
   }

TH1F* ProjectionYtoTH1F(TH2F* h2, string title, int binmin, int binmax) {
   //HistInfo(h2);
   TH1D* hd=h2->ProjectionY(title.data(), binmin, binmax);
   //HistInfo(hd);
   TH1F* hf=TH1DtoTH1F(hd);
   return hf;
   }


TH1F * Extract_Bin(TH1 * Histo, int bin,int third_dim,bool reverse=false)
{
   TH1F* Slice;
   if(!reverse){
   Slice = new TH1F("","",Histo->GetNbinsX(),Histo->GetXaxis()->GetBinLowEdge(1),Histo->GetXaxis()->GetBinLowEdge(Histo->GetNbinsX()+1));
   for(int i = 0; i< Histo->GetNbinsX(); i++)
      Slice->SetBinContent(i+1,Histo->GetBinContent(i+1,bin+1,third_dim+1));
   }
   else{
   Slice = new TH1F("","",Histo->GetNbinsY(),Histo->GetYaxis()->GetBinLowEdge(1),Histo->GetYaxis()->GetBinLowEdge(Histo->GetNbinsY()+1));
   for(int i = 0; i< Histo->GetNbinsY(); i++)
        Slice->SetBinContent(i+1,Histo->GetBinContent(bin+1,i+1,third_dim+1));
   }

   return Slice;

}


TH1 * SetErrors ( TH1 * Eff){
        TH1 * Errors = (TH1 *)Eff->Clone();
        for(int iR=0;iR<Eff->GetNbinsX();iR++)
                for(int mc_types=0;mc_types<Eff->GetNbinsY();mc_types++)
                        for(int S=0;S<Eff->GetNbinsZ();S++)
                         Errors->SetBinContent(iR+1,mc_types+1,S+1,Eff->GetBinError(iR+1,mc_types+1,S+1)/Eff->GetBinContent(iR+1,mc_types+1,S+1));
        return Errors;

}

void EvalError( TH1 * Eff, TH1 * stat, TH1 * syst){
        cout<<Eff<<" "<<stat<<" "<<syst<<endl;
        for(int iR=0;iR<Eff->GetNbinsX();iR++)
                for(int mc_types=0;mc_types<Eff->GetNbinsY();mc_types++)
                        for(int S=0;S<Eff->GetNbinsZ();S++){
                                float Stat=0;
                                float Syst=0;
                                if(stat&&syst){
                                 Stat = stat->GetBinContent(iR+1,mc_types+1,S+1);
                                 Syst = syst->GetBinContent(iR+1,mc_types+1,S+1);
                                 Eff->SetBinError(iR+1,mc_types+1,S+1,pow(pow(Stat,2)+pow(Syst,2),0.5)*Eff->GetBinContent(iR+1,mc_types+1,S+1));
                                }
                        }
                return;
}




void Disable_MCreweighting(){
Tup.mcweight =1;
return;
}
