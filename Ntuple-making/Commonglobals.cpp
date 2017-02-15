using namespace std;
int Ev_Num;
int scelta=0;
void aggiungiantupla (TTree *albero,int i,TNtuple *ntupla);
void GrandezzequalRICH (TTree *albero,int i,TNtuple *ntupla);
void Grandezzequal (TTree *albero,int i,TNtuple *ntupla);
void Grandezzesepd (TTree *albero,int i,TNtuple *ntupla);
void Grandezzesep (TTree *albero,int i,TNtuple *ntupla);
void Trigg (TTree *albero,int i,TNtuple *ntupla);

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


TF1 *protons = new TF1("f1","pow((pow(x,2)/pow(0.938,2)/(1 + pow(x,2)/pow(0.938,2))),0.5)",0.1,100);
TF1 *deutons = new TF1("f1","pow((pow(x,2)/pow(1.875,2)/(1 + pow(x,2)/pow(1.875,2))),0.5)",0.1,100);


TFile *_file1 = TFile::Open("/storage/gpfs_ams/ams/users/fdimicco/Deutons/Ntuple-making/QualityVariables.root");
TFile *_file3 = TFile::Open("/storage/gpfs_ams/ams/users/fdimicco/Deutons/Ntuple-making/QualityVariables_NaF.root");
TFile *_file3b = TFile::Open("/storage/gpfs_ams/ams/users/fdimicco/Deutons/Ntuple-making/QualityVariables_Agl.root");


TSpline3 *Bkgnd[9];
TSpline3 *Signal[9];
TSpline3 *BkgndNaF[9];
TSpline3 *SignalNaF[9];
TSpline3 *BkgndAgl[9];
TSpline3 *SignalAgl[9];





