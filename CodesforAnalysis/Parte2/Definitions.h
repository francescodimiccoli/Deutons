
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
double geomag[12]={0,0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.3};
int avanzamento=0;
int INDX=0;
int FRAC=0;
double geomagC[11]={0,0.05,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.15};
int notpassed[3]={155,139,11};
int passed[3]={187,155,139};
float B=0.4;
float B1=0;
float B2=0;
float E=0.1;
int binnum=1;
float distcut=4;
float ddiscrcut=0.15;
bool qualcut=true;
float a=(log(0.9)-log(0.1))/18;
float E2=exp(log(0.1)+1.5*a);
float Betabins[18]={0.4};
float Betacent[18]={0};
float Ekincent[18]={0};
float BetabinsR_P[18]={0};
float BetabinsR_D[18]={0};
float BetabinsNaF[18]={0.4};
float BetacentNaF[18]={0};
float EkincentNaF[18]={0};
float BetabinsNaFR_P[18]={0};
float BetabinsNaFR_D[18]={0};
float BetabinsAgl[18]={0.4};
float BetabinsAglR_P[18]={0};
float BetabinsAglR_D[18]={0};
float BetaP[18]={0};
float BetaD[18]={0};
float BetaNaFP[18]={0};
float BetaNaFD[18]={0};
float BetaAglP[18]={0};
float BetaAglD[18]={0};
float Unbias=0;
float BetacentAgl[18]={0};
float EkincentAgl[18]={0};
float bin[44];
double R_cent[43];
float encinprot[43];
float encindeut[43];
float deltaencinprot[43];
float deltaencindeut[43];
float deltaencinTOF[18];
float deltaencinNaF[18];
float deltaencinAgl[18];
float Var=0;
float Var2=0;
float Var3=0;
bool Likcut=false;
bool Distcut=false;
bool Herejcut=false;
TH1F * Esposizione[10];
TH1F * Tempi;
TH2F * esposizionegeo;
TH2F * esposizionepgeo;
TH2F * esposizionepgeoNaF;
TH2F * esposizionepgeoAgl;
TH2F * esposizionedgeo;
TH2F * esposizionedgeoNaF;
TH2F * esposizionedgeoAgl;
