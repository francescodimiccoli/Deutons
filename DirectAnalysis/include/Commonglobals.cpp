using namespace std;
int Ev_Num;
int Timebeg;


TSpline3 *Rig_p		;	
TF1 *	  beta_p	;	
TF1 *	  betaNaF_p	;	
TF1 *	  betaAgl_p	;	
TSpline3 *Rigmean_p	;	
TF1 *	  betamean_p	;	
TF1 *	  betaNaFmean_p	;	
TF1 *	  betaAglmean_p	;	
TSpline3 *etofu_p	;	
TSpline3 *etrack_p	;	
TSpline3 *etofd_p	;	
TSpline3 *EdepTOFbeta_p	;	
TSpline3 *EdepTrackbeta_p;	
TSpline3 *EdepTOFDbeta_p;	
                                
TSpline3 *Rig_d		;	
TF1 *	  beta_d	;	
TF1 *	  betaNaF_d	;	
TF1 *	  betaAgl_d	;	
TSpline3 *Rigmean_d	;	
TF1 *	  betamean_d	;	
TF1 *	  betaNaFmean_d	;	
TF1 *	  betaAglmean_d	;	
TSpline3 *etofu_d	;	
TSpline3 *etrack_d	;	
TSpline3 *etofd_d	;	
TSpline3 *EdepTOFbeta_d	;	
TSpline3 *EdepTrackbeta_d;	
TSpline3 *EdepTOFDbeta_d;	


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





