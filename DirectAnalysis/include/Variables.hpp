#ifndef VARIABLES_H
#define VARIABLES_H

#include "reweight.h"
#include <TMVA/Reader.h>
#include <TMVA/Tools.h>
#include "Globals.h"
#include "TF1.h"
#include "TSpline.h"

using namespace std;



struct Variables{

    Reweighter reweighter;
    Reweighter reweighterHe;

   static Double_t ChiXcut_X[];
   static Double_t ChiXcut_Y[];
   static Double_t ChiYcut_X[];
   static Double_t ChiYcut_Y[];

    TSpline3 * Chi2Xcut; 	
    TSpline3 * Chi2Ycut; 	



    // Event info
    int        Run;
    int        Event;
    int        NEvent;
    int        U_time;
    int        NTracks;
    float      ThetaS;
    float      PhiS;
    float      Livetime;
    float      Latitude;
    float      PrescaleFactor;

    // Cutoffs
    float      Rcutoff;
    float      IGRFRcutoff;

    //RTI
    int good_RTI;
    float       Rcutoff_RTI;
    float 	Rcutoff_IGRFRTI;
    int isinsaa;
    float	Livetime_RTI;

    // Bit fields
    int        JMembPatt;
    int        PhysBPatt;
    int        CUTMASK;
    int        RICHmask_new;
    int	       P_standard_sel;	

    // Counts
    float        NTofClusters;
    float        NTofClustersusati;
    float      NTofUsed=0;  // ? 
    float        NTRDclusters;
    float        NAnticluster;
    float        NTRDSegments;
    float        NTrackHits;                             
    float        clustertrack;
    float        clustertottrack;

    //Tracking Efficiency
    float theta_track;
    float phi_track;
    float entrypointcoo[3];	
    float beta_SA;
    float betapatt_SA;	
    float qUtof_SA;
    float qLtof_SA;			
    float qTrd_SA;	
    float EdepECAL;
    
    //L1 pick-up Efficiency
    float exthit_closest_q;     
    int   exthit_closest_status;
    float hitdistfromint; 

    // Track
    float      R_pre;  // ?
    float      Rup;
    float      Rdown;
    float      R;
    float      R_L1;
    float      R_noMS;
    float      Chisquare;
    float      Chisquare_L1;
    float      Chisquare_y;
    float      Chisquare_L1_y;
    int        hitbits;
    float        FiducialVolume;	 	    
    float      R_sec;		

    // Tracker Charge
    float      qL1;
    float      qL2;
    float      qL1InnerNoL2;	   
    float      qL1Status;
    float      qL2Status;
    float      qInner;

    // TOF
    float      Beta;
    float      BetaR;
    float      Beta_pre; // ?
    float      qUtof;
    float      qLtof;
    float      TOFchisq_s;
    float      TOFchisq_t;	
    float      NBadTOF;

    // TRD
    float       TRDEdepovPath;		
    float 	TRDLikP;
    float 	TRDLikD;
    float 	TRDLike;
    float	EdepTRD;
	
    // RICH 
    float      BetaRICH_new;
    float        Richtotused;
    float      RichPhEl;
    float      RICHprob;
    float        RICHPmts;
    float      RICHcollovertotal;
    float        RICHgetExpected;

   float RICHLipBetaConsistency =0;
   float RICHTOFBetaConsistency =0;	
   float RICHChargeConsistency  =0;	
   float tot_hyp_p_uncorr	=0;
   float Bad_ClusteringRICH     =0;
   float NSecondariesRICHrich   =0;	    
   float HitHValldir =0;
   float HitHVallrefl=0;
   float HitHVoutdir=0;
   float HitHVoutrefl=0;   
   float HVBranchCheck=0;


    //MC vars
    float      Momento_gen;
    float      Momento_gen_UTOF;
    float      Momento_gen_LTOF;
    float      Momento_gen_RICH;
    float      Massa_gen;
    float      mcweight=0;
    float      GenX, GenY, GenZ;
    float      GenPX, GenPY, GenPZ;
    UInt_t     MCClusterGeantPids; 	
    float      Charge_gen;

    

    //Other
    std::vector<float>  *       trtrack_edep=0;
    std::vector<float>  *       trtot_edep=0;
    std::vector<float>  *       Endep=0;

    float      joinCutmask=0;
    float      diffR=0;
    float	   TOF_Up_Down=0;
    float	   Layernonusati=0;

    //Discriminants
    float DistP=0;
    float DistD=0;
    float Likelihood=0;
    float BDTDiscr=0;

    //Summed E. deps.
    float EdepTOFU=0;
    float EdepTOFD=0;
    float EdepTrack=0;	
    float EdepL1=0;	

    //Checks on Variables
    float  beta_ncl;
    float  chisqcn; 
    float  chisqtn; 
    float  nTrTracks;
    float  sumclsn; 

    TMVA::Reader *readerTOF;
    TMVA::Reader *readerNaF;
    TMVA::Reader *readerAgl;

    Variables();

    void ResetVariables();
    void ReadBranches(TTree * tree);
    void RegisterBranches(TTree * tree);
    void RegisterTemplatesBranches(TTree * tree);
    void Update();
    void PrintCurrentState();
    void BDTreader();
    void Eval_Discriminants();
    inline bool IsFromNaF     (){ return (((int)joinCutmask>>11)==1024&&BetaRICH_new>0);}
    inline bool IsFromAgl     (){ return (((int)joinCutmask>>11)==0&&BetaRICH_new>0);}
    inline bool IsFromNaF_nosel     (){ return ((((int)joinCutmask>>11)&1024/*1921*/)==1024&&BetaRICH_new>0);}
    inline bool IsFromAgl_nosel     (){ return ((((int)joinCutmask>>11)&1024)==0&&BetaRICH_new>0);}

};

float GetInverseRigidity (Variables * vars); 
float GetGenMomentum     (Variables * vars); 

float GetInverseEdepUToF (Variables * vars); 
float GetInverseEdepLToF (Variables * vars); 
float GetInverseEdepTrack(Variables * vars); 
float GetInverseEdepTRD  (Variables * vars); 

float GetBetaGen         (Variables * vars);
float GetBetaSlowTOF     (Variables * vars);
float GetBetaSlowRICH    (Variables * vars);
float GetRigSlow         (Variables * vars);


float GetInverseBetaTOF  (Variables * vars);
float GetInverseBetaRICH (Variables * vars);
float GetBetaTOF         (Variables * vars);
float GetBetaRICH        (Variables * vars);
float GetRecMassTOF	     (Variables * vars);
float GetRecMassRICH     (Variables * vars);
float GetNegRecMassTOF	 (Variables * vars);
float GetNegRecMassRICH  (Variables * vars);

float GetRigidity (Variables * vars);
float GetRigiditySecondTrack (Variables * vars);


float GetUtofQ	(Variables * vars);
float GetLtofQ	(Variables * vars);
float GetInnerQ	(Variables * vars);
float GetL1Q    (Variables * vars);
float GetL2Q    (Variables * vars);
float GetInnerL1NoL2Q    (Variables * vars);

float GetTRDEdepovPath    (Variables * vars);
float GetTRDePLikRatio    (Variables * vars); 
float GetTRDDPLikRatio    (Variables * vars); 

float GetTOFSpatialChi (Variables * vars); 
float GetTOFTimeChi (Variables * vars); 
float GetRupdown (Variables * vars);
float GetChisquareX(Variables * vars);
float GetChisquareY(Variables * vars);

float GetRICHBDT(Variables * vars);

int   GetPIDatL1 (Variables * vars);
int   GetPIDatL2 (Variables * vars);
int   GetPIDatL3 (Variables * vars);

float GetLoweredBetaTOF  (Variables * vars);
float GetLoweredBetaNaF  (Variables * vars);
float GetLoweredBetaAgl  (Variables * vars);
float GetRICHBDT(Variables* vars);

float GetEdepECAL(Variables * vars);
float GetMomentumProxy(Variables *vars); 

float GetNToFClusters(Variables * vars);; 
float GetUtofQ	(Variables * vars) ;
float GetLtofQ	(Variables * vars) ;
float GetTofChisqcn (Variables * vars) ;
float GetTofChisqtn  (Variables * vars); 
float GetNTracks (Variables * vars) ;
float GetTofOnTime (Variables * vars);  


#endif
