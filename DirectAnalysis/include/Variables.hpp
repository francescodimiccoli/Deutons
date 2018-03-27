#ifndef VARIABLES_H
#define VARIABLES_H

#include "reweight.h"
#include <TMVA/Reader.h>
#include <TMVA/Tools.h>
#include "Globals.h"
#include "TF1.h"

using namespace std;

struct Variables{

    Reweighter reweighter;
    Reweighter reweighterHe;

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

    // Bit fields
    int        JMembPatt;
    int        PhysBPatt;
    int        CUTMASK;
    int        RICHmask_new;

    // Counts
    float      NAnticluster_float;
    int        NTofClusters;
    int        NTofClustersusati;
    float      NTofUsed=0;  // ? 
    int        NTRDclusters;
    int        NAnticluster;
    int        NTRDSegments;
    int        NTrackHits;                             
    int        clustertrack;
    int        clustertottrack;


    // Track
    float      R_pre;  // ?
    float      Rup;
    float      Rdown;
    float      R;
    float      R_L1;
    float      Chisquare;
    float      Chisquare_L1;
    float      Chisquare_y;
    float      Chisquare_L1_y;
    int        hitbits;
    int        FiducialVolume;	 	    

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
    int        Richtotused;
    float      Richtotused_float=0;
    float      RichPhEl;
    float      RICHprob;
    int        RICHPmts;
    float      RICHPmts_float;	
    float      RICHcollovertotal;
    int        RICHgetExpected;
    float      RICHgetExpected_float; 

   float RICHLipBetaConsistency =0;
   float RICHTOFBetaConsistency =0;	
   float RICHChargeConsistency  =0;	
   float tot_hyp_p_uncorr	=0;
   float Bad_ClusteringRICH     =0;
   float NSecondariesRICHrich   =0;	    

    //MC vars
    float      Momento_gen;
    float      Massa_gen;
    float      mcweight=0;
    float      GenX, GenY, GenZ;
    float      GenPX, GenPY, GenPZ;
    UInt_t     MCClusterGeantPids; 	
    float      Charge_gen;

    

    //Other
    float      EdepECAL;
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
    inline bool IsFromNaF     (){ return (((int)joinCutmask>>11)==512&&BetaRICH_new>0);}
    inline bool IsFromAgl     (){ return (((int)joinCutmask>>11)==0&&BetaRICH_new>0);}
    inline bool IsFromNaF_nosel     (){ return ((((int)joinCutmask>>11)&512)==512&&BetaRICH_new>0);}
    inline bool IsFromAgl_nosel     (){ return ((((int)joinCutmask>>11)&512)==0&&BetaRICH_new>0);}

};

float GetInverseRigidity (Variables * vars); 
float GetGenMomentum     (Variables * vars); 

float GetInverseEdepUToF (Variables * vars); 
float GetInverseEdepLToF (Variables * vars); 
float GetInverseEdepTrack(Variables * vars); 
float GetInverseEdepTRD  (Variables * vars); 

float GetBetaGen         (Variables * vars);
float GetInverseBetaTOF  (Variables * vars);
float GetInverseBetaRICH (Variables * vars);
float GetBetaTOF         (Variables * vars);
float GetBetaRICH        (Variables * vars);
float GetRecMassTOF	     (Variables * vars);
float GetRecMassRICH     (Variables * vars);
float GetNegRecMassTOF	 (Variables * vars);
float GetNegRecMassRICH  (Variables * vars);

float GetRigidity (Variables * vars);

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

float GetRICHBDT(Variables * vars);

int   GetPIDatL1 (Variables * vars);
int   GetPIDatL2 (Variables * vars);
int   GetPIDatL3 (Variables * vars);

float GetLoweredBetaTOF  (Variables * vars);
float GetLoweredBetaNaF  (Variables * vars);
float GetLoweredBetaAgl  (Variables * vars);
float GetRICHBDT(Variables* vars);
#endif
