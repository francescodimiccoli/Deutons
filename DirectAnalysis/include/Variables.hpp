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
    float exthit_int[3];          
    float exthit_closest_coo[3];
    float exthit_closest_q;     
    int   exthit_closest_status;
    float exthit_largest_coo[3];
    float exthit_largest_q;     
    int   exthit_largest_status;


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

    int   nparticle;                ///< Number of crossed PMTs (= number of traversing particles)
    int   tot_hit_uncorr;           ///< Uncorrected total number of hits
    int   tot_pmt_uncorr;           ///< Uncorrected total number of PMTs
    float tot_p_uncorr;             ///< Uncorrected total number of p.e. discarding those on PMTs crossed by charged particles.
    float max_p_uncorr[2];          ///< Uncorrected maximum number of p.e. in a PMT (all, excluding crossed PMTs)

    int   pmt_nhit_uncorr[5];       ///< First 5 PMTs by number of p.e., uncorrected number of hits
    float pmt_np_uncorr[5];         ///< First 5 PMTs by number of p.e., uncorrected number of p.e.
    float pmt_dist[5];              ///< First 5 PMTs by number of p.e., CoG distance from reference track [cm]
    int   pmt_pmt[5];               ///< First 5 PMTs by number of p.e., PMT number

    // hits

    int          hit_hit;                  ///< Total number of hits in the beta-hit container (too far beta hyp. are dropped)
    //! First 30 beta-hit status
    /** Status bits:
      - 1: Hit used in the ring number 1
      - 2: Hit used in the ring number 2
      - 3: Hit used in the ring number 3
      - ...
      - 10: Hit used in the ring number 10 (no more than 10 created)
      - 29: Channel taggeg as good in calibration
      - 30: Gain mode chosen for the hit 0=x1(low) 1=x5(high)
      - 31: Hit belongs to a PMT apparently crossed by a charged particle
      */
    unsigned int hit_stat[30];
    int          hit_used[30];             ///< First 30 beta-hit used-word (0: used in ring as direct, 1: used in ring as reflected, otherwise not used in the ring)
    int          hit_chan[30];             ///< First 30 beta-hit channel number (16*PMT+pixel)
    float        hit_np_uncorr[30];        ///< First 30 beta-hit uncorrected number of p.e.
    float        hit_beta[30][2];          ///< First 30 beta-hit raw beta (0: direct, 1: reflected, -1 for bad evaluation)
    int          tot_hyp_hit_uncorr[2][3]; ///< Total number of hit compatible with beta=1 hypothesys (all, out of ring|direct, reflected, reflected)
    float        tot_hyp_p_uncorr[2];      ///< Total number of p.e. for beta=1 hypothesys (all, out of ring)

    // ring

    int   selection;                ///< Javier selection (see Tools::RichQC)
    int   is_naf;                   ///< Used radiatior is NaF
    int   status;                   ///< Status word
    int   correction_status;        ///< Charge corrections fail flag (-1/0/1 : Not/Done/Failed)
    int   nhit;                     ///< Number of hits in the ring
    int   nhit_uncorr;              ///< Uncorrected number of hits in the ring
    int   nhit_refl;                ///< Number of hits which are consistent with reflected photons
    int   npmt;                     ///< Number of PMTs in the ring
    int   npmt_uncorr;              ///< Uncorrected number of PMTs in the ring
    float np;                       ///< Number of p.e. collected in the ring
    float np_uncorr;                ///< Uncorrected number of p.e. collected in the ring
    float np_w[10];                 ///< Photoelectrons associated to the ring for different windows sizes (1,2,3...10)
    float npmt_exp;                 ///< Expected number of PMTs
    float np_exp;                   ///< Expected number of p.e. for a Z=1 ring with reconstruction input pars of the current event
    float np_exp_uncorr;            ///< Uncorrected expected number of p.e.
    float np_exp_elec;              ///< Expected number of p.e. for a Z=1 ring with reconstruction input pars of the current event and beta = 1
    float prob;                     ///< Kolmogorov test to the distribution of charge along the ring
    float width;                    ///< Width of the distribution of charge around the ring over the expected one
    float udist;                    ///< (1/^2) for unused hits which do not belong to PMTs crossed by a charged particle (?)
    int   tile_id;                  ///< Tile id for the tile crossed by the particle
    float rad_coo[2];               ///< Track position interpolated to the radiator entrance (as used in the reconstruction) [cm]
    float rad_theta;                ///< Track theta interpolated to the radiator entrance (as used in the reconstruction) [rad]
    float rad_phi;                  ///< Track phi interpolated to the radiator entrance (as used in the reconstruction) [rad]
    float distance_tile_border;     ///< Distance of the track impact point in the radiator to the border of the radiator tile
    float beta;                     ///< Beta (tileCorrection)
    float beta_res;                 ///< Expected resolution
    float beta_rms;                 ///< Expected resolution RMS
    float beta_raw;                 ///< Raw beta
    float beta_refit;               ///< Refit beta
    float beta_corrected;           ///< Beta corrected for impact point and direction (best beta estimator)
    float beta_unifcorr;            ///< Beta corrected for impact point and direction (with previous correction, only for comparison)
    float q;                        ///< Charge
    float q_consistency;            ///< Statistical test to check if the hit by hit charge is consistent PMT-by-PMT
    float q_res;                    ///< Expected charge resolution
    float q_rms;                    ///< Expected charge resolution RMS
    int   nclus;                    ///< Number of beta clusters
    int   clus_size[10];            ///< Size of first 10 clusters (ordered by size)
    float clus_mean[10];            ///< Average beta of first 10 clusters (ordered by size)
    float clus_rms[10];             ///< RMS of first 10 clusters (ordered by size)
    float lip_beta;                 ///< LIP beta
    float lip_q;                    ///< LIP charge

    // veto

    float veto_np_exp_elec;         ///< Estimate number of expected p.e. in the case of electron (beta = 1) and no ring present
    float veto_np_exp_prot;         ///< Estimate number of expected p.e. in the case of proton (beta from rigidity) and no ring present
    float veto_beta_trk_prot;       ///< Beta estimated from rigidity



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
    void RegisterTemplatesBranches_RICH(TTree * tree);
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
