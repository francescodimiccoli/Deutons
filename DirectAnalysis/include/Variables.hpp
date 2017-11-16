#ifndef VARIABLES_H
#define VARIABLES_H

#include "reweight.h"
#include <TMVA/Reader.h>
#include <TMVA/Tools.h>

using namespace std;

struct Variables{

    Reweighter reweighter;
    Reweighter reweighterHe;

    int     Run;
    int     Event;
    int     NEvent;
    int 	   U_time;
    int 	NTracks;
    float 	   Latitude;
    float 	   Rcutoff;
    float   ThetaS;
    float   PhiS;
    float 	   IGRFRcutoff;
    float      Livetime;
    int        JMembPatt;
    int        PhysBPatt;
    float      R_pre;
    float      Beta_pre;
    int        CUTMASK;
    std::vector<float>  *       trtrack_edep=0;
    std::vector<float>  *       trtot_edep=0;
    std::vector<float>  *       Endep=0;
    float      BetaRICH_new;
    int        RICHmask_new;
    float      EdepECAL;
    int        NAnticluster;
    float	   NAnticluster_float;
    int        NTofClusters;
    int        NTofClustersusati;
    float	   NTofUsed=0;
    float      Rup;
    float      Rdown;
    float      R;
    float      Chisquare;
    float   Chisquare_L1;
    float   Chisquare_y;
    float   Chisquare_L1_y;
    int     hitbits;

    int     NTRDclusters;
    int     NTRDSegments;
    int     clustertrack;
    int     clustertottrack;


    float      Beta;
    float      BetaR;
    int        NTrackHits;                             
    int        Richtotused;
    float      Richtotused_float=0;
    float      RichPhEl;
    float      R_L1;
    float	   qL1;
    float   qL1Status;
    float	   qInner;
    float	   qUtof;
    float	   qLtof;
    float 	   joinCutmask=0;
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
    float EdepTRD=0;
    float EdepL1=0;	

    //MC vars
    float	   Momento_gen;
    float	   Massa_gen;		
    float	   mcweight=0;
    float GenX, GenY, GenZ;
    float GenPX, GenPY, GenPZ;
    UInt_t 	   MCClusterGeantPids; 	
    float	   Charge_gen;

    //RICH Variables
    float RICHprob;
    int RICHPmts;
    float RICHcollovertotal;
    int RICHgetExpected;

    TMVA::Reader *readerTOF;
    TMVA::Reader *readerNaF;
    TMVA::Reader *readerAgl;

    Variables();

    void ResetVariables();
    void ReadBranches(TTree * tree);
    void Update();
    void PrintCurrentState();
    void BDTreader();
    void Eval_Discriminants();
    inline bool IsFromNaF     (){ return (((int)joinCutmask>>11)==512&&BetaRICH_new>0);}
    inline bool IsFromAgl     (){ return (((int)joinCutmask>>11)==0&&BetaRICH_new>0);}
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
int   GetPIDatL1 (Variables * vars);
int   GetPIDatL2 (Variables * vars);
int   GetPIDatL3 (Variables * vars);


#endif
