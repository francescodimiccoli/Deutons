#ifndef VARIABLES_HPP
#define VARIABLES_HPP

struct Variables 
{
    int     Run;
    int     Event;
    int     NEvent;
    int     U_time;
    float   Latitude;
    float   Zenith;
    float   Rcutoff35;
    float   Rcutoff40;
    float   Livetime;
    float   ThetaS;
    float   PhiS;

    int     PhysBPatt;
    int     CUTMASK;
    int     RICHmask;

    float   R;
    float   Rup;
    float   Rdown;
    float   R_L1;
    float   Chisquare;
    float   Chisquare_L1;
    int     hitbits;

    float   BetaOld;
    float   BetaRaw;
    float   BetaHRS;
    float   BetaHR;
    float   BetaRICH;

    std::vector<float> trtrack_edep;
    std::vector<float> trtot_edep;

    std::vector<float>   TOFEndep;
    std::vector<float>   TOFEndepR;

    float   EnergyECAL;
    int     NAnticluster;
    int     NTofClusters;
    int     NTofClustersusati;
    int     NTrackHits;
    int     NTRDclusters;
    int     NTRDSegments;
    int     clustertrack;
    int     clustertottrack;
    int     Richtotused;
    float   RichPhEl;
    float   EdepTRD;
    float   qL1;
    float   qUtof;
    float   qLtof;
    float   qInner;

    // Monte-Carlo Data
	float GenMomentum;
	float GenMass;
    float GenX, GenY, GenZ;
    float GenPX, GenPY, GenPZ;

    //Constructor
    Variables():
        trtrack_edep(9,0),
          trtot_edep(9,0),
            TOFEndep(4,0),
           TOFEndepR(4,0) { }
    
    void RegisterBranches(TTree * tree);
    void ResetVariables();

};

void Variables::RegisterBranches(TTree * tree) 
{
    tree->Branch("Run",&Run);
    tree->Branch("Event",&Event);
    tree->Branch("NEvent",&NEvent);
    tree->Branch("U_time",&U_time);
    tree->Branch("Zenith",&Zenith);
    tree->Branch("Latitude",&Latitude);
    tree->Branch("Rcutoff35",&Rcutoff35);
    tree->Branch("Rcutoff40",&Rcutoff40);
    tree->Branch("Livetime",&Livetime);
    tree->Branch("ThetaS",&ThetaS);
    tree->Branch("PhiS",&PhiS);

    tree->Branch("PhysBPatt",&PhysBPatt);
    tree->Branch("CUTMASK",&CUTMASK);
    tree->Branch("RICHmask",&RICHmask);

    tree->Branch("R;",&R);
    tree->Branch("Rup",&Rup);
    tree->Branch("Rdown",&Rdown);
    tree->Branch("R_L1",&R_L1);
    tree->Branch("Chisquare",&Chisquare);
    tree->Branch("Chisquare_L1",&Chisquare_L1);
    tree->Branch("hitbits",&hitbits);

    tree->Branch("BetaOld"  , &BetaOld);
    tree->Branch("BetaRaw"  , &BetaRaw);
    tree->Branch("BetaHRS"  , &BetaHRS);
    tree->Branch("BetaHR"   , &BetaHR);
    tree->Branch("BetaRICH" , &BetaRICH);

    tree->Branch("trtrack_edep",&trtrack_edep);
    tree->Branch("trtot_edep",&trtot_edep);

    tree->Branch("TOFEndep",&TOFEndep);
    tree->Branch("TOFEndepR",&TOFEndepR);

    tree->Branch("EnergyECAL",&EnergyECAL);
    tree->Branch("NAnticluster",&NAnticluster);
    tree->Branch("NTofClusters",&NTofClusters);
    tree->Branch("NTofClustersusati",&NTofClustersusati);
    tree->Branch("NTrackHits",&NTrackHits);
    tree->Branch("NTRDclusters",&NTRDclusters);
    tree->Branch("NTRDSegments",&NTRDSegments);
    tree->Branch("clustertrack",&clustertrack);
    tree->Branch("clustertottrack",&clustertottrack);
    tree->Branch("Richtotused",&Richtotused);
    tree->Branch("RichPhEl",&RichPhEl);
    tree->Branch("EdepTRD",&EdepTRD);
    tree->Branch("qL1",&qL1);
    tree->Branch("qUtof",&qUtof);
    tree->Branch("qLtof",&qLtof);
    tree->Branch("qInner",&qInner);

    // Monte-Carlo Data
	tree->Branch("GenMomentum",&GenMomentum);
	tree->Branch("GenMass",&GenMass);
    tree->Branch("GenX",&GenX); 
    tree->Branch("GenY",&GenY); 
    tree->Branch("GenZ",&GenZ);
    tree->Branch("GenPX",&GenPX); 
    tree->Branch("GenPY",&GenPY); 
    tree->Branch("GenPZ",&GenPZ);
}

void Variables::ResetVariables() {
    Run      = 0;
    Event    = 0;
    NEvent   = 0;
    U_time   = 0;
    Latitude = 0;
    Zenith   = 0;
    Rcutoff35= 0;
    Rcutoff40= 0;
    Livetime = 0;
    ThetaS   = 0;
    PhiS     = 0;
    
    PhysBPatt = 0;
    CUTMASK   = 0;
    RICHmask  = 0;
    
    R            = 0;
    Rup          = 0;
    Rdown        = 0;
    R_L1         = 0;
    Chisquare    = 0;
    Chisquare_L1 = 0;
    hitbits      = 0;
      
    BetaOld = 0;
    BetaRaw = 0;
    BetaHRS = 0;
    BetaHR  = 0;
    BetaRICH= -1;

    EnergyECAL   = 0;
    NAnticluster = 0;
    NTofClusters = 0;
    NTofClustersusati = 0;
    NTrackHits   = 0;
    NTRDclusters = 0;
    NTRDSegments = 0;
    clustertrack = 0;
    clustertottrack = 0;
    Richtotused  = 0;
    RichPhEl     = 0;
    EdepTRD      = 0;
    qL1          = 0;
    qUtof        = 0;
    qLtof        = 0;
    qInner       = 0;

    std::fill(trtrack_edep.begin(), trtrack_edep.end(), 0);
    std::fill(trtot_edep.begin(),   trtot_edep.end(),   0);
    std::fill(TOFEndep.begin(),     TOFEndep.end(),     0);
    std::fill(TOFEndepR.begin(),    TOFEndepR.end(),    0);

    // Monte-Carlo Data
	GenMomentum = -9999;
	GenMass = 0;
    GenX  = 0; GenY  = 0; GenZ  = 0;
    GenPX = 0; GenPY = 0; GenPZ = 0;
}






#endif 
