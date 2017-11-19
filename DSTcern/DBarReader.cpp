// g++ -c DBarReader.cpp -std=c++14 -I/afs/cern.ch/user/k/kostams/data_reduction/include/ -I`root-config --incdir`


#include "TTree.h"
#include "Variables.hpp"
#include "DBarReader.h"

UInt_t DBarReader::getPackedLayers_1to4() {
    UInt_t id1 = ntpMCHeader.hitGeantID[0];
    UInt_t id2 = ntpMCHeader.hitGeantID[1];
    UInt_t id3 = ntpMCHeader.hitGeantID[2];
    UInt_t id4 = ntpMCHeader.hitGeantID[3];

    id1 = id1 < 127 ? id1 : 0;
    id2 = id2 < 127 ? id2 : 0;
    id3 = id3 < 127 ? id3 : 0;
    id4 = id4 < 127 ? id4 : 0;

    return id1 + (id2 << 8) + (id3 << 16) + (id4 << 24);
}

int countBits(int n) {
    int counter = 0;
    while(n) {
        counter += n % 2;
        n >>= 1;
    }
    return counter;
}

bool DBarReader::minTOF(){
    return ( ntpTof.flagp[0] == 0 ) &&
           ( ntpTof.flagp[1] == 0 ) &&
           ( ntpTof.flagp[2] == 0 ) &&
           ( ntpTof.flagp[3] == 0 ) ;
}

bool DBarReader::goldenTOF(){
    if ( (ntpTof.flag & 0x2) == 0 ) return false;
    if (  ntpTof.chisqcn > 10     ) return false;
    if (  ntpTof.chisqtn > 10     ) return false;
}

int DBarReader::RICHmaskConverter(){
    if(ntpRich.selection<0) return -ntpRich.selection;
    else if(ntpRich.selection==1) return 512;
    else return 0;
}


void DBarReader::FillVariables(int NEvent, Variables * vars){

    vars->ResetVariables();

    ////////////////////// EVENT INFORMATION ///////////////////////////////////////
    vars->Run       = ntpHeader.run;
    vars->Event     = ntpHeader.event;
    vars->NEvent    = NEvent;
    vars->U_time    = ntpHeader.utc_time; 
    vars->Livetime  = ntpHeader.livetime;
    vars->Latitude = ntpHeader.thetam;
    vars->ThetaS    = ntpHeader.thetas;
    vars->PhiS      = ntpHeader.phis;
    vars->Zenith    = 0; // not in dbar events 
    vars->Rcutoff35 = 0; // not in dbar event
    vars->Rcutoff40 = 0; // not in dbar event

    // Stroemer cutoff is in the tracker data
    vars->StoermerRcutoff = ntpTracker.stoermer_cutoff[0];

    vars->NAnticluster      = ntpHeader.nanti;
    vars->NTofClusters      = ntpHeader.ntof; // NOTE: That sees to be an error
    vars->NTofClustersusati = 0;
    vars->NTrackHits        = ntpHeader.ntrrechit;
    vars->NTracks           = ntpHeader.ntrtrack;
    vars->NTRDclusters      = ntpHeader.ntrdcluster;
    vars->NTRDSegments      = ntpHeader.ntrdsegment;

    ////////////////////// MONTE CARLO INFO ////////////////////////////////////////
    if(isMC) {
        vars->GenCharge   = ntpMCHeader.charge;
        vars->GenMass     = ntpMCHeader.mass;
        vars->GenMomentum = ntpMCHeader.momentum[0];
        vars->GenX = ntpMCHeader.coo[0][0];  vars->GenPX = ntpMCHeader.dir[0];;
        vars->GenY = ntpMCHeader.coo[0][1];  vars->GenPY = ntpMCHeader.dir[1];;
        vars->GenZ = ntpMCHeader.coo[0][2];  vars->GenPZ = ntpMCHeader.dir[2];;
        vars->MCClusterGeantPids = getPackedLayers_1to4();
    }    
   
    /////////////////////////////////// UNBIAS ////////////////////////////////////////
    vars->PhysBPatt = ntpHeader.sublvl1;
    vars->JMembPatt = ntpHeader.trigpatt;

    bool goodChi2 = (ntpTracker.chisqn[1][0] < 10) && 
                    (ntpTracker.chisqn[1][1] < 10);
    
    if( (ntpHeader.trigpatt & 0x2) != 0   )  vars->CUTMASK |= 1 << 0;
    if( minTOF()                          )  vars->CUTMASK |= 1 << 1;
    if( ntpTrd.nseg ==2                   )  vars->CUTMASK |= 1 << 2;
    if( ntpTracker.rig[0] != 0.0          )  vars->CUTMASK |= 1 << 3;
    if( goodChi2                          )  vars->CUTMASK |= 1 << 4;  
    if( goldenTOF()                       )  vars->CUTMASK |= 1 << 5;  
                                                                // 6
    if( ntpHeader.nparticle == 1          )  vars->CUTMASK |= 1 << 7;
    if( ntpTracker.rig[4] != 0.0          )  vars->CUTMASK |= 1 << 8;

    /////////////////////////////// TRACKER ////////////////////////////////////
    
    vars->R     = ntpTracker.rig[1]; // 1 -- Inner tracker
    vars->Rup   = ntpTracker.rig[2]; // 2 -- Upper inner tracker
    vars->Rdown = ntpTracker.rig[3]; // 3 -- Lower inner tracker
    vars->R_L1  = ntpTracker.rig[4]; // 4 -- Inner + L1

    vars->Chisquare         = ntpTracker.chisqn[1][0]; // 1 = Inner      , 0 = X side
    vars->Chisquare_L1      = ntpTracker.chisqn[4][0]; // 4 = L1 + Inner , 0 = X side
    vars->Chisquare_y       = ntpTracker.chisqn[1][1]; // 1 = Inner      , 1 = Y side
    vars->Chisquare_L1_y    = ntpTracker.chisqn[4][1]; // 4 = L1 + Inner , 1 = Y side
    vars->hitbits           = ntpTracker.pattxy; 

    vars->qL1               = ntpTracker.q_lay[1][0];
    vars->qL1Status         = ntpTracker.q_lay_status[1][0];
    vars->qInner            = ntpTracker.q_inn;
    vars->clustertottrack   = ntpHeader.ntrrechit;
    vars->clustertrack      = countBits(vars->hitbits);

    /////////////////////////////// TOF ////////////////////////////////////
    
    // TODO: proper averaging inclding errors?
    vars->qUtof             = ( ntpTof.q_lay[0] + ntpTof.q_lay[1] ) / 2.0;
    vars->qLtof             = ( ntpTof.q_lay[2] + ntpTof.q_lay[3] ) / 2.0;

    vars->BetaOld           = ntpTof.evgeni_beta;
    vars->BetaRaw           = ntpTof.beta;
    vars->BetaHRS           = ntpTof.beta;
    vars->BetaHR            = ntpTof.beta;

    /////////////////////////////// RICH ////////////////////////////////////

    vars->BetaRICH          = ntpRich.beta;
    vars->RICHmask          = RICHmaskConverter();
    vars->Richtotused       = ntpRich.nhit_used;
    vars->RichPhEl          = ntpRich.np_exp_uncorr/ntpRich.np;
    vars->RICHprob          = ntpRich.prob;
    vars->RICHPmts          = ntpRich.npmt;
    vars->RICHcollovertotal = ntpRich.np_exp_uncorr/ntpRich.tot_p;
    vars->RICHgetExpected   = ntpRich.np_exp_uncorr;

    /////////////////////////////// TRD ////////////////////////////////////

    vars->trds_like_ep      = ntpTrd.trds_like_ep ;  
    vars->trds_like_hep     = ntpTrd.trds_like_hep; 
    vars->trds_like_ehe     = ntpTrd.trds_like_ehe;

}

DBarReader::DBarReader(TFile * tfile, bool _isMC) {
    TTree * tree = (TTree *)tfile->Get("Event");
    tree->Branch( "Header"  , "Header" , &ntpHeader     );
    tree->Branch( "Trd"     , "Trd"    , &ntpTrd        );
    tree->Branch( "Tof"     , "Tof"    , &ntpTof        );
    tree->Branch( "Tracker" , "Tracker", &ntpTracker    );
    tree->Branch( "Rich"    , "Rich"   , &ntpRich       );
//  tree->Branch( "Ecal"    , "Ecal"   , &ntpEcal       );
//  tree->Branch( "Anti"    , "Anti"   , &ntpAnti       );
//  tree->Branch( "SA"      , "SA"     , &ntpStandAlone );

    isMC = _isMC;
    if (isMC) tree->Branch("MCHeader","MCHeader",&ntpMCHeader);
}
