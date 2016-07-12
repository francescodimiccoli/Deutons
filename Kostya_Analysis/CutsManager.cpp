#include <cmath>
#include "CutsManager.h"
#include "cutmask.h"


 
void CutsManager::SelectionCuts(const TupleVars & tup) {
    Cutmask cmask(tup.Cutmask);
    PassPreselection = false;
    if(    tup.Unbias == 0       
        && cmask.isPreselected() 
        && tup.Beta_pre > 0      
        && tup.R_pre    > 0   ) PassPreselection = true;

    MinBias = cmask.isMinimumBiasTrigger();
    TOF3or4 = MinBias && cmask.isMinimumBiasToF3or4Layers() && tup.Beta_pre>0;

    TOFEdepGood  = TOF3or4 && cal.TOFBetaInside1(tup.Beta_pre, tup.EdepTOFU);
    MinBiasTrack = TOFEdepGood && cmask.isMinimumBiasTracker();
}



void CutsManager::TestDiscriminats(const TupleVars & tup) {
    PassLikelihood = false;
    PassDistance   = false;

    if(isnan(tup.LDiscriminant)) return;

    float logL = - log(1 - tup.LDiscriminant);
    if(  (NaF||Agl) && logL > 2.6  ) PassLikelihood = true;
    if( !(NaF||Agl) && logL > 0.55 ) PassLikelihood = true;

    PassDistance = ( tup.Dist5D < 4 || tup.Dist5D_P < 4);
}



void CutsManager::ProcessEvent(const TupleVars & tup) {
    HeliumRejection = false;
    float EdepTOFud = (tup.EdepTOFU + tup.EdepTOFD)/2;
    if(cal.OffsetTrackBeta( tup.Beta_pre, tup.EdepTrack ) < 4  ) HeliumRejection = true;
    if(cal.OffsetTOFBeta  ( tup.Beta_pre, EdepTOFud     ) < 10 ) HeliumRejection = true;

    Cutmask cmask(tup.Cutmask);
    NaF = cmask.isFromNaF();
    Agl = cmask.isFromAgl();

    BetaControlRegion = false;
    if( NaF && tup.BetaRICH < 0.97  ) BetaControlRegion = true;
    if( Agl && tup.BetaRICH < 0.985 ) BetaControlRegion = true;
    if( !( NaF || Agl ) && tup.Beta_pre < 0.8 ) BetaControlRegion = true;

    TestDiscriminats(tup);
    SelectionCuts(tup);
}

