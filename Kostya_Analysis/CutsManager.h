#ifndef CUTSMANAGER__H__
#define CUTSMANAGER__H__

#include "TupleVars.hpp"
#include "Calibrations.h"

class CutsManager {
    const Calibrations & cal;

    bool NaF, Agl;

    bool HeliumRejection;
    bool BetaControlRegion;
    bool PassLikelihood;
    bool PassDistance;
    bool PassPreselection;
    bool MinBias;
    bool TOF3or4;
    bool TOFEdepGood;
    bool MinBiasTrack;
    void SelectionCuts(const TupleVars & tup);
    void TestDiscriminats(const TupleVars & tup);

public:
    bool isEventTOF() { return true; }
    bool isEventNaF() { return  NaF; }
    bool isEventAgl() { return  Agl; }
    bool isHelumRejected()     { return HeliumRejection;   } 
    bool isBetaControlRegion() { return BetaControlRegion; }
    bool isMinimumBiasTrigger(){ return MinBias;           }
    bool isTOF3or4()           { return TOF3or4;           }
    bool isTOFEdepGood()       { return TOFEdepGood;       }
    bool isMinBiasTrack()      { return MinBiasTrack;      }

    bool passesPreselection()  { return PassPreselection;  }
    bool passesLikelihood()    { return PassLikelihood;    } 
    bool passesDistance()      { return PassDistance;      }

    CutsManager(const Calibrations & c):cal(c) {}
    void ProcessEvent(const TupleVars & tup);
};

#endif
