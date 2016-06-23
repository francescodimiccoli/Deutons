#include "TEfficiency.h"

// Efficiencies of the  top-level selections for each subdetector
class EfficienciesQualTrue {
    TEfficiency * Distance;
    TEfficiency * DistanceLikelihood;

    EfficienciesQualMC( std::vector<float> bins ){
        Distance           = new TEfficiency("Distance",           "Distance",           bins.size()-1, &bins[0]);
        DistanceLikelihood = new TEfficiency("DistanceLikelihood", "DistanceLikelihood", bins.size()-1, &bins[0]);

    void Fill(const TupleVars & tup,, const CutsManager & cuts) {
        bool dist = cuts.passesDistance();
        bool lklh = cuts.passesLikelihood();
        Distance          ->Fill(dist,         tup.Momento_gen );
        DistanceLikelihood->Fill(dist && lklh, tup.Momento_gen );
    }
};

// The distributions (for each MC particle/type) 
// that we want to build when processing "Qual" ntuples
struct DistributionsQualMC {
    EfficienciesQualTrue * Rigidity,  * TOF, * NaF, * Agl;
    TH2F * ECAL_vs_Rgen;

    EfficienciesMC(ParticleBins b){
        Rigidity = new EfficienciesQualMC(b.Rigidity().GetRBins());
        TOF = new EfficienciesQualMC(b.TOF().GetRBins());
        NaF = new EfficienciesQualMC(b.NaF().GetRBins());
        Agl = new EfficienciesQualMC(b.Agl().GetRBins());

        ECAL_vs_Rgen = new TH2F("ECALvsR_MC","ECALvsR_MC",1000,0,100,1000,0,100);
    }

    void Fill(const TupleVars & tup, const CutsManager & cuts)
    {
        Rigidity->Fill(tup, cuts);
        TOF->Fill(tup, cuts);
        if(cuts.isEventNaF() ) NaF->Fill(tup, cuts);
        if(cuts.isEventAgl() ) Agl->Fill(tup, cuts);

		ECAL_vs_Rgen->Fill(tup.EdepECAL, tup.Momento_gen);	
    }
};


