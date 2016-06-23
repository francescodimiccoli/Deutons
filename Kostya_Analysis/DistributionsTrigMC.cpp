#include "Distributions.h"

// This reperesents efficiencies as a function of the true MC momentum
// for each subdetector
class EfficienciesTrue {

    TEfficiency * Preselection;
    TEfficiency * Unbias;
    TEfficiency * Trigger;
    TEfficiency * TOF;
    TEfficiency * Tracker;

    EfficienciesTrigMC( std::vector<float> bins ){
        Preselection = new TEfficiency("Preselection", "Preselection", bins.size()-1, &bins[0]);
        Unbias       = new TEfficiency("Unbias",       "Unbias",       bins.size()-1, &bins[0]);
        Trigger      = new TEfficiency("Trigger",      "Trigger",      bins.size()-1, &bins[0]);
        TOF          = new TEfficiency("TOF",          "TOF",          bins.size()-1, &bins[0]);
        Tracker      = new TEfficiency("Tracker",      "Tracker",      bins.size()-1, &bins[0]);
    }

    void Fill(const TupleVars & tup, const CutsManager & cuts) {
        Unbias      ->Fill(tup.Unbias == 0,           tup.Momento_gen );
        Preselection->Fill(cuts.passesPreselection(), tup.Momento_gen );

        bool minbias = cuts.isMinimumBiasTrigger();
        Trigger->Fill(minbias, tup.Momento_gen); 

        bool tof34 = cuts.isTOF3or4();
        if( minbias ) TOF->Fill(tof34, tup.Momento_gen); 

        // NOte -- that insert an extra  cut
        bool edep    = cuts.isTOFEdepGood();
        bool mbtrack = cuts.isMinimumBiasTracker();
        if(minbias && tof34 && edep ) Tracker->Fill(mbtrack, tup.Momento_gen)
    }
};

class ParticleDistributions {
    EfficienciesTrue * Rigidity,  * TOF, * NaF, * Agl;
    TH2F * MigrationMatrix;
public:
    DistributionsTrigMC(ParticleBins b){
        Rigidity = new EfficienciesTrigMC(b.Rigidity().GetRBins());
        TOF = new EfficienciesTrigMC(b.TOF().GetRBins());
        NaF = new EfficienciesTrigMC(b.NaF().GetRBins());
        Agl = new EfficienciesTrigMC(b.Agl().GetRBins());
    }

    void Fill(const TupleVars & tup, const CutsManager & cuts)
    {
        Rigidity->Fill(tup, cuts);
        TOF->Fill(tup, cuts);
        if(cuts.isEventNaF() ) NaF->Fill(tup, cuts);
        if(cuts.isEventAgl() ) Agl->Fill(tup, cuts);
    }
};

class PerParticle {
    ParticleDistributions * Rigidity,  * TOF, * NaF, * Agl;
    TH2F * MigrationMatrix;
public:
    PerParticle(ParticleBins b){
        Rigidity = new EfficienciesTrigMC(b.Rigidity().GetRBins());
        TOF = new EfficienciesTrigMC(b.TOF().GetRBins());
        NaF = new EfficienciesTrigMC(b.NaF().GetRBins());
        Agl = new EfficienciesTrigMC(b.Agl().GetRBins());
    }

    void Fill(const TupleVars & tup, const CutsManager & cuts)
    {
        Rigidity->Fill(tup, cuts);
        TOF->Fill(tup, cuts);
        if(cuts.isEventNaF() ) NaF->Fill(tup, cuts);
        if(cuts.isEventAgl() ) Agl->Fill(tup, cuts);
    }
};


class TriggerMC  : public Distributions {
    CutsManager cuts;
    AnalysisBins bins;
    PerParticle * Proton;
    std::map<int, PerParticle *> Deuteron;
public:
    TriggerMC(AnalysisBins b, const Calibrations & calibrations):
        bins(b), cuts(calibrations)
    {
        Proton = new PerParticle(b.Proton());
        Deuteron[0] = new PerParticle(b.Deuteron());
        Deuteron[1] = new PerParticle(b.Deuteron());
        Deuteron[2] = new PerParticle(b.Deuteron());
        Deuteron[3] = new PerParticle(b.Deuteron());
        Deuteron[4] = new PerParticle(b.Deuteron());
        Deuteron[5] = new PerParticle(b.Deuteron());
    }

    virtual void Fill(const TupleVars & tup)
    {
        cuts.ProcessEvent(tup);

        if( tup.GetMCType() == TupleVars::MCTYPE::Proton) proton->Fill(tup, cuts);
        if( tup.GetMCType() == TupleVars::MCTYPE::Deuteron ) {
            int mctype = tup.MCsubtype();
            deuteron[mctype]->Fill(tup, cuts);
        }
    }
};


Distributions * Distributions::TriggerMC(AnalysisBins b, const Calibrations & calibrations){
    return new TriggerMC(b, calibrations):
}

