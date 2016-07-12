#include "Distributions.h"
#include "TH1F.h"
#include "binning.h"

class EfficienciesHists {
    TH1F * TotalTriggered; 
    TH1F * UnbiasOnly;     
public:

    EfficienciesHists( std::vector<float> bins ){
        TotalTriggered = new TH1F("TotalTriggered", "TotalTriggered", bins.size()-1, &bins[0]);
        UnbiasOnly     = new TH1F("UnbiasOnly",     "UnbiasOnly",     bins.size()-1, &bins[0]);
    }

    void Fill(const TupleVars & tup, const CutsManager & cuts) {
        TotalTriggered->Fill(tup.Momento_gen);
        if(tup.Unbias == 0) UnbiasOnly->Fill(tup.Momento_gen );
    }
};

class PerParticle {
    EfficienciesHists * Rigidity, * TOF, * NaF, * Agl;
    TH2F * MigrationMatrix;
public:
    PerParticle(ParticleBins b){
        Rigidity = new EfficienciesHists(b.Rigidity().GetRBins());
        TOF = new EfficienciesHists(b.TOF().GetRBins());
        NaF = new EfficienciesHists(b.NaF().GetRBins());
        Agl = new EfficienciesHists(b.Agl().GetRBins());
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

