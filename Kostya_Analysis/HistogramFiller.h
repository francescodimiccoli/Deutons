#include "TEfficiency.h"

class NTupleLooper {
    CutsManager cuts;
    AnalysisBins bins;

    EfficienciesMC * proton;
    std::map<int, EfficienciesMC *> deuteron;

public:
    TrigMCFiller(Binning b, const Calibrations & calibrations);
    void Loop(TNtuple * ntup);
};

TrigMCFiller::TrigMCFiller(AnalysisBins b, const Calibrations & calibrations): 
    bins(b), cuts(calibrations)
{
    proton = new EfficienciesMC(b.Proton());
    deuteron[0] = new EfficienciesMC(b.Deuteron());
    deuteron[1] = new EfficienciesMC(b.Deuteron());
    deuteron[2] = new EfficienciesMC(b.Deuteron());
    deuteron[3] = new EfficienciesMC(b.Deuteron());
    deuteron[4] = new EfficienciesMC(b.Deuteron());
    deuteron[5] = new EfficienciesMC(b.Deuteron());
}

void Loop(TNtuple *  ntup)
{
    TupleVars tup(ntup);
    
    int nentries = ntup->GetEntries();
    for(int i=0; i<ntup->GetEntries(); i++) {
        ntup->GetEvent(i);
        cuts.ProcessEvent(tup);


        


        // First I create a vector of efficiencies to fill
        EfficienciesMC * effsToFill = NULL;
        if( tup.GetMCType() == TupleVars::MCTYPE::Proton) effsToFill = proton;
        if( tup.GetMCType() == TupleVars::MCTYPE::Deuteron ) {
            int mctype = tup.MCsubtype();
            effsToFill = deuteron[mctype];
        }

        effsToFill->Fill(tup);
        
/*        
      MigrationMatrix_Fill();
      Correlazione_Preselezioni();
      FluxFactorizationtest_Pre_Fill();

      DVSMCTrackeff_Fill();
      DVSMCPreSeleff_Fill();
      DVSMCPreSeleffD_Fill();
*/

   }
   return;
}


