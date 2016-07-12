



void LoopOnMCTrig(TNtuple*  ntup)
{
   int nentries=ntupMCTrig->GetEntries();
   EffUnbiasMCP = new Efficiency("EffUnbiasMCP");
   EffUnbiasMCD = new Efficiency("EffUnbiasMCD");
   for(int i=0; i<ntupMCTrig->GetEntries(); i++) {
      ntupMCTrig->GetEvent(i);
      cmask.setMask(Tup.Cutmask);
      Cuts_Pre();
      Massa_gen = ReturnMass_Gen();
      RUsed=Tup.R_pre;
      UpdateProgressBar(i, nentries);

      MCpreseff_Fill();
      MCUnbiaseff_Fill();
      MCTrackeff_Fill();
      MigrationMatrix_Fill();
      Correlazione_Preselezioni();
      FluxFactorizationtest_Pre_Fill();
      DVSMCTrackeff_Fill();
      DVSMCPreSeleff_Fill();
      DVSMCPreSeleffD_Fill();
   }
   cout << endl;
   return;
}


