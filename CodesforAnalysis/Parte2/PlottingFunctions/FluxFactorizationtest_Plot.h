


void FluxFactorizationtest_Plot(
        TH1 * Eff_FullSETMCP_R_TH1F,
        TH1 *FactorizedEffMCP_R
        ){

   TCanvas *c9 = new TCanvas("MC Protons Factorization Test");
   c9->cd();
   gPad->SetLogx();
   gPad->SetGridx();
   gPad->SetGridy();

   TGraphErrors * FullsetEfficiency    = new TGraphErrors();
   TGraphErrors * FactorizedEfficiency = new TGraphErrors();

   for(int i=0; i<Eff_FullSETMCP_R_TH1F->GetNbinsX(); i++) {
      FullsetEfficiency    ->SetPoint(i,PRB.RigBinCent(i),Eff_FullSETMCP_R_TH1F->GetBinContent(i+1));
      FactorizedEfficiency ->SetPoint(i,PRB.RigBinCent(i),FactorizedEffMCP_R->GetBinContent(i+1));
   }

   FullsetEfficiency->SetMarkerColor(2);
   FullsetEfficiency->SetMarkerStyle(8);
   FullsetEfficiency->SetLineColor(2);
   FullsetEfficiency->SetLineWidth(2);

   FactorizedEfficiency->SetMarkerColor(2);
   FactorizedEfficiency->SetMarkerStyle(4);
   FactorizedEfficiency->SetLineColor(2);
   FactorizedEfficiency->SetLineWidth(2);

   FullsetEfficiency->Draw("APC");
   FactorizedEfficiency->Draw("PCsame");

   finalPlots.Add(c9);
   finalPlots.writeObjsInFolder("MC Results/Eff. Factorization Test");

   return;
}
