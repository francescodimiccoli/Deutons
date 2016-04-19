using namespace std;

Efficiency * EffUnbiasDATA = new Efficiency("EffUnbiasDATA");

void DATAUnbiaseff_Fill(TNtuple *ntupla, int l) {
   int k = ntupla->GetEvent(l);
   if(((int)Cutmask&187)!=187||R_pre<=0||Beta_pre<=0||R_pre<1.2*Rcutoff) return;
   
   int Kbin=GetRBin(fabs(R_pre));
   if(EdepTrack<EdepTrackbeta->Eval(Beta_pre)+0.2&&EdepTrack>EdepTrackbeta->Eval(Beta_pre)-0.2) {
      EffUnbiasDATA->beforeR->Fill(Kbin);
      if(Unbias==1) EffUnbiasDATA->afterR->Fill(Kbin);
   }

   if(EdepTrack<EdepTrackbeta->Eval(Beta_pre)+0.2&&EdepTrack>EdepTrackbeta->Eval(Beta_pre)-0.2) {
      Kbin=GetArrayBin(Var, BetaP, nbinsToF);
      EffUnbiasDATA->beforeTOF->Fill(Kbin);
      if(Unbias==1) EffUnbiasDATA->afterTOF->Fill(Kbin);
   }
   
   return;
}


void DATAUnbiaseff_Write() {
   EffUnbiasDATA -> Write();
   return;
}


void DATAUnbiaseff(TFile * file1) {
   Efficiency * EffUnbiasDATA = new Efficiency(file1,"EffUnbiasDATA");

   cout<<"********** DATA Unbias TRIGG. EFFICIENCY ******************************"<<endl;

   EffUnbiasDATA -> Eval_Efficiency();

   TH1F *EffUnbDATA_R_TH1F = (TH1F *) EffUnbiasDATA -> effR   ->Clone();
   TH1F *EffUnbDATA_TH1F   = (TH1F *) EffUnbiasDATA -> effTOF ->Clone();

   cout<<"*** Updating P1 file ****"<<endl;
   string nomefile="../Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
   file1 =TFile::Open(nomefile.c_str(),"UPDATE");

   file1->cd("Results");
   EffUnbDATA_R_TH1F ->Write();
   EffUnbDATA_TH1F	  ->Write();
   file1-> Write();
   file1-> Close();

   TCanvas *c12=new TCanvas("DATA: Unb. Trigger Efficiency");

   c12->Divide(2,1);
   c12->cd(1);
   gPad->SetLogx();
   gPad->SetGridx();
   gPad->SetGridy();
   string MCLegend[2]= {"protons","deutons"};
   TGraph * EffUnbDATA_R = new TGraph();
   EffUnbDATA_R->SetTitle(MCLegend[0].c_str());
   for(int i=0; i<nbinsr; i++) EffUnbDATA_R->SetPoint(i,R_cent[i],EffUnbDATA_R_TH1F->GetBinContent(i+1));
   TGraph * EffUnbMCD_R[6];
   EffUnbDATA_R->SetMarkerColor(2);
   EffUnbDATA_R->SetMarkerStyle(8);
   EffUnbDATA_R->SetLineColor(2);
   EffUnbDATA_R->SetLineWidth(2);
   EffUnbDATA_R->SetTitle("Physical Trigg. Efficiency  (R bins)");
   EffUnbDATA_R->GetXaxis()->SetTitle("R [GV]");
   EffUnbDATA_R->GetYaxis()->SetTitle("Efficiency");
   EffUnbDATA_R->GetXaxis()->SetTitleSize(0.045);
   EffUnbDATA_R->GetYaxis()->SetTitleSize(0.045);
   {
      EffUnbDATA_R->Draw("ACP");
      TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
      leg->AddEntry(EffUnbDATA_R,MCLegend[0].c_str(), "ep");

   }

   c12->cd(2);
   gPad->SetLogx();
   gPad->SetGridx();
   gPad->SetGridy();
   TGraph * EffUnbDATA = new TGraph();
   for(int i=0; i<nbinsbeta; i++) EffUnbDATA->SetPoint(i,Ekincent[i],EffUnbDATA_TH1F->GetBinContent(i+1));
   TGraph * EffUnbMCD[6];
   EffUnbDATA->SetMarkerColor(2);
   EffUnbDATA->SetMarkerStyle(8);
   EffUnbDATA->SetLineColor(2);
   EffUnbDATA->SetLineWidth(2);
   EffUnbDATA->SetTitle("Physical Trigg. Efficiency  (Beta bins)");
   EffUnbDATA->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
   EffUnbDATA->GetYaxis()->SetTitle("Efficiency");
   EffUnbDATA->GetXaxis()->SetTitleSize(0.045);
   EffUnbDATA->GetYaxis()->SetTitleSize(0.045);
   {
      EffUnbDATA->Draw("ACP");
      TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
      leg->AddEntry(EffUnbDATA,MCLegend[0].c_str(), "ep");

   }

   cout<<"*** Updating Results file ***"<<endl;
   nomefile="./Final_plots/"+mese+".root";
   TFile *f_out=new TFile(nomefile.c_str(), "UPDATE");
   f_out->mkdir("DATA-driven Results");
   f_out->cd("DATA-driven Results");
   c12->Write();
   f_out->Write();
   f_out->Close();


}

