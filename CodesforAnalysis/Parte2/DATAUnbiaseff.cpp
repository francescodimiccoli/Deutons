using namespace std;

Efficiency * EffUnbiasDATA = new Efficiency ("EffUnbiasDATA");

void DATAUnbiaseff_Fill () {
   if (!cmask.isPreselected() ||Tup.R_pre<=0||Tup.R_pre<1.2*Tup.Rcutoff) return;
   if (!(Tup.EdepTrack<EdepTrackbeta->Eval (Tup.Beta_pre)+0.2&&Tup.EdepTrack>EdepTrackbeta->Eval (Tup.Beta_pre)-0.2)) return;	

   int Kbin=PRB.GetRBin (fabs (Tup.R_pre) );
   if(Tup.Unbias==0) EffUnbiasDATA->beforeR->Fill(Kbin);
   if(Tup.Unbias==1) EffUnbiasDATA->beforeR->Fill(Kbin,100);
   if(Tup.Unbias==0) EffUnbiasDATA->afterR->Fill(Kbin);

   return;
}


void DATAUnbiaseff_Write() {
   EffUnbiasDATA -> Write();
   return;
}


void DATAUnbiaseff (TFile * inputHistoFile) {
   Efficiency * EffUnbiasDATA = new Efficiency (inputHistoFile,"EffUnbiasDATA");

   cout<<"********** DATA Tup.Unbias TRIGG. EFFICIENCY ******************************"<<endl;

   EffUnbiasDATA -> Eval_Efficiency();

   TH1F *EffUnbDATA_R_TH1F = (TH1F *) EffUnbiasDATA -> effR   ->Clone();
   TH1F *EffUnbDATA_TH1F   = (TH1F *) EffUnbiasDATA -> effTOF ->Clone();

  float EffMean = 0;
   for(int n = 0; n< EffUnbDATA_R_TH1F->GetNbinsX();n++){
		EffMean+=EffUnbDATA_R_TH1F->GetBinContent(n+1);
	}
   EffMean=EffMean/(float)EffUnbDATA_R_TH1F->GetNbinsX();	

   TH1F * TriggerGlobalFactor = new TH1F("TriggerGlobalFactor","TriggerGlobalFactor",1,0,1);
   TriggerGlobalFactor -> SetBinContent(1,EffMean);
   TriggerGlobalFactor -> SetBinError(1,0.01); 

   cout<<"*** Updating P1 file ****"<<endl;
   inputHistoFile->ReOpen("UPDATE");

   inputHistoFile->cd ("Results");
   EffUnbDATA_R_TH1F ->Write();
   TriggerGlobalFactor -> Write();
   inputHistoFile-> Write();
   inputHistoFile->Close();

   TCanvas *c12=new TCanvas ("DATA: Unb. Trigger Efficiency");

   c12->cd ();
   gPad->SetLogx();
   gPad->SetGridx();
   gPad->SetGridy();
   string MCLegend[2]= {"protons","deutons"};
   TGraph * EffUnbDATA_R = new TGraph();
   EffUnbDATA_R->SetTitle (MCLegend[0].c_str() );
   for (int i=0; i<nbinsr; i++) EffUnbDATA_R->SetPoint (i,PRB.RigBinCent (i),EffUnbDATA_R_TH1F->GetBinContent (i+1) );
   EffUnbDATA_R->SetMarkerColor (2);
   EffUnbDATA_R->SetMarkerStyle (8);
   EffUnbDATA_R->SetLineColor (2);
   EffUnbDATA_R->SetLineWidth (2);
   EffUnbDATA_R->SetTitle ("Physical Trigg. Efficiency  (R bins)");
   EffUnbDATA_R->GetYaxis()->SetRangeUser(0,1);
   EffUnbDATA_R->GetXaxis()->SetTitle ("R [GV]");
   EffUnbDATA_R->GetYaxis()->SetTitle ("Efficiency");
   EffUnbDATA_R->GetXaxis()->SetTitleSize (0.045);
   EffUnbDATA_R->GetYaxis()->SetTitleSize (0.045);
   {
      EffUnbDATA_R->Draw ("ACP");
      TLegend* leg =new TLegend (0.4, 0.7,0.95,0.95);
      leg->AddEntry (EffUnbDATA_R,MCLegend[0].c_str(), "ep");

   }

   finalPlots.Add(c12);
   finalPlots.writeObjsInFolder("DATA-driven Results");

}
