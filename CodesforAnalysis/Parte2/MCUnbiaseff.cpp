using namespace std;

Efficiency * EffUnbiasMCP = new Efficiency("EffUnbiasMCP");
Efficiency * EffUnbiasMCD = new Efficiency("EffUnbiasMCD");

void MCUnbiaseff_Fill(TNtuple *ntupla, int l) {
   int k = ntupla->GetEvent(l);
   if((Cutmask&187)!=187||Beta_pre<=0||R_pre<=0) return;
   if(!(EdepTrack<EdepTrackbeta->Eval(Beta_pre)+0.2&&EdepTrack>EdepTrackbeta->Eval(Beta_pre)-0.2)) return;

   int Kbin;

   if(Massa_gen<1&&Massa_gen>0.5) {
      //R bins
      Kbin=GetRBin(fabs(Momento_gen));
      EffUnbiasMCP->beforeR->Fill(Kbin);
      if(Unbias==0) EffUnbiasMCP->afterR->Fill(Kbin);

      //Beta bins
      Kbin=GetArrayBin(Var3, BetaP, nbinsToF);
      EffUnbiasMCP->beforeTOF->Fill(Kbin);
      if(Unbias==0) EffUnbiasMCP->afterTOF->Fill(Kbin);
      
   }

   if(Massa_gen>1&&Massa_gen<2) {
      //R bins
      Kbin=GetRBin(fabs(Momento_gen));
      FillBinMGen(EffUnbiasMCD->beforeR, Kbin);
      if(Unbias==0) FillBinMGen(EffUnbiasMCD->afterR , Kbin);
      
      //Beta bins
      Kbin=GetArrayBin(Var3, BetaD, nbinsToF);
      FillBinMGen(EffUnbiasMCD->beforeTOF, Kbin);
      if(Unbias==0) FillBinMGen(EffUnbiasMCD->afterTOF , Kbin);
      
   }

   return;
}


void MCUnbiaseff_Write() {
   EffUnbiasMCP->Write();
   EffUnbiasMCD->Write();
   return;
}



void MCUnbiaseff(TFile * file1) {

   Efficiency * EffUnbiasMCP = new Efficiency(file1,"EffUnbiasMCP" );
   Efficiency * EffUnbiasMCD = new Efficiency(file1,"EffUnbiasMCD" );

   string nome;
   Tempi = (TH1F *)file1->Get("Tempi");

   cout<<"**** MC Unbias TRIGGER EFF. ****"<<endl;

   EffUnbiasMCP -> Eval_Efficiency();
   EffUnbiasMCD -> Eval_Efficiency();

   TH1F *EffUnbMCP_R_TH1F = (TH1F*)  EffUnbiasMCP->effR   ->Clone();
   TH1F *EffUnbMCP_TH1F   = (TH1F*)  EffUnbiasMCP->effTOF->Clone();
   TH2F *EffUnbMCD_R_TH2F = (TH2F*)  EffUnbiasMCD->effR   ->Clone();
   TH2F *EffUnbMCD_TH2F   = (TH2F*)  EffUnbiasMCD->effTOF->Clone();


   cout<<"*** Updating P1 file ****"<<endl;
   string nomefile="../Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
   file1 =TFile::Open(nomefile.c_str(),"UPDATE");

   file1->cd("Results");
   EffUnbMCP_R_TH1F  -> Write();
   EffUnbMCP_TH1F   -> Write();
   EffUnbMCD_R_TH2F  -> Write();
   EffUnbMCD_TH2F    -> Write();
   file1->Write();
   file1->Close();

   TCanvas *c11=new TCanvas("Unbias Trigger Efficiency");
   c11->Divide(2,1);
   c11->cd(1);
   gPad->SetLogx();
   gPad->SetGridx();
   gPad->SetGridy();
   string MCLegend[7]= {"protons.B800","d.pl1.0_520_GG_Blic","d.pl1.0_520_GG_BlicDPMJet","d.pl1.0_520_GG_QMD","d.pl1.0_520_Shen_Blic","d.pl1.0_520_Shen_BlicDPMJet","d.pl1.0_520_Shen_QMD"};
   TGraph * EffUnbMCP_R = new TGraph();
   EffUnbMCP_R->SetTitle(MCLegend[0].c_str());
   for(int i=0; i<nbinsr; i++) EffUnbMCP_R->SetPoint(i,R_cent[i],1);
   TGraph * EffUnbMCD_R[6];
   EffUnbMCP_R->SetMarkerColor(2);
   EffUnbMCP_R->SetMarkerStyle(8);
   EffUnbMCP_R->SetLineColor(2);
   EffUnbMCP_R->SetLineWidth(2);
   EffUnbMCP_R->SetTitle("Unbias Trigger Efficiency (R bins)");
   EffUnbMCP_R->GetXaxis()->SetTitle("R [GV]");
   EffUnbMCP_R->GetYaxis()->SetTitle("Efficiency");
   EffUnbMCP_R->GetXaxis()->SetTitleSize(0.045);
   EffUnbMCP_R->GetYaxis()->SetTitleSize(0.045);
   {
      EffUnbMCP_R->Draw("ACP");
      TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
      leg->AddEntry(EffUnbMCP_R,MCLegend[0].c_str(), "ep");

      for(int h=0; h<6; h++) {
         EffUnbMCD_R[h]= new TGraph();
         EffUnbMCD_R[h]->SetTitle(MCLegend[h+1].c_str());
         for(int i=1; i<nbinsr; i++) EffUnbMCD_R[h]->SetPoint(i,R_cent[i],EffUnbMCD_R_TH2F->GetBinContent(i+1,h+1));
         leg->AddEntry(EffUnbMCD_R[h],MCLegend[h+1].c_str(), "ep");
         EffUnbMCD_R[h]->SetMarkerColor(4);
         EffUnbMCD_R[h]->SetMarkerStyle(h+3);
         EffUnbMCD_R[h]->SetMarkerSize(2);
         EffUnbMCD_R[h]->SetLineColor(4);
         EffUnbMCD_R[h]->SetLineWidth(2);
         // EffUnbMCD_R[h]->Draw("Psame");
         leg->Draw();
      }
   }

   c11->cd(2);
   gPad->SetLogx();
   gPad->SetGridx();
   gPad->SetGridy();
   TGraph * EffUnbMCP = new TGraph();
   for(int i=0; i<17; i++) EffUnbMCP->SetPoint(i,Ekincent[i],1);
   TGraph * EffUnbMCD[6];
   EffUnbMCP->SetMarkerColor(2);
   EffUnbMCP->SetMarkerStyle(8);
   EffUnbMCP->SetLineColor(2);
   EffUnbMCP->SetLineWidth(2);
   EffUnbMCP->SetTitle("Unbias Trigger Efficiency (Beta bins)");
   EffUnbMCP->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
   EffUnbMCP->GetYaxis()->SetTitle("Efficiency");
   EffUnbMCP->GetXaxis()->SetTitleSize(0.045);
   EffUnbMCP->GetYaxis()->SetTitleSize(0.045);
   {
      EffUnbMCP->Draw("ACP");
      TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
      leg->AddEntry(EffUnbMCP,MCLegend[0].c_str(), "ep");

      for(int h=0; h<6; h++) {
         EffUnbMCD[h]= new TGraph();
         for(int i=0; i<17; i++) EffUnbMCD[h]->SetPoint(i,Ekincent[i],EffUnbMCD_TH2F->GetBinContent(i+1,h+1));
         EffUnbMCD[h]->SetMarkerColor(4);
         EffUnbMCD[h]->SetMarkerStyle(h+3);
         leg->AddEntry(EffUnbMCD[h],MCLegend[h+1].c_str(), "ep");
         EffUnbMCD[h]->SetMarkerSize(2);
         EffUnbMCD[h]->SetLineColor(4);
         EffUnbMCD[h]->SetLineWidth(2);
         //    EffUnbMCD[h]->Draw("Psame");
         leg->Draw();
      }
   }
   cout<<"*** Updating Results file ***"<<endl;
   nomefile="./Final_plots/"+mese+".root";
   TFile *f_out=new TFile(nomefile.c_str(), "UPDATE");
   f_out->cd("MC Results/Preselections");
   c11->Write();
   f_out->Write();
   f_out->Close();
}

