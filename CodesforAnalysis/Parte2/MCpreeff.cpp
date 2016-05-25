using namespace std;



Efficiency * EffpreselMCP = new Efficiency("EffpreselMCP");
Efficiency * EffpreselMCD = new Efficiency("EffpreselMCD", 6);


void MCpreseff_Fill(TNtuple *ntupla, int l) {
    ntupla->GetEvent(l);

   if(Massa_gen<1&&Massa_gen>0.5) {
      //R bins
      EffpreselMCP->beforeR->Fill(RB.GetRBin(fabs(Tup.Momento_gen)));
      if(Tup.Unbias==0&&((int)Tup.Cutmask&187)==187&&Tup.Beta_pre>0&&Tup.R_pre>0) EffpreselMCP->afterR->Fill(RB.GetRBin(fabs(Tup.R_pre)));
         
      // Beta bins
 ToFPB.wRFill(EffpreselMCP->beforeTOF, Tup.Momento_gen);
 NaFPB.wRFill(EffpreselMCP->beforeNaF, Tup.Momento_gen);
 AglPB.wRFill(EffpreselMCP->beforeAgl, Tup.Momento_gen);

      if(Tup.Unbias==0 && ((int)Tup.Cutmask&187)==187 && Tup.Beta_pre>0 && Tup.R_pre>0)
      {
                                     ToFPB.wRFill(EffpreselMCP->afterTOF, RUsed);
         if(((int)Tup.Cutmask)>>11==512) NaFPB.wRFill(EffpreselMCP->afterNaF, RUsed);
         if(((int)Tup.Cutmask)>>11==0  ) AglPB.wRFill(EffpreselMCP->afterAgl, RUsed);
      }

   }

   if(Massa_gen>1&&Massa_gen<2) {
      // R bins      
      FillBinMGen(EffpreselMCD->beforeR, RB.GetRBin(fabs(Tup.Momento_gen)));
      if(((int)Tup.Cutmask&187)==187&&Tup.Beta_pre>0&&Tup.Unbias==0&&Tup.R_pre>0)
         FillBinMGen(EffpreselMCD->afterR, RB.GetRBin(fabs(Tup.R_pre)));

      // Beta bins

         FillBinMGen(EffpreselMCD->beforeTOF, ToFDB.GetRBin(Tup.Momento_gen) );
         FillBinMGen(EffpreselMCD->beforeNaF, NaFDB.GetRBin(Tup.Momento_gen) );
         FillBinMGen(EffpreselMCD->beforeAgl, AglDB.GetRBin(Tup.Momento_gen) );

         if(((int)Tup.Cutmask&187)==187 && Tup.Beta_pre>0 && Tup.Unbias==0 && Tup.R_pre>0)
         {
                                           FillBinMGen(EffpreselMCD->afterTOF, ToFDB.GetRBin(RUsed));
            if(((int)Tup.Cutmask)>>11 == 512 ) FillBinMGen(EffpreselMCD->afterNaF, NaFDB.GetRBin(RUsed));
            if(((int)Tup.Cutmask)>>11 == 0   ) FillBinMGen(EffpreselMCD->afterAgl, AglDB.GetRBin(RUsed));
         }
      
   }
   return;
}


void MCpreeff_Write() {
   EffpreselMCP->Write();
   EffpreselMCD->Write();
   return;
}



void MCpreeff(TFile * file1) {
   Efficiency * EffpreselMCP = new Efficiency(file1, "EffpreselMCP");
   Efficiency * EffpreselMCD = new Efficiency(file1, "EffpreselMCD");

   string tagli[10]= {"Trigger","3of4 TOF","TRD Segments","Rigidity exists","Chi^2 R","Matching TOF","Matching TRD","In TRD Accept.","1 Particle","1 Tr. Track"};
   string nome;

   cout<<"**** MC PRESELECTIONS EFFICIENCY (FULL SET) ****"<<endl;

   EffpreselMCP -> Eval_Efficiency();
   EffpreselMCD -> Eval_Efficiency();

   TH1F * EffPreMCP_R_TH1F  =  (TH1F *)EffpreselMCP->effR	->Clone();
   TH1F * EffPreMCP_TH1F    =  (TH1F *)EffpreselMCP->effTOF->Clone();
   TH1F * EffPreMCPNaF_TH1F =  (TH1F *)EffpreselMCP->effNaF->Clone();
   TH1F * EffPreMCPAgl_TH1F =  (TH1F *)EffpreselMCP->effAgl->Clone();
   TH2F * EffPreMCD_R_TH2F  =  (TH2F *)EffpreselMCD->effR  ->Clone();
   TH2F * EffPreMCD_TH2F    =  (TH2F *)EffpreselMCD->effTOF->Clone();
   TH2F * EffPreMCDNaF_TH2F =  (TH2F *)EffpreselMCD->effNaF->Clone();
   TH2F * EffPreMCDAgl_TH2F =  (TH2F *)EffpreselMCD->effAgl->Clone();


   cout<<"*** Updating P1 file ****"<<endl;
   string nomefile="../Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
   file1 =TFile::Open(nomefile.c_str(),"UPDATE");

   file1->mkdir("Results");
   file1->cd("Results");
   EffPreMCP_R_TH1F->Write();
   EffPreMCP_TH1F ->Write();
   EffPreMCPNaF_TH1F->Write();
   EffPreMCPAgl_TH1F ->Write();
   EffPreMCD_R_TH2F->Write();
   EffPreMCD_TH2F->Write();
   EffPreMCDNaF_TH2F ->Write();
   EffPreMCDAgl_TH2F->Write();
   file1->Write();
   file1->Close();


   TCanvas *c4=new TCanvas("Preselections Efficiency (R bins)");
   TCanvas *c4_bis=new TCanvas("Preselections Efficiency (Beta bins)");
   c4_bis->Divide(3,1);

   c4->cd();
   gPad->SetLogx();
   gPad->SetGridx();
   gPad->SetGridy();
   string MCLegend[7]= {"protons.B800","d.pl1.0_520_GG_Blic","d.pl1.0_520_GG_BlicDPMJet","d.pl1.0_520_GG_QMD","d.pl1.0_520_Shen_Blic","d.pl1.0_520_Shen_BlicDPMJet","d.pl1.0_520_Shen_QMD"};
   TGraph * EffPreMCP_R = new TGraph();
   EffPreMCP_R->SetTitle(MCLegend[0].c_str());
   for(int i=0; i<nbinsr; i++) EffPreMCP_R->SetPoint(i,R_cent[i],EffPreMCP_R_TH1F->GetBinContent(i+1));
   TGraph * EffPreMCD_R[6];
   EffPreMCP_R->SetMarkerColor(2);
   EffPreMCP_R->SetMarkerStyle(8);
   EffPreMCP_R->SetLineColor(2);
   EffPreMCP_R->SetLineWidth(2);
   EffPreMCP_R->SetTitle("Preselections Efficiency MC (R bins)");
   EffPreMCP_R->GetXaxis()->SetTitle("R [GV]");
   EffPreMCP_R->GetYaxis()->SetTitle("Pres. Efficiency");
   EffPreMCP_R->GetXaxis()->SetTitleSize(0.045);
   EffPreMCP_R->GetYaxis()->SetTitleSize(0.045);
   {
      EffPreMCP_R->Draw("ACP");
      TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
      leg->AddEntry(EffPreMCP_R,MCLegend[0].c_str(), "ep");

      for(int h=0; h<6; h++) {
         EffPreMCD_R[h]= new TGraph();
         EffPreMCD_R[h]->SetTitle(MCLegend[h+1].c_str());
         for(int i=1; i<nbinsr; i++) EffPreMCD_R[h]->SetPoint(i,R_cent[i],EffPreMCD_R_TH2F->GetBinContent(i+1,h+1));
         leg->AddEntry(EffPreMCD_R[h],MCLegend[h+1].c_str(), "ep");
         EffPreMCD_R[h]->SetMarkerColor(4);
         EffPreMCD_R[h]->SetMarkerStyle(h+3);
         EffPreMCD_R[h]->SetMarkerSize(2);
         EffPreMCD_R[h]->SetLineColor(4);
         EffPreMCD_R[h]->SetLineWidth(2);
         EffPreMCD_R[h]->Draw("Psame");
         leg->Draw();
      }
   }

   c4_bis->cd(1);
   gPad->SetLogx();
   gPad->SetGridx();
   gPad->SetGridy();
   TGraph * EffPreMCP = new TGraph();
   for(int i=0; i<nbinsToF; i++) EffPreMCP->SetPoint(i,Ekincent[i],EffPreMCP_TH1F->GetBinContent(i+1));
   TGraph * EffPreMCD[6];
   EffPreMCP->SetMarkerColor(2);
   EffPreMCP->SetMarkerStyle(8);
   EffPreMCP->SetLineColor(2);
   EffPreMCP->SetLineWidth(2);
   EffPreMCP->SetTitle("Preselections Efficiency MC (Beta bins)");
   EffPreMCP->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
   EffPreMCP->GetYaxis()->SetTitle("Pres. Efficiency");
   EffPreMCP->GetXaxis()->SetTitleSize(0.045);
   EffPreMCP->GetYaxis()->SetTitleSize(0.045);
   {
      EffPreMCP->Draw("ACP");
      TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
      leg->AddEntry(EffPreMCP,MCLegend[0].c_str(), "ep");

      for(int h=0; h<6; h++) {
         EffPreMCD[h]= new TGraph();
         for(int i=0; i<nbinsToF; i++) EffPreMCD[h]->SetPoint(i,Ekincent[i], EffPreMCD_TH2F->GetBinContent(i+1,h+1));
         EffPreMCD[h]->SetMarkerColor(4);
         EffPreMCD[h]->SetMarkerStyle(h+3);
         leg->AddEntry(EffPreMCD[h],MCLegend[h+1].c_str(), "ep");
         EffPreMCD[h]->SetMarkerSize(2);
         EffPreMCD[h]->SetLineColor(4);
         EffPreMCD[h]->SetLineWidth(2);
         EffPreMCD[h]->Draw("Psame");
         leg->Draw();
      }
   }

   c4_bis->cd(2);
   gPad->SetLogx();
   gPad->SetGridx();
   gPad->SetGridy();
   TGraph * EffPreMCPNaF = new TGraph();
   for(int i=0; i<nbinsNaF; i++) EffPreMCPNaF->SetPoint(i,EkincentNaF[i],EffPreMCPNaF_TH1F->GetBinContent(i+1));
   TGraph * EffPreMCDNaF[6];
   EffPreMCPNaF->SetMarkerColor(2);
   EffPreMCPNaF->SetMarkerStyle(8);
   EffPreMCPNaF->SetLineColor(2);
   EffPreMCPNaF->SetLineWidth(2);
   EffPreMCPNaF->SetTitle("Preselections Efficiency MC (Beta bins NaF)");
   EffPreMCPNaF->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
   EffPreMCPNaF->GetYaxis()->SetTitle("Pres. Efficiency");
   EffPreMCPNaF->GetXaxis()->SetTitleSize(0.045);
   EffPreMCPNaF->GetYaxis()->SetTitleSize(0.045);
   {
      EffPreMCPNaF->Draw("ACP");
      TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
      leg->AddEntry(EffPreMCPNaF,MCLegend[0].c_str(), "ep");

      for(int h=0; h<6; h++) {
         EffPreMCDNaF[h]= new TGraph();
         for(int i=0; i<nbinsNaF; i++) EffPreMCDNaF[h]->SetPoint(i,EkincentNaF[i],
                  EffPreMCDNaF_TH2F->GetBinContent(i+1,h+1));
         EffPreMCDNaF[h]->SetMarkerColor(4);
         EffPreMCDNaF[h]->SetMarkerStyle(h+3);
         leg->AddEntry(EffPreMCDNaF[h],MCLegend[h+1].c_str(), "ep");
         EffPreMCDNaF[h]->SetMarkerSize(2);
         EffPreMCDNaF[h]->SetLineColor(4);
         EffPreMCDNaF[h]->SetLineWidth(2);
         EffPreMCDNaF[h]->Draw("Psame");
         leg->Draw();
      }
   }

   c4_bis->cd(3);
   gPad->SetLogx();
   gPad->SetGridx();
   gPad->SetGridy();
   TGraph * EffPreMCPAgl = new TGraph();
   for(int i=0; i<nbinsAgl; i++) EffPreMCPAgl->SetPoint(i,EkincentAgl[i],EffPreMCPAgl_TH1F->GetBinContent(i+1));
   TGraph * EffPreMCDAgl[6];
   EffPreMCPAgl->SetMarkerColor(2);
   EffPreMCPAgl->SetMarkerStyle(8);
   EffPreMCPAgl->SetLineColor(2);
   EffPreMCPAgl->SetLineWidth(2);
   EffPreMCPAgl->SetTitle("Preselections Efficiency MC (Beta bins Agl)");
   EffPreMCPAgl->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
   EffPreMCPAgl->GetYaxis()->SetTitle("Pres. Efficiency");
   EffPreMCPAgl->GetXaxis()->SetTitleSize(0.045);
   EffPreMCPAgl->GetYaxis()->SetTitleSize(0.045);
   {
      EffPreMCPAgl->Draw("ACP");
      TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
      leg->AddEntry(EffPreMCPAgl,MCLegend[0].c_str(), "ep");

      for(int h=0; h<6; h++) {
         EffPreMCDAgl[h]= new TGraph();
         for(int i=0; i<nbinsAgl; i++) EffPreMCDAgl[h]->SetPoint(i,EkincentAgl[i], EffPreMCDAgl_TH2F->GetBinContent(i+1,h+1));
         EffPreMCDAgl[h]->SetMarkerColor(4);
         EffPreMCDAgl[h]->SetMarkerStyle(h+3);
         leg->AddEntry(EffPreMCDAgl[h],MCLegend[h+1].c_str(), "ep");
         EffPreMCDAgl[h]->SetMarkerSize(2);
         EffPreMCDAgl[h]->SetLineColor(4);
         EffPreMCDAgl[h]->SetLineWidth(2);
         EffPreMCDAgl[h]->Draw("Psame");
         leg->Draw();
      }
   }


   cout<<"*** Updating Results file ***"<<endl;
   nomefile="./Final_plots/"+mese+".root";
   TFile *f_out=new TFile(nomefile.c_str(), "UPDATE");
   f_out->mkdir("MC Results/Preselections");
   f_out->cd("MC Results/Preselections");
   c4->Write();
   c4_bis->Write();
   f_out->Close();

}

