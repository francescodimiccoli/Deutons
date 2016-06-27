
void 	MCTrackeff_Plot(TH1* EffTriggerMCP_R_TH1F,
                        TH1* EffTriggerMCP_TH1F 	,
                        TH1* EffTriggerMCD_R_TH2F,
                        TH1* EffTriggerMCD_TH2F,	
                        TH1* EffTrackerMCP_R_TH1F,
                        TH1* EffTrackerMCP_TH1F,	
                        TH1* EffTrackerMCD_R_TH2F,
                        TH1* EffTrackerMCD_TH2F,	
                        TH1* EffTOF_MCP_R_TH1F,	
                        TH1* EffTOF_MCP_TH1F	,
                        TH1* EffTOF_MCD_R_TH2F,	
                        TH1* EffTOF_MCD_TH2F	
	){

   TCanvas *c_7 	=new TCanvas ("Trigger sel. efficiency");
   TCanvas *c7	=new TCanvas ("Tracker rec. efficiency");
   TCanvas *c8	=new TCanvas ("TOF rec. efficiency");
   c_7->Divide (2,1);
   c7->Divide (2,1);
   c8->Divide (2,1);
   c_7->cd (1);
   gPad->SetLogx();
   gPad->SetGridx();
   gPad->SetGridy();
   string MCLegend[7]= {"protons.B800","d.pl1.0_520_GG_Blic","d.pl1.0_520_GG_BlicDPMJet","d.pl1.0_520_GG_QMD","d.pl1.0_520_Shen_Blic","d.pl1.0_520_Shen_BlicDPMJet","d.pl1.0_520_Shen_QMD"};
   TGraph * EffTriggerMCP_R = new TGraph();
   EffTriggerMCP_R->SetTitle (MCLegend[0].c_str() );
   for (int i=0; i<nbinsr; i++) EffTriggerMCP_R->SetPoint (i,PRB.RigBinCent (i),EffTriggerMCP_R_TH1F->GetBinContent (i+1) );
   TGraph * EffTriggerMCD_R[6];
   EffTriggerMCP_R->SetMarkerColor (2);
   EffTriggerMCP_R->SetMarkerStyle (8);
   EffTriggerMCP_R->SetLineColor (2);
   EffTriggerMCP_R->SetLineWidth (2);
   EffTriggerMCP_R->SetTitle ("Trigger rec. Efficiency MC (R bins)");
   EffTriggerMCP_R->GetXaxis()->SetTitle ("R [GV]");
   EffTriggerMCP_R->GetYaxis()->SetTitle ("Pres. Efficiency");
   EffTriggerMCP_R->GetXaxis()->SetTitleSize (0.045);
   EffTriggerMCP_R->GetYaxis()->SetTitleSize (0.045);
   {
      EffTriggerMCP_R->Draw ("ACP");
      TLegend* leg =new TLegend (0.4, 0.7,0.95,0.95);
      leg->AddEntry (EffTriggerMCP_R,MCLegend[0].c_str(), "ep");

      for (int h=0; h<6; h++) {
         EffTriggerMCD_R[h]= new TGraph();
         EffTriggerMCD_R[h]->SetTitle (MCLegend[h+1].c_str() );
         for (int i=1; i<nbinsr; i++) EffTriggerMCD_R[h]->SetPoint (i,PRB.RigBinCent (i),EffTriggerMCD_R_TH2F->GetBinContent (i+1,h+1) );
         leg->AddEntry (EffTriggerMCD_R[h],MCLegend[h+1].c_str(), "ep");
         EffTriggerMCD_R[h]->SetMarkerColor (4);
         EffTriggerMCD_R[h]->SetMarkerStyle (h+3);
         EffTriggerMCD_R[h]->SetMarkerSize (2);
         EffTriggerMCD_R[h]->SetLineColor (4);
         EffTriggerMCD_R[h]->SetLineWidth (2);
         EffTriggerMCD_R[h]->Draw ("Psame");
         leg->Draw();
      }
   }
   c_7->cd (2);
   gPad->SetLogx();
   gPad->SetGridx();
   gPad->SetGridy();
   TGraph * EffTriggerMCP = new TGraph();
   for (int i=0; i<17; i++) EffTriggerMCP->SetPoint (i,ToFPB.EkBinCent (i),EffTriggerMCP_TH1F->GetBinContent (i+1) );
   TGraph * EffTriggerMCD[6];
   EffTriggerMCP->SetMarkerColor (2);
   EffTriggerMCP->SetMarkerStyle (8);
   EffTriggerMCP->SetLineColor (2);
   EffTriggerMCP->SetLineWidth (2);
   EffTriggerMCP->SetTitle ("Trigger rec. Efficiency MC (Beta bins)");
   EffTriggerMCP->GetXaxis()->SetTitle ("Kin. En. / nucl. [GeV/nucl.]");
   EffTriggerMCP->GetYaxis()->SetTitle ("Trigger rec. Efficiency");
   EffTriggerMCP->GetXaxis()->SetTitleSize (0.045);
   EffTriggerMCP->GetYaxis()->SetTitleSize (0.045);
   {
      EffTriggerMCP->Draw ("ACP");
      TLegend* leg =new TLegend (0.4, 0.7,0.95,0.95);
      leg->AddEntry (EffTriggerMCP,MCLegend[0].c_str(), "ep");

      for (int h=0; h<6; h++) {
         EffTriggerMCD[h]= new TGraph();
         for (int i=0; i<17; i++) EffTriggerMCD[h]->SetPoint (i,ToFPB.EkBinCent (i),EffTriggerMCD_TH2F->GetBinContent (i+1,h+1) );
         EffTriggerMCD[h]->SetMarkerColor (4);
         EffTriggerMCD[h]->SetMarkerStyle (h+3);
         leg->AddEntry (EffTriggerMCD[h],MCLegend[h+1].c_str(), "ep");
         EffTriggerMCD[h]->SetMarkerSize (2);
         EffTriggerMCD[h]->SetLineColor (4);
         EffTriggerMCD[h]->SetLineWidth (2);
         EffTriggerMCD[h]->Draw ("Psame");
         leg->Draw();
      }
   }

   c7->cd (1);
   gPad->SetLogx();
   gPad->SetGridx();
   gPad->SetGridy();
   TGraph * EffTrackerMCP_R = new TGraph();
   EffTrackerMCP_R->SetTitle (MCLegend[0].c_str() );
   for (int i=0; i<nbinsr; i++) EffTrackerMCP_R->SetPoint (i,PRB.RigBinCent (i),EffTrackerMCP_R_TH1F->GetBinContent (i+1) );
   TGraph * EffTrackerMCD_R[6];
   EffTrackerMCP_R->SetMarkerColor (2);
   EffTrackerMCP_R->SetMarkerStyle (8);
   EffTrackerMCP_R->SetLineColor (2);
   EffTrackerMCP_R->SetLineWidth (2);
   EffTrackerMCP_R->SetTitle ("Tracker rec. Efficiency MC (R bins)");
   EffTrackerMCP_R->GetXaxis()->SetTitle ("R [GV]");
   EffTrackerMCP_R->GetYaxis()->SetTitle ("Pres. Efficiency");
   EffTrackerMCP_R->GetXaxis()->SetTitleSize (0.045);
   EffTrackerMCP_R->GetYaxis()->SetTitleSize (0.045);
   {
      EffTrackerMCP_R->Draw ("ACP");
      TLegend* leg =new TLegend (0.4, 0.7,0.95,0.95);
      leg->AddEntry (EffTrackerMCP_R,MCLegend[0].c_str(), "ep");

      for (int h=0; h<6; h++) {
         EffTrackerMCD_R[h]= new TGraph();
         EffTrackerMCD_R[h]->SetTitle (MCLegend[h+1].c_str() );
         for (int i=1; i<nbinsr; i++) EffTrackerMCD_R[h]->SetPoint (i,PRB.RigBinCent (i),EffTrackerMCD_R_TH2F->GetBinContent (i+1,h+1) );
         leg->AddEntry (EffTrackerMCD_R[h],MCLegend[h+1].c_str(), "ep");
         EffTrackerMCD_R[h]->SetMarkerColor (4);
         EffTrackerMCD_R[h]->SetMarkerStyle (h+3);
         EffTrackerMCD_R[h]->SetMarkerSize (2);
         EffTrackerMCD_R[h]->SetLineColor (4);
         EffTrackerMCD_R[h]->SetLineWidth (2);
         EffTrackerMCD_R[h]->Draw ("Psame");
         leg->Draw();
      }
   }
   c7->cd (2);
   gPad->SetLogx();
   gPad->SetGridx();
   gPad->SetGridy();
   TGraph * EffTrackerMCP = new TGraph();
   for (int i=0; i<17; i++) EffTrackerMCP->SetPoint (i,ToFPB.EkBinCent (i),EffTrackerMCP_TH1F->GetBinContent (i+1) );
   TGraph * EffTrackerMCD[6];
   EffTrackerMCP->SetMarkerColor (2);
   EffTrackerMCP->SetMarkerStyle (8);
   EffTrackerMCP->SetLineColor (2);
   EffTrackerMCP->SetLineWidth (2);
   EffTrackerMCP->SetTitle ("Tracker rec. Efficiency MC (Beta bins)");
   EffTrackerMCP->GetXaxis()->SetTitle ("Kin. En. / nucl. [GeV/nucl.]");
   EffTrackerMCP->GetYaxis()->SetTitle ("Tracker rec. Efficiency");
   EffTrackerMCP->GetXaxis()->SetTitleSize (0.045);
   EffTrackerMCP->GetYaxis()->SetTitleSize (0.045);
   {
      EffTrackerMCP->Draw ("ACP");
      TLegend* leg =new TLegend (0.4, 0.7,0.95,0.95);
      leg->AddEntry (EffTrackerMCP,MCLegend[0].c_str(), "ep");

      for (int h=0; h<6; h++) {
         EffTrackerMCD[h]= new TGraph();
         for (int i=0; i<17; i++) EffTrackerMCD[h]->SetPoint (i,ToFPB.EkBinCent (i),EffTrackerMCD_TH2F->GetBinContent (i+1,h+1) );
         EffTrackerMCD[h]->SetMarkerColor (4);
         EffTrackerMCD[h]->SetMarkerStyle (h+3);
         leg->AddEntry (EffTrackerMCD[h],MCLegend[h+1].c_str(), "ep");
         EffTrackerMCD[h]->SetMarkerSize (2);
         EffTrackerMCD[h]->SetLineColor (4);
         EffTrackerMCD[h]->SetLineWidth (2);
         EffTrackerMCD[h]->Draw ("Psame");
         leg->Draw();
      }
   }

   c8->cd (1);
   gPad->SetLogx();
   gPad->SetGridx();
   gPad->SetGridy();
   TGraph * EffTOF_MCP_R = new TGraph();
   EffTOF_MCP_R->SetTitle (MCLegend[0].c_str() );
   for (int i=0; i<nbinsr; i++) EffTOF_MCP_R->SetPoint (i,PRB.RigBinCent (i),EffTOF_MCP_R_TH1F->GetBinContent (i+1) );
   TGraph * EffTOF_MCD_R[6];
   EffTOF_MCP_R->SetMarkerColor (2);
   EffTOF_MCP_R->SetMarkerStyle (8);
   EffTOF_MCP_R->SetLineColor (2);
   EffTOF_MCP_R->SetLineWidth (2);
   EffTOF_MCP_R->SetTitle ("TOF rec. Efficiency MC (R bins)");
   EffTOF_MCP_R->GetXaxis()->SetTitle ("R [GV]");
   EffTOF_MCP_R->GetYaxis()->SetTitle ("Pres. Efficiency");
   EffTOF_MCP_R->GetXaxis()->SetTitleSize (0.045);
   EffTOF_MCP_R->GetYaxis()->SetTitleSize (0.045);
   {
      EffTOF_MCP_R->Draw ("ACP");
      TLegend* leg =new TLegend (0.4, 0.7,0.95,0.95);
      leg->AddEntry (EffTOF_MCP_R,MCLegend[0].c_str(), "ep");

      for (int h=0; h<6; h++) {
         EffTOF_MCD_R[h]= new TGraph();
         EffTOF_MCD_R[h]->SetTitle (MCLegend[h+1].c_str() );
         for (int i=1; i<nbinsr; i++) EffTOF_MCD_R[h]->SetPoint (i,PRB.RigBinCent (i),EffTOF_MCD_R_TH2F->GetBinContent (i+1,h+1) );
         leg->AddEntry (EffTOF_MCD_R[h],MCLegend[h+1].c_str(), "ep");
         EffTOF_MCD_R[h]->SetMarkerColor (4);
         EffTOF_MCD_R[h]->SetMarkerStyle (h+3);
         EffTOF_MCD_R[h]->SetMarkerSize (2);
         EffTOF_MCD_R[h]->SetLineColor (4);
         EffTOF_MCD_R[h]->SetLineWidth (2);
         EffTOF_MCD_R[h]->Draw ("Psame");
         leg->Draw();
      }
   }

   c8->cd (2);
   gPad->SetLogx();
   gPad->SetGridx();
   gPad->SetGridy();
   TGraph * EffTOF_MCP = new TGraph();
   for (int i=0; i<17; i++) EffTOF_MCP->SetPoint (i,ToFPB.EkBinCent (i),EffTOF_MCP_TH1F->GetBinContent (i+1) );
   TGraph * EffTOF_MCD[6];
   EffTOF_MCP->SetMarkerColor (2);
   EffTOF_MCP->SetMarkerStyle (8);
   EffTOF_MCP->SetLineColor (2);
   EffTOF_MCP->SetLineWidth (2);
   EffTOF_MCP->SetTitle ("TOF_ rec. Efficiency MC (Beta bins)");
   EffTOF_MCP->GetXaxis()->SetTitle ("Kin. En. / nucl. [GeV/nucl.]");
   EffTOF_MCP->GetYaxis()->SetTitle ("TOF_ rec. Efficiency");
   EffTOF_MCP->GetXaxis()->SetTitleSize (0.045);
   EffTOF_MCP->GetYaxis()->SetTitleSize (0.045);
   {
      EffTOF_MCP->Draw ("ACP");
      TLegend* leg =new TLegend (0.4, 0.7,0.95,0.95);
      leg->AddEntry (EffTOF_MCP,MCLegend[0].c_str(), "ep");

      for (int h=0; h<6; h++) {
         EffTOF_MCD[h]= new TGraph();
         for (int i=0; i<17; i++) EffTOF_MCD[h]->SetPoint (i,ToFPB.EkBinCent (i),EffTOF_MCD_TH2F->GetBinContent (i+1,h+1) );
         EffTOF_MCD[h]->SetMarkerColor (4);
         EffTOF_MCD[h]->SetMarkerStyle (h+3);
         leg->AddEntry (EffTOF_MCD[h],MCLegend[h+1].c_str(), "ep");
         EffTOF_MCD[h]->SetMarkerSize (2);
         EffTOF_MCD[h]->SetLineColor (4);
         EffTOF_MCD[h]->SetLineWidth (2);
         EffTOF_MCD[h]->Draw ("Psame");
         leg->Draw();
      }
   }

   finalPlots.Add(c_7);
   finalPlots.Add(c7);
   finalPlots.Add(c8);
   finalPlots.writeObjsInFolder("MC Results/Preselections/Basic Selections");

}
