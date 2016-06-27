
void SlidesforPlot (string filename)
{
    
   cout<<"************** Plots for slides *******************"<<endl;

        cout<<"*** Reading  P1 file ****"<<endl;
        TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");

   TH2F * RvsBetaTOF_P= (TH2F*) inputHistoFile->Get ("RvsBetaTOF_P");
   TH2F * RvsBetaNaF_P= (TH2F*) inputHistoFile->Get ("RvsBetaNaF_P");
   TH2F * RvsBetaAgl_P= (TH2F*) inputHistoFile->Get ("RvsBetaAgl_P");
   TH2F * RvsBetaTOF_D= (TH2F*) inputHistoFile->Get ("RvsBetaTOF_D");
   TH2F * RvsBetaNaF_D= (TH2F*) inputHistoFile->Get ("RvsBetaNaF_D");
   TH2F * RvsBetaAgl_D= (TH2F*) inputHistoFile->Get ("RvsBetaAgl_D");
   TH2F * RvsBetaTOF_He= (TH2F*) inputHistoFile->Get ("RvsBetaTOF_He");
   TH2F * RvsBetaNaF_He= (TH2F*) inputHistoFile->Get ("RvsBetaNaF_He");
   TH2F * RvsBetaAgl_He= (TH2F*) inputHistoFile->Get ("RvsBetaAgl_He");
   TH2F * EdepUTOFvsR_P= (TH2F*) inputHistoFile->Get ("EdepUTOFvsR_P");
   TH2F * EdepUTOFvsR_D= (TH2F*) inputHistoFile->Get ("EdepUTOFvsR_D");
   TH2F * EdepUTOFvsR_He= (TH2F*) inputHistoFile->Get ("EdepUTOFvsR_He");
   TH2F * EdepLTOFvsR_P= (TH2F*) inputHistoFile->Get ("EdepLTOFvsR_P");
   TH2F * EdepLTOFvsR_D= (TH2F*) inputHistoFile->Get ("EdepLTOFvsR_D");
   TH2F * EdepLTOFvsR_He= (TH2F*) inputHistoFile->Get ("EdepLTOFvsR_He");
   TH2F * EdepTrackvsR_P= (TH2F*) inputHistoFile->Get ("EdepTrackvsR_P");
   TH2F * EdepTrackvsR_D= (TH2F*) inputHistoFile->Get ("EdepTrackvsR_D");
   TH2F * EdepTrackvsR_He= (TH2F*) inputHistoFile->Get ("EdepTrackvsR_He");
   TH1F * MassTOF_P= (TH1F*) inputHistoFile->Get ("MassTOF_P");
   TH1F * MassTOF_D= (TH1F*) inputHistoFile->Get ("MassTOF_D");
   TH1F * MassNaF_P= (TH1F*) inputHistoFile->Get ("MassNaF_P");
   TH1F * MassNaF_D= (TH1F*) inputHistoFile->Get ("MassNaF_D");
   TH1F * MassAgl_P= (TH1F*) inputHistoFile->Get ("MassAgl_P");
   TH1F * MassAgl_D= (TH1F*) inputHistoFile->Get ("MassAgl_D");
   TH1F * MassTOF_PQ= (TH1F*) inputHistoFile->Get ("MassTOF_PQ");
   TH1F * MassTOF_DQ= (TH1F*) inputHistoFile->Get ("MassTOF_DQ");
   TH1F * MassNaF_PQ= (TH1F*) inputHistoFile->Get ("MassNaF_PQ");
   TH1F * MassNaF_DQ= (TH1F*) inputHistoFile->Get ("MassNaF_DQ");
   TH1F * MassAgl_PQ= (TH1F*) inputHistoFile->Get ("MassAgl_PQ");
   TH1F * MassAgl_DQ= (TH1F*) inputHistoFile->Get ("MassAgl_DQ");
   TH2F * RvsBetaTOF= (TH2F*) inputHistoFile->Get ("RvsBetaTOF");
   TH2F * RvsBetaNaF= (TH2F*) inputHistoFile->Get ("RvsBetaNaF");
   TH2F * RvsBetaAgl= (TH2F*) inputHistoFile->Get ("RvsBetaAgl");
   TH1F * MassTOF= (TH1F*) inputHistoFile->Get ("MassTOF");
   TH1F * MassNaF= (TH1F*) inputHistoFile->Get ("MassNaF");
   TH1F * MassAgl= (TH1F*) inputHistoFile->Get ("MassAgl");
   TH1F * MassTOFQ= (TH1F*) inputHistoFile->Get ("MassTOFQ");
   TH1F * MassNaFQ= (TH1F*) inputHistoFile->Get ("MassNaFQ");
   TH1F * MassAglQ= (TH1F*) inputHistoFile->Get ("MassAglQ");
   TH1F * DistTOF_P= (TH1F*) inputHistoFile->Get ("DistTOF_P");
   TH1F * DistNaF_P= (TH1F*) inputHistoFile->Get ("DistNaF_P");
   TH1F * DistAgl_P= (TH1F*) inputHistoFile->Get ("DistAgl_P");
   TH1F * DistTOF_D= (TH1F*) inputHistoFile->Get ("DistTOF_D");
   TH1F * DistNaF_D= (TH1F*) inputHistoFile->Get ("DistNaF_D");
   TH1F * DistAgl_D= (TH1F*) inputHistoFile->Get ("DistAgl_D");
   TH1F * DistTOF_He= (TH1F*) inputHistoFile->Get ("DistTOF_He");
   TH1F * DistNaF_He= (TH1F*) inputHistoFile->Get ("DistNaF_He");
   TH1F * DistAgl_He= (TH1F*) inputHistoFile->Get ("DistAgl_He");
   TH2F * RvsDistTOF_P= (TH2F*) inputHistoFile->Get ("RvsDistTOF_P");
   TH2F * RvsDistNaF_P= (TH2F*) inputHistoFile->Get ("RvsDistNaF_P");
   TH2F * RvsDistAgl_P= (TH2F*) inputHistoFile->Get ("RvsDistAgl_P");
   TH2F * RvsDistTOF_D= (TH2F*) inputHistoFile->Get ("RvsDistTOF_D");
   TH2F * RvsDistNaF_D= (TH2F*) inputHistoFile->Get ("RvsDistNaF_D");
   TH2F * RvsDistAgl_D= (TH2F*) inputHistoFile->Get ("RvsDistAgl_D");

   TCanvas *p1 =new TCanvas ("RvsBeta TOF MC");
   TCanvas *p2 =new TCanvas ("RvsBeta NaF MC");
   TCanvas *p3 =new TCanvas ("RvsBeta Agl MC");
   TCanvas *p4 =new TCanvas ("RvsBeta TOF H.L. data");
   TCanvas *p5 =new TCanvas ("RvsBeta NaF H.L. data");
   TCanvas *p6 =new TCanvas ("RvsBeta Agl H.L. data");
   TCanvas *p7 =new TCanvas ("RvsEdep Upper TOF MC");
   TCanvas *p8 =new TCanvas ("RvsEdep Lower TOF MC");
   TCanvas *p9 =new TCanvas ("RvsEdep Tracker MC");
   TCanvas *p10=new TCanvas ("Mass TOF H.L. data");
   TCanvas *p11=new TCanvas ("Mass NaF H.L. data");
   TCanvas *p12=new TCanvas ("Mass Agl H.L. data");
   TCanvas *p13=new TCanvas ("Mass TOF MC");
   TCanvas *p14=new TCanvas ("Mass NaF MC");
   TCanvas *p15=new TCanvas ("Mass Agl MC");
   TCanvas *p10Q=new TCanvas ("Qual. Mass TOF H.L. data");
   TCanvas *p11Q=new TCanvas ("Qual. Mass NaF H.L. data");
   TCanvas *p12Q=new TCanvas ("Qual. Mass Agl H.L. data");
   TCanvas *p13Q=new TCanvas ("Qual. Mass TOF MC");
   TCanvas *p14Q=new TCanvas ("Qual. Mass NaF MC");
   TCanvas *p15Q=new TCanvas ("Qual. Mass Agl MC");
   TCanvas *p16=new TCanvas ("Distance discr. TOF MC");
   TCanvas *p17=new TCanvas ("Distance discr. NaF MC");
   TCanvas *p18=new TCanvas ("Distance discr. Agl MC");
   TCanvas *p19=new TCanvas ("DistvsR  TOF MC");
   TCanvas *p20=new TCanvas ("DistvsR  NaF MC");
   TCanvas *p21=new TCanvas ("DistvsR  Agl MC");



   cout<<"******************* R vs Beta plots ******************"<<endl;

   {
      p1->cd();
      gPad->SetGridx();
      gPad->SetGridy();
      RvsBetaTOF_P->SetMarkerColor (2);
      RvsBetaTOF_D->SetMarkerColor (4);
      RvsBetaTOF_He->SetMarkerColor (3);
      RvsBetaTOF_P->SetFillColor (2);
      RvsBetaTOF_D->SetFillColor (4);
      RvsBetaTOF_He->SetFillColor (3);

      RvsBetaTOF_D->SetTitle ("R vs Beta TOF (MC)");
      RvsBetaTOF_D->GetXaxis()->SetTitle ("R [GV]");
      RvsBetaTOF_D->GetYaxis()->SetTitle ("Beta TOF");
      RvsBetaTOF_He->GetZaxis()->SetRangeUser (1,100);
      RvsBetaTOF_D->Draw();
      RvsBetaTOF_P->Draw ("same");
      RvsBetaTOF_He->Draw ("same");
      TLegend* leg =new TLegend (0.4, 0.7,0.95,0.95);
      leg->AddEntry (RvsBetaTOF_P,"Protons MC", "ep");
      leg->AddEntry (RvsBetaTOF_D,"Deutons MC", "ep");
      leg->AddEntry (RvsBetaTOF_He,"Q=1 from He fragm.", "ep");
      leg->Draw ("same");
   }

   {
      p2->cd();
      gPad->SetGridx();
      gPad->SetGridy();
      RvsBetaNaF_P->SetMarkerColor (2);
      RvsBetaNaF_D->SetMarkerColor (4);
      RvsBetaNaF_He->SetMarkerColor (3);
      RvsBetaNaF_P->SetFillColor (2);
      RvsBetaNaF_D->SetFillColor (4);
      RvsBetaNaF_He->SetFillColor (3);
      RvsBetaNaF_D->SetTitle ("R vs Beta NaF (MC)");
      RvsBetaNaF_D->GetXaxis()->SetTitle ("R [GV]");
      RvsBetaNaF_D->GetYaxis()->SetTitle ("Beta NaF");
      RvsBetaNaF_D->Draw();
      RvsBetaNaF_P->Draw ("same");
      RvsBetaNaF_He->Draw ("same");
      TLegend* leg =new TLegend (0.4, 0.7,0.95,0.95);
      leg->AddEntry (RvsBetaTOF_P,"Protons MC", "ep");
      leg->AddEntry (RvsBetaTOF_D,"Deutons MC", "ep");
      leg->AddEntry (RvsBetaTOF_He,"Q=1 from He fragm.", "ep");
      leg->Draw ("same");
   }

   {
      p3->cd();
      gPad->SetGridx();
      gPad->SetGridy();
      RvsBetaAgl_P->SetMarkerColor (2);
      RvsBetaAgl_D->SetMarkerColor (4);
      RvsBetaAgl_He->SetMarkerColor (3);
      RvsBetaAgl_P->SetFillColor (2);
      RvsBetaAgl_D->SetFillColor (4);
      RvsBetaAgl_He->SetFillColor (3);
      RvsBetaAgl_D->SetTitle ("R vs Beta Agl (MC)");
      RvsBetaAgl_D->GetXaxis()->SetTitle ("R [GV]");
      RvsBetaAgl_D->GetYaxis()->SetTitle ("Beta Agl");
      RvsBetaAgl_D->Draw();
      RvsBetaAgl_P->Draw ("same");
      RvsBetaAgl_He->Draw ("same");
      TLegend* leg =new TLegend (0.4, 0.7,0.95,0.95);
      leg->AddEntry (RvsBetaTOF_P,"Protons MC", "ep");
      leg->AddEntry (RvsBetaTOF_D,"Deutons MC", "ep");
      leg->AddEntry (RvsBetaTOF_He,"Q=1 from He fragm.", "ep");
      leg->Draw ("same");
   }

   {
      p4->cd();
      gPad->SetGridx();
      gPad->SetGridy();
      gPad->SetLogz();
      RvsBetaTOF->SetTitle ("R vs Beta TOF (H.L. ISS data)");
      RvsBetaTOF->GetXaxis()->SetTitle ("R [GV]");
      RvsBetaTOF->GetYaxis()->SetTitle ("Beta TOF");
      RvsBetaTOF->Draw ("col");
      protons->SetLineColor (2);
      deutons->SetLineColor (4);
      protons->SetLineWidth (3);
      deutons->SetLineWidth (3);
      protons->Draw ("same");
      deutons->Draw ("same");
      TLegend* leg =new TLegend (0.4, 0.7,0.95,0.95);
      leg->AddEntry (protons,"Protons Mass curve", "ep");
      leg->AddEntry (deutons,"Deutons Mass curve", "ep");
      leg->Draw ("same");
   }

   {
      p5->cd();
      gPad->SetGridx();
      gPad->SetGridy();
      gPad->SetLogz();
      RvsBetaNaF->SetTitle ("R vs Beta NaF (H.L. ISS data)");
      RvsBetaNaF->GetXaxis()->SetTitle ("R [GV]");
      RvsBetaNaF->GetYaxis()->SetTitle ("Beta NaF");
      RvsBetaNaF->Draw ("col");
      protons->SetLineColor (2);
      deutons->SetLineColor (4);
      protons->SetLineWidth (3);
      deutons->SetLineWidth (3);
      protons->Draw ("same");
      deutons->Draw ("same");
      TLegend* leg =new TLegend (0.4, 0.7,0.95,0.95);
      leg->AddEntry (protons,"Protons Mass curve", "ep");
      leg->AddEntry (deutons,"Deutons Mass curve", "ep");
      leg->Draw ("same");
   }

   {
      p6->cd();
      gPad->SetGridx();
      gPad->SetGridy();
      gPad->SetLogz();
      RvsBetaAgl->SetTitle ("R vs Beta Agl (H.L. ISS data)");
      RvsBetaAgl->GetXaxis()->SetTitle ("R [GV]");
      RvsBetaAgl->GetYaxis()->SetTitle ("Beta Agl");
      RvsBetaAgl->Draw ("col");
      protons->SetLineColor (2);
      deutons->SetLineColor (4);
      protons->SetLineWidth (3);
      deutons->SetLineWidth (3);
      protons->Draw ("same");
      deutons->Draw ("same");
      TLegend* leg =new TLegend (0.4, 0.7,0.95,0.95);
      leg->AddEntry (protons,"Protons Mass curve", "ep");
      leg->AddEntry (deutons,"Deutons Mass curve", "ep");
      leg->Draw ("same");
   }

   cout<<"******************* R vs Edep plots ******************"<<endl;
   p7->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   EdepUTOFvsR_P->SetMarkerColor (2);
   EdepUTOFvsR_D->SetMarkerColor (4);
   EdepUTOFvsR_He->SetMarkerColor (3);
   EdepUTOFvsR_D->SetTitle ("R vs E.dep. Upper TOF (MC)");
   EdepUTOFvsR_D->GetXaxis()->SetTitle ("R [GV]");
   EdepUTOFvsR_D->GetYaxis()->SetTitle ("E. dep. Upper TOF [MeV]");
   EdepUTOFvsR_He->GetZaxis()->SetRangeUser (100,10000);
   EdepUTOFvsR_D->GetZaxis()->SetRangeUser (10,1000);
   EdepUTOFvsR_P->GetZaxis()->SetRangeUser (150,10000);
   EdepUTOFvsR_D->GetXaxis()->SetRangeUser (0,4);
   EdepUTOFvsR_D->GetYaxis()->SetRangeUser (0,20);
   EdepUTOFvsR_D->Draw();
   EdepUTOFvsR_P->Draw ("same");
   EdepUTOFvsR_He->Draw ("same");

   p8->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   EdepLTOFvsR_P->SetMarkerColor (2);
   EdepLTOFvsR_D->SetMarkerColor (4);
   EdepLTOFvsR_He->SetMarkerColor (3);
   EdepLTOFvsR_D->SetTitle ("R vs E.dep. Lower TOF (MC)");
   EdepLTOFvsR_D->GetXaxis()->SetTitle ("R [GV]");
   EdepLTOFvsR_D->GetYaxis()->SetTitle ("E. dep. Upper TOF [MeV]");
   EdepLTOFvsR_He->GetZaxis()->SetRangeUser (100,10000);
   EdepLTOFvsR_D->GetZaxis()->SetRangeUser (10,1000);
   EdepLTOFvsR_P->GetZaxis()->SetRangeUser (150,10000);
   EdepLTOFvsR_D->GetXaxis()->SetRangeUser (0,4);
   EdepLTOFvsR_D->GetYaxis()->SetRangeUser (0,20);
   EdepLTOFvsR_D->Draw();
   EdepLTOFvsR_P->Draw ("same");
   EdepLTOFvsR_He->Draw ("same");

   p9->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   EdepTrackvsR_P->SetMarkerColor (2);
   EdepTrackvsR_D->SetMarkerColor (4);
   EdepTrackvsR_He->SetMarkerColor (3);
   EdepTrackvsR_D->SetTitle ("R vs E.dep. Tracker (MC)");
   EdepTrackvsR_D->GetXaxis()->SetTitle ("R [GV]");
   EdepTrackvsR_D->GetYaxis()->SetTitle ("E. dep. Upper TOF [MeV]");
   EdepTrackvsR_He->GetZaxis()->SetRangeUser (100,10000);
   EdepTrackvsR_D->GetZaxis()->SetRangeUser (10,1000);
   EdepTrackvsR_P->GetZaxis()->SetRangeUser (150,10000);
   EdepTrackvsR_D->GetXaxis()->SetRangeUser (0,4);
   EdepTrackvsR_D->GetYaxis()->SetRangeUser (0,1.22);
   EdepTrackvsR_D->Draw();
   EdepTrackvsR_P->Draw ("same");
   EdepTrackvsR_He->Draw ("same");


   cout<<"******************* Mass plots ******************"<<endl;

   p10->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   gPad->SetLogy();
   MassTOF->SetLineColor (1);
   MassTOF->SetLineWidth (2);
   MassTOF->GetXaxis()->SetTitle ("Mass [GeV/c^2]");
   MassTOF->GetXaxis()->SetTitleSize (0.045);
   MassTOF->SetTitle ("Mass TOF H.L. data");
   MassTOF->Draw();

   p11->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   gPad->SetLogy();
   MassNaF->SetLineColor (1);
   MassNaF->SetLineWidth (2);
   MassNaF->GetXaxis()->SetTitle ("Mass [GeV/c^2]");
   MassNaF->SetTitle ("Mass NaF H.L. data");
   MassNaF->GetXaxis()->SetTitleSize (0.045);
   MassNaF->Draw();

   p12->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   gPad->SetLogy();
   MassAgl->SetLineColor (1);
   MassAgl->SetLineWidth (2);
   MassAgl->GetXaxis()->SetTitle ("Mass [GeV/c^2]");
   MassAgl->SetTitle ("Mass Agl H.L. data");
   MassAgl->GetXaxis()->SetTitleSize (0.045);
   MassAgl->Draw();

   p13->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   gPad->SetLogy();
   MassTOF_P->SetLineColor (2);
   MassTOF_P->SetLineWidth (2);
   MassTOF_D->SetLineColor (4);
   MassTOF_D->SetLineWidth (2);
   MassTOF_D->SetFillColor (4);
   MassTOF_P->SetFillColor (2);
   MassTOF_P->SetFillStyle (3001);
   MassTOF_D->SetFillStyle (3002);
   MassTOF_P->GetXaxis()->SetTitle ("Mass [GeV/c^2]");
   MassTOF_P->GetXaxis()->SetTitleSize (0.045);
   MassTOF_P->SetTitle ("Mass TOF MC");
   MassTOF_P->Draw();
   MassTOF_D->Draw ("same");

   p14->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   gPad->SetLogy();
   MassNaF_P->SetLineColor (2);
   MassNaF_P->SetLineWidth (2);
   MassNaF_D->SetLineColor (4);
   MassNaF_D->SetLineWidth (2);
   MassNaF_D->SetFillColor (4);
   MassNaF_P->SetFillColor (2);
   MassNaF_P->SetFillStyle (3001);
   MassNaF_D->SetFillStyle (3002);
   MassNaF_P->GetXaxis()->SetTitle ("Mass [GeV/c^2]");
   MassNaF_P->SetTitle ("Mass NaF MC");
   MassNaF_P->GetXaxis()->SetTitleSize (0.045);
   MassNaF_P->Draw();
   MassNaF_D->Draw ("same");

   p15->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   gPad->SetLogy();
   MassAgl_P->SetLineColor (2);
   MassAgl_P->SetLineWidth (2);
   MassAgl_D->SetLineColor (4);
   MassAgl_D->SetLineWidth (2);
   MassAgl_D->SetFillColor (4);
   MassAgl_P->SetFillColor (2);
   MassAgl_P->SetFillStyle (3001);
   MassAgl_D->SetFillStyle (3002);
   MassAgl_P->GetXaxis()->SetTitle ("Mass [GeV/c^2]");
   MassAgl_P->SetTitle ("Mass Agl MC");
   MassAgl_P->GetXaxis()->SetTitleSize (0.045);
   MassAgl_P->Draw();
   MassAgl_D->Draw ("same");


   cout<<"******************* Quality Mass plots ******************"<<endl;

   p10Q->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   gPad->SetLogy();
   MassTOFQ->SetLineColor (1);
   MassTOFQ->SetLineWidth (2);
   MassTOFQ->GetXaxis()->SetTitle ("Mass [GeV/c^2]");
   MassTOFQ->SetTitle ("Mass TOF H.L. data");
   MassTOFQ->GetXaxis()->SetTitleSize (0.045);
   MassTOFQ->Draw();

   p11Q->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   gPad->SetLogy();
   MassNaFQ->SetLineColor (1);
   MassNaFQ->SetLineWidth (2);
   MassNaFQ->GetXaxis()->SetTitle ("Mass [GeV/c^2]");
   MassNaFQ->SetTitle ("Mass NaF H.L. data");
   MassNaFQ->GetXaxis()->SetTitleSize (0.045);
   MassNaFQ->Draw();

   p12Q->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   gPad->SetLogy();
   MassAglQ->SetLineColor (1);
   MassAglQ->SetLineWidth (2);
   MassAglQ->GetXaxis()->SetTitle ("Mass [GeV/c^2]");
   MassAglQ->SetTitle ("Mass Agl H.L. data");
   MassAglQ->GetXaxis()->SetTitleSize (0.045);
   MassAglQ->Draw();

   p13Q->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   gPad->SetLogy();
   MassTOF_PQ->SetLineColor (2);
   MassTOF_PQ->SetLineWidth (2);
   MassTOF_DQ->SetLineColor (4);
   MassTOF_DQ->SetLineWidth (2);
   MassTOF_DQ->SetFillColor (4);
   MassTOF_PQ->SetFillColor (2);
   MassTOF_PQ->SetFillStyle (3001);
   MassTOF_DQ->SetFillStyle (3002);
   MassTOF_PQ->GetXaxis()->SetTitle ("Mass [GeV/c^2]");
   MassTOF_PQ->GetXaxis()->SetTitleSize (0.045);
   MassTOF_PQ->SetTitle ("Mass TOF MC");
   MassTOF_PQ->Draw();
   MassTOF_DQ->Draw ("same");

   p14Q->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   gPad->SetLogy();
   MassNaF_PQ->SetLineColor (2);
   MassNaF_PQ->SetLineWidth (2);
   MassNaF_DQ->SetLineColor (4);
   MassNaF_DQ->SetLineWidth (2);
   MassNaF_DQ->SetFillColor (4);
   MassNaF_PQ->SetFillColor (2);
   MassNaF_PQ->SetFillStyle (3001);
   MassNaF_DQ->SetFillStyle (3002);
   MassNaF_PQ->GetXaxis()->SetTitle ("Mass [GeV/c^2]");
   MassNaF_PQ->SetTitle ("Mass NaF MC");
   MassNaF_PQ->GetXaxis()->SetTitleSize (0.045);
   MassNaF_PQ->Draw();
   MassNaF_DQ->Draw ("same");

   p15Q->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   gPad->SetLogy();
   MassAgl_PQ->SetLineColor (2);
   MassAgl_PQ->SetLineWidth (2);
   MassAgl_DQ->SetLineColor (4);
   MassAgl_DQ->SetLineWidth (2);
   MassAgl_DQ->SetFillColor (4);
   MassAgl_PQ->SetFillColor (2);
   MassAgl_PQ->SetFillStyle (3001);
   MassAgl_DQ->SetFillStyle (3002);
   MassAgl_PQ->GetXaxis()->SetTitle ("Mass [GeV/c^2]");
   MassAgl_PQ->SetTitle ("Mass Agl MC");
   MassAgl_PQ->GetXaxis()->SetTitleSize (0.045);
   MassAgl_PQ->Draw();
   MassAgl_DQ->Draw ("same");

   cout<<"******************* Distance discr. plots ******************"<<endl;


   p16->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   gPad->SetLogy();
   DistTOF_P->SetLineColor (2);
   DistTOF_P->SetLineWidth (2);
   DistTOF_He->SetLineWidth (2);
   DistTOF_D->SetLineColor (4);
   DistTOF_He->SetLineColor (3);
   DistTOF_D->SetLineWidth (2);
   DistTOF_D->SetFillColor (4);
   DistTOF_P->SetFillColor (2);
   DistTOF_He->SetFillColor (3);
   DistTOF_P->SetFillStyle (3001);
   DistTOF_D->SetFillStyle (3002);
   DistTOF_He->SetFillStyle (3002);
   DistTOF_He->GetXaxis()->SetTitle ("Distance discriminant");
   DistTOF_He->GetXaxis()->SetTitleSize (0.045);
   DistTOF_He->SetTitle ("Distance discr. distribution  TOF MC");
   DistTOF_He->Draw();
   DistTOF_D->Draw ("same");
   DistTOF_P->Draw ("same");

   p17->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   gPad->SetLogy();
   DistNaF_P->SetLineColor (2);
   DistNaF_P->SetLineWidth (2);
   DistNaF_He->SetLineWidth (2);
   DistNaF_D->SetLineColor (4);
   DistNaF_He->SetLineColor (3);
   DistNaF_D->SetLineWidth (2);
   DistNaF_D->SetFillColor (4);
   DistNaF_P->SetFillColor (2);
   DistNaF_He->SetFillColor (3);
   DistNaF_P->SetFillStyle (3001);
   DistNaF_D->SetFillStyle (3002);
   DistNaF_He->SetFillStyle (3002);
   DistNaF_He->GetXaxis()->SetTitle ("Distance discriminant");
   DistNaF_He->GetXaxis()->SetTitleSize (0.045);
   DistNaF_He->SetTitle ("Distance discr. distribution NaF MC");
   DistNaF_He->Draw();
   DistNaF_D->Draw ("same");
   DistNaF_P->Draw ("same");

   p18->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   gPad->SetLogy();
   DistAgl_P->SetLineColor (2);
   DistAgl_P->SetLineWidth (2);
   DistAgl_He->SetLineWidth (2);
   DistAgl_D->SetLineColor (4);
   DistAgl_He->SetLineColor (3);
   DistAgl_D->SetLineWidth (2);
   DistAgl_D->SetFillColor (4);
   DistAgl_P->SetFillColor (2);
   DistAgl_He->SetFillColor (3);
   DistAgl_P->SetFillStyle (3001);
   DistAgl_D->SetFillStyle (3002);
   DistAgl_He->SetFillStyle (3002);
   DistAgl_He->GetXaxis()->SetTitle ("Distance discriminant");
   DistAgl_He->GetXaxis()->SetTitleSize (0.045);
   DistAgl_He->SetTitle ("Distance discr. distribution Agl MC");
   DistAgl_He->Draw();
   DistAgl_D->Draw ("same");
   DistAgl_P->Draw ("same");

   cout<<"******************* R vs Dist plots ******************"<<endl;
   p19->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   RvsDistTOF_P->SetMarkerColor (2);
   RvsDistTOF_D->SetMarkerColor (4);
   RvsDistTOF_D->SetTitle ("R vs Distance discr. TOF (MC)");
   RvsDistTOF_D->GetXaxis()->SetTitle ("R [GV]");
   RvsDistTOF_D->GetYaxis()->SetTitle ("Distance discr. TOF");
   RvsDistTOF_P->GetZaxis()->SetRangeUser (10,4000);
   RvsDistTOF_D->GetZaxis()->SetRangeUser (1,4000);
   RvsDistTOF_D->Draw();
   RvsDistTOF_P->Draw ("same");


   p20->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   RvsDistNaF_P->SetMarkerColor (2);
   RvsDistNaF_D->SetMarkerColor (4);
   RvsDistNaF_D->SetTitle ("R vs Distance discr. NaF (MC)");
   RvsDistNaF_D->GetXaxis()->SetTitle ("R [GV]");
   RvsDistNaF_D->GetYaxis()->SetTitle ("Distance discr. NaF");
   RvsDistNaF_P->GetZaxis()->SetRangeUser (5,20);
   RvsDistNaF_D->GetZaxis()->SetRangeUser (1,20);
   RvsDistNaF_D->Draw();
   RvsDistNaF_P->Draw ("same");

   p21->cd();
   gPad->SetGridx();
   gPad->SetGridy();
   RvsDistAgl_P->SetMarkerColor (2);
   RvsDistAgl_D->SetMarkerColor (4);
   RvsDistAgl_D->SetTitle ("R vs Distance discr. Agl (MC)");
   RvsDistAgl_D->GetXaxis()->SetTitle ("R [GV]");
   RvsDistAgl_D->GetYaxis()->SetTitle ("Distance discr. Agl");
   RvsDistAgl_P->GetZaxis()->SetRangeUser (40,400);
   RvsDistAgl_D->GetZaxis()->SetRangeUser (1,400);
   RvsDistAgl_D->Draw();
   RvsDistAgl_P->Draw ("same");


   //fileFinalPlots->Flush();
   //fileFinalPlots->Close();
   
   
   finalPlots.Add(p1  );
   finalPlots.Add(p2  );
   finalPlots.Add(p3  );
   finalPlots.Add(p4  );
   finalPlots.Add(p5  );
   finalPlots.Add(p6  );
   finalPlots.Add(p7  );
   finalPlots.Add(p8  );
   finalPlots.Add(p9  );
   finalPlots.Add(p10 );
   finalPlots.Add(p11 );
   finalPlots.Add(p12 );
   finalPlots.Add(p13 );
   finalPlots.Add(p14 );
   finalPlots.Add(p15 );
   finalPlots.Add(p10Q);
   finalPlots.Add(p11Q);
   finalPlots.Add(p12Q);
   finalPlots.Add(p13Q);
   finalPlots.Add(p14Q);
   finalPlots.Add(p15Q);
   finalPlots.Add(p16 );
   finalPlots.Add(p17 );
   finalPlots.Add(p18 );
   finalPlots.Add(p19 );
   finalPlots.Add(p20 );
   finalPlots.Add(p21 );
   finalPlots.writeObjsInFolder("Common plots for slides");
   
}
