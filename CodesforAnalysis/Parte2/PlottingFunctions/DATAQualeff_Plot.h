

void 	DATAQualeff_Plot( TH1 *LATDistDATATOF      ,
                          TH1 *LATLikDATATOF       ,
                          TH1 *LATDistDATANaF      ,
                          TH1 *LATLikDATANaF       ,
                          TH1 *LATDistDATAAgl      ,
                          TH1 *LATLikDATAAgl       ,
                            
                          TH1 *LikLATcorr_TOF      ,
                          TH1 *DistLATcorr_TOF     ,
                          TH1 *LikLATcorr_NaF      ,
                          TH1 *DistLATcorr_NaF     ,
                          TH1 *LikLATcorr_Agl      ,
                          TH1 *DistLATcorr_Agl     ,
                            
                          TH1 *LikLATcorr_TOF_fit  ,
                          TH1 *DistLATcorr_TOF_fit ,
                          TH1 *LikLATcorr_NaF_fit  ,
                          TH1 *DistLATcorr_NaF_fit ,
                          TH1 *LikLATcorr_Agl_fit  ,
                          TH1 *DistLATcorr_Agl_fit ){ 


   TCanvas *c15=new TCanvas ("Latitude Likelihood Efficiency");
   TCanvas *c16=new TCanvas ("Latitude Distance Efficiency");
   string Legend[11]={"Lat. Zone 0","Lat. Zone 1","Lat. Zone 2","Lat. Zone 3","Lat. Zone 4","Lat. Zone 5","Lat. Zone 6","Lat. Zone 7","Lat. Zone 8","Lat. Zone 9","Lat. Zone 10" };	

   c15->Divide (2,3);
   c16->Divide (2,3);

   c15->cd (1);
   gPad->SetLogx();
   gPad->SetGridx();
   gPad->SetGridy();
   TGraphErrors *EffDATALikP[11];
   for (int l=0; l<11; l++) {
      EffDATALikP[l]=new TGraphErrors();
      int j=0;
      for (int i=1; i<nbinsr; i++) {
         EffDATALikP[l]->SetPoint (j,PRB.RigBinCent (i),LATLikDATATOF->GetBinContent (i+1,l+1) );
         EffDATALikP[l]->SetPointError (j,0,LATLikDATATOF->GetBinError (i+1,l+1) );
         j++;
      }
   }
   cout<<endl;
   EffDATALikP[10]->SetMarkerColor (1);
   EffDATALikP[10]->SetMarkerStyle (8);
   EffDATALikP[10]->SetLineColor (1);
   EffDATALikP[10]->SetTitle ("Likelihood Latitude Efficiency ");
   EffDATALikP[10]->GetXaxis()->SetTitle ("R [GV]");
   EffDATALikP[10]->GetYaxis()->SetTitle ("Efficiency");
   EffDATALikP[10]->GetXaxis()->SetTitleSize (0.045);
   EffDATALikP[10]->GetYaxis()->SetTitleSize (0.045);
   EffDATALikP[10]->GetYaxis()->SetRangeUser (0.1,1.1);
   EffDATALikP[10]->Draw ("AP");
   for (int l=0; l<10; l++) {
      EffDATALikP[l]->SetMarkerColor (55+2*l);
      EffDATALikP[l]->SetMarkerStyle (8);
      EffDATALikP[l]->SetLineColor (55+2*l);
      EffDATALikP[l]->Draw ("Psame");
   }
   {
   TLegend* leg =new TLegend(0.8, 0.1,0.98,0.95);
  	for (int l=0; l<11; l++) leg->AddEntry(EffDATALikP[l],Legend[l].c_str(), "p");
   	leg->SetLineWidth(2);
	leg->Draw("same");
   }
	
   c15->cd (3);
   gPad->SetLogx();
   gPad->SetGridx();
   gPad->SetGridy();
   TGraphErrors *EffDATALikNaFP[11];
   for (int l=0; l<11; l++) {
      EffDATALikNaFP[l]=new TGraphErrors();
      int j=0;
      for (int i=1; i<nbinsr; i++) {
         EffDATALikNaFP[l]->SetPoint (j,PRB.RigBinCent (i),LATLikDATANaF ->GetBinContent (i+1,l+1) );
         EffDATALikNaFP[l]->SetPointError (j,0,LATLikDATANaF ->GetBinError (i+1,l+1) );
         j++;
      }
   }
   EffDATALikNaFP[10]->SetMarkerColor (1);
   EffDATALikNaFP[10]->SetMarkerStyle (8);
   EffDATALikNaFP[10]->SetLineColor (1);
   EffDATALikNaFP[10]->SetTitle ("Likelihood Latitude Efficiency (NaF)");
   EffDATALikNaFP[10]->GetXaxis()->SetTitle ("R [GV]");
   EffDATALikNaFP[10]->GetYaxis()->SetTitle ("Efficiency");
   EffDATALikNaFP[10]->GetXaxis()->SetTitleSize (0.045);
   EffDATALikNaFP[10]->GetYaxis()->SetTitleSize (0.045);
   EffDATALikNaFP[10]->GetYaxis()->SetRangeUser (0.1,1.1);
   EffDATALikNaFP[10]->Draw ("AP");
   for (int l=0; l<10; l++) {
      EffDATALikNaFP[l]->SetMarkerColor (55+2*l);
      EffDATALikNaFP[l]->SetMarkerStyle (8);
      EffDATALikNaFP[l]->SetLineColor (55+2*l);
      EffDATALikNaFP[l]->Draw ("Psame");
   }
   {
   TLegend* leg =new TLegend(0.8, 0.1,0.98,0.95);
        for (int l=0; l<11; l++) leg->AddEntry(EffDATALikNaFP[l],Legend[l].c_str(), "p");
   	leg->SetLineWidth(2);
	leg->Draw("same");
   }



   c15->cd (5);
   gPad->SetLogx();
   gPad->SetGridx();
   gPad->SetGridy();
   TGraphErrors *EffDATALikAglP[11];
   for (int l=0; l<11; l++) {
      EffDATALikAglP[l]=new TGraphErrors();
      int j=0;
      for (int i=1; i<nbinsr; i++) {
         EffDATALikAglP[l]->SetPoint (j,PRB.RigBinCent (i),LATLikDATAAgl ->GetBinContent (i+1,l+1) );
         EffDATALikAglP[l]->SetPointError (j,0,LATLikDATAAgl ->GetBinError (i+1,l+1) );
         j++;
      }
   }
   EffDATALikAglP[10]->SetMarkerColor (1);
   EffDATALikAglP[10]->SetMarkerStyle (8);
   EffDATALikAglP[10]->SetLineColor (1);
   EffDATALikAglP[10]->SetTitle ("Likelihood Latitude Efficiency (Agl)");
   EffDATALikAglP[10]->GetXaxis()->SetTitle ("R [GV]");
   EffDATALikAglP[10]->GetYaxis()->SetTitle ("Efficiency");
   EffDATALikAglP[10]->GetXaxis()->SetTitleSize (0.045);
   EffDATALikAglP[10]->GetYaxis()->SetTitleSize (0.045);
   EffDATALikAglP[10]->GetYaxis()->SetRangeUser (0.1,1.1);
   EffDATALikAglP[10]->Draw ("AP");
   for (int l=0; l<10; l++) {
      EffDATALikAglP[l]->SetMarkerColor (55+2*l);
      EffDATALikAglP[l]->SetMarkerStyle (8);
      EffDATALikAglP[l]->SetLineColor (55+2*l);
      EffDATALikAglP[l]->Draw ("Psame");
   }
   {
   TLegend* leg =new TLegend(0.8, 0.1,0.98,0.95);
        for (int l=0; l<11; l++) leg->AddEntry(EffDATALikAglP[l],Legend[l].c_str(), "p");
   	leg->SetLineWidth(2);
	leg->Draw("same");
   }
	


   TGraphErrors *CorrLATLik;
   c15->cd (2);
   gPad->SetGridy();
   gPad->SetGridx();
   CorrLATLik=new TGraphErrors();
   CorrLATLik->SetTitle ("Latitude Efficiency Corr.");
   CorrLATLik->GetXaxis()->SetTitle ("Latitude");
   CorrLATLik->GetYaxis()->SetTitle ("Eff. Corr. Factor");
   CorrLATLik->GetYaxis()->SetRangeUser (0.96,1.04);
   CorrLATLik->SetMarkerStyle (8);
   for (int i=1; i<11; i++) {
      CorrLATLik->SetPoint (i-1,geomagC[i],LikLATcorr_TOF->GetBinContent (i+1,1) );
      CorrLATLik->SetPointError (i-1,0,LikLATcorr_TOF->GetBinError (i+1,1) );
   }
   CorrLATLik->Draw ("AP");
   TGraphErrors *CorrLAT_Lik_Spl=new TGraphErrors ();
   CorrLAT_Lik_Spl->SetName ("CorrLAT_Lik_Spl");
   for (int i=1; i<11; i++) {
      CorrLAT_Lik_Spl->SetPoint (i-1,geomagC[i],LikLATcorr_TOF_fit->GetBinContent (i+1) );
      CorrLAT_Lik_Spl->SetPointError (i-1,0,LikLATcorr_TOF_fit->GetBinError (i+1) );
   }
   CorrLAT_Lik_Spl->SetLineColor (2);
   CorrLAT_Lik_Spl->SetMarkerColor (2);
   CorrLAT_Lik_Spl->SetFillColor (2);
   CorrLAT_Lik_Spl->SetFillStyle (3001);
   CorrLAT_Lik_Spl->SetLineWidth (2);
   CorrLAT_Lik_Spl->Draw ("Csame");


   TGraphErrors *CorrLATLikNaF;
   c15->cd (4);
   gPad->SetGridy();
   gPad->SetGridx();
   CorrLATLikNaF=new TGraphErrors();
   CorrLATLikNaF->SetTitle ("Latitude Efficiency Corr.");
   CorrLATLikNaF->GetXaxis()->SetTitle ("Latitude");
   CorrLATLikNaF->GetYaxis()->SetTitle ("Eff. Corr. Factor");
   CorrLATLikNaF->GetYaxis()->SetRangeUser (0.96,1.04);
   CorrLATLikNaF->SetMarkerStyle (8);
   for (int i=1; i<11; i++) {
      CorrLATLikNaF->SetPoint (i-1,geomagC[i],LikLATcorr_NaF->GetBinContent (i+1,1) );
      CorrLATLikNaF->SetPointError (i-1,0,LikLATcorr_NaF->GetBinError (i+1,1) );
   }
   CorrLATLikNaF->Draw ("AP");
   TGraphErrors *CorrLAT_LikNaF_Spl=new TGraphErrors ();
   CorrLAT_LikNaF_Spl->SetName ("CorrLAT_LikNaF_Spl");
   for (int i=1; i<11; i++) {
      CorrLAT_LikNaF_Spl->SetPoint (i-1,geomagC[i],LikLATcorr_NaF_fit->GetBinContent (i+1) );
      CorrLAT_LikNaF_Spl->SetPointError (i-1,0,LikLATcorr_NaF_fit->GetBinError (i+1) );
   }
   CorrLAT_LikNaF_Spl->SetLineColor (2);
   CorrLAT_LikNaF_Spl->SetMarkerColor (2);
   CorrLAT_LikNaF_Spl->SetFillColor (2);
   CorrLAT_LikNaF_Spl->SetFillStyle (3001);
   CorrLAT_LikNaF_Spl->SetLineWidth (2);
   CorrLAT_LikNaF_Spl->Draw ("Csame");

   TGraphErrors *CorrLATLikAgl;
   c15->cd (6);
   gPad->SetGridy();
   gPad->SetGridx();
   CorrLATLikAgl=new TGraphErrors();
   CorrLATLikAgl->SetTitle ("Latitude Efficiency Corr.");
   CorrLATLikAgl->GetXaxis()->SetTitle ("Latitude");
   CorrLATLikAgl->GetYaxis()->SetTitle ("Eff. Corr. Factor");
   CorrLATLikAgl->GetYaxis()->SetRangeUser (0.96,1.04);
   CorrLATLikAgl->SetMarkerStyle (8);
   for (int i=1; i<11; i++) {
      CorrLATLikAgl->SetPoint (i-1,geomagC[i],LikLATcorr_Agl->GetBinContent (i+1,1) );
      CorrLATLikAgl->SetPointError (i-1,0,LikLATcorr_Agl->GetBinError (i+1,1) );
   }
   CorrLATLikAgl->Draw ("AP");
   TGraphErrors *CorrLAT_LikAgl_Spl=new TGraphErrors ();
   CorrLAT_LikAgl_Spl->SetName ("CorrLAT_LikAgl_Spl");
   for (int i=1; i<11; i++) {
      CorrLAT_LikAgl_Spl->SetPoint (i-1,geomagC[i],LikLATcorr_Agl_fit->GetBinContent (i+1) );
      CorrLAT_LikAgl_Spl->SetPointError (i-1,0,LikLATcorr_Agl_fit->GetBinError (i+1) );
   }
   CorrLAT_LikAgl_Spl->SetLineColor (2);
   CorrLAT_LikAgl_Spl->SetMarkerColor (2);
   CorrLAT_LikAgl_Spl->SetFillColor (2);
   CorrLAT_LikAgl_Spl->SetFillStyle (3001);
   CorrLAT_LikAgl_Spl->SetLineWidth (2);
   CorrLAT_LikAgl_Spl->Draw ("Csame");



   c16->cd (1);
   gPad->SetLogx();
   gPad->SetGridx();
   gPad->SetGridy();
   TGraphErrors *EffDATADistP[11];
   for (int l=0; l<11; l++) {
      EffDATADistP[l]=new TGraphErrors();
      int j=0;
      for (int i=1; i<nbinsr; i++) {
         EffDATADistP[l]->SetPoint (j,PRB.RigBinCent (i),LATDistDATATOF ->GetBinContent (i+1,l+1) );
         EffDATADistP[l]->SetPointError (j,0, LATDistDATATOF ->GetBinError (i+1,l+1) );
         j++;
      }
   }
   EffDATADistP[10]->SetMarkerColor (1);
   EffDATADistP[10]->SetMarkerStyle (8);
   EffDATADistP[10]->SetLineColor (1);
   EffDATADistP[10]->SetTitle ("Distance Latitude Efficiency (TOF)");
   EffDATADistP[10]->GetXaxis()->SetTitle ("R [GV]");
   EffDATADistP[10]->GetYaxis()->SetTitle ("Efficiency");
   EffDATADistP[10]->GetXaxis()->SetTitleSize (0.045);
   EffDATADistP[10]->GetYaxis()->SetTitleSize (0.045);
   EffDATADistP[10]->GetYaxis()->SetRangeUser (0.1,1.1);
   EffDATADistP[10]->Draw ("AP");
   for (int l=0; l<10; l++) {
      EffDATADistP[l]->SetMarkerColor (55+2*l);
      EffDATADistP[l]->SetMarkerStyle (8);
      EffDATADistP[l]->SetLineColor (55+2*l);
      EffDATADistP[l]->Draw ("Psame");
   }
   {
   TLegend* leg =new TLegend(0.8, 0.1,0.98,0.95);
        for (int l=0; l<11; l++) leg->AddEntry(EffDATADistP[l],Legend[l].c_str(), "p");
   	leg->SetLineWidth(2);
	leg->Draw("same");
   }


   c16->cd (3);
   gPad->SetLogx();
   gPad->SetGridx();
   gPad->SetGridy();
   TGraphErrors *EffDATADistNaFP[11];
   for (int l=0; l<11; l++) {
      EffDATADistNaFP[l]=new TGraphErrors();
      int j=0;
      for (int i=1; i<nbinsr; i++) {
         EffDATADistNaFP[l]->SetPoint (j,PRB.RigBinCent (i),LATDistDATANaF ->GetBinContent (i+1,l+1) );
         EffDATADistNaFP[l]->SetPointError (j,0, LATDistDATANaF ->GetBinError (i+1,l+1) );
         j++;
      }
   }
   EffDATADistNaFP[10]->SetMarkerColor (1);
   EffDATADistNaFP[10]->SetMarkerStyle (8);
   EffDATADistNaFP[10]->SetLineColor (1);
   EffDATADistNaFP[10]->SetTitle ("Distance Latitude Efficiency (NaF)");
   EffDATADistNaFP[10]->GetXaxis()->SetTitle ("R [GV]");
   EffDATADistNaFP[10]->GetYaxis()->SetTitle ("Efficiency");
   EffDATADistNaFP[10]->GetXaxis()->SetTitleSize (0.045);
   EffDATADistNaFP[10]->GetYaxis()->SetTitleSize (0.045);
   EffDATADistNaFP[10]->GetYaxis()->SetRangeUser (0.1,1.1);
   EffDATADistNaFP[10]->Draw ("AP");
   for (int l=0; l<10; l++) {
      EffDATADistNaFP[l]->SetMarkerColor (55+2*l);
      EffDATADistNaFP[l]->SetMarkerStyle (8);
      EffDATADistNaFP[l]->SetLineColor (55+2*l);
      EffDATADistNaFP[l]->Draw ("Psame");
   }
   {
   TLegend* leg =new TLegend(0.8, 0.1,0.98,0.95);
        for (int l=0; l<11; l++) leg->AddEntry(EffDATADistNaFP[l],Legend[l].c_str(), "p");
   	leg->SetLineWidth(2);
	leg->Draw("same");
   }


   c16->cd (5);
   gPad->SetLogx();
   gPad->SetGridx();
   gPad->SetGridy();
   TGraphErrors *EffDATADistAglP[11];
   for (int l=0; l<11; l++) {
      EffDATADistAglP[l]=new TGraphErrors();
      int j=0;
      for (int i=1; i<nbinsr; i++) {
         EffDATADistAglP[l]->SetPoint (j,PRB.RigBinCent (i),LATDistDATAAgl ->GetBinContent (i+1,l+1) );
         EffDATADistAglP[l]->SetPointError (j,0, LATDistDATAAgl ->GetBinError (i+1,l+1) );
         j++;
      }
   }
   EffDATADistAglP[10]->SetMarkerColor (1);
   EffDATADistAglP[10]->SetMarkerStyle (8);
   EffDATADistAglP[10]->SetLineColor (1);
   EffDATADistAglP[10]->SetTitle ("Distance Latitude Efficiency (Agl)");
   EffDATADistAglP[10]->GetXaxis()->SetTitle ("R [GV]");
   EffDATADistAglP[10]->GetYaxis()->SetTitle ("Efficiency");
   EffDATADistAglP[10]->GetXaxis()->SetTitleSize (0.045);
   EffDATADistAglP[10]->GetYaxis()->SetTitleSize (0.045);
   EffDATADistAglP[10]->GetYaxis()->SetRangeUser (0.1,1.1);
   EffDATADistAglP[10]->Draw ("AP");
   for (int l=0; l<10; l++) {
      EffDATADistAglP[l]->SetMarkerColor (55+2*l);
      EffDATADistAglP[l]->SetMarkerStyle (8);
      EffDATADistAglP[l]->SetLineColor (55+2*l);
      EffDATADistAglP[l]->Draw ("Psame");
   }
   {
   TLegend* leg =new TLegend(0.8, 0.1,0.98,0.95);
        for (int l=0; l<11; l++) leg->AddEntry(EffDATADistAglP[l],Legend[l].c_str(), "p");
   	leg->SetLineWidth(2);
	leg->Draw("same");
   }


   TGraphErrors *CorrLATDist;
   c16->cd (2);
   gPad->SetGridy();
   gPad->SetGridx();
   CorrLATDist=new TGraphErrors();
   CorrLATDist->SetTitle ("Latitude Efficiency Corr.");
   CorrLATDist->GetXaxis()->SetTitle ("Latitude");
   CorrLATDist->GetYaxis()->SetTitle ("Eff. Corr. Factor");
   CorrLATDist->GetYaxis()->SetRangeUser (0.96,1.04);
   CorrLATDist->SetMarkerStyle (8);
   for (int i=1; i<11; i++) {
      CorrLATDist->SetPoint (i-1,geomagC[i],DistLATcorr_TOF->GetBinContent (i+1,1) );
      CorrLATDist->SetPointError (i-1,0,DistLATcorr_TOF->GetBinError (i+1,1) );
   }
   CorrLATDist->Draw ("AP");
   TGraphErrors *CorrLAT_Dist_Spl=new TGraphErrors ();
   CorrLAT_Dist_Spl->SetName ("CorrLAT_Dist_Spl");
   for (int i=1; i<11; i++) {
      CorrLAT_Dist_Spl->SetPoint (i-1,geomagC[i],DistLATcorr_TOF_fit->GetBinContent (i+1,1) );
      CorrLAT_Dist_Spl->SetPointError (i-1,0,DistLATcorr_TOF_fit->GetBinError (i+1,1) );
   }
   CorrLAT_Dist_Spl->SetLineColor (2);
   CorrLAT_Dist_Spl->SetMarkerColor (2);
   CorrLAT_Dist_Spl->SetFillColor (2);
   CorrLAT_Dist_Spl->SetFillStyle (3001);
   CorrLAT_Dist_Spl->SetLineWidth (2);
   CorrLAT_Dist_Spl->Draw ("Csame");


   TGraphErrors *CorrLATDistNaF;
   c16->cd (4);
   gPad->SetGridy();
   gPad->SetGridx();
   CorrLATDistNaF=new TGraphErrors();
   CorrLATDistNaF->SetTitle ("Latitude Efficiency Corr.");
   CorrLATDistNaF->GetXaxis()->SetTitle ("Latitude");
   CorrLATDistNaF->GetYaxis()->SetTitle ("Eff. Corr. Factor");
   CorrLATDistNaF->GetYaxis()->SetRangeUser (0.96,1.04);
   CorrLATDistNaF->SetMarkerStyle (8);
   for (int i=1; i<11; i++) {
      CorrLATDistNaF->SetPoint (i-1,geomagC[i],DistLATcorr_NaF->GetBinContent (i+1,1) );
      CorrLATDistNaF->SetPointError (i-1,0,DistLATcorr_NaF->GetBinError (i+1,1) );
   }
   CorrLATDistNaF->Draw ("AP");
   TGraphErrors *CorrLAT_DistNaF_Spl=new TGraphErrors ();
   CorrLAT_DistNaF_Spl->SetName ("CorrLAT_DistNaF_Spl");
   for (int i=1; i<11; i++) {
      CorrLAT_DistNaF_Spl->SetPoint (i-1,geomagC[i],DistLATcorr_NaF_fit->GetBinContent (i+1) );
      CorrLAT_DistNaF_Spl->SetPointError (i-1,0,DistLATcorr_NaF_fit->GetBinError (i+1) );
   }
   CorrLAT_DistNaF_Spl->SetLineColor (2);
   CorrLAT_DistNaF_Spl->SetMarkerColor (2);
   CorrLAT_DistNaF_Spl->SetFillColor (2);
   CorrLAT_DistNaF_Spl->SetFillStyle (3001);
   CorrLAT_DistNaF_Spl->SetLineWidth (2);
   CorrLAT_DistNaF_Spl->Draw ("Csame");

   TGraphErrors *CorrLATDistAgl;
   c16->cd (6);
   gPad->SetGridy();
   gPad->SetGridx();
   CorrLATDistAgl=new TGraphErrors();
   CorrLATDistAgl->SetTitle ("Latitude Efficiency Corr.");
   CorrLATDistAgl->GetXaxis()->SetTitle ("Latitude");
   CorrLATDistAgl->GetYaxis()->SetTitle ("Eff. Corr. Factor");
   CorrLATDistAgl->GetYaxis()->SetRangeUser (0.96,1.04);
   CorrLATDistAgl->SetMarkerStyle (8);
   for (int i=1; i<11; i++) {
      CorrLATDistAgl->SetPoint (i-1,geomagC[i],DistLATcorr_Agl->GetBinContent (i+1,1) );
      CorrLATDistAgl->SetPointError (i-1,0,DistLATcorr_Agl->GetBinError (i+1,1) );
   }
   CorrLATDistAgl->Draw ("AP");
   TGraphErrors *CorrLAT_DistAgl_Spl=new TGraphErrors ();
   CorrLAT_DistAgl_Spl->SetName ("CorrLAT_DistAgl_Spl");
   for (int i=1; i<11; i++) {
      CorrLAT_DistAgl_Spl->SetPoint (i-1,geomagC[i],DistLATcorr_Agl_fit->GetBinContent (i+1) );
      CorrLAT_DistAgl_Spl->SetPointError (i-1,0,DistLATcorr_Agl_fit->GetBinError (i+1) );
   }
   CorrLAT_DistAgl_Spl->SetLineColor (2);
   CorrLAT_DistAgl_Spl->SetMarkerColor (2);
   CorrLAT_DistAgl_Spl->SetFillColor (2);
   CorrLAT_DistAgl_Spl->SetFillStyle (3001);
   CorrLAT_DistAgl_Spl->SetLineWidth (2);
   CorrLAT_DistAgl_Spl->Draw ("Csame");

finalPlots.Add(c15);
	
	finalPlots.Add(c15);
	finalPlots.Add(c16);
   finalPlots.writeObjsInFolder("DATA-driven Results/Latitude effect/Quality");
	finalPlots.Add(CorrLAT_Lik_Spl    );
	finalPlots.Add(CorrLAT_LikNaF_Spl );
	finalPlots.Add(CorrLAT_LikAgl_Spl );
	finalPlots.Add(CorrLAT_Dist_Spl   );
	finalPlots.Add(CorrLAT_DistNaF_Spl);
	finalPlots.Add(CorrLAT_DistAgl_Spl);
	finalPlots.writeObjsInFolder("Export");



}
