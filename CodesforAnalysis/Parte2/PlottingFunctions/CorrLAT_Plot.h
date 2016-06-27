
void CorrLAT_Plot(	TH1 *TOTLATCorrTOF,
                        TH1 *PreLATCorr,
                       TH1 *CorrezioneLATpre_pR,
			TH1 * CorrezioneLAT_pR,
                       TH1 * CorrezioneLAT_pTOF,
                       TH1 * CorrezioneLAT_dTOF,
                       TH1 * CorrezioneLAT_pNaF,
                       TH1 * CorrezioneLAT_dNaF,
                       TH1 * CorrezioneLAT_pAgl,
                       TH1 * CorrezioneLAT_dAgl
        ){


   TCanvas * c26 = new TCanvas("Latitude pile-up correction (R bins)");
   TCanvas * c26_bis = new TCanvas("Latitude pile-up correction (Beta bins)");



   c26->Divide(1,2);
   c26->cd(1);
   gPad->SetGridy();
   gPad->SetGridx();
   TGraphErrors *CorrLAT_tot_Spl=new TGraphErrors();
   for(int m=1; m<11; m++) {
      CorrLAT_tot_Spl->SetPoint(m-1,geomagC[m],TOTLATCorrTOF -> GetBinContent(m+1));
      CorrLAT_tot_Spl->SetPointError(m-1,0,TOTLATCorrTOF -> GetBinError(m+1));
   }
   TGraphErrors *CorrLAT_pre_Spl=new TGraphErrors();
   for(int m=1; m<11; m++) {
      CorrLAT_pre_Spl->SetPoint(m-1,geomagC[m],PreLATCorr -> GetBinContent(m+1));
      CorrLAT_pre_Spl->SetPointError(m-1,0,PreLATCorr -> GetBinError(m+1));
   }
   CorrLAT_tot_Spl->SetLineColor(2);
   CorrLAT_tot_Spl->SetMarkerColor(2);
   CorrLAT_tot_Spl->SetLineWidth(2);
   CorrLAT_tot_Spl->SetMarkerStyle(8);
   CorrLAT_pre_Spl->SetLineColor(4);
   CorrLAT_pre_Spl->SetMarkerColor(4);
   CorrLAT_pre_Spl->SetLineWidth(2);
   CorrLAT_pre_Spl->SetMarkerStyle(8);

   CorrLAT_tot_Spl->SetFillStyle(3002);
   CorrLAT_tot_Spl->GetXaxis()->SetTitle("Latitude");
   CorrLAT_tot_Spl->GetYaxis()->SetTitle("Eff. Corr. Factor");
   CorrLAT_tot_Spl->Draw("APC");
   CorrLAT_pre_Spl->Draw("PCsame");

   c26->cd(2);
   gPad->SetGridy();
   gPad->SetGridx();
   gPad->SetLogx();
   TGraphErrors *CorrLAT_totM1_Spl=new TGraphErrors();
   for(int i=0; i<nbinsr; i++) {
      CorrLAT_totM1_Spl->SetPoint(i,PRB.RigBinCent(i),CorrezioneLAT_pR->GetBinContent(i+1));
      CorrLAT_totM1_Spl->SetPointError(i,0,CorrezioneLAT_pR->GetBinError(i+1));
   }
   CorrLAT_totM1_Spl->SetLineColor(2);
   CorrLAT_totM1_Spl->SetMarkerColor(2);
   CorrLAT_totM1_Spl->SetLineWidth(2);
   CorrLAT_totM1_Spl->SetMarkerStyle(8);
   CorrLAT_totM1_Spl->SetFillStyle(3002);
   CorrLAT_totM1_Spl->GetXaxis()->SetTitle("R [GV]");
   CorrLAT_totM1_Spl->GetYaxis()->SetTitle("Eff. Corr. Factor");
   CorrLAT_totM1_Spl->Draw("APC");
   TGraphErrors *CorrLAT_totM2_Spl=new TGraphErrors();
   for(int i=0; i<nbinsr; i++) {
      CorrLAT_totM2_Spl->SetPoint(i,PRB.RigBinCent(i),CorrezioneLATpre_pR->GetBinContent(i+1));
      CorrLAT_totM2_Spl->SetPointError(i,0,CorrezioneLATpre_pR->GetBinError(i+1));	
   }
   CorrLAT_totM2_Spl->SetLineColor(4);
   CorrLAT_totM2_Spl->SetMarkerColor(4);
   CorrLAT_totM2_Spl->SetLineWidth(2);
   CorrLAT_totM2_Spl->SetMarkerStyle(8);
   CorrLAT_totM2_Spl->SetFillStyle(3002);
   CorrLAT_totM2_Spl->Draw("PCsame");

   c26_bis->Divide(1,3);
   c26_bis->cd(1);
   gPad->SetGridy();
   gPad->SetGridx();
   TGraphErrors * CorrLATp_TOF_Spl=new TGraphErrors();
   int point =0;
   for(int m=0; m<nbinsToF; m++) {
      if(CorrezioneLAT_pTOF->GetBinContent(m+1)>0)
         CorrLATp_TOF_Spl->SetPoint(point,ToFPB.EkBinCent(m),CorrezioneLAT_pTOF->GetBinContent(m+1));
      	 CorrLATp_TOF_Spl->SetPointError(point,0,CorrezioneLAT_pTOF->GetBinError(m+1));
	 point++;
   }
   CorrLATp_TOF_Spl->SetLineColor(2);
   CorrLATp_TOF_Spl->SetMarkerColor(2);
   CorrLATp_TOF_Spl->SetLineWidth(2);
   CorrLATp_TOF_Spl->SetLineStyle(2);
   CorrLATp_TOF_Spl->SetMarkerStyle(8);
   CorrLATp_TOF_Spl->SetFillStyle(3002);
   CorrLATp_TOF_Spl->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
   CorrLATp_TOF_Spl->GetYaxis()->SetTitle("Eff. Corr. Factor");
   CorrLATp_TOF_Spl->GetYaxis()->SetRangeUser(1.14,1.27);
   CorrLATp_TOF_Spl->Draw("APC");
   TGraphErrors * CorrLATd_TOF_Spl=new TGraphErrors();
   point=0;
   for(int m=0; m<nbinsToF; m++) {
      if(CorrezioneLAT_dTOF->GetBinContent(m+1)>0)
         CorrLATd_TOF_Spl->SetPoint(point,ToFPB.EkBinCent(m),CorrezioneLAT_dTOF->GetBinContent(m+1));
         CorrLATd_TOF_Spl->SetPointError(point,0,CorrezioneLAT_dTOF->GetBinError(m+1)); 
	point++;
   }
   CorrLATd_TOF_Spl->SetLineColor(4);
   CorrLATd_TOF_Spl->SetMarkerColor(4);
   CorrLATd_TOF_Spl->SetLineWidth(2);
   CorrLATd_TOF_Spl->SetLineStyle(2);
   CorrLATd_TOF_Spl->SetMarkerStyle(8);
   CorrLATd_TOF_Spl->SetFillStyle(3002);
   CorrLATd_TOF_Spl->Draw("PCsame");

   c26_bis->cd(2);
   gPad->SetGridy();
   gPad->SetGridx();
   TGraphErrors * CorrLATp_NaF_Spl=new TGraphErrors();
   for(int m=0; m<nbinsNaF; m++) {
      CorrLATp_NaF_Spl->SetPoint(m,NaFPB.EkBinCent(m),CorrezioneLAT_pNaF->GetBinContent(m+1));
      CorrLATp_NaF_Spl->SetPointError(m,0,CorrezioneLAT_pNaF->GetBinError(m+1));	
   }
   CorrLATp_NaF_Spl->SetLineColor(2);
   CorrLATp_NaF_Spl->SetMarkerColor(2);
   CorrLATp_NaF_Spl->SetLineWidth(2);
   CorrLATp_NaF_Spl->SetLineStyle(2);
   CorrLATp_NaF_Spl->SetMarkerStyle(8);
   CorrLATp_NaF_Spl->SetFillStyle(3002);
   CorrLATp_NaF_Spl->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
   CorrLATp_NaF_Spl->GetYaxis()->SetTitle("Eff. Corr. Factor");
   CorrLATp_NaF_Spl->GetYaxis()->SetRangeUser(1.08,1.22);
   CorrLATp_NaF_Spl->Draw("APC");
   TGraphErrors * CorrLATd_NaF_Spl=new TGraphErrors();
   for(int m=0; m<nbinsNaF; m++) {
      CorrLATd_NaF_Spl->SetPoint(m,NaFPB.EkBinCent(m),CorrezioneLAT_dNaF->GetBinContent(m+1));
      CorrLATd_NaF_Spl->SetPointError(m,0,CorrezioneLAT_dNaF->GetBinError(m+1));	
   }
   CorrLATd_NaF_Spl->SetLineColor(4);
   CorrLATd_NaF_Spl->SetMarkerColor(4);
   CorrLATd_NaF_Spl->SetLineWidth(2);
   CorrLATd_NaF_Spl->SetLineStyle(2);
   CorrLATd_NaF_Spl->SetMarkerStyle(8);
   CorrLATd_NaF_Spl->SetFillStyle(3002);
   CorrLATd_NaF_Spl->Draw("PCsame");

   c26_bis->cd(3);
   gPad->SetGridy();
   gPad->SetGridx();
   TGraphErrors * CorrLATp_Agl_Spl=new TGraphErrors();
   for(int m=0; m<nbinsAgl; m++) {
      CorrLATp_Agl_Spl->SetPoint(m,AglPB.EkBinCent(m),CorrezioneLAT_pAgl->GetBinContent(m+1));
      CorrLATp_Agl_Spl->SetPointError(m,0,CorrezioneLAT_pAgl->GetBinError(m+1));
   }
   CorrLATp_Agl_Spl->SetLineColor(2);
   CorrLATp_Agl_Spl->SetMarkerColor(2);
   CorrLATp_Agl_Spl->SetLineWidth(2);
   CorrLATp_Agl_Spl->SetLineStyle(2);
   CorrLATp_Agl_Spl->SetMarkerStyle(8);
   CorrLATp_Agl_Spl->SetFillStyle(3002);
   CorrLATp_Agl_Spl->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
   CorrLATp_Agl_Spl->GetYaxis()->SetTitle("Eff. Corr. Factor");
   CorrLATp_Agl_Spl->GetYaxis()->SetRangeUser(1.04,1.25);
   CorrLATp_Agl_Spl->Draw("APC");
   TGraphErrors * CorrLATd_Agl_Spl=new TGraphErrors();
   for(int m=0; m<nbinsAgl; m++) {
      CorrLATd_Agl_Spl->SetPoint(m,AglPB.EkBinCent(m),CorrezioneLAT_dAgl->GetBinContent(m+1));
      CorrLATd_Agl_Spl->SetPointError(m,0,CorrezioneLAT_dAgl->GetBinError(m+1));
  }
   CorrLATd_Agl_Spl->SetLineColor(4);
   CorrLATd_Agl_Spl->SetMarkerColor(4);
   CorrLATd_Agl_Spl->SetLineWidth(2);
   CorrLATd_Agl_Spl->SetLineStyle(2);
   CorrLATd_Agl_Spl->SetMarkerStyle(8);
   CorrLATd_Agl_Spl->Draw("CPsame");

   
   finalPlots.Add(c26);
   finalPlots.Add(c26_bis);
   finalPlots.writeObjsInFolder("DATA-driven Results/Latitude effect/Correction");

   return;


}

