TH1 * Weighted_CorrLAT(TH2F * esposizionegeo, TH1 * LATcorr);

void CorrLAT(string histoName) {
  inputHistoFile=TFile::Open(histoName.data(), "READ");

   TH2F * esposizionegeo_R    = (TH2F*)inputHistoFile->Get( "esposizionegeo"       );
   TH2F * esposizionepgeoTOF  = (TH2F*)inputHistoFile->Get(	"esposizionepgeo"	);
   TH2F * esposizionepgeoNaF  = (TH2F*)inputHistoFile->Get(	"esposizionepgeoNaF"	);
   TH2F * esposizionepgeoAgl  = (TH2F*)inputHistoFile->Get(	"esposizionepgeoAgl"	);
   TH2F * esposizionedgeoTOF  = (TH2F*)inputHistoFile->Get(	"esposizionedgeo"	);
   TH2F * esposizionedgeoNaF  = (TH2F*)inputHistoFile->Get(	"esposizionedgeoNaF"	);
   TH2F * esposizionedgeoAgl  = (TH2F*)inputHistoFile->Get(	"esposizionedgeoAgl"	);

   LATcorr * LATpreSelDATA      	= new LATcorr(inputHistoFile,"LATpreSelDATA"  	 ,"Results");

   LATcorr * LATLikelihoodDATA_TOF = new LATcorr(inputHistoFile,"LATLikDATA_TOF"   	 ,"Results");
   LATcorr * LATDistanceDATA_TOF   = new LATcorr(inputHistoFile,"LATDistDATA_TOF" 	 ,"Results");

   LATcorr * LATLikelihoodDATA_NaF = new LATcorr(inputHistoFile,"LATLikDATA_NaF"  	 ,"Results");
   LATcorr * LATDistanceDATA_NaF   = new LATcorr(inputHistoFile,"LATDistDATA_NaF" 	 ,"Results");

   LATcorr * LATLikelihoodDATA_Agl = new LATcorr(inputHistoFile,"LATLikDATA_Agl"  	 ,"Results");
   LATcorr * LATDistanceDATA_Agl   = new LATcorr(inputHistoFile,"LATDistDATA_Agl" 	 ,"Results");

   LATcorr * LATrichDATA_NaF       = new LATcorr(inputHistoFile,"LATrichDATA_NaF" 	 ,"Results");
   LATcorr * LATrichDATA_Agl       = new LATcorr(inputHistoFile,"LATrichDATA_Agl" 	 ,"Results");	

   cout<<"******* TOTAL LAT. CORRECTION *************"<<endl;


	TH2F*  PreLATCorr(static_cast<TH2F *>(LATpreSelDATA   -> LATcorrR_fit));
cout << PreLATCorr->GetEntries() << " " << PreLATCorr->GetName() << " " << PreLATCorr->ClassName ()  << endl;
   TH1F * LATpreSelDATA1 =      ProjectionXtoTH1F(PreLATCorr , "LATpreSelDATA1"    ,1,1) ;
   TH1F * LATpreSelDATA2 =  ProjectionXtoTH1F(PreLATCorr , "LATpreSelDATA2",2,2) ;
   TH1F * LATpreSelDATA3 =  ProjectionXtoTH1F(PreLATCorr , "LATpreSelDATA3",3,3) ;

   LATpreSelDATA1 -> Multiply( LATpreSelDATA2	);
   LATpreSelDATA1 -> Multiply( LATpreSelDATA3	);

   TH1F *  TOTLATCorrTOF( LATpreSelDATA1);
   TH1F *  TOTLATCorrNaF( LATpreSelDATA1);
   TH1F *  TOTLATCorrAgl( LATpreSelDATA1);

   TOTLATCorrTOF  -> Multiply (LATLikelihoodDATA_TOF ->  LATcorrR_fit);
   TOTLATCorrTOF  -> Multiply (LATDistanceDATA_TOF   ->  LATcorrR_fit);

   TOTLATCorrNaF  -> Multiply (LATLikelihoodDATA_NaF ->  LATcorrR_fit);
   TOTLATCorrNaF  -> Multiply (LATDistanceDATA_NaF   ->  LATcorrR_fit);
   TOTLATCorrNaF  -> Multiply (LATrichDATA_NaF       ->  LATcorrR_fit);   

   TOTLATCorrAgl  -> Multiply (LATLikelihoodDATA_Agl ->  LATcorrR_fit);
   TOTLATCorrAgl  -> Multiply (LATDistanceDATA_Agl   ->  LATcorrR_fit);
   TOTLATCorrNaF  -> Multiply (LATrichDATA_Agl       ->  LATcorrR_fit);

   //Only pres.
   TH1F * CorrezioneLATpre_pR = (TH1F *) Weighted_CorrLAT ( esposizionegeo_R , LATpreSelDATA1);

   TH1F * CorrezioneLATpre_dR = (TH1F *) Weighted_CorrLAT ( esposizionegeo_R , LATpreSelDATA1);

   //Full set
   TH1F * CorrezioneLAT_pR   = (TH1F *) Weighted_CorrLAT ( esposizionegeo_R   , TOTLATCorrTOF 	);
   TH1F * CorrezioneLAT_pTOF = (TH1F *) Weighted_CorrLAT ( esposizionepgeoTOF , TOTLATCorrTOF 	);
   TH1F * CorrezioneLAT_pNaF = (TH1F *) Weighted_CorrLAT ( esposizionepgeoNaF , TOTLATCorrNaF 	);
   TH1F * CorrezioneLAT_pAgl = (TH1F *) Weighted_CorrLAT ( esposizionepgeoAgl , TOTLATCorrAgl 	);

   TH1F * CorrezioneLAT_dR   = (TH1F *) Weighted_CorrLAT ( esposizionegeo_R   , TOTLATCorrTOF 	);
   TH1F * CorrezioneLAT_dTOF = (TH1F *) Weighted_CorrLAT ( esposizionedgeoTOF , TOTLATCorrTOF	);
   TH1F * CorrezioneLAT_dNaF = (TH1F *) Weighted_CorrLAT ( esposizionedgeoNaF , TOTLATCorrNaF	);
   TH1F * CorrezioneLAT_dAgl = (TH1F *) Weighted_CorrLAT ( esposizionedgeoAgl , TOTLATCorrAgl	);

   cout<<"*** Updating P1 file ****"<<endl;
   inputHistoFile->ReOpen("UPDATE");
   inputHistoFile->cd("Results");

    LATpreSelDATA1    -> Write(     "PreLATCorr_LATcorrR_fit" 	);
        LATpreSelDATA1-> Write(     "PreLATCorr_LATcorrTOF_fit");
        LATpreSelDATA1-> Write(     "PreLATCorr_LATcorrNaF_fit");
        LATpreSelDATA1-> Write(     "PreLATCorr_LATcorrAgl_fit");

   TOTLATCorrTOF  -> Write(     "TOTLATCorr_LATcorrR_fit"  );
   TOTLATCorrTOF  -> Write(     "TOTLATCorr_LATcorrTOF_fit");
   TOTLATCorrNaF  -> Write(     "TOTLATCorr_LATcorrNaF_fit");
   TOTLATCorrAgl  -> Write(     "TOTLATCorr_LATcorrAgl_fit");

   CorrezioneLATpre_pR -> Write(  "CorrezioneLATPrep_R"	);
   CorrezioneLATpre_dR -> Write(  "CorrezioneLATPred_R"	);

   CorrezioneLAT_pR  -> Write(  "CorrezioneLATp_R"  	);
   CorrezioneLAT_pTOF-> Write(  "CorrezioneLATp_TOF"	);
   CorrezioneLAT_pNaF-> Write(  "CorrezioneLATp_NaF"       );
   CorrezioneLAT_pAgl-> Write(  "CorrezioneLATp_Agl"	);

   CorrezioneLAT_dR  -> Write(  "CorrezioneLATd_R"  	);
   CorrezioneLAT_dTOF-> Write(  "CorrezioneLATd_TOF"	);
   CorrezioneLAT_dNaF-> Write(  "CorrezioneLATd_NaF"	);
   CorrezioneLAT_dAgl-> Write(  "CorrezioneLATd_Agl"       );

   inputHistoFile->Write();



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



TH1 * Weighted_CorrLAT(TH2F * esposizionegeo, TH1 * LATcorr) {
   TH2F * temp    = (TH2F *)esposizionegeo -> Clone();
   TH2F * temp_err= (TH2F *)esposizionegeo -> Clone(); 
  for(int m=0; m<11; m++) {
      for(int i=0; i< temp -> GetNbinsX(); i++) {
         temp    ->SetBinContent(i+1,m,esposizionegeo->GetBinContent(i+1,m)*LATcorr -> GetBinContent(m+1));
	 temp_err->SetBinContent(i+1,m,pow(esposizionegeo->GetBinContent(i+1,m)*LATcorr -> GetBinError(m+1),2));
	}
   }
   //summ all over latitudes	
   TH1F * temp2     = ProjectionXtoTH1F(temp     , "",0,10);
   TH1F * temp2_err = ProjectionXtoTH1F(temp_err , "",0,10);
  
   for(int i=0; i< temp2_err -> GetNbinsX(); i++) temp2_err->SetBinContent(i+1,pow(temp2_err->GetBinContent(i+1),0.5));
   //	
   
   //Divide by total exposure time
   TH1F * Exptime =ProjectionXtoTH1F( esposizionegeo , "",0,10);
   temp2     -> Divide ( Exptime );
   temp2_err -> Divide ( Exptime );
   //

   //setting errors
   for(int i=0; i< temp2 -> GetNbinsX(); i++)	temp2 -> SetBinError(i+1,temp2_err -> GetBinContent(i+1));

   return (TH1 *)temp2 -> Clone();
}
