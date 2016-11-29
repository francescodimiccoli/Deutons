void	Acceptance_Plot (TH1* PGen_Acceptance_R   	,
                        TH1 * PGen_Acceptance_TOF	,
                        TH1 * PGen_Acceptance_NaF 	,
                        TH1 * PGen_Acceptance_Agl	,
                        TH1 * PrePMCAcceptance_R  	,
                        TH1 * PrePGeomag_Acceptance_R 	,
                        TH1 * PMCAcceptance_R        	,
                        TH1 * PMCAcceptance_TOF	,
                        TH1 * PMCAcceptance_NaF	,
                        TH1 * PMCAcceptance_Agl	,
                        TH1 * PGeomag_Acceptance_R 	,
                        TH1 * PGeomag_Acceptance_TOF	,
                        TH1 * PGeomag_Acceptance_NaF	,
                        TH1 * PGeomag_Acceptance_Agl	,
                        TH1 * PCorrectedAcceptance_R   ,
                        TH1 * PCorrectedAcceptance_TOF,	
                        TH1 * PCorrectedAcceptance_NaF,	
                        TH1 * PCorrectedAcceptance_Agl,	
                       	TH1 * DGen_Acceptance_R  	,
                        TH1 * DGen_Acceptance_TOF	,
                        TH1 * DGen_Acceptance_NaF     ,	
                        TH1 * DGen_Acceptance_Agl	,
                        TH1 * DMCAcceptance_R         ,	
                        TH1 * DMCAcceptance_TOF	,
                        TH1 * DMCAcceptance_NaF	,
                        TH1 * DMCAcceptance_Agl	,
                        TH1 * DGeomag_Acceptance_R  	,
                        TH1 * DGeomag_Acceptance_TOF	,
                        TH1 * DGeomag_Acceptance_NaF	,
                        TH1 * DGeomag_Acceptance_Agl	,
                        TH1 * DCorrectedAcceptance_R  ,	
                        TH1 * DCorrectedAcceptance_TOF,	
                        TH1 * DCorrectedAcceptance_NaF,	
	                TH1 * DCorrectedAcceptance_Agl,

			TH1F * LatCorrErrR_P   	,	
                        TH1F * LatCorrErrTOF_P ,
                        TH1F * LatCorrErrNaF_P, 
                        TH1F * LatCorrErrAgl_P, 
                        
                        TH1F * EffStatErrR_P   ,
                        TH1F * EffStatErrTOF_P, 
                        TH1F * EffStatErrNaF_P, 
                        TH1F * EffStatErrAgl_P, 
                                               
                        TH1F * EffSystErrR_P   ,
                        TH1F * EffSystErrTOF_P, 
                        TH1F * EffSystErrNaF_P, 
                        TH1F * EffSystErrAgl_P, 
                                               
                        TH1F * CorrStatErrR_P  ,
                        TH1F * CorrStatErrTOF_P,
                        TH1F * CorrStatErrNaF_P,
                        TH1F * CorrStatErrAgl_P,
                                               
                        TH1F * CorrSystErrR_P  ,
                        TH1F * CorrSystErrTOF_P,
                        TH1F * CorrSystErrNaF_P,
                        TH1F * CorrSystErrAgl_P,
			
			TH1F * RICHStatErrNaF_P,
                        TH1F * RICHStatErrAgl_P,
                                               
                        TH1F * RICHSystErrNaF_P,
                        TH1F * RICHSystErrAgl_P,

			TH1F * EffStatErrR_D  , 
                        TH1F * EffStatErrTOF_D, 
                        TH1F * EffStatErrNaF_D, 
                        TH1F * EffStatErrAgl_D, 
                                               
                        TH1F * EffSystErrR_D  , 
                        TH1F * EffSystErrTOF_D, 
                        TH1F * EffSystErrNaF_D, 
                        TH1F * EffSystErrAgl_D 


){


	int plottingstyles[6]={3,4,20,25,29,26};

	TCanvas * c31_tris = new TCanvas("Gen. Efficiency");

	string MCLegend[7]= {"Protons MC B800","Deuteorons MC \"GG_Blic\"","Deuterons MC \"GG_BlicDPMJet\"","Deuterons MC \"GG_QMD\"","Deuterons MC \"Shen_Blic\"","Deuterons MC \"Shen_BlicDPMJet\"","Deuterons MC \"Shen_QMD\""};

	//c31_tris->Divide(2,1);
	c31_tris->cd();
	gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetGridx();
        gPad->SetGridy();
	TGraphErrors * EffgenbetaD[6];
	int p =0;
	for(int h=0;h<6;h++){
		EffgenbetaD[h]=new TGraphErrors();
		p=0;
		for(int i=0;i<nbinsr;i++) {EffgenbetaD[h]->SetPoint(p,DRB.EkPerMassBinCent(i),DGen_Acceptance_R  ->GetBinContent(i+1,h+1));p++;}
					
	
	EffgenbetaD[h]->SetMarkerStyle(8);
        EffgenbetaD[h]->SetMarkerColor(0);
        EffgenbetaD[h]->SetMarkerSize(1);
        EffgenbetaD[h]->SetLineColor(0);
        EffgenbetaD[h]->SetLineWidth(0);
        EffgenbetaD[h]->SetMarkerStyle(8);
        EffgenbetaD[h]->SetTitle("");
        EffgenbetaD[h]->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
        EffgenbetaD[h]->GetYaxis()->SetTitle("Efficiency");
        EffgenbetaD[h]->GetXaxis()->SetTitleSize(0.045);
        EffgenbetaD[h]->GetYaxis()->SetTitleSize(0.045);
	}
	EffgenbetaD[0]->Draw("ACP");
	


	
	TGraphErrors * EffgenbetaDTOF[6];
        p=0;
        for(int h=0;h<6;h++){
                EffgenbetaDTOF[h]=new TGraphErrors();
                p=0;
                for(int i=0;i<nbinsToF;i++) {EffgenbetaDTOF[h]->SetPoint(p,ToFDB.EkPerMassBinCent(i),DGen_Acceptance_TOF  ->GetBinContent(i+1,h+1));p++;}
        EffgenbetaDTOF[h]->SetMarkerStyle(8);
        EffgenbetaDTOF[h]->SetMarkerColor(4);
        EffgenbetaDTOF[h]->SetMarkerSize(2);
        EffgenbetaDTOF[h]->SetLineColor(4);
        EffgenbetaDTOF[h]->SetLineWidth(1);
        EffgenbetaDTOF[h]->SetMarkerStyle(plottingstyles[h]);
        EffgenbetaDTOF[h]->SetTitle("");
        EffgenbetaDTOF[h]->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
        EffgenbetaDTOF[h]->GetYaxis()->SetTitle("Efficiency");
        EffgenbetaDTOF[h]->GetXaxis()->SetTitleSize(0.045);
        EffgenbetaDTOF[h]->GetYaxis()->SetTitleSize(0.045);
        }
        EffgenbetaDTOF[0]->Draw("CPsame");
        for(int h=1;h<6;h++){
                EffgenbetaDTOF[h]->Draw("CPsame");
        }


	
	TGraphErrors * EffgenbetaDNaF[6];
        p=0;
        for(int h=0;h<6;h++){
                EffgenbetaDNaF[h]=new TGraphErrors();
                p=0;
                for(int i=0;i<nbinsNaF;i++) {EffgenbetaDNaF[h]->SetPoint(p,NaFDB.EkPerMassBinCent(i),DGen_Acceptance_NaF  ->GetBinContent(i+1,h+1));p++;}

        EffgenbetaDNaF[h]->SetMarkerStyle(8);
        EffgenbetaDNaF[h]->SetMarkerColor(4);
        EffgenbetaDNaF[h]->SetMarkerSize(2);
        EffgenbetaDNaF[h]->SetLineColor(4);
        EffgenbetaDNaF[h]->SetLineWidth(1);
        EffgenbetaDNaF[h]->SetMarkerStyle(plottingstyles[h]);
        EffgenbetaDNaF[h]->SetTitle("");
        EffgenbetaDNaF[h]->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
        EffgenbetaDNaF[h]->GetYaxis()->SetTitle("Efficiency");
        EffgenbetaDNaF[h]->GetXaxis()->SetTitleSize(0.045);
        EffgenbetaDNaF[h]->GetYaxis()->SetTitleSize(0.045);
        }
        EffgenbetaDNaF[0]->Draw("CPsame");
        for(int h=1;h<6;h++){
                EffgenbetaDNaF[h]->Draw("CPsame");
        }
	
	TGraphErrors * EffgenbetaDAgl[6];
        p=0;
        for(int h=0;h<6;h++){
                EffgenbetaDAgl[h]=new TGraphErrors();
                p=0;
                for(int i=0;i<nbinsAgl;i++) {EffgenbetaDAgl[h]->SetPoint(p,AglDB.EkPerMassBinCent(i),DGen_Acceptance_Agl ->GetBinContent(i+1,h+1));p++;}

        EffgenbetaDAgl[h]->SetMarkerStyle(8);
        EffgenbetaDAgl[h]->SetMarkerColor(4);
        EffgenbetaDAgl[h]->SetMarkerSize(2);
        EffgenbetaDAgl[h]->SetLineColor(4);
        EffgenbetaDAgl[h]->SetLineWidth(1);
        EffgenbetaDAgl[h]->SetMarkerStyle(plottingstyles[h]);
        EffgenbetaDAgl[h]->SetTitle("");
        EffgenbetaDAgl[h]->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
        EffgenbetaDAgl[h]->GetYaxis()->SetTitle("Efficiency");
        EffgenbetaDAgl[h]->GetXaxis()->SetTitleSize(0.045);
        EffgenbetaDAgl[h]->GetYaxis()->SetTitleSize(0.045);
        }
        EffgenbetaDAgl[0]->Draw("CPsame");
        for(int h=1;h<6;h++){
                EffgenbetaDAgl[h]->Draw("CPsame");
        }
	
	TGraphErrors * EffgenbetaP;
	p =0;
	EffgenbetaP=new TGraphErrors();
	p=0;
	for(int i=0;i<nbinsr;i++) {EffgenbetaP->SetPoint(p,PRB.EkPerMassBinCent(i),PGen_Acceptance_R  ->GetBinContent(i+1));p++;}

	EffgenbetaP->SetMarkerStyle(8);
	EffgenbetaP->SetMarkerColor(0);
	EffgenbetaP->SetMarkerSize(1);
	EffgenbetaP->SetLineColor(0);
	EffgenbetaP->SetLineWidth(1);
	EffgenbetaP->SetTitle("");
	EffgenbetaP->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
	EffgenbetaP->GetYaxis()->SetTitle("Efficiency");
	EffgenbetaP->GetXaxis()->SetRangeUser(0.1,20);
	EffgenbetaP->GetXaxis()->SetTitleSize(0.045);
	EffgenbetaP->GetYaxis()->SetTitleSize(0.045);
	//EffgenbetaP->Draw("ACP");

	TGraphErrors * EffgenbetaPTOF;
	p=0;
	EffgenbetaPTOF=new TGraphErrors();
	p=0;
	for(int i=0;i<nbinsToF;i++) {EffgenbetaPTOF->SetPoint(p,ToFPB.EkPerMassBinCent(i),PGen_Acceptance_TOF  ->GetBinContent(i+1));p++;}
	EffgenbetaPTOF->SetMarkerStyle(8);
	EffgenbetaPTOF->SetMarkerColor(2);
	EffgenbetaPTOF->SetMarkerSize(2);
	EffgenbetaPTOF->SetLineColor(2);
	EffgenbetaPTOF->SetLineWidth(1);
	EffgenbetaPTOF->SetTitle("");
	EffgenbetaPTOF->GetXaxis()->SetTitle("Kin. En. / nucl.  [GeV/nucl.]");
	EffgenbetaPTOF->GetYaxis()->SetTitle("Efficiency");
	EffgenbetaPTOF->GetXaxis()->SetTitleSize(0.045);
	EffgenbetaPTOF->GetYaxis()->SetTitleSize(0.045);
	EffgenbetaPTOF->GetXaxis()->SetRangeUser(0.1,20);

	EffgenbetaPTOF->Draw("CPsame");
	{ TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(EffgenbetaDTOF[0],MCLegend[1].c_str(), "p");
        for(int h=1;h<6;h++){
                leg->AddEntry(EffgenbetaDTOF[h],MCLegend[h+1].c_str(), "p");
        }
               leg->AddEntry(EffgenbetaPTOF,MCLegend[0].c_str(), "p");	       
               leg->SetLineWidth(2);
                leg->Draw("same");
        }	




	TGraphErrors * EffgenbetaPNaF;
	p=0;
	EffgenbetaPNaF=new TGraphErrors();
	p=0;
	for(int i=0;i<nbinsNaF;i++) {EffgenbetaPNaF->SetPoint(p,NaFPB.EkPerMassBinCent(i),PGen_Acceptance_NaF  ->GetBinContent(i+1));p++;}
	EffgenbetaPNaF->SetMarkerStyle(8);
	EffgenbetaPNaF->SetMarkerColor(2);
	EffgenbetaPNaF->SetMarkerSize(2);
	EffgenbetaPNaF->SetLineColor(2);
	EffgenbetaPNaF->SetLineWidth(1);
	EffgenbetaPNaF->SetTitle("");
	EffgenbetaPNaF->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
	EffgenbetaPNaF->GetYaxis()->SetTitle("Efficiency");
	EffgenbetaPNaF->GetXaxis()->SetTitleSize(0.045);
	EffgenbetaPNaF->GetYaxis()->SetTitleSize(0.045);
	EffgenbetaPNaF->Draw("CPsame");

	TGraphErrors * EffgenbetaPAgl;
	p=0;
	EffgenbetaPAgl=new TGraphErrors();
	p=0;
	for(int i=0;i<nbinsAgl;i++) {EffgenbetaPAgl->SetPoint(p,AglPB.EkPerMassBinCent(i),PGen_Acceptance_Agl  ->GetBinContent(i+1));p++;}
	EffgenbetaPAgl->SetMarkerStyle(8);
	EffgenbetaPAgl->SetMarkerColor(2);
	EffgenbetaPAgl->SetMarkerSize(2);
	EffgenbetaPAgl->SetLineColor(2);
	EffgenbetaPAgl->SetLineWidth(1);
	EffgenbetaPAgl->SetTitle("");
	EffgenbetaPAgl->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
	EffgenbetaPAgl->GetYaxis()->SetTitle("Efficiency");
	EffgenbetaPAgl->GetXaxis()->SetTitleSize(0.045);
	EffgenbetaPAgl->GetYaxis()->SetTitleSize(0.045);
	EffgenbetaPAgl->Draw("CPsame");

	TCanvas * c22 = new TCanvas("Protons Acceptance");
	c22->cd();
	gPad->SetLogx();
        gPad->SetLogy();
	gPad->SetGridx();
        gPad->SetGridy();
	TGraphErrors * AccgeoP= new TGraphErrors();
	TGraphErrors * AccPreMCP= new TGraphErrors();
	TGraphErrors * AccSelMCP= new TGraphErrors();
	TGraphErrors * AccSelP[11];
	TGraphErrors * AccpreP[11];
	p=0;
	for(int i=0;i<nbinsr;i++) {AccgeoP->SetPoint(p,PRB.EkPerMassBinCent(i),PGen_Acceptance_R->GetBinContent(i+1));
				   AccgeoP->SetPointError(p,0,PGen_Acceptance_R->GetBinError(i+1));
				   p++;}	
	p=0;
        for(int i=0;i<nbinsr;i++) {AccPreMCP->SetPoint(p,PRB.EkPerMassBinCent(i),PrePMCAcceptance_R->GetBinContent(i+1));
				   AccPreMCP->SetPointError(p,0,PrePMCAcceptance_R->GetBinError(i+1));
				   p++;}
	p=0;
	for(int i=0;i<nbinsr;i++) {AccSelMCP->SetPoint(p,PRB.EkPerMassBinCent(i),PMCAcceptance_R->GetBinContent(i+1));
				   AccSelMCP->SetPointError(p,0,PMCAcceptance_R->GetBinError(i+1));
				   p++;}
	

	for(int j=0;j<11;j++) {
		AccSelP[j]=new TGraphErrors();
		p=0;
		for(int i=0;i<nbinsr;i++) {AccSelP[j]->SetPoint(p,PRB.EkPerMassBinCent(i),PGeomag_Acceptance_R -> GetBinContent(i+1,j+1) );
					  AccSelP[j]->SetPointError(p,0,PGeomag_Acceptance_R -> GetBinError(i+1,j+1));
					  p++;}
		AccSelP[j]->SetMarkerStyle(8);
        	AccSelP[j]->SetMarkerColor(j-1);
        	AccSelP[j]->SetLineColor(j-1);
        	AccSelP[j]->SetLineWidth(2);
	}
	for(int j=0;j<11;j++) {
                AccpreP[j]=new TGraphErrors();
                p=0;
 		for(int i=0;i<nbinsr;i++) {AccpreP[j]->SetPoint(p,PRB.EkPerMassBinCent(i),PrePGeomag_Acceptance_R -> GetBinContent(i+1,j+1));
					   AccpreP[j]->SetPointError(p,0,PrePGeomag_Acceptance_R -> GetBinError(i+1,j+1));
					   p++;}
                AccpreP[j]->SetMarkerStyle(8);
                AccpreP[j]->SetMarkerColor(j-1);
                AccpreP[j]->SetLineColor(j-1);
                AccpreP[j]->SetLineWidth(2);
        }

	AccgeoP->SetMarkerStyle(8);
	AccgeoP->SetMarkerColor(2);
	AccgeoP->SetLineColor(2);
	AccgeoP->SetLineWidth(4);
	AccPreMCP->SetMarkerStyle(8);
        AccPreMCP->SetMarkerColor(1);
        AccPreMCP->SetLineColor(2);
        AccPreMCP->SetLineWidth(4);
	AccSelMCP->SetMarkerStyle(8);
        AccSelMCP->SetMarkerColor(1);
        AccSelMCP->SetLineColor(2);
	AccgeoP->SetTitle("Protons Acceptance");
        AccgeoP->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
        AccgeoP->GetYaxis()->SetTitle("Acceptance [m^2 sr]");
        AccgeoP->GetXaxis()->SetTitleSize(0.045);
        AccgeoP->GetYaxis()->SetTitleSize(0.045);
	AccgeoP->GetYaxis()->SetRangeUser(1e-2,1.3);
	AccgeoP->Draw("AC");
	for(int j=0;j<11;j++) AccSelP[j]->Draw("PCsame");
	for(int j=0;j<11;j++) AccpreP[j]->Draw("PCsame");
	AccPreMCP->Draw("Csame");
	AccSelMCP->Draw("Csame");



	TCanvas * c31_bis = new TCanvas("MC Final Acceptance (Beta bins) (Beta bins)");
	
	c31_bis->cd();

        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetGridx();
        gPad->SetGridy();
	TH2F * Frame = new TH2F("Frame","Frame",1000,0,50,1000,0.0001,1);
      Frame->SetTitle("Final effective Acceptance");
       Frame->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
      Frame->GetYaxis()->SetTitle("Acceptance [m^2 sr]");
       Frame->GetXaxis()->SetTitleSize(0.045);
       Frame->GetYaxis()->SetTitleSize(0.045);
       Frame->GetYaxis()->SetRangeUser(1e-3,1.3);
        	

	  Frame->Draw();

        TGraphErrors * AccSelMCDbeta[6];
        TGraphErrors * AccSelMCPbeta=new TGraphErrors();
	for(int h=0;h<6;h++){
                AccSelMCDbeta[h]= new TGraphErrors();
                int p=0;
                for(int m=0;m<nbinsToF;m++) {AccSelMCDbeta[h]->SetPoint(p,ToFDB.EkPerMassBinCent(m),DMCAcceptance_TOF -> GetBinContent(m+1,h+1));
					    AccSelMCDbeta[h]->SetPointError(p,0,DMCAcceptance_TOF -> GetBinError(m+1,h+1));
					     p++;}
		if(h==1) 
		{
			p=0;
			for(int m=0;m<nbinsToF;m++) {AccSelMCPbeta->SetPoint(p,ToFPB.EkPerMassBinCent(m),PMCAcceptance_TOF -> GetBinContent(m+1));
						     AccSelMCPbeta->SetPointError(p,0,PMCAcceptance_TOF -> GetBinError(m+1));
						     p++;}
			AccSelMCPbeta->SetMarkerStyle(8);
                	AccSelMCPbeta->SetFillColor(2);
                	AccSelMCPbeta->SetMarkerSize(2);
                	AccSelMCPbeta->SetLineColor(2);
                	AccSelMCPbeta->SetLineWidth(4);
			AccSelMCPbeta->SetFillStyle(3001);		

		AccSelMCDbeta[h]->SetMarkerStyle(8);
                AccSelMCDbeta[h]->SetFillColor(4);
                AccSelMCDbeta[h]->SetMarkerSize(2);
                AccSelMCDbeta[h]->SetLineColor(4);
		AccSelMCDbeta[h]->SetLineStyle(1);
		AccSelMCDbeta[h]->SetFillStyle(3001);
                AccSelMCDbeta[h]->SetLineWidth(4);
                AccSelMCDbeta[h]->SetMarkerStyle(h+3);
                AccSelMCDbeta[1]->SetTitle("Final effective Acceptance");
                AccSelMCDbeta[1]->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
                AccSelMCDbeta[1]->GetYaxis()->SetTitle("Acceptance [m^2 sr]");
                AccSelMCDbeta[1]->GetXaxis()->SetTitleSize(0.045);
                AccSelMCDbeta[1]->GetYaxis()->SetTitleSize(0.045);
                AccSelMCDbeta[1]->GetYaxis()->SetRangeUser(1e-2,1.3);
        	}
	}
	{	TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
		AccSelMCDbeta[1]->Draw("C4same");
		AccSelMCPbeta->SetLineStyle(1);
		AccSelMCPbeta->Draw("C4same");
        	leg->AddEntry(AccSelMCPbeta,"Proton Acceptance","fl");		
		leg->AddEntry(AccSelMCDbeta[1],"Deuteron Acceptance","fl");
		leg->SetLineWidth(2);
		leg->Draw("same");
	}


        TGraphErrors * AccSelMCDbetaNaF[6];
        TGraphErrors * AccSelMCPbetaNaF=new TGraphErrors();
	for(int h=0;h<6;h++){
                AccSelMCDbetaNaF[h]= new TGraphErrors();
                int p=0;
                for(int m=0;m<nbinsNaF;m++) {AccSelMCDbetaNaF[h]->SetPoint(p,NaFDB.EkPerMassBinCent(m),DMCAcceptance_NaF -> GetBinContent(m+1,h+1));
					    AccSelMCDbetaNaF[h]->SetPointError(p,0,DMCAcceptance_NaF -> GetBinError(m+1,h+1));
					     p++;}
		if(h==0)
                {
                        p=0;
			for(int m=0;m<nbinsNaF;m++) {AccSelMCPbetaNaF->SetPoint(p,NaFPB.EkPerMassBinCent(m),PMCAcceptance_NaF -> GetBinContent(m+1));
					 	     AccSelMCPbetaNaF->SetPointError(p,0,PMCAcceptance_NaF -> GetBinError(m+1));
						     p++;}
                        AccSelMCPbetaNaF->SetMarkerStyle(8);
                        AccSelMCPbetaNaF->SetFillColor(2);
                        AccSelMCPbetaNaF->SetMarkerSize(2);
                        AccSelMCPbetaNaF->SetLineColor(2);
                        AccSelMCPbetaNaF->SetLineWidth(4);
			AccSelMCPbetaNaF->SetFillStyle(3001);
                }
                AccSelMCDbetaNaF[h]->SetMarkerStyle(8);
                AccSelMCDbetaNaF[h]->SetFillColor(4);
               AccSelMCDbetaNaF[h]->SetFillStyle(3001);
		 AccSelMCDbetaNaF[h]->SetMarkerSize(2);
                AccSelMCDbetaNaF[h]->SetLineColor(4);
                AccSelMCDbetaNaF[h]->SetLineStyle(1);
		AccSelMCDbetaNaF[h]->SetLineWidth(4);
                AccSelMCDbetaNaF[h]->SetMarkerStyle(h+3);
        }
                AccSelMCDbetaNaF[1]->Draw("C4same");
		AccSelMCPbetaNaF->Draw("C4same");
	

        TGraphErrors * AccSelMCDbetaAgl[6];
        TGraphErrors * AccSelMCPbetaAgl=new TGraphErrors();
	for(int h=0;h<6;h++){
                AccSelMCDbetaAgl[h]= new TGraphErrors();
                int p=0;
                for(int m=0;m<nbinsAgl;m++)  {AccSelMCDbetaAgl[h]->SetPoint(p,AglDB.EkPerMassBinCent(m),DMCAcceptance_Agl -> GetBinContent(m+1,h+1));
					      AccSelMCDbetaAgl[h]->SetPointError(p,0,DMCAcceptance_Agl -> GetBinError(m+1,h+1));
					      p++;}
		if(h==0)
                {
                        p=0;
			for(int m=0;m<nbinsAgl;m++) {AccSelMCPbetaAgl->SetPoint(p,AglPB.EkPerMassBinCent(m),PMCAcceptance_Agl -> GetBinContent(m+1));
						     AccSelMCPbetaAgl->SetPointError(p,0,PMCAcceptance_Agl -> GetBinError(m+1));
						     p++;}
                        AccSelMCPbetaAgl->SetMarkerStyle(8);
                        AccSelMCPbetaAgl->SetFillColor(2);
                        AccSelMCPbetaAgl->SetMarkerSize(2);
                        AccSelMCPbetaAgl->SetLineColor(2);
                        AccSelMCPbetaAgl->SetLineWidth(4);
			AccSelMCPbetaAgl->SetFillStyle(3001);
                }
                AccSelMCDbetaAgl[h]->SetMarkerStyle(8);
                AccSelMCDbetaAgl[h]->SetFillColor(4);
                AccSelMCDbetaAgl[h]->SetFillStyle(3001);
		AccSelMCDbetaAgl[h]->SetMarkerSize(2);
                AccSelMCDbetaAgl[h]->SetLineColor(4);
                AccSelMCDbetaAgl[h]->SetLineStyle(1);
		AccSelMCDbetaAgl[h]->SetLineWidth(4);
                AccSelMCDbetaAgl[h]->SetMarkerStyle(h+3);
        }
                AccSelMCDbetaAgl[1]->Draw("C4same");
		AccSelMCPbetaAgl->SetLineStyle(1);
		AccSelMCPbetaAgl->Draw("C4same");

	TCanvas *e = new TCanvas("Errors Breakdown (P)");
	e-> Divide(2,2);

	e->cd(1);
	gPad->SetGridx();
	gPad->SetGridy();
	gPad->SetLogy();
	TH1F * ErrTOTP_R=(TH1F *)EffStatErrR_P->Clone();
	for(int iR=0;iR<ErrTOTP_R->GetNbinsX();iR++) 
		ErrTOTP_R->SetBinContent(iR+1,PCorrectedAcceptance_R ->GetBinError(iR+1)/PCorrectedAcceptance_R ->GetBinContent(iR+1));
		
	ErrTOTP_R->GetXaxis()->SetTitle("Bin nr.");
	ErrTOTP_R->GetYaxis()->SetTitle("Relative Error");
	ErrTOTP_R->GetXaxis()->SetTitleSize(0.045);
	ErrTOTP_R->GetYaxis()->SetTitleSize(0.045);
	ErrTOTP_R->SetTitle("Acceptance TOTAL Error breakdown (R bins)");
	ErrTOTP_R->GetYaxis()->SetRangeUser(1e-4,1);

	ErrTOTP_R->SetLineColor(1);
	ErrTOTP_R->SetLineWidth(7);	
	LatCorrErrR_P->SetLineColor(2);
	LatCorrErrR_P->SetLineWidth(4);
	EffStatErrR_P->SetLineColor(3);   	
	EffStatErrR_P->SetLineWidth(4);  
	EffSystErrR_P->SetLineColor(4);   	
	EffSystErrR_P->SetLineWidth(4);  
	CorrStatErrR_P->SetLineColor(5);   	
	CorrStatErrR_P->SetLineWidth(4);  
	CorrSystErrR_P->SetLineColor(6);   	
	CorrSystErrR_P->SetLineWidth(4);  
	
	CorrSystErrR_P->Smooth(2);
	ErrTOTP_R->Smooth(2);
	
	ErrTOTP_R->Draw("hist");	
	LatCorrErrR_P  ->Draw("hist same"); 
	EffStatErrR_P  ->Draw("hist same"); 
	EffSystErrR_P  ->Draw("hist same"); 
	CorrStatErrR_P ->Draw("hist same"); 
	CorrSystErrR_P ->Draw("hist same"); 
	{       TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(  ErrTOTP_R       ,"TOTAL Error", "l");
                leg->AddEntry(  LatCorrErrR_P   ,"LAT. dep.  syst. Error", "l");
                leg->AddEntry(  EffStatErrR_P   ,"Efficiency stat. Error", "l");
                leg->AddEntry(  EffSystErrR_P   ,"Efficiency Fit   Error", "l");
                leg->AddEntry(  CorrStatErrR_P  ,"Eff. Corr. Fit   Error", "l");
                leg->AddEntry(  CorrSystErrR_P  ,"Eff. Corr. syst. Error", "l");
		leg->Draw("same");

        }


	e->cd(2);
	gPad->SetGridx();
	gPad->SetGridy();
	gPad->SetLogy();
	TH1F * ErrTOTP_TOF=(TH1F *)EffStatErrTOF_P->Clone();
	for(int iR=0;iR<ErrTOTP_TOF->GetNbinsX();iR++) 
		ErrTOTP_TOF->SetBinContent(iR+1,PCorrectedAcceptance_TOF ->GetBinError(iR+1)/PCorrectedAcceptance_TOF ->GetBinContent(iR+1));
	ErrTOTP_TOF->GetXaxis()->SetTitle("Bin nr.");
        ErrTOTP_TOF->GetYaxis()->SetTitle("Relative Error");
        ErrTOTP_TOF->GetXaxis()->SetTitleSize(0.045);
        ErrTOTP_TOF->GetYaxis()->SetTitleSize(0.045);
        ErrTOTP_TOF->SetTitle("Acceptance TOTAL Error breakdown (TOF range)");
	ErrTOTP_TOF->GetYaxis()->SetRangeUser(1e-4,1);


	ErrTOTP_TOF->SetLineColor(1);
	ErrTOTP_TOF->SetLineWidth(7);	
	LatCorrErrTOF_P->SetLineColor(2);
	LatCorrErrTOF_P->SetLineWidth(4);
	EffStatErrTOF_P->SetLineColor(3);   	
	EffStatErrTOF_P->SetLineWidth(4);  
	EffSystErrTOF_P->SetLineColor(4);   	
	EffSystErrTOF_P->SetLineWidth(4);  
	CorrStatErrTOF_P->SetLineColor(5);   	
	CorrStatErrTOF_P->SetLineWidth(4);  
	CorrSystErrTOF_P->SetLineColor(6);   	
	CorrSystErrTOF_P->SetLineWidth(4);  
	
	CorrSystErrTOF_P->Smooth(2);
	ErrTOTP_TOF->Smooth(2);
	
	ErrTOTP_TOF->Draw("hist");	
	LatCorrErrTOF_P->Draw("hist same"); 
	EffStatErrTOF_P->Draw("hist same"); 
	EffSystErrTOF_P->Draw("hist same"); 
	CorrStatErrTOF_P->Draw("hist same"); 
	CorrSystErrTOF_P->Draw("hist same"); 
	{       TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(  ErrTOTP_TOF       ,"TOTAL Error", "l");
                leg->AddEntry(  LatCorrErrTOF_P   ,"LAT. dep.  syst. Error", "l");
                leg->AddEntry(  EffStatErrTOF_P   ,"Efficiency stat. Error", "l");
                leg->AddEntry(  EffSystErrTOF_P   ,"Efficiency Fit   Error", "l");
                leg->AddEntry(  CorrStatErrTOF_P  ,"Eff. Corr. Fit   Error", "l");
                leg->AddEntry(  CorrSystErrTOF_P  ,"Eff. Corr. syst. Error", "l");
		leg->Draw("same");

        }


	
	e->cd(3);
	gPad->SetGridx();
	gPad->SetGridy();
	gPad->SetLogy();
	TH1F * ErrTOTP_NaF=(TH1F *)EffStatErrNaF_P->Clone();
	for(int iR=0;iR<ErrTOTP_NaF->GetNbinsX();iR++) 
		ErrTOTP_NaF->SetBinContent(iR+1,PCorrectedAcceptance_NaF ->GetBinError(iR+1)/PCorrectedAcceptance_NaF ->GetBinContent(iR+1));

	ErrTOTP_NaF->GetXaxis()->SetTitle("Bin nr.");
        ErrTOTP_NaF->GetYaxis()->SetTitle("Relative Error");
        ErrTOTP_NaF->GetXaxis()->SetTitleSize(0.045);
        ErrTOTP_NaF->GetYaxis()->SetTitleSize(0.045);
        ErrTOTP_NaF->SetTitle("Acceptance TOTAL Error breakdown (NaF range)");
	ErrTOTP_NaF->GetYaxis()->SetRangeUser(1e-4,1);



	ErrTOTP_NaF->SetLineColor(1);
	ErrTOTP_NaF->SetLineWidth(7);	
	LatCorrErrNaF_P->SetLineColor(2);
	LatCorrErrNaF_P->SetLineWidth(4);
	EffStatErrNaF_P->SetLineColor(3);   	
	EffStatErrNaF_P->SetLineWidth(4);  
	EffSystErrNaF_P->SetLineColor(4);   	
	EffSystErrNaF_P->SetLineWidth(4);  
	CorrStatErrNaF_P->SetLineColor(5);   	
	CorrStatErrNaF_P->SetLineWidth(4);  
	CorrSystErrNaF_P->SetLineColor(6);   	
	CorrSystErrNaF_P->SetLineWidth(4);  

	CorrSystErrNaF_P->Smooth(2);	
	ErrTOTP_NaF->Smooth(2);
	
	RICHStatErrNaF_P->SetLineColor(7);   	
	RICHStatErrNaF_P->SetLineWidth(4);  
	RICHSystErrNaF_P->SetLineColor(8);   	
	RICHSystErrNaF_P->SetLineWidth(4);  
	

	ErrTOTP_NaF->Draw("hist");	
	LatCorrErrNaF_P->Draw("hist same"); 
	EffStatErrNaF_P->Draw("hist same"); 
	EffSystErrNaF_P->Draw("hist same"); 
	CorrStatErrNaF_P->Draw("hist same"); 
	CorrSystErrNaF_P->Draw("hist same"); 
	RICHStatErrNaF_P->Draw("hist same"); 
	RICHSystErrNaF_P->Draw("hist same"); 
	{       TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(  ErrTOTP_NaF       ,"TOTAL Error", "l");
                leg->AddEntry(  LatCorrErrNaF_P   ,"LAT. dep.  syst. Error", "l");
                leg->AddEntry(  EffStatErrNaF_P   ,"Efficiency stat. Error", "l");
                leg->AddEntry(  EffSystErrNaF_P   ,"Efficiency Fit   Error", "l");
                leg->AddEntry(  CorrStatErrNaF_P  ,"Eff. Corr. Fit   Error", "l");
                leg->AddEntry(  CorrSystErrNaF_P  ,"Eff. Corr. syst. Error", "l");
		leg->AddEntry(  RICHStatErrNaF_P  ,"RICH eff.  Fit   Error", "l");
                leg->AddEntry(  RICHSystErrNaF_P  ,"RICH eff.  syst. Error", "l");
        


	leg->Draw("same");
        }






	e->cd(4);
	gPad->SetGridx();
	gPad->SetGridy();
	gPad->SetLogy();
	TH1F * ErrTOTP_Agl=(TH1F *)EffStatErrAgl_P->Clone();
	for(int iR=0;iR<ErrTOTP_Agl->GetNbinsX();iR++) 
		ErrTOTP_Agl->SetBinContent(iR+1,PCorrectedAcceptance_Agl ->GetBinError(iR+1)/PCorrectedAcceptance_Agl ->GetBinContent(iR+1));
	
	ErrTOTP_Agl->GetXaxis()->SetTitle("Bin nr.");
        ErrTOTP_Agl->GetYaxis()->SetTitle("Relative Error");
        ErrTOTP_Agl->GetXaxis()->SetTitleSize(0.045);
        ErrTOTP_Agl->GetYaxis()->SetTitleSize(0.045);
        ErrTOTP_Agl->SetTitle("Acceptance TOTAL Error breakdown (Agl range)");
	ErrTOTP_NaF->GetYaxis()->SetRangeUser(1e-4,1);

	ErrTOTP_Agl->SetLineColor(1);
	ErrTOTP_Agl->SetLineWidth(4);	
	LatCorrErrAgl_P->SetLineColor(2);
	LatCorrErrAgl_P->SetLineWidth(4);
	EffStatErrAgl_P->SetLineColor(3);   	
	EffStatErrAgl_P->SetLineWidth(4);  
	EffSystErrAgl_P->SetLineColor(4);   	
	EffSystErrAgl_P->SetLineWidth(4);  
	CorrStatErrAgl_P->SetLineColor(5);   	
	CorrStatErrAgl_P->SetLineWidth(4);  
	CorrSystErrAgl_P->SetLineColor(6);   	
	CorrSystErrAgl_P->SetLineWidth(4);  
	RICHStatErrAgl_P->SetLineColor(7);   	
	RICHStatErrAgl_P->SetLineWidth(4);  
	RICHSystErrAgl_P->SetLineColor(8);   	
	RICHSystErrAgl_P->SetLineWidth(4);  
	
	CorrSystErrAgl_P->Smooth(2);
	ErrTOTP_Agl->Smooth(2);

	ErrTOTP_Agl->Draw("hist");	
	LatCorrErrAgl_P->Draw("hist same"); 
	EffStatErrAgl_P->Draw("hist same"); 
	EffSystErrAgl_P->Draw("hist same"); 
	CorrStatErrAgl_P->Draw("hist same"); 
	CorrSystErrAgl_P->Draw("hist same"); 
	RICHStatErrAgl_P->Draw("hist same"); 
	RICHSystErrAgl_P->Draw("hist same"); 



	{       TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(  ErrTOTP_Agl       ,"TOTAL Error", "l");
                leg->AddEntry(  LatCorrErrAgl_P   ,"LAT. dep.  syst. Error", "l");
                leg->AddEntry(  EffStatErrAgl_P   ,"Efficiency stat. Error", "l");
                leg->AddEntry(  EffSystErrAgl_P   ,"Efficiency Fit   Error", "l");
                leg->AddEntry(  CorrStatErrAgl_P  ,"Eff. Corr. Fit   Error", "l");
                leg->AddEntry(  CorrSystErrAgl_P  ,"Eff. Corr. syst. Error", "l");
        	leg->AddEntry(  RICHStatErrAgl_P  ,"RICH eff.  Fit   Error", "l");
                leg->AddEntry(  RICHSystErrAgl_P  ,"RICH eff.  syst. Error", "l");
        


	        leg->Draw("same");

        }



	TCanvas *e1 = new TCanvas("Errors Breakdown (D)");
	e1-> Divide(3,1);

	e1->cd(1);
	gPad->SetGridx();
	gPad->SetGridy();
	gPad->SetLogy();
	TH1F * ErrTOTD_TOF=(TH1F *)EffStatErrTOF_P->Clone();
	for(int iR=0;iR<ErrTOTD_TOF->GetNbinsX();iR++) 
		ErrTOTD_TOF->SetBinContent(iR+1, DCorrectedAcceptance_TOF->GetBinError(iR+1,1)/DCorrectedAcceptance_TOF->GetBinContent(iR+1,1));
	ErrTOTD_TOF->GetXaxis()->SetTitle("Bin nr.");
        ErrTOTD_TOF->GetYaxis()->SetTitle("Relative Error");
        ErrTOTD_TOF->GetXaxis()->SetTitleSize(0.045);
        ErrTOTD_TOF->GetYaxis()->SetTitleSize(0.045);
        ErrTOTD_TOF->SetTitle("Acceptance TOTAL Error breakdown (TOF bins)");
        ErrTOTD_TOF->GetYaxis()->SetRangeUser(1e-4,1);
	ErrTOTD_TOF->SetLineColor(1);
	ErrTOTD_TOF->SetLineWidth(7);
	LatCorrErrTOF_P->SetLineColor(2);
        LatCorrErrTOF_P->SetLineWidth(4);
        EffStatErrTOF_D->SetLineColor(3);
        EffStatErrTOF_D->SetLineWidth(4);
        EffSystErrTOF_D->SetLineColor(4);
        EffSystErrTOF_D->SetLineWidth(4);
        CorrStatErrTOF_P->SetLineColor(5);
        CorrStatErrTOF_P->SetLineWidth(4);
        CorrSystErrTOF_P->SetLineColor(6);
        CorrSystErrTOF_P->SetLineWidth(4);

	ErrTOTD_TOF->Smooth(2);
	

	ErrTOTD_TOF->Draw("hist");
        LatCorrErrTOF_P->Draw("hist same");
        EffStatErrTOF_D->Draw("hist same");
        EffSystErrTOF_D->Draw("hist same");
        CorrStatErrTOF_P->Draw("hist same");
        CorrSystErrTOF_P->Draw("hist same");
	{       TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(  ErrTOTD_TOF       ,"TOTAL Error", "l");
                leg->AddEntry(  LatCorrErrTOF_P   ,"LAT. dep.  syst. Error", "l");
                leg->AddEntry(  EffStatErrTOF_D   ,"Efficiency stat. Error", "l");
                leg->AddEntry(  EffSystErrTOF_D   ,"Efficiency Fit   Error", "l");
                leg->AddEntry(  CorrStatErrTOF_P  ,"Eff. Corr. Fit   Error", "l");
                leg->AddEntry(  CorrSystErrTOF_P  ,"Eff. Corr. syst. Error", "l");
                leg->Draw("same");

        }


	e1->cd(2);
	gPad->SetGridx();
	gPad->SetGridy();
	gPad->SetLogy();
	TH1F * ErrTOTD_NaF=(TH1F *)EffStatErrNaF_P->Clone();
	for(int iR=0;iR<ErrTOTD_NaF->GetNbinsX();iR++) 
		ErrTOTD_NaF->SetBinContent(iR+1, DCorrectedAcceptance_NaF->GetBinError(iR+1,1)/DCorrectedAcceptance_NaF->GetBinContent(iR+1,1));
	ErrTOTD_NaF->GetXaxis()->SetTitle("Bin nr.");
        ErrTOTD_NaF->GetYaxis()->SetTitle("Relative Error");
        ErrTOTD_NaF->GetXaxis()->SetTitleSize(0.045);
        ErrTOTD_NaF->GetYaxis()->SetTitleSize(0.045);
        ErrTOTD_NaF->SetTitle("Acceptance TOTAL Error breakdown (NaF bins)");
        ErrTOTD_NaF->GetYaxis()->SetRangeUser(1e-4,1);
	ErrTOTD_NaF->SetLineColor(1);
	ErrTOTD_NaF->SetLineWidth(7);
	LatCorrErrNaF_P->SetLineColor(2);
        LatCorrErrNaF_P->SetLineWidth(4);
        EffStatErrNaF_D->SetLineColor(3);
        EffStatErrNaF_D->SetLineWidth(4);
        EffSystErrNaF_D->SetLineColor(4);
        EffSystErrNaF_D->SetLineWidth(4);
        CorrStatErrNaF_P->SetLineColor(5);
        CorrStatErrNaF_P->SetLineWidth(4);
        CorrSystErrNaF_P->SetLineColor(6);
        CorrSystErrNaF_P->SetLineWidth(4);
	RICHStatErrNaF_P->SetLineColor(7);
        RICHStatErrNaF_P->SetLineWidth(4);
        RICHSystErrNaF_P->SetLineColor(8);
        RICHSystErrNaF_P->SetLineWidth(4);

	ErrTOTD_NaF->Smooth(2);
	


	ErrTOTD_NaF->Draw("hist");
        LatCorrErrNaF_P->Draw("hist same");
        EffStatErrNaF_D->Draw("hist same");
        EffSystErrNaF_D->Draw("hist same");
        CorrStatErrNaF_P->Draw("hist same");
        CorrSystErrNaF_P->Draw("hist same");
	RICHStatErrNaF_P->Draw("hist same");
        RICHSystErrNaF_P->Draw("hist same");
	
	{       TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(  ErrTOTD_NaF       ,"TOTAL Error", "l");
                leg->AddEntry(  LatCorrErrNaF_P   ,"LAT. dep.  syst. Error", "l");
                leg->AddEntry(  EffStatErrNaF_D   ,"Efficiency stat. Error", "l");
                leg->AddEntry(  EffSystErrNaF_D   ,"Efficiency Fit   Error", "l");
                leg->AddEntry(  CorrStatErrNaF_P  ,"Eff. Corr. Fit   Error", "l");
                leg->AddEntry(  CorrSystErrNaF_P  ,"Eff. Corr. syst. Error", "l");
                leg->AddEntry(  RICHStatErrNaF_P  ,"RICH eff.  Fit   Error", "l");
                leg->AddEntry(  RICHSystErrNaF_P  ,"RICH eff.  syst. Error", "l");

		leg->Draw("same");

        }

	e1->cd(3);
        gPad->SetGridx();
        gPad->SetGridy();
        gPad->SetLogy();
        TH1F * ErrTOTD_Agl=(TH1F *)EffStatErrAgl_P->Clone();
        for(int iR=0;iR<ErrTOTD_Agl->GetNbinsX();iR++)
                ErrTOTD_Agl->SetBinContent(iR+1, DCorrectedAcceptance_Agl->GetBinError(iR+1,1)/DCorrectedAcceptance_Agl->GetBinContent(iR+1,1));
        ErrTOTD_Agl->GetXaxis()->SetTitle("Bin nr.");
        ErrTOTD_Agl->GetYaxis()->SetTitle("Relative Error");
        ErrTOTD_Agl->GetXaxis()->SetTitleSize(0.045);
        ErrTOTD_Agl->GetYaxis()->SetTitleSize(0.045);
        ErrTOTD_Agl->SetTitle("Acceptance TOTAL Error breakdown (Agl bins)");
        ErrTOTD_Agl->GetYaxis()->SetRangeUser(1e-4,1);
        ErrTOTD_Agl->SetLineColor(1);
        ErrTOTD_Agl->SetLineWidth(7);
        LatCorrErrAgl_P->SetLineColor(2);
        LatCorrErrAgl_P->SetLineWidth(4);
        EffStatErrAgl_D->SetLineColor(3);
        EffStatErrAgl_D->SetLineWidth(4);
        EffSystErrAgl_D->SetLineColor(4);
        EffSystErrAgl_D->SetLineWidth(4);
        CorrStatErrAgl_P->SetLineColor(5);
        CorrStatErrAgl_P->SetLineWidth(4);
        CorrSystErrAgl_P->SetLineColor(6);
        CorrSystErrAgl_P->SetLineWidth(4);
        RICHStatErrAgl_P->SetLineColor(7);
        RICHStatErrAgl_P->SetLineWidth(4);
        RICHSystErrAgl_P->SetLineColor(8);
        RICHSystErrAgl_P->SetLineWidth(4);

	ErrTOTD_Agl->Smooth(2);

        ErrTOTD_Agl->Draw("hist");
        LatCorrErrAgl_P->Draw("hist same");
        EffStatErrAgl_D->Draw("hist same");
        EffSystErrAgl_D->Draw("hist same");
        CorrStatErrAgl_P->Draw("hist same");
        CorrSystErrAgl_P->Draw("hist same");
        RICHStatErrNaF_P->Draw("hist same");
	RICHSystErrNaF_P->Draw("hist same");
	{       TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(  ErrTOTD_Agl       ,"TOTAL Error", "l");
                leg->AddEntry(  LatCorrErrAgl_P   ,"LAT. dep.  syst. Error", "l");
                leg->AddEntry(  EffStatErrAgl_D   ,"Efficiency stat. Error", "l");
                leg->AddEntry(  EffSystErrAgl_D   ,"Efficiency Fit   Error", "l");
                leg->AddEntry(  CorrStatErrAgl_P  ,"Eff. Corr. Fit   Error", "l");
                leg->AddEntry(  CorrSystErrAgl_P  ,"Eff. Corr. syst. Error", "l");
                leg->AddEntry(  RICHStatErrAgl_P  ,"RICH eff.  Fit   Error", "l");
                leg->AddEntry(  RICHSystErrAgl_P  ,"RICH eff.  syst. Error", "l");

                leg->Draw("same");

        }




	
	finalPlots.Add(c31_tris );
	finalPlots.Add(c22 );
	finalPlots.Add(c31_bis  );
	finalPlots.Add(e  );
	finalPlots.Add(e1 );
	finalPlots.writeObjsInFolder("Acceptance");

	return;
}
