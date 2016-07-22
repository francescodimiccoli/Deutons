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
                        TH1 * DCorrectedAcceptance_Agl){


	TCanvas * c31_tris = new TCanvas("Gen. Efficiency");

	c31_tris->Divide(2,1);
	c31_tris->cd(1);
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
        EffgenbetaD[h]->SetMarkerColor(4);
        EffgenbetaD[h]->SetMarkerSize(1);
        EffgenbetaD[h]->SetLineColor(4);
        EffgenbetaD[h]->SetLineWidth(1);
        EffgenbetaD[h]->SetMarkerStyle(h+3);
        EffgenbetaD[h]->SetTitle("");
        EffgenbetaD[h]->GetXaxis()->SetTitle("R [GV]");
        EffgenbetaD[h]->GetYaxis()->SetTitle("Gen. Eff.");
        EffgenbetaD[h]->GetXaxis()->SetTitleSize(0.045);
        EffgenbetaD[h]->GetYaxis()->SetTitleSize(0.045);
	}
	EffgenbetaD[0]->Draw("ACP");
        for(int h=1;h<6;h++){
                EffgenbetaD[h]->Draw("CPsame");
        }
	
	TGraphErrors * EffgenbetaDTOF[6];
        p=0;
        for(int h=0;h<6;h++){
                EffgenbetaDTOF[h]=new TGraphErrors();
                p=0;
                for(int i=0;i<nbinsToF;i++) {EffgenbetaDTOF[h]->SetPoint(p,ToFDB.EkPerMassBinCent(i),DGen_Acceptance_TOF  ->GetBinContent(i+1,h+1));p++;}
        EffgenbetaDTOF[h]->SetMarkerStyle(8);
        EffgenbetaDTOF[h]->SetMarkerColor(4);
        EffgenbetaDTOF[h]->SetMarkerSize(1.4);
        EffgenbetaDTOF[h]->SetLineColor(4);
        EffgenbetaDTOF[h]->SetLineWidth(1);
        EffgenbetaDTOF[h]->SetMarkerStyle(h+3);
        EffgenbetaDTOF[h]->SetTitle("");
        EffgenbetaDTOF[h]->GetXaxis()->SetTitle("R [GV]");
        EffgenbetaDTOF[h]->GetYaxis()->SetTitle("Gen. Eff.");
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
        EffgenbetaDNaF[h]->SetMarkerSize(1.4);
        EffgenbetaDNaF[h]->SetLineColor(4);
        EffgenbetaDNaF[h]->SetLineWidth(1);
        EffgenbetaDNaF[h]->SetMarkerStyle(h+3);
        EffgenbetaDNaF[h]->SetTitle("");
        EffgenbetaDNaF[h]->GetXaxis()->SetTitle("R [GV]");
        EffgenbetaDNaF[h]->GetYaxis()->SetTitle("Gen. Eff.");
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
        EffgenbetaDAgl[h]->SetMarkerSize(1.4);
        EffgenbetaDAgl[h]->SetLineColor(4);
        EffgenbetaDAgl[h]->SetLineWidth(1);
        EffgenbetaDAgl[h]->SetMarkerStyle(h+3);
        EffgenbetaDAgl[h]->SetTitle("");
        EffgenbetaDAgl[h]->GetXaxis()->SetTitle("R [GV]");
        EffgenbetaDAgl[h]->GetYaxis()->SetTitle("Gen. Eff.");
        EffgenbetaDAgl[h]->GetXaxis()->SetTitleSize(0.045);
        EffgenbetaDAgl[h]->GetYaxis()->SetTitleSize(0.045);
        }
        EffgenbetaDAgl[0]->Draw("CPsame");
        for(int h=1;h<6;h++){
                EffgenbetaDAgl[h]->Draw("CPsame");
        }
	
	c31_tris->cd(2);
	gPad->SetLogx();
	gPad->SetLogy();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraphErrors * EffgenbetaP;
	p =0;
	EffgenbetaP=new TGraphErrors();
	p=0;
	for(int i=0;i<nbinsr;i++) {EffgenbetaP->SetPoint(p,PRB.EkPerMassBinCent(i),PGen_Acceptance_R  ->GetBinContent(i+1));p++;}

	EffgenbetaP->SetMarkerStyle(8);
	EffgenbetaP->SetMarkerColor(2);
	EffgenbetaP->SetMarkerSize(1);
	EffgenbetaP->SetLineColor(2);
	EffgenbetaP->SetLineWidth(1);
	EffgenbetaP->SetTitle("");
	EffgenbetaP->GetXaxis()->SetTitle("R [GV]");
	EffgenbetaP->GetYaxis()->SetTitle("Gen. Eff.");
	EffgenbetaP->GetXaxis()->SetTitleSize(0.045);
	EffgenbetaP->GetYaxis()->SetTitleSize(0.045);
	EffgenbetaP->Draw("ACP");

	TGraphErrors * EffgenbetaPTOF;
	p=0;
	EffgenbetaPTOF=new TGraphErrors();
	p=0;
	for(int i=0;i<nbinsToF;i++) {EffgenbetaPTOF->SetPoint(p,ToFPB.EkPerMassBinCent(i),PGen_Acceptance_TOF  ->GetBinContent(i+1));p++;}
	EffgenbetaPTOF->SetMarkerStyle(8);
	EffgenbetaPTOF->SetMarkerColor(2);
	EffgenbetaPTOF->SetMarkerSize(1.4);
	EffgenbetaPTOF->SetLineColor(2);
	EffgenbetaPTOF->SetLineWidth(1);
	EffgenbetaPTOF->SetTitle("");
	EffgenbetaPTOF->GetXaxis()->SetTitle("R [GV]");
	EffgenbetaPTOF->GetYaxis()->SetTitle("Gen. Eff.");
	EffgenbetaPTOF->GetXaxis()->SetTitleSize(0.045);
	EffgenbetaPTOF->GetYaxis()->SetTitleSize(0.045);
	EffgenbetaPTOF->Draw("CPsame");

	TGraphErrors * EffgenbetaPNaF;
	p=0;
	EffgenbetaPNaF=new TGraphErrors();
	p=0;
	for(int i=0;i<nbinsNaF;i++) {EffgenbetaPNaF->SetPoint(p,NaFPB.EkPerMassBinCent(i),PGen_Acceptance_NaF  ->GetBinContent(i+1));p++;}
	EffgenbetaPNaF->SetMarkerStyle(8);
	EffgenbetaPNaF->SetMarkerColor(2);
	EffgenbetaPNaF->SetMarkerSize(1.4);
	EffgenbetaPNaF->SetLineColor(2);
	EffgenbetaPNaF->SetLineWidth(1);
	EffgenbetaPNaF->SetTitle("");
	EffgenbetaPNaF->GetXaxis()->SetTitle("R [GV]");
	EffgenbetaPNaF->GetYaxis()->SetTitle("Gen. Eff.");
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
	EffgenbetaPAgl->SetMarkerSize(1.4);
	EffgenbetaPAgl->SetLineColor(2);
	EffgenbetaPAgl->SetLineWidth(1);
	EffgenbetaPAgl->SetTitle("");
	EffgenbetaPAgl->GetXaxis()->SetTitle("R [GV]");
	EffgenbetaPAgl->GetYaxis()->SetTitle("Gen. Eff.");
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
        AccgeoP->GetXaxis()->SetTitle("R [GV]");
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
	string MCLegend[7]={"protons.B800","d.pl1.0_520_GG_Blic","d.pl1.0_520_GG_BlicDPMJet","d.pl1.0_520_GG_QMD","d.pl1.0_520_Shen_Blic","d.pl1.0_520_Shen_BlicDPMJet","d.pl1.0_520_Shen_QMD"};
	
	c31_bis->cd();
        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetGridx();
        gPad->SetGridy();
        TGraphErrors * AccSelMCDbeta[6];
        TGraphErrors * AccSelMCPbeta=new TGraphErrors();
	for(int h=0;h<6;h++){
                AccSelMCDbeta[h]= new TGraphErrors();
                int p=0;
                for(int m=0;m<nbinsToF;m++) {AccSelMCDbeta[h]->SetPoint(p,ToFDB.EkPerMassBinCent(m),DMCAcceptance_TOF -> GetBinContent(m+1,h+1));
					    AccSelMCDbeta[h]->SetPointError(p,0,DMCAcceptance_TOF -> GetBinError(m+1,h+1));
					     p++;}
		if(h==0) 
		{
			p=0;
			for(int m=0;m<nbinsToF;m++) {AccSelMCPbeta->SetPoint(p,ToFPB.EkPerMassBinCent(m),PMCAcceptance_TOF -> GetBinContent(m+1));
						     AccSelMCPbeta->SetPointError(p,0,PMCAcceptance_TOF -> GetBinError(m+1));
						     p++;}
			AccSelMCDbeta[0]->SetPoint(p,50,0.001);
			AccSelMCPbeta->SetMarkerStyle(8);
                	AccSelMCPbeta->SetMarkerColor(2);
                	AccSelMCPbeta->SetMarkerSize(1.4);
                	AccSelMCPbeta->SetLineColor(2);
                	AccSelMCPbeta->SetLineWidth(1);
		}
		AccSelMCDbeta[h]->SetMarkerStyle(8);
                AccSelMCDbeta[h]->SetMarkerColor(4);
                AccSelMCDbeta[h]->SetMarkerSize(1.4);
                AccSelMCDbeta[h]->SetLineColor(4);
		AccSelMCDbeta[h]->SetLineStyle(2);
                AccSelMCDbeta[h]->SetLineWidth(1);
                AccSelMCDbeta[h]->SetMarkerStyle(h+3);
                AccSelMCDbeta[0]->SetTitle("Final effective Acceptance");
                AccSelMCDbeta[0]->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
                AccSelMCDbeta[0]->GetYaxis()->SetTitle("Acceptance [m^2 sr]");
                AccSelMCDbeta[0]->GetXaxis()->SetTitleSize(0.045);
                AccSelMCDbeta[0]->GetYaxis()->SetTitleSize(0.045);
                AccSelMCDbeta[0]->GetYaxis()->SetRangeUser(1e-2,1.3);
        }
	{	TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(AccSelMCDbeta[0],MCLegend[1].c_str(), "ep");
		AccSelMCDbeta[0]->Draw("AP");
		AccSelMCPbeta->SetLineStyle(2);
		AccSelMCPbeta->Draw("Psame");
        for(int h=1;h<6;h++){
                leg->AddEntry(AccSelMCDbeta[h],MCLegend[h+1].c_str(), "ep");
		AccSelMCDbeta[h]->Draw("Psame");
        }
        	leg->AddEntry(AccSelMCPbeta,"Protons B800","ep");		
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
                        AccSelMCPbetaNaF->SetMarkerColor(2);
                        AccSelMCPbetaNaF->SetMarkerSize(2);
                        AccSelMCPbetaNaF->SetLineColor(2);
                        AccSelMCPbetaNaF->SetLineWidth(1);
                }
                AccSelMCDbetaNaF[h]->SetMarkerStyle(8);
                AccSelMCDbetaNaF[h]->SetMarkerColor(4);
                AccSelMCDbetaNaF[h]->SetMarkerSize(2);
                AccSelMCDbetaNaF[h]->SetLineColor(4);
                AccSelMCDbetaNaF[h]->SetLineStyle(2);
		AccSelMCDbetaNaF[h]->SetLineWidth(1);
                AccSelMCDbetaNaF[h]->SetMarkerStyle(h+3);
                AccSelMCDbetaNaF[h]->SetTitle("Deutons Acceptance");
                AccSelMCDbetaNaF[h]->GetXaxis()->SetTitle("R [GV]");
                AccSelMCDbetaNaF[h]->GetYaxis()->SetTitle("Acceptance [m^2 sr]");
                AccSelMCDbetaNaF[h]->GetXaxis()->SetTitleSize(0.045);
                AccSelMCDbetaNaF[h]->GetYaxis()->SetTitleSize(0.045);
                AccSelMCDbetaNaF[h]->GetYaxis()->SetRangeUser(1e-2,1.3);
        }
                AccSelMCDbetaNaF[0]->Draw("Psame");
		AccSelMCPbetaNaF->SetLineStyle(2);
		AccSelMCPbetaNaF->Draw("Psame");
        for(int h=1;h<6;h++){
                AccSelMCDbetaNaF[h]->Draw("Psame");
        }
	

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
                        AccSelMCPbetaAgl->SetMarkerColor(2);
                        AccSelMCPbetaAgl->SetMarkerSize(2);
                        AccSelMCPbetaAgl->SetLineColor(2);
                        AccSelMCPbetaAgl->SetLineWidth(1);
                }
                AccSelMCDbetaAgl[h]->SetMarkerStyle(8);
                AccSelMCDbetaAgl[h]->SetMarkerColor(4);
                AccSelMCDbetaAgl[h]->SetMarkerSize(2);
                AccSelMCDbetaAgl[h]->SetLineColor(4);
                AccSelMCDbetaAgl[h]->SetLineStyle(2);
		AccSelMCDbetaAgl[h]->SetLineWidth(1);
                AccSelMCDbetaAgl[h]->SetMarkerStyle(h+3);
                AccSelMCDbetaAgl[h]->SetTitle("Deutons Acceptance");
                AccSelMCDbetaAgl[h]->GetXaxis()->SetTitle("R [GV]");
                AccSelMCDbetaAgl[h]->GetYaxis()->SetTitle("Acceptance [m^2 sr]");
                AccSelMCDbetaAgl[h]->GetXaxis()->SetTitleSize(0.045);
                AccSelMCDbetaAgl[h]->GetYaxis()->SetTitleSize(0.045);
                AccSelMCDbetaAgl[h]->GetYaxis()->SetRangeUser(1e-2,1.3);
        }
                AccSelMCDbetaAgl[0]->Draw("Psame");
		AccSelMCPbetaAgl->SetLineStyle(2);
		AccSelMCPbetaAgl->Draw("Psame");
        for(int h=1;h<6;h++){
                AccSelMCDbetaAgl[h]->Draw("Psame");
        }

	
	finalPlots.Add(c31_tris );
	finalPlots.Add(c22 );
	finalPlots.Add(c31_bis  );
	finalPlots.writeObjsInFolder("Acceptance");

	return;
}
