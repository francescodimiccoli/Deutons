
void 	MCFullseteff_Plot(
			  
			  TH1 * EffFullsetMCP_R_Fit  , 
                          TH1 * EffFullsetMCP_Fit,
                          TH1 * EffFullsetMCPNaF_Fit,
                          TH1 * EffFullsetMCPAgl_Fit,
                          TH1 * EffFullsetMCD_R_Fit  ,
                          TH1 * EffFullsetMCD_Fit,
                          TH1 * EffFullsetMCDNaF_Fit,
                          TH1 * EffFullsetMCDAgl_Fit,

			  
			  TH1 * EffFullsetMCP_R_TH1F, 
                          TH1 * EffFullsetMCP_TH1F   ,
                          TH1 * EffFullsetMCPNaF_TH1F,
                          TH1 * EffFullsetMCPAgl_TH1F,
                          TH1 * EffFullsetMCD_R_TH2F ,
                          TH1 * EffFullsetMCD_TH2F   ,
                          TH1 * EffFullsetMCDNaF_TH2F,
                          TH1 * EffFullsetMCDAgl_TH2F){


	TCanvas *c41=new TCanvas("FULL-SET Efficiency (R bins)");
	TCanvas *c41_bis=new TCanvas("FULL-SET Efficiency (Beta bins)");	
	c41_bis->Divide(3,1);

	c41->cd();
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	string MCLegend[7]={"protons.B800","d.pl1.0_520_GG_Blic","d.pl1.0_520_GG_BlicDPMJet","d.pl1.0_520_GG_QMD","d.pl1.0_520_Shen_Blic","d.pl1.0_520_Shen_BlicDPMJet","d.pl1.0_520_Shen_QMD"};
	
	TGraphErrors * EffFullsetMCP_R = new TGraphErrors();
	TGraphErrors * EffFullsetMCP_R_fit = new TGraphErrors();	

	EffFullsetMCP_R->SetTitle(MCLegend[0].c_str());
	for(int i=0;i<nbinsr;i++) EffFullsetMCP_R->SetPoint(i,PRB.RigBinCent(i),EffFullsetMCP_R_TH1F->GetBinContent(i+1));
	for(int i=0;i<nbinsr;i++) EffFullsetMCP_R->SetPointError(i,0,EffFullsetMCP_R_TH1F->GetBinError(i+1));

	EffFullsetMCP_R_fit->SetTitle("Protons Fit");
	for(int i=0;i<nbinsr;i++) EffFullsetMCP_R_fit->SetPoint(i,PRB.RigBinCent(i),EffFullsetMCP_R_Fit->GetBinContent(i+1));
	for(int i=0;i<nbinsr;i++) EffFullsetMCP_R_fit->SetPointError(i,0,EffFullsetMCP_R_Fit->GetBinError(i+1));
	
	TGraphErrors * EffFullsetMCD_R[6];
	TGraphErrors * EffFullsetMCD_R_fit[6];	

	EffFullsetMCP_R->SetMarkerColor(2);
	EffFullsetMCP_R->SetMarkerStyle(8);
	EffFullsetMCP_R->SetMarkerSize(2);
	EffFullsetMCP_R->SetLineColor(2);
	EffFullsetMCP_R->SetLineWidth(2);
	EffFullsetMCP_R_fit->SetLineColor(2);
        EffFullsetMCP_R_fit->SetLineWidth(4);
	EffFullsetMCP_R_fit->SetFillColor(2);
	EffFullsetMCP_R_fit->SetFillStyle(3001);

	EffFullsetMCP_R->SetTitle("Full-set selections Efficiency MC (R bins)");
	EffFullsetMCP_R->GetXaxis()->SetTitle("R [GV]");
	EffFullsetMCP_R->GetYaxis()->SetTitle("Full-set Efficiency");
	EffFullsetMCP_R->GetXaxis()->SetTitleSize(0.045);
	EffFullsetMCP_R->GetYaxis()->SetTitleSize(0.045);
	
	{
		EffFullsetMCP_R->Draw("AP");
		EffFullsetMCP_R_fit->Draw("C4same");
		TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
		leg->AddEntry(EffFullsetMCP_R,MCLegend[0].c_str(), "ep");
		leg->AddEntry(EffFullsetMCP_R_fit,"Protons Fit", "l");

		for(int h=0;h<1;h++){
			EffFullsetMCD_R[h]= new TGraphErrors();
			EffFullsetMCD_R_fit[h]= new TGraphErrors();

			EffFullsetMCD_R[h]->SetTitle(MCLegend[h+1].c_str());
			for(int i=0;i<nbinsr;i++) EffFullsetMCD_R[h]->SetPoint(i,PRB.RigBinCent(i),EffFullsetMCD_R_TH2F->GetBinContent(i+1,h+1));
			for(int i=0;i<nbinsr;i++) EffFullsetMCD_R[h]->SetPointError(i,0,EffFullsetMCD_R_TH2F->GetBinError(i+1,h+1));
	
			EffFullsetMCD_R[h]->SetTitle("Deuterons Fit");
                        for(int i=0;i<nbinsr;i++) EffFullsetMCD_R_fit[h]->SetPoint(i,PRB.RigBinCent(i),EffFullsetMCD_R_Fit->GetBinContent(i+1,h+1));
                        for(int i=0;i<nbinsr;i++) EffFullsetMCD_R_fit[h]->SetPointError(i,0,EffFullsetMCD_R_Fit->GetBinError(i+1,h+1));

			leg->AddEntry(EffFullsetMCD_R[h],MCLegend[h+1].c_str(), "ep");
      		        leg->AddEntry(EffFullsetMCD_R_fit[h],"Deuterons Fit", "l");

			EffFullsetMCD_R[h]->SetMarkerColor(4);
			EffFullsetMCD_R[h]->SetMarkerStyle(8);
			EffFullsetMCD_R[h]->SetMarkerSize(2);
			EffFullsetMCD_R[h]->SetLineColor(4);
			EffFullsetMCD_R[h]->SetLineWidth(2);
			EffFullsetMCD_R_fit[h]->SetLineColor(4);
			EffFullsetMCD_R_fit[h]->SetLineWidth(4);
			EffFullsetMCD_R_fit[h]->SetFillColor(4);
			EffFullsetMCD_R_fit[h]->SetFillStyle(3001);
			EffFullsetMCD_R[h]->Draw("Psame");
			EffFullsetMCD_R_fit[h]->Draw("c4same");
		
			leg->Draw();
		}
	}

	c41_bis->cd(1);
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();

	TGraphErrors * EffFullsetMCP = new TGraphErrors();
	TGraphErrors * EffFullsetMCP_fit = new TGraphErrors();	

	EffFullsetMCP->SetTitle(MCLegend[0].c_str());
	for(int i=0;i<nbinsToF;i++) EffFullsetMCP->SetPoint(i,ToFPB.EkPerMassBinCent(i),EffFullsetMCP_TH1F->GetBinContent(i+1));
	for(int i=0;i<nbinsToF;i++) EffFullsetMCP->SetPointError(i,0,EffFullsetMCP_TH1F->GetBinError(i+1));

	EffFullsetMCP_fit->SetTitle("Protons Fit");
	for(int i=0;i<nbinsToF;i++) EffFullsetMCP_fit->SetPoint(i,ToFPB.EkPerMassBinCent(i),EffFullsetMCP_Fit->GetBinContent(i+1));
	for(int i=0;i<nbinsToF;i++) EffFullsetMCP_fit->SetPointError(i,0,EffFullsetMCP_Fit->GetBinError(i+1));
	
	TGraphErrors * EffFullsetMCD[6];
	TGraphErrors * EffFullsetMCD_fit[6];	

	EffFullsetMCP->SetMarkerColor(2);
	EffFullsetMCP->SetMarkerStyle(8);
	EffFullsetMCP->SetMarkerSize(2);
	EffFullsetMCP->SetLineColor(2);
	EffFullsetMCP->SetLineWidth(2);
	EffFullsetMCP_fit->SetLineColor(2);
        EffFullsetMCP_fit->SetLineWidth(4);
	EffFullsetMCP_fit->SetFillColor(2);
	EffFullsetMCP_fit->SetFillStyle(3001);

	EffFullsetMCP->SetTitle("Full-set selections Efficiency MC (TOF range)");
	EffFullsetMCP->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
	EffFullsetMCP->GetYaxis()->SetTitle("Full-set Efficiency");
	EffFullsetMCP->GetXaxis()->SetTitleSize(0.045);
	EffFullsetMCP->GetYaxis()->SetTitleSize(0.045);
	
	{
		EffFullsetMCP->Draw("AP");
		EffFullsetMCP_fit->Draw("C4same");
		TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
		leg->AddEntry(EffFullsetMCP,MCLegend[0].c_str(), "ep");
		leg->AddEntry(EffFullsetMCP_fit,"Protons Fit", "l");

		for(int h=0;h<1;h++){
			EffFullsetMCD[h]= new TGraphErrors();
			EffFullsetMCD_fit[h]= new TGraphErrors();

			EffFullsetMCD[h]->SetTitle(MCLegend[h+1].c_str());
			for(int i=0;i<nbinsToF;i++) EffFullsetMCD[h]->SetPoint(i,ToFPB.EkPerMassBinCent(i),EffFullsetMCD_TH2F->GetBinContent(i+1,h+1));
			for(int i=0;i<nbinsToF;i++) EffFullsetMCD[h]->SetPointError(i,0,EffFullsetMCD_TH2F->GetBinError(i+1,h+1));
	
			EffFullsetMCD[h]->SetTitle("Deuterons Fit");
                        for(int i=0;i<nbinsToF;i++) EffFullsetMCD_fit[h]->SetPoint(i,ToFPB.EkPerMassBinCent(i),EffFullsetMCD_Fit->GetBinContent(i+1,h+1));
                        for(int i=0;i<nbinsToF;i++) EffFullsetMCD_fit[h]->SetPointError(i,0,EffFullsetMCD_Fit->GetBinError(i+1,h+1));

			leg->AddEntry(EffFullsetMCD[h],MCLegend[h+1].c_str(), "ep");
      		        leg->AddEntry(EffFullsetMCD_fit[h],"Deuterons Fit", "l");

			EffFullsetMCD[h]->SetMarkerColor(4);
			EffFullsetMCD[h]->SetMarkerStyle(8);
			EffFullsetMCD[h]->SetMarkerSize(2);
			EffFullsetMCD[h]->SetLineColor(4);
			EffFullsetMCD[h]->SetLineWidth(2);
			EffFullsetMCD_fit[h]->SetLineColor(4);
			EffFullsetMCD_fit[h]->SetLineWidth(4);
			EffFullsetMCD_fit[h]->SetFillColor(4);
			EffFullsetMCD_fit[h]->SetFillStyle(3001);
			EffFullsetMCD[h]->Draw("Psame");
			EffFullsetMCD_fit[h]->Draw("c4same");
		
			leg->Draw();
		}
	}


	c41_bis->cd(2);
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();

	TGraphErrors * EffFullsetMCPNaF = new TGraphErrors();
	TGraphErrors * EffFullsetMCPNaF_fit = new TGraphErrors();	

	EffFullsetMCPNaF->SetTitle(MCLegend[0].c_str());
	for(int i=0;i<nbinsNaF;i++) EffFullsetMCPNaF->SetPoint(i,NaFPB.EkPerMassBinCent(i),EffFullsetMCPNaF_TH1F->GetBinContent(i+1));
	for(int i=0;i<nbinsNaF;i++) EffFullsetMCPNaF->SetPointError(i,0,EffFullsetMCPNaF_TH1F->GetBinError(i+1));

	EffFullsetMCPNaF_fit->SetTitle("Protons Fit");
	for(int i=0;i<nbinsNaF;i++) EffFullsetMCPNaF_fit->SetPoint(i,NaFPB.EkPerMassBinCent(i),EffFullsetMCPNaF_Fit->GetBinContent(i+1));
	for(int i=0;i<nbinsNaF;i++) EffFullsetMCPNaF_fit->SetPointError(i,0,EffFullsetMCPNaF_Fit->GetBinError(i+1));
	
	TGraphErrors * EffFullsetMCDNaF[6];
	TGraphErrors * EffFullsetMCDNaF_fit[6];	

	EffFullsetMCPNaF->SetMarkerColor(2);
	EffFullsetMCPNaF->SetMarkerStyle(8);
	EffFullsetMCPNaF->SetMarkerSize(2);
	EffFullsetMCPNaF->SetLineColor(2);
	EffFullsetMCPNaF->SetLineWidth(2);
	EffFullsetMCPNaF_fit->SetLineColor(2);
        EffFullsetMCPNaF_fit->SetLineWidth(4);
	EffFullsetMCPNaF_fit->SetFillColor(2);
	EffFullsetMCPNaF_fit->SetFillStyle(3001);

	EffFullsetMCPNaF->SetTitle("Full-set selections Efficiency MC (NaF range)");
	EffFullsetMCPNaF->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
	EffFullsetMCPNaF->GetYaxis()->SetTitle("Full-set Efficiency");
	EffFullsetMCPNaF->GetXaxis()->SetTitleSize(0.045);
	EffFullsetMCPNaF->GetYaxis()->SetTitleSize(0.045);
	
	{
		EffFullsetMCPNaF->Draw("AP");
		EffFullsetMCPNaF_fit->Draw("C4same");
		TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
		leg->AddEntry(EffFullsetMCPNaF,MCLegend[0].c_str(), "ep");
		leg->AddEntry(EffFullsetMCPNaF_fit,"Protons Fit", "l");

		for(int h=0;h<1;h++){
			EffFullsetMCDNaF[h]= new TGraphErrors();
			EffFullsetMCDNaF_fit[h]= new TGraphErrors();

			EffFullsetMCDNaF[h]->SetTitle(MCLegend[h+1].c_str());
			for(int i=0;i<nbinsNaF;i++) EffFullsetMCDNaF[h]->SetPoint(i,NaFPB.EkPerMassBinCent(i),EffFullsetMCDNaF_TH2F->GetBinContent(i+1,h+1));
			for(int i=0;i<nbinsNaF;i++) EffFullsetMCDNaF[h]->SetPointError(i,0,EffFullsetMCDNaF_TH2F->GetBinError(i+1,h+1));
	
			EffFullsetMCDNaF[h]->SetTitle("Deuterons Fit");
                        for(int i=0;i<nbinsNaF;i++) EffFullsetMCDNaF_fit[h]->SetPoint(i,NaFPB.EkPerMassBinCent(i),EffFullsetMCDNaF_Fit->GetBinContent(i+1,h+1));
                        for(int i=0;i<nbinsNaF;i++) EffFullsetMCDNaF_fit[h]->SetPointError(i,0,EffFullsetMCDNaF_Fit->GetBinError(i+1,h+1));

			leg->AddEntry(EffFullsetMCDNaF[h],MCLegend[h+1].c_str(), "ep");
      		        leg->AddEntry(EffFullsetMCDNaF_fit[h],"Deuterons Fit", "l");

			EffFullsetMCDNaF[h]->SetMarkerColor(4);
			EffFullsetMCDNaF[h]->SetMarkerStyle(8);
			EffFullsetMCDNaF[h]->SetMarkerSize(2);
			EffFullsetMCDNaF[h]->SetLineColor(4);
			EffFullsetMCDNaF[h]->SetLineWidth(2);
			EffFullsetMCDNaF_fit[h]->SetLineColor(4);
			EffFullsetMCDNaF_fit[h]->SetLineWidth(4);
			EffFullsetMCDNaF_fit[h]->SetFillColor(4);
			EffFullsetMCDNaF_fit[h]->SetFillStyle(3001);
			EffFullsetMCDNaF[h]->Draw("Psame");
			EffFullsetMCDNaF_fit[h]->Draw("c4same");
		
			leg->Draw();
		}
	}


	 c41_bis->cd(3);
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();

        TGraphErrors * EffFullsetMCPAgl = new TGraphErrors();
        TGraphErrors * EffFullsetMCPAgl_fit = new TGraphErrors();

        EffFullsetMCPAgl->SetTitle(MCLegend[0].c_str());
        for(int i=0;i<nbinsNaF;i++) EffFullsetMCPAgl->SetPoint(i,AglPB.EkPerMassBinCent(i),EffFullsetMCPAgl_TH1F->GetBinContent(i+1));
        for(int i=0;i<nbinsNaF;i++) EffFullsetMCPAgl->SetPointError(i,0,EffFullsetMCPAgl_TH1F->GetBinError(i+1));

        EffFullsetMCPAgl_fit->SetTitle("Protons Fit");
        for(int i=0;i<nbinsNaF;i++) EffFullsetMCPAgl_fit->SetPoint(i,AglPB.EkPerMassBinCent(i),EffFullsetMCPAgl_Fit->GetBinContent(i+1));
        for(int i=0;i<nbinsNaF;i++) EffFullsetMCPAgl_fit->SetPointError(i,0,EffFullsetMCPAgl_Fit->GetBinError(i+1));

        TGraphErrors * EffFullsetMCDAgl[6];
        TGraphErrors * EffFullsetMCDAgl_fit[6];

        EffFullsetMCPAgl->SetMarkerColor(2);
        EffFullsetMCPAgl->SetMarkerStyle(8);
        EffFullsetMCPAgl->SetMarkerSize(2);
        EffFullsetMCPAgl->SetLineColor(2);
        EffFullsetMCPAgl->SetLineWidth(2);
        EffFullsetMCPAgl_fit->SetLineColor(2);
        EffFullsetMCPAgl_fit->SetLineWidth(4);
        EffFullsetMCPAgl_fit->SetFillColor(2);
        EffFullsetMCPAgl_fit->SetFillStyle(3001);

        EffFullsetMCPAgl->SetTitle("Full-set selections Efficiency MC (Agl range)");
        EffFullsetMCPAgl->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
        EffFullsetMCPAgl->GetYaxis()->SetTitle("Full-set Efficiency");
        EffFullsetMCPAgl->GetXaxis()->SetTitleSize(0.045);
        EffFullsetMCPAgl->GetYaxis()->SetTitleSize(0.045);

        {
                EffFullsetMCPAgl->Draw("AP");
                EffFullsetMCPAgl_fit->Draw("C4same");
                TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(EffFullsetMCPAgl,MCLegend[0].c_str(), "ep");
                leg->AddEntry(EffFullsetMCPAgl_fit,"Protons Fit", "l");

                for(int h=0;h<1;h++){
                        EffFullsetMCDAgl[h]= new TGraphErrors();
                        EffFullsetMCDAgl_fit[h]= new TGraphErrors();

                        EffFullsetMCDAgl[h]->SetTitle(MCLegend[h+1].c_str());
                        for(int i=0;i<nbinsAgl;i++) EffFullsetMCDAgl[h]->SetPoint(i,AglPB.EkPerMassBinCent(i),EffFullsetMCDAgl_TH2F->GetBinContent(i+1,h+1));
                        for(int i=0;i<nbinsAgl;i++) EffFullsetMCDAgl[h]->SetPointError(i,0,EffFullsetMCDAgl_TH2F->GetBinError(i+1,h+1));

                        EffFullsetMCDAgl[h]->SetTitle("Deuterons Fit");
                        for(int i=0;i<nbinsAgl;i++) EffFullsetMCDAgl_fit[h]->SetPoint(i,AglPB.EkPerMassBinCent(i),EffFullsetMCDAgl_Fit->GetBinContent(i+1,h+1));
                        for(int i=0;i<nbinsAgl;i++) EffFullsetMCDAgl_fit[h]->SetPointError(i,0,EffFullsetMCDAgl_Fit->GetBinError(i+1,h+1));

                        leg->AddEntry(EffFullsetMCDAgl[h],MCLegend[h+1].c_str(), "ep");
                        leg->AddEntry(EffFullsetMCDAgl_fit[h],"Deuterons Fit", "l");

                        EffFullsetMCDAgl[h]->SetMarkerColor(4);
                        EffFullsetMCDAgl[h]->SetMarkerStyle(8);
                        EffFullsetMCDAgl[h]->SetMarkerSize(2);
                        EffFullsetMCDAgl[h]->SetLineColor(4);
                        EffFullsetMCDAgl[h]->SetLineWidth(2);
                        EffFullsetMCDAgl_fit[h]->SetLineColor(4);
                        EffFullsetMCDAgl_fit[h]->SetLineWidth(4);
                        EffFullsetMCDAgl_fit[h]->SetFillColor(4);
                        EffFullsetMCDAgl_fit[h]->SetFillStyle(3001);
                        EffFullsetMCDAgl[h]->Draw("Psame");
                        EffFullsetMCDAgl_fit[h]->Draw("c4same");

                        leg->Draw();
                }
        }

	
   finalPlots.Add(c41);
   finalPlots.Add(c41_bis);
   
   finalPlots.writeObjsInFolder("MC Results/Full-set selections");

}
