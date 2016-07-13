
void 	MCFullseteff_Plot(TH1 * EffFullsetMCP_R_TH1F, 
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
	TGraph * EffFullsetMCP_R = new TGraph();
	EffFullsetMCP_R->SetTitle(MCLegend[0].c_str());
	for(int i=0;i<nbinsr;i++) EffFullsetMCP_R->SetPoint(i,PRB.RigBinCent(i),EffFullsetMCP_R_TH1F->GetBinContent(i+1));
	TGraph * EffFullsetMCD_R[6];
	EffFullsetMCP_R->SetMarkerColor(2);
	EffFullsetMCP_R->SetMarkerStyle(8);
	EffFullsetMCP_R->SetLineColor(2);
	EffFullsetMCP_R->SetLineWidth(2);
	EffFullsetMCP_R->SetTitle("Full-set selections Efficiency MC (R bins)");
	EffFullsetMCP_R->GetXaxis()->SetTitle("R [GV]");
	EffFullsetMCP_R->GetYaxis()->SetTitle("Full-set Efficiency");
	EffFullsetMCP_R->GetXaxis()->SetTitleSize(0.045);
	EffFullsetMCP_R->GetYaxis()->SetTitleSize(0.045);
	{
		EffFullsetMCP_R->Draw("ACP");
		TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
		leg->AddEntry(EffFullsetMCP_R,MCLegend[0].c_str(), "ep");

		for(int h=0;h<6;h++){
			EffFullsetMCD_R[h]= new TGraph();
			EffFullsetMCD_R[h]->SetTitle(MCLegend[h+1].c_str());
			for(int i=1;i<nbinsr;i++) EffFullsetMCD_R[h]->SetPoint(i,PRB.RigBinCent(i),EffFullsetMCD_R_TH2F->GetBinContent(i+1,h+1));
			leg->AddEntry(EffFullsetMCD_R[h],MCLegend[h+1].c_str(), "ep");
			EffFullsetMCD_R[h]->SetMarkerColor(4);
			EffFullsetMCD_R[h]->SetMarkerStyle(h+3);
			EffFullsetMCD_R[h]->SetMarkerSize(2);
			EffFullsetMCD_R[h]->SetLineColor(4);
			EffFullsetMCD_R[h]->SetLineWidth(2);
			EffFullsetMCD_R[h]->Draw("Psame");
			leg->Draw();
		}
	}

	c41_bis->cd(1);
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraph * EffFullsetMCP = new TGraph();
	for(int i=0;i<nbinsToF;i++) EffFullsetMCP->SetPoint(i,ToFPB.EkPerMassBinCent(i),EffFullsetMCP_TH1F->GetBinContent(i+1));
	TGraph * EffFullsetMCD[6];
	EffFullsetMCP->SetMarkerColor(2);
	EffFullsetMCP->SetMarkerStyle(8);
	EffFullsetMCP->SetLineColor(2);
	EffFullsetMCP->SetLineWidth(2);
	EffFullsetMCP->SetTitle("Full-set selections Efficiency MC (Beta bins)");
	EffFullsetMCP->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
	EffFullsetMCP->GetYaxis()->SetTitle("Full-set Efficiency");
	EffFullsetMCP->GetXaxis()->SetTitleSize(0.045);
	EffFullsetMCP->GetYaxis()->SetTitleSize(0.045);
	{
		EffFullsetMCP->Draw("ACP");
		TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
		leg->AddEntry(EffFullsetMCP,MCLegend[0].c_str(), "ep");

		for(int h=0;h<6;h++){
			EffFullsetMCD[h]= new TGraph();
			for(int i=0;i<nbinsToF;i++) EffFullsetMCD[h]->SetPoint(i,ToFPB.EkPerMassBinCent(i), EffFullsetMCD_TH2F->GetBinContent(i+1,h+1));
			EffFullsetMCD[h]->SetMarkerColor(4);
			EffFullsetMCD[h]->SetMarkerStyle(h+3);
			leg->AddEntry(EffFullsetMCD[h],MCLegend[h+1].c_str(), "ep");
			EffFullsetMCD[h]->SetMarkerSize(2);
			EffFullsetMCD[h]->SetLineColor(4);
			EffFullsetMCD[h]->SetLineWidth(2);
			EffFullsetMCD[h]->Draw("Psame");
			leg->Draw();
		}
	}
	
	c41_bis->cd(2);
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
        TGraph * EffFullsetMCPNaF = new TGraph();
        for(int i=0;i<nbinsNaF;i++) EffFullsetMCPNaF->SetPoint(i,NaFPB.EkPerMassBinCent(i),EffFullsetMCPNaF_TH1F->GetBinContent(i+1));
        TGraph * EffFullsetMCDNaF[6];
        EffFullsetMCPNaF->SetMarkerColor(2);
        EffFullsetMCPNaF->SetMarkerStyle(8);
        EffFullsetMCPNaF->SetLineColor(2);
        EffFullsetMCPNaF->SetLineWidth(2);
        EffFullsetMCPNaF->SetTitle("Full-set selections Efficiency MC (Beta bins NaF)");
        EffFullsetMCPNaF->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
        EffFullsetMCPNaF->GetYaxis()->SetTitle("Full-set Efficiency");
        EffFullsetMCPNaF->GetXaxis()->SetTitleSize(0.045);
        EffFullsetMCPNaF->GetYaxis()->SetTitleSize(0.045);
        {
                EffFullsetMCPNaF->Draw("ACP");
                TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(EffFullsetMCPNaF,MCLegend[0].c_str(), "ep");

                for(int h=0;h<6;h++){
                        EffFullsetMCDNaF[h]= new TGraph();
                        for(int i=0;i<nbinsNaF;i++) EffFullsetMCDNaF[h]->SetPoint(i,NaFPB.EkPerMassBinCent(i), 
                            EffFullsetMCDNaF_TH2F->GetBinContent(i+1,h+1));
                        EffFullsetMCDNaF[h]->SetMarkerColor(4);
                        EffFullsetMCDNaF[h]->SetMarkerStyle(h+3);
                        leg->AddEntry(EffFullsetMCDNaF[h],MCLegend[h+1].c_str(), "ep");
                        EffFullsetMCDNaF[h]->SetMarkerSize(2);
                        EffFullsetMCDNaF[h]->SetLineColor(4);
                        EffFullsetMCDNaF[h]->SetLineWidth(2);
                        EffFullsetMCDNaF[h]->Draw("Psame");
                        leg->Draw();
                }
        }

	c41_bis->cd(3);
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
        TGraph * EffFullsetMCPAgl = new TGraph();
        for(int i=0;i<nbinsAgl;i++) EffFullsetMCPAgl->SetPoint(i,AglPB.EkPerMassBinCent(i),EffFullsetMCPAgl_TH1F->GetBinContent(i+1));
        TGraph * EffFullsetMCDAgl[6];
        EffFullsetMCPAgl->SetMarkerColor(2);
        EffFullsetMCPAgl->SetMarkerStyle(8);
        EffFullsetMCPAgl->SetLineColor(2);
        EffFullsetMCPAgl->SetLineWidth(2);
        EffFullsetMCPAgl->SetTitle("Full-set selections Efficiency MC (Beta bins Agl)");
        EffFullsetMCPAgl->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
        EffFullsetMCPAgl->GetYaxis()->SetTitle("Full-set Efficiency");
        EffFullsetMCPAgl->GetXaxis()->SetTitleSize(0.045);
        EffFullsetMCPAgl->GetYaxis()->SetTitleSize(0.045);
        {
                EffFullsetMCPAgl->Draw("ACP");
                TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(EffFullsetMCPAgl,MCLegend[0].c_str(), "ep");

                for(int h=0;h<6;h++){
                        EffFullsetMCDAgl[h]= new TGraph();
                        for(int i=0;i<nbinsAgl;i++) EffFullsetMCDAgl[h]->SetPoint(i,AglPB.EkPerMassBinCent(i), EffFullsetMCDAgl_TH2F->GetBinContent(i+1,h+1));
                        EffFullsetMCDAgl[h]->SetMarkerColor(4);
                        EffFullsetMCDAgl[h]->SetMarkerStyle(h+3);
                        leg->AddEntry(EffFullsetMCDAgl[h],MCLegend[h+1].c_str(), "ep");
                        EffFullsetMCDAgl[h]->SetMarkerSize(2);
                        EffFullsetMCDAgl[h]->SetLineColor(4);
                        EffFullsetMCDAgl[h]->SetLineWidth(2);
                        EffFullsetMCDAgl[h]->Draw("Psame");
                        leg->Draw();
                }
        }
	
   finalPlots.Add(c41);
   finalPlots.Add(c41_bis);
   
   finalPlots.writeObjsInFolder("MC Results/Full-set selections");

}
