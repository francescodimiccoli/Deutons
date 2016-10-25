
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

		for(int h=3;h<4;h++){
			EffFullsetMCD_R[h]= new TGraphErrors();
			EffFullsetMCD_R_fit[h]= new TGraphErrors();

			EffFullsetMCD_R[h]->SetTitle(MCLegend[h+1].c_str());
			for(int i=1;i<nbinsr;i++) EffFullsetMCD_R[h]->SetPoint(i,PRB.RigBinCent(i),EffFullsetMCD_R_TH2F->GetBinContent(i+1,h+1));
			for(int i=1;i<nbinsr;i++) EffFullsetMCD_R[h]->SetPointError(i,0,EffFullsetMCD_R_TH2F->GetBinError(i+1,h+1));
	
			EffFullsetMCD_R[h]->SetTitle("Deuterons Fit");
                        for(int i=1;i<nbinsr;i++) EffFullsetMCD_R_fit[h]->SetPoint(i,PRB.RigBinCent(i),EffFullsetMCD_R_Fit->GetBinContent(i+1,h+1));
                        for(int i=1;i<nbinsr;i++) EffFullsetMCD_R_fit[h]->SetPointError(i,0,EffFullsetMCD_R_Fit->GetBinError(i+1,h+1));

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

	
   finalPlots.Add(c41);
   finalPlots.Add(c41_bis);
   
   finalPlots.writeObjsInFolder("MC Results/Full-set selections");

}
