
void MCpreeff_Plot(

		TH1 * EffPreMCP_R_TH1F , 
		TH1 * EffPreMCP_TH1F   , 
		TH1 * EffPreMCPNaF_TH1F, 
		TH1 * EffPreMCPAgl_TH1F, 
		TH1 * EffPreMCD_R_TH2F , 
		TH1 * EffPreMCD_TH2F   , 
		TH1 * EffPreMCDNaF_TH2F, 
		TH1 * EffPreMCDAgl_TH2F, 
		TH1F * EffCascade[],
		TH1F * EffCascadeRICH[],
		TH1F * Effpre ,
		TH1F * Effagl ,
		TH1F * Effnaf, 
		TH1F * EffRInner, 
		TH1F * EffR_L1   

		){

	TCanvas *c4=new TCanvas("Preselections Efficiency (R bins)");
	TCanvas *c4_bis=new TCanvas("Preselections Efficiency (Beta bins)");
	TCanvas *c5 = new TCanvas("Preselection Cascade");
	TCanvas *c6 = new TCanvas("RICH Cascade");
	TCanvas *c7 = new TCanvas("RICH Fullset");
	TCanvas *c8 = new TCanvas("L1+Inner vs Inner");	

	c4_bis->Divide(3,1);

	c4->cd();
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();

	string MCLegend[7]= {"Protons MC B800","Deuteorons MC \"GG_Blic\"","Deuterons MC \"GG_BlicDPMJet\"","Deuterons MC \"GG_QMD\"","Deuterons MC \"Shen_Blic\"","Deuterons MC \"Shen_BlicDPMJet\"","Deuterons MC \"Shen_QMD\""};
	
	TGraph * EffPreMCP_R = new TGraph();
	EffPreMCP_R->SetTitle(MCLegend[0].c_str());
	for(int i=0; i<nbinsr; i++) EffPreMCP_R->SetPoint(i,PRB.RigBinCent(i),EffPreMCP_R_TH1F->GetBinContent(i+1));
	TGraph * EffPreMCD_R[6];
	EffPreMCP_R->SetMarkerColor(2);
	EffPreMCP_R->SetMarkerStyle(8);
	EffPreMCP_R->SetLineColor(2);
	EffPreMCP_R->SetLineWidth(2);
	EffPreMCP_R->SetTitle("Preselections Efficiency MC (R bins)");
	EffPreMCP_R->GetXaxis()->SetTitle("R [GV]");
	EffPreMCP_R->GetYaxis()->SetTitle("Pres. Efficiency");
	EffPreMCP_R->GetXaxis()->SetTitleSize(0.045);
	EffPreMCP_R->GetYaxis()->SetTitleSize(0.045);
	{
		EffPreMCP_R->Draw("ACP");
		TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
		leg->AddEntry(EffPreMCP_R,MCLegend[0].c_str(), "ep");

		for(int h=0; h<6; h++) {
			EffPreMCD_R[h]= new TGraph();
			EffPreMCD_R[h]->SetTitle(MCLegend[h+1].c_str());
			for(int i=0; i<nbinsr; i++) EffPreMCD_R[h]->SetPoint(i,PRB.RigBinCent(i),EffPreMCD_R_TH2F->GetBinContent(i+1,h+1));
			leg->AddEntry(EffPreMCD_R[h],MCLegend[h+1].c_str(), "ep");
			EffPreMCD_R[h]->SetMarkerColor(4);
			EffPreMCD_R[h]->SetMarkerStyle(h+3);
			EffPreMCD_R[h]->SetMarkerSize(2);
			EffPreMCD_R[h]->SetLineColor(4);
			EffPreMCD_R[h]->SetLineWidth(2);
			EffPreMCD_R[h]->Draw("Psame");
			leg->Draw();
		}
	}

        int plottingstyles[6]={3,4,20,25,29,26};
	
	c4_bis->cd(1);
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraph * EffPreMCP = new TGraph();
	for(int i=0; i<nbinsToF; i++) EffPreMCP->SetPoint(i,ToFPB.EkPerMassBinCent(i),EffPreMCP_TH1F->GetBinContent(i+1));
	TGraph * EffPreMCD[6];
	EffPreMCP->SetMarkerColor(2);
	EffPreMCP->SetMarkerStyle(8);
	EffPreMCP->SetMarkerSize(2);
	EffPreMCP->SetLineColor(2);
	EffPreMCP->SetLineWidth(2);
	EffPreMCP->SetTitle("Preselections Efficiency MC (Beta bins)");
	EffPreMCP->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
	EffPreMCP->GetYaxis()->SetTitle("Pres. Efficiency");
	EffPreMCP->GetXaxis()->SetTitleSize(0.045);
	EffPreMCP->GetYaxis()->SetTitleSize(0.045);
	{
		EffPreMCP->Draw("AP");
		TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
		leg->AddEntry(EffPreMCP,MCLegend[0].c_str(), "ep");

		for(int h=0; h<6; h++) {
			EffPreMCD[h]= new TGraph();
			for(int i=0; i<nbinsToF; i++) EffPreMCD[h]->SetPoint(i,ToFPB.EkPerMassBinCent(i), EffPreMCD_TH2F->GetBinContent(i+1,h+1));
			EffPreMCD[h]->SetMarkerColor(4);
			EffPreMCD[h]->SetMarkerStyle(plottingstyles[h]);
			leg->AddEntry(EffPreMCD[h],MCLegend[h+1].c_str(), "p");
			EffPreMCD[h]->SetMarkerSize(2);
			EffPreMCD[h]->SetLineColor(4);
			EffPreMCD[h]->SetLineWidth(2);
			EffPreMCD[h]->Draw("Psame");
			leg->SetLineWidth(2);
			leg->Draw();
		}
	}

	c4_bis->cd(2);
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraph * EffPreMCPNaF = new TGraph();
	for(int i=0; i<nbinsNaF; i++) EffPreMCPNaF->SetPoint(i,NaFPB.EkPerMassBinCent(i),EffPreMCPNaF_TH1F->GetBinContent(i+1));
	TGraph * EffPreMCDNaF[6];
	EffPreMCPNaF->SetMarkerColor(2);
	EffPreMCPNaF->SetMarkerStyle(8);
	EffPreMCPNaF->SetMarkerSize(2);
	EffPreMCPNaF->SetLineColor(2);
	EffPreMCPNaF->SetLineWidth(2);
	EffPreMCPNaF->SetTitle("Preselections Efficiency MC (Beta bins NaF)");
	EffPreMCPNaF->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
	EffPreMCPNaF->GetYaxis()->SetTitle("Pres. Efficiency");
	EffPreMCPNaF->GetXaxis()->SetTitleSize(0.045);
	EffPreMCPNaF->GetYaxis()->SetTitleSize(0.045);
	{
		EffPreMCPNaF->Draw("AP");
		TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
		leg->AddEntry(EffPreMCPNaF,MCLegend[0].c_str(), "p");

		for(int h=0; h<6; h++) {
			EffPreMCDNaF[h]= new TGraph();
			for(int i=0; i<nbinsNaF; i++) EffPreMCDNaF[h]->SetPoint(i,NaFPB.EkPerMassBinCent(i),
					EffPreMCDNaF_TH2F->GetBinContent(i+1,h+1));
			EffPreMCDNaF[h]->SetMarkerColor(4);
			EffPreMCDNaF[h]->SetMarkerStyle(plottingstyles[h]);
			leg->AddEntry(EffPreMCDNaF[h],MCLegend[h+1].c_str(), "p");
			EffPreMCDNaF[h]->SetMarkerSize(2);
			EffPreMCDNaF[h]->SetLineColor(4);
			EffPreMCDNaF[h]->SetLineWidth(2);
			EffPreMCDNaF[h]->Draw("Psame");
			leg->SetLineWidth(2);
			leg->Draw();
		}
	}

	c4_bis->cd(3);
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraph * EffPreMCPAgl = new TGraph();
	for(int i=0; i<nbinsAgl; i++) EffPreMCPAgl->SetPoint(i,AglPB.EkPerMassBinCent(i),EffPreMCPAgl_TH1F->GetBinContent(i+1));
	TGraph * EffPreMCDAgl[6];
	EffPreMCPAgl->SetMarkerColor(2);
	EffPreMCPAgl->SetMarkerStyle(8);
	EffPreMCPAgl->SetMarkerSize(2);
	EffPreMCPAgl->SetLineColor(2);
	EffPreMCPAgl->SetLineWidth(2);
	EffPreMCPAgl->SetTitle("Preselections Efficiency MC (Beta bins Agl)");
	EffPreMCPAgl->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
	EffPreMCPAgl->GetYaxis()->SetTitle("Pres. Efficiency");
	EffPreMCPAgl->GetXaxis()->SetTitleSize(0.045);
	EffPreMCPAgl->GetYaxis()->SetTitleSize(0.045);
	{
		EffPreMCPAgl->Draw("AP");
		TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
		leg->AddEntry(EffPreMCPAgl,MCLegend[0].c_str(), "p");

		for(int h=0; h<6; h++) {
			EffPreMCDAgl[h]= new TGraph();
			for(int i=0; i<nbinsAgl; i++) EffPreMCDAgl[h]->SetPoint(i,AglPB.EkPerMassBinCent(i), EffPreMCDAgl_TH2F->GetBinContent(i+1,h+1));
			EffPreMCDAgl[h]->SetMarkerColor(4);
			EffPreMCDAgl[h]->SetMarkerStyle(plottingstyles[h]);
			leg->AddEntry(EffPreMCDAgl[h],MCLegend[h+1].c_str(), "p");
			EffPreMCDAgl[h]->SetMarkerSize(2);
			EffPreMCDAgl[h]->SetLineColor(4);
			EffPreMCDAgl[h]->SetLineWidth(2);
			EffPreMCDAgl[h]->Draw("Psame");
			leg->SetLineWidth(2);
			leg->Draw();
		}
	}

	c5->cd();
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraphErrors * Cascade[5];
	for (int i=0;i<5;i++){
		Cascade[i]=new TGraphErrors;
	}			
	for (int i=0;i<5;i++){
		for(int j=0; j<nbinsr; j++){
			Cascade[i]->SetPoint(j,PRB.RigBinCent(j),EffCascade[i]->GetBinContent(j+1));
			Cascade[i]->SetPointError(j,0,EffCascade[i]->GetBinError(j+1));	
		}
		Cascade[i]->SetLineColor(i+2);	
		Cascade[i]->SetLineWidth(3);
		Cascade[i]->SetMarkerStyle(8);
		Cascade[i]->SetMarkerColor(i+2);
	}	
	Cascade[0]->Draw("APC");
	Cascade[0]->GetXaxis()->SetTitle("Momentum [GeV/c]");
	Cascade[0]->GetYaxis()->SetTitle("Efficiency");	   
	Cascade[0]->GetXaxis()->SetTitleSize(0.045);	  
	Cascade[0]->GetYaxis()->SetTitleSize(0.045);

	Cascade[0]->GetYaxis()->SetRangeUser(0,1);
	for (int i=1;i<5;i++) Cascade[i]->Draw("PCsame");

	{
		TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
		leg->AddEntry(Cascade[0],"Minimum Bias TOF", "ep");
		leg->AddEntry(Cascade[1],"R Exists", "ep");
		leg->AddEntry(Cascade[2],"TOF-Tracker Match", "ep"); 
		leg->AddEntry(Cascade[3],"Good R Chisquare", "ep");
		leg->AddEntry(Cascade[4],"Only One Track", "ep");
		leg->Draw("same");
	}

	c6->cd();
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraphErrors * RICHCascade[4];
	for (int i=0;i<4;i++){
		RICHCascade[i]=new TGraphErrors;
	}
	for (int i=0;i<4;i++){
		for(int j=0; j<nbinsr; j++){
			RICHCascade[i]->SetPoint(j,PRB.RigBinCent(j),EffCascadeRICH[i]->GetBinContent(j+1));
			RICHCascade[i]->SetPointError(j,0,EffCascadeRICH[i]->GetBinError(j+1));
		}
		RICHCascade[i]->SetLineColor(i+2);
		RICHCascade[i]->SetLineWidth(3);
		RICHCascade[i]->SetMarkerStyle(8);
		RICHCascade[i]->SetMarkerColor(i+2);
	}
	RICHCascade[0]->Draw("APC");
	RICHCascade[0]->GetXaxis()->SetTitle("Momentum [GeV/c]");
	RICHCascade[0]->GetYaxis()->SetTitle("Efficiency");
	RICHCascade[0]->GetXaxis()->SetTitleSize(0.045);
	RICHCascade[0]->GetYaxis()->SetTitleSize(0.045);

	RICHCascade[0]->GetYaxis()->SetRangeUser(0,1);
	for (int i=1;i<4;i++) RICHCascade[i]->Draw("PCsame");

	{
		TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
		leg->AddEntry(RICHCascade[0],"Good RICH Ring", "ep");
		leg->AddEntry(RICHCascade[1],"PMTs > 3", "ep");
		leg->AddEntry(RICHCascade[2],"Photoelectrons", "ep");
		leg->AddEntry(RICHCascade[3],"Beta Consistency", "ep");
		leg->Draw("same");
	}

	c7->cd();
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();

	TGraphErrors * preEff=new TGraphErrors();
	TGraphErrors * nafEff=new TGraphErrors();
	TGraphErrors * aglEff=new TGraphErrors();		

	for(int j=0; j<nbinsr; j++){
		preEff->SetPoint(j,PRB.RigBinCent(j),Effpre->GetBinContent(j+1));
		preEff->SetPointError(j,0,Effpre->GetBinError(j+1));
		nafEff->SetPoint(j,PRB.RigBinCent(j),Effnaf->GetBinContent(j+1));
		nafEff->SetPointError(j,0,Effnaf->GetBinError(j+1));
		aglEff->SetPoint(j,PRB.RigBinCent(j),Effagl->GetBinContent(j+1));
		aglEff->SetPointError(j,0,Effagl->GetBinError(j+1));

	}

	preEff->SetMarkerColor(2);
	preEff->SetLineColor(2);
	preEff->SetMarkerStyle(8);
	preEff->SetLineWidth(2);
	nafEff->SetMarkerColor(3);
	nafEff->SetLineColor(3);
	nafEff->SetMarkerStyle(8);
	nafEff->SetLineWidth(2);
	aglEff->SetMarkerColor(4);
	aglEff->SetLineColor(4);
	aglEff->SetMarkerStyle(8);
	aglEff->SetLineWidth(2);


	preEff->GetXaxis()->SetTitle("Momentum [GeV/c]");
	preEff->GetYaxis()->SetTitle("Efficiency");


	preEff->Draw("APC");
	nafEff->Draw("PCsame");		
	aglEff->Draw("PCsame");

	{
		TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95); 
		leg->AddEntry(preEff,"Clean-Event sel.", "ep");
		leg->AddEntry(nafEff,"Clean-Event sel. + RICH NaF", "ep");
		leg->AddEntry(aglEff,"Clean-Event sel. + RICH Agl", "ep");
		leg->Draw("same");
	}

	c8->cd();
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	gPad->SetLogy();
	TGraphErrors * RInner=new TGraphErrors();
	TGraphErrors * RL1=new TGraphErrors();
	for(int j=0; j<nbinsr; j++){
                RInner->SetPoint(j,PRB.RigBinCent(j),EffRInner->GetBinContent(j+1));
                RInner->SetPointError(j,0,EffRInner->GetBinError(j+1));
                RL1->SetPoint(j,PRB.RigBinCent(j),EffR_L1->GetBinContent(j+1));
                RL1->SetPointError(j,0,EffR_L1->GetBinError(j+1));
	}	

	RInner->SetMarkerColor(2);
	RInner->SetLineColor(2);	
	RInner->SetLineWidth(2);
	RInner->SetMarkerStyle(8);
	RL1->SetMarkerColor(2);
        RL1->SetLineColor(2);
        RL1->SetLineWidth(2);
        RL1->SetMarkerStyle(4);
	RInner->GetYaxis()->SetRangeUser(0.001,1.2);	
	RInner->GetXaxis()->SetTitle("Momentum [GeV/c]");
	RInner->GetYaxis()->SetTitle("Efficiency");

	RInner->Draw("APC");
	RL1   ->Draw("PCsame");

	{
		TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95); 
		leg->AddEntry(RInner,"Inner Tracker Fit (Protons MC)", "p");
		leg->AddEntry(RL1   ,"Inner Tracker + L1 Fit (Protons MC)", "p");
		leg->Draw("same");
	}



	finalPlots.Add(c4);
	finalPlots.Add(c4_bis);
	finalPlots.Add(c5);
	finalPlots.Add(c6);   
	finalPlots.Add(c7);
	finalPlots.Add(c8);

	finalPlots.writeObjsInFolder("MC Results/Preselections");



}

