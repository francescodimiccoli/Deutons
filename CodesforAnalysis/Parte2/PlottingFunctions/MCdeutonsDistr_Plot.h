void MCdeutonsDistr_Plot(
				
			TH1F * SliceMassDMCTOF[],
			TH1F * SliceMassDMCNaF[],
			TH1F * SliceMassDMCAgl[],

			TH1F * SliceLikDMCTOF[],
			TH1F * SliceLikDMCNaF[],
			TH1F * SliceLikDMCAgl[],

			TH1F * SliceDistDMCTOF[],
			TH1F * SliceDistDMCNaF[],
			TH1F * SliceDistDMCAgl[],
			
			TH2F * EffFragmMCD_TOF, 			
                       	TH2F * EffFragmMCD_NaF, 
                       	TH2F * EffFragmMCD_Agl 


				){

	TCanvas * c  = new TCanvas("Mass Distributions");	
	TCanvas * c1 = new TCanvas("Lik Distributions");	
	TCanvas * c2 = new TCanvas("Dist Distributions");	
	TCanvas * c3 = new TCanvas("Fragm. Efficiency");


	c->Divide(3,1);

	for (int mc_type=0;mc_type<6;mc_type++){		
	
		SliceMassDMCTOF[mc_type]->SetLineColor(mc_type+2);
                SliceMassDMCNaF[mc_type]->SetLineColor(mc_type+2);
                SliceMassDMCAgl[mc_type]->SetLineColor(mc_type+2);
	}

	c->cd(1);
	gPad->SetLogy();
        gPad->SetGridx();
        gPad->SetGridy();

	SliceMassDMCTOF[0]->Draw();
	 for (int mc_type=1;mc_type<6;mc_type++)  SliceMassDMCTOF[mc_type]->Draw("same");	

	c->cd(2);
	gPad->SetLogy();
        gPad->SetGridx();
        gPad->SetGridy();

	SliceMassDMCNaF[0]->Draw();
	 for (int mc_type=1;mc_type<6;mc_type++)  SliceMassDMCNaF[mc_type]->Draw("same");	

	c->cd(3);
	gPad->SetLogy();
        gPad->SetGridx();
        gPad->SetGridy();

	SliceMassDMCAgl[0]->Draw();
	for (int mc_type=1;mc_type<6;mc_type++)  SliceMassDMCAgl[mc_type]->Draw("same");	

	
	c1->Divide(3,1);

	for (int mc_type=0;mc_type<6;mc_type++){		
	
		SliceLikDMCTOF[mc_type]->SetLineColor(mc_type+2);
                SliceLikDMCNaF[mc_type]->SetLineColor(mc_type+2);
                SliceLikDMCAgl[mc_type]->SetLineColor(mc_type+2);
	}

	c1->cd(1);
	gPad->SetLogy();
        gPad->SetGridx();
        gPad->SetGridy();

	SliceLikDMCTOF[0]->Draw();
	 for (int mc_type=1;mc_type<6;mc_type++)  SliceLikDMCTOF[mc_type]->Draw("same");	

	c1->cd(2);
	gPad->SetLogy();
        gPad->SetGridx();
        gPad->SetGridy();

	SliceLikDMCNaF[0]->Draw();
	 for (int mc_type=1;mc_type<6;mc_type++)  SliceLikDMCNaF[mc_type]->Draw("same");	

	c1->cd(3);
	gPad->SetLogy();
        gPad->SetGridx();
        gPad->SetGridy();

	SliceLikDMCAgl[0]->Draw();
	for (int mc_type=1;mc_type<6;mc_type++)  SliceLikDMCAgl[mc_type]->Draw("same");	


	c2->Divide(3,1);

	for (int mc_type=0;mc_type<6;mc_type++){		
	
		SliceDistDMCTOF[mc_type]->SetLineColor(mc_type+2);
                SliceDistDMCNaF[mc_type]->SetLineColor(mc_type+2);
                SliceDistDMCAgl[mc_type]->SetLineColor(mc_type+2);
	}

	c2->cd(1);
	gPad->SetLogy();
        gPad->SetGridx();
        gPad->SetGridy();

	SliceDistDMCTOF[0]->Draw();
	 for (int mc_type=1;mc_type<6;mc_type++)  SliceDistDMCTOF[mc_type]->Draw("same");	

	c2->cd(2);
	gPad->SetLogy();
        gPad->SetGridx();
        gPad->SetGridy();

	SliceDistDMCNaF[0]->Draw();
	 for (int mc_type=1;mc_type<6;mc_type++)  SliceDistDMCNaF[mc_type]->Draw("same");	

	c2->cd(3);
	gPad->SetLogy();
        gPad->SetGridx();
        gPad->SetGridy();

	SliceDistDMCAgl[0]->Draw();
	for (int mc_type=1;mc_type<6;mc_type++)  SliceDistDMCAgl[mc_type]->Draw("same");	

	string MCLegend[7]= {"Protons MC B800","Deuteorons MC \"GG_Blic\"","Deuterons MC \"GG_BlicDPMJet\"","Deuterons MC \"GG_QMD\"","Deuterons MC \"Shen_Blic\"","Deuterons MC \"Shen_BlicDPMJet\"","Deuterons MC \"Shen_QMD\""};

	int plottingstyles[6]={3,4,20,25,29,26};
	
	c3->Divide(3,1);
	c3->cd(1);
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraphErrors * EffFragmMCDTOF[6];
	EffFragmMCD_TOF -> Smooth();
	int p=0;
	for(int h=0;h<6;h++){
		EffFragmMCDTOF[h]=new TGraphErrors();
		p=0;
		for(int i=0;i<nbinsToF;i++) {EffFragmMCDTOF[h]->SetPoint(p,ToFDB.EkPerMassBinCent(i),EffFragmMCD_TOF ->GetBinContent(i+1,h+1));p++;}
		EffFragmMCDTOF[h]->SetMarkerStyle(8);
		EffFragmMCDTOF[h]->SetMarkerColor(4);
		EffFragmMCDTOF[h]->SetMarkerSize(2);
		EffFragmMCDTOF[h]->SetLineColor(4);
		EffFragmMCDTOF[h]->SetLineWidth(1);
		EffFragmMCDTOF[h]->SetMarkerStyle(plottingstyles[h]);
		EffFragmMCDTOF[h]->SetTitle("");
        EffFragmMCDTOF[h]->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
        EffFragmMCDTOF[h]->GetYaxis()->SetTitle("Efficiency");
        EffFragmMCDTOF[h]->GetXaxis()->SetTitleSize(0.045);
        EffFragmMCDTOF[h]->GetYaxis()->SetTitleSize(0.045);
        }
        EffFragmMCDTOF[0]->Draw("AP");
        for(int h=1;h<6;h++){
                EffFragmMCDTOF[h]->Draw("Psame");
        }
        {       TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(EffFragmMCDTOF[0],MCLegend[1].c_str(), "p");
        for(int h=1;h<6;h++){
                leg->AddEntry(EffFragmMCDTOF[h],MCLegend[h+1].c_str(), "p");
        }
                leg->SetLineWidth(2);
                leg->Draw("same");
        }

 
	c3->cd(2);
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraphErrors * EffFragmMCDNaF[6];
	EffFragmMCD_NaF -> Smooth();
	p=0;
	for(int h=0;h<6;h++){
		EffFragmMCDNaF[h]=new TGraphErrors();
		p=0;
		for(int i=0;i<nbinsNaF;i++) {EffFragmMCDNaF[h]->SetPoint(p,NaFDB.EkPerMassBinCent(i),EffFragmMCD_NaF ->GetBinContent(i+1,h+1));p++;}
		EffFragmMCDNaF[h]->SetMarkerStyle(8);
		EffFragmMCDNaF[h]->SetMarkerColor(4);
		EffFragmMCDNaF[h]->SetMarkerSize(2);
		EffFragmMCDNaF[h]->SetLineColor(4);
		EffFragmMCDNaF[h]->SetLineWidth(1);
		EffFragmMCDNaF[h]->SetMarkerStyle(plottingstyles[h]);
		EffFragmMCDNaF[h]->SetTitle("");
        EffFragmMCDNaF[h]->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
        EffFragmMCDNaF[h]->GetYaxis()->SetTitle("Efficiency");
        EffFragmMCDNaF[h]->GetXaxis()->SetTitleSize(0.045);
        EffFragmMCDNaF[h]->GetYaxis()->SetTitleSize(0.045);
        }
        EffFragmMCDNaF[0]->Draw("AP");
        for(int h=1;h<6;h++){
                EffFragmMCDNaF[h]->Draw("Psame");
        }
        {       TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(EffFragmMCDNaF[0],MCLegend[1].c_str(), "p");
        for(int h=1;h<6;h++){
                leg->AddEntry(EffFragmMCDNaF[h],MCLegend[h+1].c_str(), "p");
        }
                leg->SetLineWidth(2);
                leg->Draw("same");
        }

 
	 
	c3->cd(3);
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraphErrors * EffFragmMCDAgl[6];
	EffFragmMCD_Agl -> Smooth();
	p=0;
	for(int h=0;h<6;h++){
		EffFragmMCDAgl[h]=new TGraphErrors();
		p=0;
		for(int i=0;i<nbinsAgl;i++) {EffFragmMCDAgl[h]->SetPoint(p,NaFDB.EkPerMassBinCent(i),EffFragmMCD_Agl ->GetBinContent(i+1,h+1));p++;}
		EffFragmMCDAgl[h]->SetMarkerStyle(8);
		EffFragmMCDAgl[h]->SetMarkerColor(4);
		EffFragmMCDAgl[h]->SetMarkerSize(2);
		EffFragmMCDAgl[h]->SetLineColor(4);
		EffFragmMCDAgl[h]->SetLineWidth(1);
		EffFragmMCDAgl[h]->SetMarkerStyle(plottingstyles[h]);
		EffFragmMCDAgl[h]->SetTitle("");
        EffFragmMCDAgl[h]->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
        EffFragmMCDAgl[h]->GetYaxis()->SetTitle("Efficiency");
        EffFragmMCDAgl[h]->GetXaxis()->SetTitleSize(0.045);
        EffFragmMCDAgl[h]->GetYaxis()->SetTitleSize(0.045);
        }
        EffFragmMCDAgl[0]->Draw("AP");
        for(int h=1;h<6;h++){
                EffFragmMCDAgl[h]->Draw("Psame");
        }
        {       TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(EffFragmMCDAgl[0],MCLegend[1].c_str(), "p");
        for(int h=1;h<6;h++){
                leg->AddEntry(EffFragmMCDAgl[h],MCLegend[h+1].c_str(), "p");
        }
                leg->SetLineWidth(2);
                leg->Draw("same");
        }










	finalPlots.Add(c);
        finalPlots.Add(c1);
        finalPlots.Add(c2);
 	finalPlots.Add(c3);


        finalPlots.writeObjsInFolder("MC Results/Quality/D MC distrib");
	
	return;
}
