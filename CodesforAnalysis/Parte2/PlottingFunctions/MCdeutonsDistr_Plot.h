void MCdeutonsDistr_Plot(
				
			TH1F * SliceMassDMCTOF[],
			TH1F * SliceMassDMCNaF[],
			TH1F * SliceMassDMCAgl[],

			TH1F * SliceLikDMCTOF[],
			TH1F * SliceLikDMCNaF[],
			TH1F * SliceLikDMCAgl[],

			TH1F * SliceDistDMCTOF[],
			TH1F * SliceDistDMCNaF[],
			TH1F * SliceDistDMCAgl[]


				){

	TCanvas * c  = new TCanvas("Mass Distributions");	
	TCanvas * c1 = new TCanvas("Lik Distributions");	
	TCanvas * c2 = new TCanvas("Dist Distributions");	

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

	finalPlots.Add(c);
        finalPlots.Add(c1);
        finalPlots.Add(c2);

        finalPlots.writeObjsInFolder("MC Results/Quality/D MC distrib");
	
	return;
}
