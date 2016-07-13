
void DVSMCRICHeff_Plot (TH1 *  RICH_Correction_P_NaF,
                        TH1 *     RICH_Correction_P_Agl,
                         TH1 *    RICH_Correction_D_NaF,
                         TH1 *    RICH_Correction_D_Agl
        ){	

	TCanvas *c20_bis=new TCanvas("Data vs MC: RICH");


	c20_bis->Divide(2,1);

	c20_bis->cd(1);
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraphErrors * RICHDVSMC_P_GraphNaF=new TGraphErrors();
	RICHDVSMC_P_GraphNaF->SetName("RICHDVSMC_P_GraphNaF");
	int j=0;
	for(int i=1;i<nbinsNaF;i++) {
		if(RICH_Correction_P_NaF -> GetBinContent(i+1)>0){
			RICHDVSMC_P_GraphNaF->SetPoint(j,NaFPB.EkPerMassBinCent(i),RICH_Correction_P_NaF -> GetBinContent(i+1));
			RICHDVSMC_P_GraphNaF->SetPointError(j,0,RICH_Correction_P_NaF -> GetBinError(i+1));
			j++;
		}
	}
	RICHDVSMC_P_GraphNaF->SetLineColor(2);
	RICHDVSMC_P_GraphNaF->SetFillColor(2);
	RICHDVSMC_P_GraphNaF->SetFillStyle(3001);
	RICHDVSMC_P_GraphNaF->SetLineWidth(4);
	RICHDVSMC_P_GraphNaF->Draw("AP4C");

	c20_bis->cd(2);
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraphErrors * RICHDVSMC_P_GraphAgl=new TGraphErrors();
        RICHDVSMC_P_GraphAgl->SetName("RICHDVSMC_P_GraphAgl");
        j=0;
        for(int i=1;i<nbinsAgl;i++) {
                if(RICH_Correction_P_Agl -> GetBinContent(i+1)>0){
                        RICHDVSMC_P_GraphAgl->SetPoint(j,AglPB.EkPerMassBinCent(i),RICH_Correction_P_Agl -> GetBinContent(i+1));
                        RICHDVSMC_P_GraphAgl->SetPointError(j,0,RICH_Correction_P_Agl -> GetBinError(i+1));
                        j++;
                }
        }
        RICHDVSMC_P_GraphAgl->SetLineColor(2);
        RICHDVSMC_P_GraphAgl->SetFillColor(2);
        RICHDVSMC_P_GraphAgl->SetFillStyle(3001);
        RICHDVSMC_P_GraphAgl->SetLineWidth(4);
        RICHDVSMC_P_GraphAgl->Draw("AP4C");


	finalPlots.Add(c20_bis);
	finalPlots.writeObjsInFolder("DATA-driven Results/Data vs MC/RICH");
	

	return;
}
