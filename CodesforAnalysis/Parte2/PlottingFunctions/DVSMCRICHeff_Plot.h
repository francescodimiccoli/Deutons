
void DVSMCRICHeff_Plot (TH1 *  RICH_Correction_P_NaF,
                        TH1 *     RICH_Correction_P_Agl,
                         TH1 *    RICH_Correction_D_NaF,
                         TH1 *    RICH_Correction_D_Agl,


			TH1F* RICH_CorrectionFit_P_NaF, 	
                        TH1F* RICH_CorrectionFit_P_Agl, 
                                                       
                        TH2F* RICH_CorrectionFit_D_NaF, 
                        TH2F* RICH_CorrectionFit_D_Agl 

        ){	

	TCanvas *c20_bis=new TCanvas("Data vs MC: RICH");


	c20_bis->Divide(2,1);

	c20_bis->cd(1);
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraphErrors * RICHDVSMC_P_GraphNaF=new TGraphErrors();
	RICHDVSMC_P_GraphNaF->SetName("RICH Agl Efficiency correction");
	TGraphErrors * RICHDVSMC_P_TH1FNaF=new TGraphErrors();
	int j=0;
	for(int i=1;i<nbinsNaF;i++) {
			RICHDVSMC_P_GraphNaF->SetPoint(j,NaFPB.EkPerMassBinCent(i),RICH_CorrectionFit_P_NaF -> GetBinContent(i+1));
			RICHDVSMC_P_GraphNaF->SetPointError(j,0,RICH_CorrectionFit_P_NaF -> GetBinError(i+1));
			RICHDVSMC_P_TH1FNaF->SetPoint(j,NaFPB.EkPerMassBinCent(i),RICH_Correction_P_NaF -> GetBinContent(i+1));
			RICHDVSMC_P_TH1FNaF->SetPointError(j,0,RICH_Correction_P_NaF -> GetBinError(i+1));
			
			j++;
	}
	RICHDVSMC_P_GraphNaF->SetLineColor(2);
	RICHDVSMC_P_GraphNaF->SetFillColor(2);
	RICHDVSMC_P_GraphNaF->SetLineWidth(4);
	RICHDVSMC_P_GraphNaF->SetFillStyle(3001);
	RICHDVSMC_P_GraphNaF->SetLineWidth(4);
	RICHDVSMC_P_TH1FNaF->SetMarkerColor(2);
	RICHDVSMC_P_TH1FNaF->SetMarkerStyle(8);
	RICHDVSMC_P_TH1FNaF->SetLineColor(2);
	RICHDVSMC_P_GraphNaF->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
        RICHDVSMC_P_GraphNaF->GetYaxis()->SetTitle("Efficiency (Data/MC)");

	RICHDVSMC_P_GraphNaF->Draw("AP4C");
	RICHDVSMC_P_TH1FNaF->Draw("Psame");
	{
                TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(RICHDVSMC_P_TH1FNaF,"Efficiency correction", "ep");
                leg->AddEntry(RICHDVSMC_P_GraphNaF,"Param.", "l");
                leg->SetLineWidth(2);
                leg->Draw("same");
        }




	c20_bis->cd(2);
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraphErrors * RICHDVSMC_P_GraphAgl=new TGraphErrors();
	RICHDVSMC_P_GraphAgl->SetName("RICH Agl Efficiency correction");
	TGraphErrors * RICHDVSMC_P_TH1FAgl=new TGraphErrors();
	
	j=0;
	for(int i=1;i<nbinsAgl;i++) {
			RICHDVSMC_P_GraphAgl->SetPoint(j,AglPB.EkPerMassBinCent(i),RICH_CorrectionFit_P_Agl -> GetBinContent(i+1));
			RICHDVSMC_P_GraphAgl->SetPointError(j,0,RICH_CorrectionFit_P_Agl -> GetBinError(i+1));
			RICHDVSMC_P_TH1FAgl->SetPoint(j,AglPB.EkPerMassBinCent(i),RICH_Correction_P_Agl -> GetBinContent(i+1));
			RICHDVSMC_P_TH1FAgl->SetPointError(j,0,RICH_Correction_P_Agl -> GetBinError(i+1));
			
			j++;
	}
	RICHDVSMC_P_GraphAgl->SetLineColor(2);
	RICHDVSMC_P_GraphAgl->SetFillColor(2);
	RICHDVSMC_P_GraphAgl->SetLineWidth(4);
	RICHDVSMC_P_GraphAgl->SetFillStyle(3001);
	RICHDVSMC_P_GraphAgl->SetLineWidth(4);
	RICHDVSMC_P_TH1FAgl->SetMarkerColor(2);
	RICHDVSMC_P_TH1FAgl->SetMarkerStyle(8);
	RICHDVSMC_P_TH1FAgl->SetLineColor(2);
	RICHDVSMC_P_GraphAgl->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
        RICHDVSMC_P_GraphAgl->GetYaxis()->SetTitle("Efficiency (Data/MC)");

	RICHDVSMC_P_GraphAgl->Draw("AP4C");
	RICHDVSMC_P_TH1FAgl->Draw("Psame");

	{
                TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(RICHDVSMC_P_TH1FAgl,"Efficiency correction", "ep");
                leg->AddEntry(RICHDVSMC_P_GraphAgl,"Param.", "l");
                leg->SetLineWidth(2);
                leg->Draw("same");
        }



	finalPlots.Add(c20_bis);
	finalPlots.writeObjsInFolder("DATA-driven Results/Data vs MC/RICH");
	

	return;
}
