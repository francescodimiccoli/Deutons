
void DVSMCQualeff_Plot(      TH1 *   DistP_Correction_R      ,
                             TH1 *   DistP_Correction_TOF,
                             TH1 *   DistP_Correction_NaF,
                             TH1 *   DistP_Correction_Agl,

                              TH1 *  LikP_Correction_R   ,
                              TH1 *  LikP_Correction_TOF ,
                              TH1 *  LikP_Correction_NaF ,
                              TH1 *  LikP_Correction_Agl
        ){


	TCanvas *c20=new TCanvas("Data vs MC: Likelihood (R bins)");
	TCanvas *c21=new TCanvas("Data vs MC: Distance (R bins)");

	TCanvas *c20_bis=new TCanvas("Data vs MC: Likelihood (Beta bins)");
	TCanvas *c21_bis=new TCanvas("Data vs MC: Distance (Beta bins)");

	TGraphErrors *LikDVSMC_P_Graph;
	TGraphErrors *DistDVSMC_P_Graph;

	c20->cd();
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	LikDVSMC_P_Graph=new TGraphErrors();
	LikDVSMC_P_Graph->SetName("LikDVSMC_P_Graph");
	int j=0;
	for(int i=1;i<nbinsr;i++) {
		if(LikP_Correction_R -> GetBinContent(i+1)>0){
			LikDVSMC_P_Graph->SetPoint(j,PRB.RigBinCent(i),LikP_Correction_R -> GetBinContent(i+1));
			LikDVSMC_P_Graph->SetPointError(j,0,LikP_Correction_R -> GetBinError(i+1));
			j++;
		}
	}
	LikDVSMC_P_Graph->SetLineColor(2);
	LikDVSMC_P_Graph->SetFillColor(2);
	LikDVSMC_P_Graph->SetFillStyle(3001);
	LikDVSMC_P_Graph->SetLineWidth(4);
	LikDVSMC_P_Graph->Draw("AP4C");

	c21->cd();
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	DistDVSMC_P_Graph=new TGraphErrors();
	DistDVSMC_P_Graph->SetName("DistDVSMC_P_Graph");
	j=0;
	for(int i=1;i<nbinsr;i++) {
		if(DistP_Correction_R -> GetBinContent(i+1)>0){
			DistDVSMC_P_Graph->SetPoint(j,PRB.RigBinCent(i),DistP_Correction_R -> GetBinContent(i+1));
			DistDVSMC_P_Graph->SetPointError(j,0,DistP_Correction_R -> GetBinError(i+1));
			j++;
		}
	}
	DistDVSMC_P_Graph->SetLineColor(2);
	DistDVSMC_P_Graph->SetFillColor(2);
	DistDVSMC_P_Graph->SetFillStyle(3001);
	DistDVSMC_P_Graph->SetLineWidth(4);
	DistDVSMC_P_Graph->Draw("AP4C");

	c20_bis->Divide(3,1);

	c20_bis->cd(1);
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraphErrors * LikDVSMC_P_GraphTOF=new TGraphErrors();
	LikDVSMC_P_GraphTOF->SetName("LikDVSMC_P_GraphTOF");
	j=0;
	for(int i=1;i<nbinsToF;i++) {
		if(LikP_Correction_TOF -> GetBinContent(i+1)>0){
			LikDVSMC_P_GraphTOF->SetPoint(j,ToFPB.EkPerMassBinCent(i),LikP_Correction_TOF -> GetBinContent(i+1));
			LikDVSMC_P_GraphTOF->SetPointError(j,0,LikP_Correction_TOF -> GetBinError(i+1));
			j++;
		}
	}
	LikDVSMC_P_GraphTOF->SetLineColor(2);
	LikDVSMC_P_GraphTOF->SetFillColor(2);
	LikDVSMC_P_GraphTOF->SetFillStyle(3001);
	LikDVSMC_P_GraphTOF->SetLineWidth(4);
	LikDVSMC_P_GraphTOF->Draw("AP4C");

	c20_bis->cd(2);
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraphErrors * LikDVSMC_P_GraphNaF=new TGraphErrors();
	LikDVSMC_P_GraphNaF->SetName("LikDVSMC_P_GraphNaF");
	j=0;
	for(int i=1;i<nbinsNaF;i++) {
		if(LikP_Correction_NaF -> GetBinContent(i+1)>0){
			LikDVSMC_P_GraphNaF->SetPoint(j,NaFPB.EkPerMassBinCent(i),LikP_Correction_NaF -> GetBinContent(i+1));
			LikDVSMC_P_GraphNaF->SetPointError(j,0,LikP_Correction_NaF -> GetBinError(i+1));
			j++;
		}
	}
	LikDVSMC_P_GraphNaF->SetLineColor(2);
	LikDVSMC_P_GraphNaF->SetFillColor(2);
	LikDVSMC_P_GraphNaF->SetFillStyle(3001);
	LikDVSMC_P_GraphNaF->SetLineWidth(4);
	LikDVSMC_P_GraphNaF->Draw("AP4C");

	c20_bis->cd(3);
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraphErrors * LikDVSMC_P_GraphAgl=new TGraphErrors();
	LikDVSMC_P_GraphAgl->SetName("LikDVSMC_P_GraphAgl");
	j=0;
	for(int i=1;i<nbinsAgl;i++) {
		if(LikP_Correction_Agl -> GetBinContent(i+1)>0){
			LikDVSMC_P_GraphAgl->SetPoint(j,AglPB.EkPerMassBinCent(i),LikP_Correction_Agl -> GetBinContent(i+1));
			LikDVSMC_P_GraphAgl->SetPointError(j,0,LikP_Correction_Agl -> GetBinError(i+1));
			j++;
		}
	}
	LikDVSMC_P_GraphAgl->SetLineColor(2);
	LikDVSMC_P_GraphAgl->SetFillColor(2);
	LikDVSMC_P_GraphAgl->SetFillStyle(3001);
	LikDVSMC_P_GraphAgl->SetLineWidth(4);
	LikDVSMC_P_GraphAgl->Draw("AP4C");

	c21_bis->Divide(3,1);

	c21_bis->cd(1);
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraphErrors * DistDVSMC_P_GraphTOF=new TGraphErrors();
	DistDVSMC_P_GraphTOF->SetName("DistDVSMC_P_GraphTOF");
	j=0;
	for(int i=1;i<nbinsToF;i++) {
		if(DistP_Correction_TOF -> GetBinContent(i+1)>0){
			DistDVSMC_P_GraphTOF->SetPoint(j,ToFPB.EkPerMassBinCent(i),DistP_Correction_TOF -> GetBinContent(i+1));
			DistDVSMC_P_GraphTOF->SetPointError(j,0,DistP_Correction_TOF -> GetBinError(i+1));
			j++;
		}
	}
	DistDVSMC_P_GraphTOF->SetLineColor(2);
	DistDVSMC_P_GraphTOF->SetFillColor(2);
	DistDVSMC_P_GraphTOF->SetFillStyle(3001);
	DistDVSMC_P_GraphTOF->SetLineWidth(4);
	DistDVSMC_P_GraphTOF->Draw("AP4C");

	c21_bis->cd(2);
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraphErrors * DistDVSMC_P_GraphNaF=new TGraphErrors();
	DistDVSMC_P_GraphNaF->SetName("DistDVSMC_P_GraphNaF");
	j=0;
	for(int i=1;i<nbinsNaF;i++) {
		if(DistP_Correction_NaF -> GetBinContent(i+1)>0){
			DistDVSMC_P_GraphNaF->SetPoint(j,NaFPB.EkPerMassBinCent(i),DistP_Correction_NaF -> GetBinContent(i+1));
			DistDVSMC_P_GraphNaF->SetPointError(j,0,DistP_Correction_NaF -> GetBinError(i+1));
			j++;
		}
	}
	DistDVSMC_P_GraphNaF->SetLineColor(2);
	DistDVSMC_P_GraphNaF->SetFillColor(2);
	DistDVSMC_P_GraphNaF->SetFillStyle(3001);
	DistDVSMC_P_GraphNaF->SetLineWidth(4);
	DistDVSMC_P_GraphNaF->Draw("AP4C");

	c21_bis->cd(3);
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraphErrors * DistDVSMC_P_GraphAgl=new TGraphErrors();
	DistDVSMC_P_GraphAgl->SetName("DistDVSMC_P_GraphAgl");
	j=0;
	for(int i=1;i<nbinsAgl;i++) {
		if(DistP_Correction_Agl -> GetBinContent(i+1)>0){
			DistDVSMC_P_GraphAgl->SetPoint(j,AglPB.EkPerMassBinCent(i),DistP_Correction_Agl -> GetBinContent(i+1));
			DistDVSMC_P_GraphAgl->SetPointError(j,0,DistP_Correction_Agl -> GetBinError(i+1));
			j++;
		}
	}
	DistDVSMC_P_GraphAgl->SetLineColor(2);
	DistDVSMC_P_GraphAgl->SetFillColor(2);
	DistDVSMC_P_GraphAgl->SetFillStyle(3001);
	DistDVSMC_P_GraphAgl->SetLineWidth(4);
	DistDVSMC_P_GraphAgl->Draw("AP4C");


	finalPlots.Add(c20);
	finalPlots.Add(c21);
	finalPlots.Add(c20_bis);
	finalPlots.Add(c21_bis);
	finalPlots.writeObjsInFolder("DATA-driven Results/Data vs MC/Protons");

        
	finalPlots.Add(DistDVSMC_P_Graph   );
	finalPlots.Add(DistDVSMC_P_GraphTOF);
	finalPlots.Add(DistDVSMC_P_GraphNaF);
	finalPlots.Add(DistDVSMC_P_GraphAgl);
	finalPlots.Add(LikDVSMC_P_Graph    );
	finalPlots.Add(LikDVSMC_P_GraphTOF );
	finalPlots.Add(LikDVSMC_P_GraphNaF );
	finalPlots.Add(LikDVSMC_P_GraphAgl );
	finalPlots.writeObjsInFolder("Export/DvsMC");

	return;
}
