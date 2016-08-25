TH2F * TransposeTH2F (TH2F * Histo){
	TH2F * transposed=new TH2F("","",Histo->GetNbinsY(),0,18,Histo->GetNbinsX(),-1.2,1.2);
	for(int x=0;x<Histo->GetNbinsY();x++)
		for(int y=0;y<Histo->GetNbinsX();y++)
			transposed->SetBinContent(x+1,y+1,Histo->GetBinContent(y+1,x+1))	;
	return transposed;

}


void AntiDCutOptimization_Plot(OptimizationCut *   DiscriminantCutTOF,
			       OptimizationCut *   DiscriminantCutNaF,
			       OptimizationCut *   DiscriminantCutAgl
		){
	TCanvas *  cTOF = new TCanvas("Discriminant MC distrib. (TOF)");
	cTOF -> Divide(6,3);

	TH1F *BinTOF_P[nbinsToF];
	TH1F *BinTOF_D[nbinsToF];
	for (int i=0;i<nbinsToF;i++){
		BinTOF_P[i]= ((TH1F*)DiscriminantCutTOF -> GetBinP(i));
		BinTOF_D[i]= ((TH1F*)DiscriminantCutTOF -> GetBinD(i));
		BinTOF_P[i] -> SetFillStyle(3001);
		BinTOF_D[i] -> SetFillStyle(3001);
		BinTOF_P[i] -> SetLineWidth(2);
		BinTOF_D[i] -> SetLineWidth(2);
		BinTOF_P[i] -> SetLineColor(2);
		BinTOF_D[i] -> SetLineColor(4);
		BinTOF_P[i] -> SetFillColor(2);
		BinTOF_D[i] -> SetFillColor(4);
		cTOF->cd(i+1);
		gPad->SetGridx();
		gPad->SetGridy();
		gPad->SetLogy();
		BinTOF_P[i] -> Draw();
		BinTOF_D[i] -> Draw("same");
		BinTOF_P[i] -> Draw("same");	
	}


	TCanvas *  cNaF = new TCanvas("Discriminant MC distrib. (NaF)");
	cNaF -> Divide(6,3);

	TH1F *BinNaF_P[nbinsNaF];
	TH1F *BinNaF_D[nbinsNaF];
	for (int i=0;i<nbinsToF;i++){
		BinNaF_P[i]= ((TH1F*)DiscriminantCutNaF -> GetBinP(i));
		BinNaF_D[i]= ((TH1F*)DiscriminantCutNaF -> GetBinD(i));
		BinNaF_P[i] -> SetFillStyle(3001);
		BinNaF_D[i] -> SetFillStyle(3001);
		BinNaF_P[i] -> SetLineWidth(2);
		BinNaF_D[i] -> SetLineWidth(2);
		BinNaF_P[i] -> SetLineColor(2);
		BinNaF_D[i] -> SetLineColor(4);
		BinNaF_P[i] -> SetFillColor(2);
		BinNaF_D[i] -> SetFillColor(4);
		cNaF->cd(i+1);
		gPad->SetGridx();
		gPad->SetGridy();
		gPad->SetLogy();
		BinNaF_P[i] -> Draw();
		BinNaF_D[i] -> Draw("same");
		BinNaF_P[i] -> Draw("same");	
	}


	TCanvas *  cAgl = new TCanvas("Discriminant MC distrib. (Agl)");
	cAgl -> Divide(6,3);

	TH1F *BinAgl_P[nbinsAgl];
	TH1F *BinAgl_D[nbinsAgl];
	for (int i=0;i<nbinsToF;i++){
		BinAgl_P[i]= ((TH1F*)DiscriminantCutAgl -> GetBinP(i));
		BinAgl_D[i]= ((TH1F*)DiscriminantCutAgl -> GetBinD(i));
		BinAgl_P[i] -> SetFillStyle(3001);
		BinAgl_D[i] -> SetFillStyle(3001);
		BinAgl_P[i] -> SetLineWidth(2);
		BinAgl_D[i] -> SetLineWidth(2);
		BinAgl_P[i] -> SetLineColor(2);
		BinAgl_D[i] -> SetLineColor(4);
		BinAgl_P[i] -> SetFillColor(2);
		BinAgl_D[i] -> SetFillColor(4);
		cAgl->cd(i+1);
		gPad->SetGridx();
		gPad->SetGridy();
		gPad->SetLogy();
		BinAgl_P[i] -> Draw();
		BinAgl_D[i] -> Draw("same");
		BinAgl_P[i] -> Draw("same");	
	}


	TCanvas *  OTOF = new TCanvas("Cut Optimization (TOF)");
        OTOF -> Divide(6,3);

	TH1F *OptTOF[nbinsAgl];

	for (int i=0;i<nbinsToF;i++){
		OptTOF[i] = ((TH1F*)DiscriminantCutTOF -> GetOptimizationCurve(i));
		OptTOF[i] -> SetLineWidth(2);
		OptTOF[i] -> SetLineColor(3);
		OTOF->cd(i+1);
                gPad->SetGridx();
                gPad->SetGridy();
                OptTOF[i] -> Draw();

	}

	TCanvas *  ONaF = new TCanvas("Cut Optimization (NaF)");
        ONaF -> Divide(6,3);

	TH1F *OptNaF[nbinsToF];

	for (int i=0;i<nbinsToF;i++){
		OptNaF[i] = ((TH1F*)DiscriminantCutNaF -> GetOptimizationCurve(i));
		OptNaF[i] -> SetLineWidth(2);
		OptNaF[i] -> SetLineColor(3);
		ONaF->cd(i+1);
                gPad->SetGridx();
                gPad->SetGridy();
                OptNaF[i] -> Draw();

	}


	TCanvas *  OAgl = new TCanvas("Cut Optimization (Agl)");
        OAgl -> Divide(6,3);

	TH1F *OptAgl[nbinsAgl];

	for (int i=0;i<nbinsToF;i++){
		OptAgl[i] = ((TH1F*)DiscriminantCutAgl -> GetOptimizationCurve(i));
		OptAgl[i] -> SetLineWidth(2);
		OptAgl[i] -> SetLineColor(3);
		OAgl->cd(i+1);
                gPad->SetGridx();
                gPad->SetGridy();
                OptAgl[i] -> Draw();

	}


	TCanvas * Cuts = new TCanvas("Optimized values for cut");
	Cuts->Divide(1,3);
	
	TH1F * cutTOF;
	TH1F * cutNaF;
	TH1F * cutAgl;	

	cutTOF = DiscriminantCutTOF -> Getcuts(); 
        cutNaF = DiscriminantCutNaF -> Getcuts(); 
        cutAgl = DiscriminantCutAgl -> Getcuts(); 

	cutTOF ->SetLineWidth(5);
	cutNaF ->SetLineWidth(5);
	cutAgl ->SetLineWidth(5);
	
	cutTOF ->SetLineColor(1);
	cutNaF ->SetLineColor(1);
	cutAgl ->SetLineColor(1);

	cutTOF ->GetYaxis() -> SetRangeUser(-1,1);
        cutNaF ->GetYaxis() -> SetRangeUser(-1,1);
        cutAgl ->GetYaxis() -> SetRangeUser(-1,1);

	TH2F * DistribDTOF=TransposeTH2F((TH2F*) DiscriminantCutTOF -> Distrib_D);
	TH2F * DistribPTOF=TransposeTH2F((TH2F*) DiscriminantCutTOF -> Distrib_P);

	TH2F * DistribDNaF=TransposeTH2F((TH2F*) DiscriminantCutNaF -> Distrib_D);
	TH2F * DistribPNaF=TransposeTH2F((TH2F*) DiscriminantCutNaF -> Distrib_P);

	TH2F * DistribDAgl=TransposeTH2F((TH2F*) DiscriminantCutAgl -> Distrib_D);
	TH2F * DistribPAgl=TransposeTH2F((TH2F*) DiscriminantCutAgl -> Distrib_P);



	Cuts->cd(1);
	DistribPTOF -> SetMarkerColor(2);
	DistribDTOF -> SetMarkerColor(4);
	DistribDTOF -> SetLineWidth(3);
	DistribDTOF -> SetContour(2);
	DistribPTOF -> Draw();
	DistribDTOF -> Draw("CONT3 same");
	cutTOF ->Draw("same");
	Cuts->cd(2);
	DistribPNaF -> SetMarkerColor(2);
	DistribDNaF -> SetMarkerColor(4);
	DistribDNaF -> SetLineWidth(3);
	DistribDNaF -> SetContour(2);
	DistribPNaF -> Draw();
	DistribDNaF -> Draw("CONT3 same");
	cutNaF ->Draw("same");
	Cuts->cd(3);
	DistribPAgl -> SetMarkerColor(2);
	DistribDAgl -> SetMarkerColor(4);
	DistribDAgl -> SetLineWidth(3);
	DistribDAgl -> SetContour(2);
	DistribPAgl -> Draw();
	DistribDAgl -> Draw("CONT3 same");
	cutAgl ->Draw("same");
	
	finalPlots.Add(cTOF);
	finalPlots.Add(cNaF);
	finalPlots.Add(cAgl);
	finalPlots.Add(OTOF);
	finalPlots.Add(ONaF);
	finalPlots.Add(OAgl);
	finalPlots.Add(Cuts);
	finalPlots.writeObjsInFolder("Anti_D Predictions/DistanceCut Optimization ");

}
