void 	DVSMCPreSeleff_Plot(TH1 *PreSel_Correction_R  ,
		TH1 *PreSel_Correction_TOF,
		TH1 *PreSel_Correction_NaF,
		TH1 *PreSel_Correction_Agl,
		TH1 *	EffData_R 	,		
		TH1 *       EffMC_R   

		){



	string tagli[3]={"Matching TOF","Chi^2 R","1 Tr. Track"};
	string nome;

	TCanvas *c21[3];
	TCanvas *c20[3];

	TGraphErrors * PreSel_Correction_R_Graph[3];

	TGraphErrors * PreSel_Correction_TOF_Graph[3];
	TGraphErrors * PreSel_Correction_NaF_Graph[3];
	TGraphErrors * PreSel_Correction_Agl_Graph[3];

	TGraphErrors * DATAEff_Graph[3];
	TGraphErrors * MCEff_Graph[3];
	for(int S=0;S<3;S++){
		c20[S] = new TCanvas(("Data vs MC: "+tagli[S] +"(R Bins)").c_str());

		c20[S]->Divide(2,1);
		c20[S]->cd(1);	
		gPad->SetLogx();
		gPad->SetGridx();
		gPad->SetGridy();
		PreSel_Correction_R_Graph[S] = new TGraphErrors(nbinsr);
		PreSel_Correction_R_Graph[S]->SetName(("Data vs MC: "+tagli[S] +"_R").c_str());
		int j=0;
		for(int i=0;i<nbinsr;i++) {
			if(PreSel_Correction_R -> GetBinContent(i+1,S+1)>0){
				PreSel_Correction_R_Graph[S]->SetPoint(j,PRB.RigBinCent(i),PreSel_Correction_R -> GetBinContent(i+1,S+1));
				PreSel_Correction_R_Graph[S]->SetPointError(j,0,PreSel_Correction_R -> GetBinError(i+1,S+1));
				j++;
			}
		}
		PreSel_Correction_R_Graph[S]->SetLineColor(2);
		PreSel_Correction_R_Graph[S]->SetFillColor(2);
		PreSel_Correction_R_Graph[S]->SetFillStyle(3001);
		PreSel_Correction_R_Graph[S]->SetLineWidth(4);
		PreSel_Correction_R_Graph[S]->Draw("AP4C");


		c20[S]->cd(2);
		gPad->SetLogx();
		gPad->SetGridx();
		gPad->SetGridy();
		MCEff_Graph[S] = new TGraphErrors(nbinsr);
		MCEff_Graph[S]->SetName(("Data vs MC: "+tagli[S] +"_R").c_str());
		j=0;
		for(int i=0;i<nbinsr;i++) {
			if(EffMC_R  -> GetBinContent(i+1,S+1)>0){
				MCEff_Graph[S]->SetPoint(j,PRB.RigBinCent(i),EffMC_R   -> GetBinContent(i+1,S+1));
				MCEff_Graph[S]->SetPointError(j,0,EffMC_R   -> GetBinError(i+1,S+1));
				j++;
			}
		}
		MCEff_Graph[S]->SetLineColor(2);
		MCEff_Graph[S]->SetFillColor(2);
		MCEff_Graph[S]->SetFillStyle(3001);
		MCEff_Graph[S]->SetLineWidth(4);
		MCEff_Graph[S]->Draw("AP4C");

		DATAEff_Graph[S] = new TGraphErrors(); DATAEff_Graph[S]->SetName(("Data vs MC: "+tagli[S] +"_R").c_str());
		j=0;
		for(int i=1;i<nbinsr;i++) {
			if(EffData_R  -> GetBinContent(i+1,S+1)>0){
				DATAEff_Graph[S]->SetPoint(j,PRB.RigBinCent(i),EffData_R   -> GetBinContent(i+1,S+1));
				DATAEff_Graph[S]->SetPointError(j,0,EffData_R   -> GetBinError(i+1,S+1));
				j++;
			}
		}
		DATAEff_Graph[S]->SetLineColor(1);
		DATAEff_Graph[S]->SetFillColor(1);
		DATAEff_Graph[S]->SetFillStyle(3001);
		DATAEff_Graph[S]->SetLineWidth(4);
		DATAEff_Graph[S]->Draw("P4Csame");


		c21[S] = new TCanvas(("Data vs MC: "+tagli[S] +"(Beta Bins)").c_str());
		c21[S] -> Divide(3,1);

		c21[S]->cd(1);
		gPad->SetLogx();
		gPad->SetGridx();
		gPad->SetGridy();

		PreSel_Correction_TOF_Graph[S] = new TGraphErrors(); PreSel_Correction_TOF_Graph[S]->SetName(("Data vs MC: "+tagli[S] +"_TOF").c_str());
		PreSel_Correction_NaF_Graph[S] = new TGraphErrors(); PreSel_Correction_NaF_Graph[S]->SetName(("Data vs MC: "+tagli[S] +"_NaF").c_str());
		PreSel_Correction_Agl_Graph[S] = new TGraphErrors(); PreSel_Correction_Agl_Graph[S]->SetName(("Data vs MC: "+tagli[S] +"_Agl").c_str());

		j=0;
		for(int i=1;i<nbinsToF;i++) {
			if(PreSel_Correction_TOF -> GetBinContent(i+1,S+1)>0){
				PreSel_Correction_TOF_Graph[S]->SetPoint(j,ToFPB.EkPerMassBinCent(i),PreSel_Correction_TOF -> GetBinContent(i+1,S+1));
				PreSel_Correction_TOF_Graph[S]->SetPointError(j,0,PreSel_Correction_TOF -> GetBinError(i+1,S+1));
				j++;
			}
		}
		PreSel_Correction_TOF_Graph[S]->SetLineColor(2);
		PreSel_Correction_TOF_Graph[S]->SetFillColor(2);
		PreSel_Correction_TOF_Graph[S]->SetFillStyle(3001);
		PreSel_Correction_TOF_Graph[S]->SetLineWidth(4);
		PreSel_Correction_TOF_Graph[S]->Draw("AP4C");

		c21[S]->cd(2);
		gPad->SetLogx();
		gPad->SetGridx();
		gPad->SetGridy();
		PreSel_Correction_NaF_Graph[S] = new TGraphErrors();
		j=0;
		for(int i=1;i<nbinsNaF;i++) {
			if(PreSel_Correction_NaF -> GetBinContent(i+1,S+1)>0){
				PreSel_Correction_NaF_Graph[S]->SetPoint(j,ToFPB.EkPerMassBinCent(i),PreSel_Correction_NaF -> GetBinContent(i+1,S+1));
				PreSel_Correction_NaF_Graph[S]->SetPointError(j,0,PreSel_Correction_NaF -> GetBinError(i+1,S+1));
				j++;
			}
		}
		PreSel_Correction_NaF_Graph[S]->SetLineColor(2);
		PreSel_Correction_NaF_Graph[S]->SetFillColor(2);
		PreSel_Correction_NaF_Graph[S]->SetFillStyle(3001);
		PreSel_Correction_NaF_Graph[S]->SetLineWidth(4);
		PreSel_Correction_NaF_Graph[S]->Draw("AP4C");

		c21[S]->cd(3);
		gPad->SetLogx();
		gPad->SetGridx();
		gPad->SetGridy();
		PreSel_Correction_Agl_Graph[S] = new TGraphErrors();
		j=0;
		for(int i=1;i<nbinsToF;i++) {
			if(PreSel_Correction_Agl -> GetBinContent(i+1,S+1)>0){
				PreSel_Correction_Agl_Graph[S]->SetPoint(j,ToFPB.EkPerMassBinCent(i),PreSel_Correction_Agl -> GetBinContent(i+1,S+1));
				PreSel_Correction_Agl_Graph[S]->SetPointError(j,0,PreSel_Correction_Agl -> GetBinError(i+1,S+1));
				j++;
			}
		}
		PreSel_Correction_Agl_Graph[S]->SetLineColor(2);
		PreSel_Correction_Agl_Graph[S]->SetFillColor(2);
		PreSel_Correction_Agl_Graph[S]->SetFillStyle(3001);
		PreSel_Correction_Agl_Graph[S]->SetLineWidth(4);
		PreSel_Correction_Agl_Graph[S]->Draw("AP4C");


	}


	for(int S=0;S<3;S++){
		finalPlots.Add(c20[S]);
		finalPlots.Add(c21[S]);
	}
	finalPlots.writeObjsInFolder("DATA-driven Results/Data vs MC/Protons");

	for(int S=0;S<3;S++){

		PreSel_Correction_R_Graph[S]  ->SetName(("DvsMC: " + tagli[S] + "_R").c_str()); 
		PreSel_Correction_TOF_Graph[S]->SetName(("DvsMC: " + tagli[S] + "_TOF").c_str());
		PreSel_Correction_NaF_Graph[S]->SetName(("DvsMC: " + tagli[S] + "_NaF").c_str());
		PreSel_Correction_Agl_Graph[S]->SetName(("DvsMC: " + tagli[S] + "_Agl").c_str());
	}

	for(int S=0;S<3;S++){
		finalPlots.Add(PreSel_Correction_R_Graph[S]  );
		finalPlots.Add(PreSel_Correction_TOF_Graph[S]);
		finalPlots.Add(PreSel_Correction_NaF_Graph[S]);
		finalPlots.Add(PreSel_Correction_Agl_Graph[S]);
	}
	finalPlots.writeObjsInFolder("Export/DvsMC");	


	return;
}
