void 	DVSMCPreSeleff_Plot(
		TH1F *PreSel_Correction_R[]  ,
		TH1F *PreSel_Correction_TOF[],
		TH1F *PreSel_Correction_NaF[],
		TH1F *PreSel_Correction_Agl[],

		TH1F *PreSel_CorrectionFit_R[]  ,
		TH1F *PreSel_CorrectionFit_TOF[],
		TH1F *PreSel_CorrectionFit_NaF[],
		TH1F *PreSel_CorrectionFit_Agl[],


		TH1F* EffData_R[],  		
		TH1F* EffData_TOF[],  
		TH1F* EffData_NaF[],  
		TH1F* EffData_Agl[],  

		TH1F* EffMC_R[],		
		TH1F* EffMC_TOF[],  
		TH1F* EffMC_NaF[],  
		TH1F* EffMC_Agl[]  


){



	string tagli[3]={"Matching TOF","Chi^2 R","1 Tr. Track"};
	string nome;

	TCanvas *c21[3];
	TCanvas *c20[3];

	TGraphErrors * PreSel_Correction_R_Graph[3];
	TGraphErrors * PreSel_Correction_TOF_Graph[3];
	TGraphErrors * PreSel_Correction_NaF_Graph[3];
	TGraphErrors * PreSel_Correction_Agl_Graph[3];

	TGraphErrors * PreSel_CorrectionFit_R_Graph[3];
	TGraphErrors * PreSel_CorrectionFit_TOF_Graph[3];
	TGraphErrors * PreSel_CorrectionFit_NaF_Graph[3];
	TGraphErrors * PreSel_CorrectionFit_Agl_Graph[3];




	TGraphErrors * PreSel_EffMC_R_Graph[3];
	TGraphErrors * PreSel_EffMC_TOF_Graph[3];
	TGraphErrors * PreSel_EffMC_NaF_Graph[3];
	TGraphErrors * PreSel_EffMC_Agl_Graph[3];

	TGraphErrors * PreSel_EffData_R_Graph[3];
	TGraphErrors * PreSel_EffData_TOF_Graph[3];
	TGraphErrors * PreSel_EffData_NaF_Graph[3];
	TGraphErrors * PreSel_EffData_Agl_Graph[3];




	for(int S=0;S<3;S++){
		c20[S] = new TCanvas(( tagli[S] + " (Data&MC Eff.)").c_str());
		c20[S]->Divide(2,2);
		c20[S]->cd(1);
		gPad->SetLogx();
		gPad->SetGridx();
		gPad->SetGridy();
		PreSel_EffMC_R_Graph[S]=new TGraphErrors();
		PreSel_EffData_R_Graph[S]=new TGraphErrors();

		int j=0;
		for(int i=0;i<nbinsr;i++) {
			if(EffMC_R[S] -> GetBinContent(i+1)>0){
				PreSel_EffMC_R_Graph[S]->SetPoint(j,PRB.RigBinCent(i),EffMC_R[S] -> GetBinContent(i+1));
				PreSel_EffMC_R_Graph[S]->SetPointError(j,0,EffMC_R[S] -> GetBinError(i+1));
				PreSel_EffData_R_Graph[S]->SetPoint(j,PRB.RigBinCent(i),EffData_R[S] -> GetBinContent(i+1));
				PreSel_EffData_R_Graph[S]->SetPointError(j,0,EffData_R[S] -> GetBinError(i+1));

				j++;
			}
		}
		PreSel_EffMC_R_Graph[S]->SetMarkerColor(2);
		PreSel_EffData_R_Graph[S]->SetMarkerColor(1);
		PreSel_EffMC_R_Graph[S]->SetMarkerStyle(8);
		PreSel_EffData_R_Graph[S]->SetMarkerStyle(21);
		PreSel_EffMC_R_Graph[S]->SetLineColor(2);
		PreSel_EffData_R_Graph[S]->SetLineColor(1);
		PreSel_EffMC_R_Graph[S]->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
		PreSel_EffMC_R_Graph[S]->GetYaxis()->SetTitle("Efficiency ");
		PreSel_EffMC_R_Graph[S]->SetTitle("R binning");
		PreSel_EffMC_R_Graph[S]->GetYaxis()->SetRangeUser(0.9,1.1);

		PreSel_EffMC_R_Graph[S]->GetXaxis()->SetTitleSize(0.045);
                PreSel_EffMC_R_Graph[S]->GetYaxis()->SetTitleSize(0.045);

		PreSel_EffMC_R_Graph[S]->Draw("AP");
		PreSel_EffData_R_Graph[S]->Draw("Psame");
		{
			TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
			leg->AddEntry(PreSel_EffMC_R_Graph[S],"Protons MC", "ep");
			leg->AddEntry(PreSel_EffData_R_Graph[S],"ISS Data", "ep");
			leg->SetLineWidth(2);
			leg->Draw("same");
		}


		c20[S]->cd(2);
		gPad->SetLogx();
		gPad->SetGridx();
		gPad->SetGridy();
		PreSel_EffMC_TOF_Graph[S]=new TGraphErrors();
		PreSel_EffData_TOF_Graph[S]=new TGraphErrors();

		j=0;
		for(int i=0;i<nbinsToF;i++) {
			if(EffMC_TOF[S] -> GetBinContent(i+1)>0){
				PreSel_EffMC_TOF_Graph[S]->SetPoint(j,ToFPB.EkPerMassBinCent(i),EffMC_TOF[S] -> GetBinContent(i+1));
				PreSel_EffMC_TOF_Graph[S]->SetPointError(j,0,EffMC_TOF[S] -> GetBinError(i+1));
				PreSel_EffData_TOF_Graph[S]->SetPoint(j,ToFPB.EkPerMassBinCent(i),EffData_TOF[S] -> GetBinContent(i+1));
				PreSel_EffData_TOF_Graph[S]->SetPointError(j,0,EffData_TOF[S] -> GetBinError(i+1));

				j++;
			}
		}
		PreSel_EffMC_TOF_Graph[S]->SetMarkerColor(2);
		PreSel_EffData_TOF_Graph[S]->SetMarkerColor(1);
		PreSel_EffMC_TOF_Graph[S]->SetMarkerStyle(8);
		PreSel_EffData_TOF_Graph[S]->SetMarkerStyle(21);
		PreSel_EffMC_TOF_Graph[S]->SetLineColor(2);
		PreSel_EffData_TOF_Graph[S]->SetLineColor(1);
		PreSel_EffMC_TOF_Graph[S]->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
		PreSel_EffMC_TOF_Graph[S]->GetYaxis()->SetTitle("Efficiency ");
		PreSel_EffMC_TOF_Graph[S]->SetTitle("R binning");
		PreSel_EffMC_TOF_Graph[S]->GetYaxis()->SetRangeUser(0.9,1.1);
		
		PreSel_EffMC_TOF_Graph[S]->GetXaxis()->SetTitleSize(0.045);
                PreSel_EffMC_TOF_Graph[S]->GetYaxis()->SetTitleSize(0.045);


		PreSel_EffMC_TOF_Graph[S]->Draw("AP");
		PreSel_EffData_TOF_Graph[S]->Draw("Psame");
		{
			TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
			leg->AddEntry(PreSel_EffMC_TOF_Graph[S],"Protons MC", "ep");
			leg->AddEntry(PreSel_EffData_TOF_Graph[S],"ISS Data", "ep");
			leg->SetLineWidth(2);
			leg->Draw("same");
		}

		c20[S]->cd(3);
		gPad->SetLogx();
		gPad->SetGridx();
		gPad->SetGridy();
		PreSel_EffMC_NaF_Graph[S]=new TGraphErrors();
		PreSel_EffData_NaF_Graph[S]=new TGraphErrors();

		j=0;
		for(int i=0;i<nbinsNaF;i++) {
			if(EffMC_NaF[S] -> GetBinContent(i+1)>0){
				PreSel_EffMC_NaF_Graph[S]->SetPoint(j,NaFPB.EkPerMassBinCent(i),EffMC_NaF[S] -> GetBinContent(i+1));
				PreSel_EffMC_NaF_Graph[S]->SetPointError(j,0,EffMC_NaF[S] -> GetBinError(i+1));
				PreSel_EffData_NaF_Graph[S]->SetPoint(j,NaFPB.EkPerMassBinCent(i),EffData_NaF[S] -> GetBinContent(i+1));
				PreSel_EffData_NaF_Graph[S]->SetPointError(j,0,EffData_NaF[S] -> GetBinError(i+1));

				j++;
			}
		}
		PreSel_EffMC_NaF_Graph[S]->SetMarkerColor(2);
		PreSel_EffData_NaF_Graph[S]->SetMarkerColor(1);
		PreSel_EffMC_NaF_Graph[S]->SetMarkerStyle(8);
		PreSel_EffData_NaF_Graph[S]->SetMarkerStyle(21);
		PreSel_EffMC_NaF_Graph[S]->SetLineColor(2);
		PreSel_EffData_NaF_Graph[S]->SetLineColor(1);
		PreSel_EffMC_NaF_Graph[S]->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
		PreSel_EffMC_NaF_Graph[S]->GetYaxis()->SetTitle("Efficiency ");
		PreSel_EffMC_NaF_Graph[S]->SetTitle("R binning");
		PreSel_EffMC_NaF_Graph[S]->GetYaxis()->SetRangeUser(0.9,1.1);
		
		PreSel_EffMC_NaF_Graph[S]->GetXaxis()->SetTitleSize(0.045);
                PreSel_EffMC_NaF_Graph[S]->GetYaxis()->SetTitleSize(0.045);

		PreSel_EffMC_NaF_Graph[S]->Draw("AP");
		PreSel_EffData_NaF_Graph[S]->Draw("Psame");
		{
			TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
			leg->AddEntry(PreSel_EffMC_NaF_Graph[S],"Protons MC", "ep");
			leg->AddEntry(PreSel_EffData_NaF_Graph[S],"ISS Data", "ep");
			leg->SetLineWidth(2);
			leg->Draw("same");
		}

		
		c20[S]->cd(4);
		gPad->SetLogx();
		gPad->SetGridx();
		gPad->SetGridy();
		PreSel_EffMC_Agl_Graph[S]=new TGraphErrors();
		PreSel_EffData_Agl_Graph[S]=new TGraphErrors();

		j=0;
		for(int i=0;i<nbinsNaF;i++) {
			if(EffMC_Agl[S] -> GetBinContent(i+1)>0){
				PreSel_EffMC_Agl_Graph[S]->SetPoint(j,AglPB.EkPerMassBinCent(i),EffMC_Agl[S] -> GetBinContent(i+1));
				PreSel_EffMC_Agl_Graph[S]->SetPointError(j,0,EffMC_Agl[S] -> GetBinError(i+1));
				PreSel_EffData_Agl_Graph[S]->SetPoint(j,AglPB.EkPerMassBinCent(i),EffData_Agl[S] -> GetBinContent(i+1));
				PreSel_EffData_Agl_Graph[S]->SetPointError(j,0,EffData_Agl[S] -> GetBinError(i+1));

				j++;
			}
		}
		PreSel_EffMC_Agl_Graph[S]->SetMarkerColor(2);
		PreSel_EffData_Agl_Graph[S]->SetMarkerColor(1);
		PreSel_EffMC_Agl_Graph[S]->SetMarkerStyle(8);
		PreSel_EffData_Agl_Graph[S]->SetMarkerStyle(21);
		PreSel_EffMC_Agl_Graph[S]->SetLineColor(2);
		PreSel_EffData_Agl_Graph[S]->SetLineColor(1);
		PreSel_EffMC_Agl_Graph[S]->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
		PreSel_EffMC_Agl_Graph[S]->GetYaxis()->SetTitle("Efficiency ");
		PreSel_EffMC_Agl_Graph[S]->SetTitle("R binning");
		PreSel_EffMC_Agl_Graph[S]->GetYaxis()->SetRangeUser(0.9,1.1);
		
		PreSel_EffMC_Agl_Graph[S]->GetXaxis()->SetTitleSize(0.045);
                PreSel_EffMC_Agl_Graph[S]->GetYaxis()->SetTitleSize(0.045);

		PreSel_EffMC_Agl_Graph[S]->Draw("AP");
		PreSel_EffData_Agl_Graph[S]->Draw("Psame");
		{
			TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
			leg->AddEntry(PreSel_EffMC_Agl_Graph[S],"Protons MC", "ep");
			leg->AddEntry(PreSel_EffData_Agl_Graph[S],"ISS Data", "ep");
			leg->SetLineWidth(2);
			leg->Draw("same");
		}





	}


	for(int S=0;S<3;S++) finalPlots.Add(c20[S]);
	finalPlots.writeObjsInFolder("DATA-driven Results/Data vs MC/Selections/Eff.");

	for(int S=0;S<3;S++){

		c21[S] = new TCanvas(( tagli[S] + " (Corrections)").c_str());
		c21[S]->Divide(2,2);

		c21[S]->cd(1);
		gPad->SetLogx();
		gPad->SetGridx();
		gPad->SetGridy();
		PreSel_Correction_R_Graph[S]=new TGraphErrors();
		PreSel_CorrectionFit_R_Graph[S]=new TGraphErrors();

		int j=0;
		for(int i=0;i<nbinsr;i++) {
			if(PreSel_CorrectionFit_R[S] -> GetBinContent(i+1)>0){
				PreSel_Correction_R_Graph[S]->SetPoint(j,PRB.RigBinCent(i),PreSel_Correction_R[S] -> GetBinContent(i+1));
				PreSel_Correction_R_Graph[S]->SetPointError(j,0,PreSel_Correction_R[S] -> GetBinError(i+1));
				PreSel_CorrectionFit_R_Graph[S]->SetPoint(j,PRB.RigBinCent(i),PreSel_CorrectionFit_R[S] -> GetBinContent(i+1));
				PreSel_CorrectionFit_R_Graph[S]->SetPointError(j,0,PreSel_CorrectionFit_R[S] -> GetBinError(i+1));

				j++;
			}
		}
		PreSel_Correction_R_Graph[S]->SetLineColor(2);
		PreSel_CorrectionFit_R_Graph[S]->SetFillColor(2);
		PreSel_CorrectionFit_R_Graph[S]->SetFillStyle(3001);
		PreSel_CorrectionFit_R_Graph[S]->SetLineWidth(4);
		PreSel_Correction_R_Graph[S]->SetMarkerColor(2);
		PreSel_Correction_R_Graph[S]->SetMarkerStyle(8);
		PreSel_Correction_R_Graph[S]->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
		PreSel_Correction_R_Graph[S]->GetYaxis()->SetTitle("Efficiency (Data/MC)");
		PreSel_Correction_R_Graph[S]->SetTitle("R binning");
		PreSel_Correction_R_Graph[S]->GetYaxis()->SetRangeUser(0.9,1.2);

		PreSel_Correction_R_Graph[S]->GetXaxis()->SetTitleSize(0.045);
                PreSel_Correction_R_Graph[S]->GetYaxis()->SetTitleSize(0.045);

		PreSel_Correction_R_Graph[S]->Draw("AP");
		PreSel_CorrectionFit_R_Graph[S]->Draw("P4Csame");
		{
			TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
			leg->AddEntry(PreSel_Correction_R_Graph[S],"Efficiency correction", "ep");
			leg->AddEntry(PreSel_CorrectionFit_R_Graph[S],"Param.", "l");
			leg->SetLineWidth(2);
			leg->Draw("same");
		}



		c21[S]->cd(2);
		gPad->SetLogx();
		gPad->SetGridx();
		gPad->SetGridy();
		PreSel_Correction_TOF_Graph[S]=new TGraphErrors();
		PreSel_CorrectionFit_TOF_Graph[S]=new TGraphErrors();

		j=0;
		for(int i=0;i<nbinsToF;i++) {
			if(PreSel_CorrectionFit_TOF[S] -> GetBinContent(i+1)>0){
				PreSel_Correction_TOF_Graph[S]->SetPoint(j,ToFPB.EkPerMassBinCent(i),PreSel_Correction_TOF[S] -> GetBinContent(i+1));
				PreSel_Correction_TOF_Graph[S]->SetPointError(j,0,PreSel_Correction_TOF[S] -> GetBinError(i+1));
				PreSel_CorrectionFit_TOF_Graph[S]->SetPoint(j,ToFPB.EkPerMassBinCent(i),PreSel_CorrectionFit_TOF[S] -> GetBinContent(i+1));
				PreSel_CorrectionFit_TOF_Graph[S]->SetPointError(j,0,PreSel_CorrectionFit_TOF[S] -> GetBinError(i+1));

				j++;
			}
		}
		PreSel_Correction_TOF_Graph[S]->SetLineColor(2);
		PreSel_CorrectionFit_TOF_Graph[S]->SetFillColor(2);
		PreSel_CorrectionFit_TOF_Graph[S]->SetFillStyle(3001);
		PreSel_CorrectionFit_TOF_Graph[S]->SetLineWidth(4);
		PreSel_Correction_TOF_Graph[S]->SetMarkerColor(2);
		PreSel_Correction_TOF_Graph[S]->SetMarkerStyle(8);
		PreSel_Correction_TOF_Graph[S]->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
		PreSel_Correction_TOF_Graph[S]->GetYaxis()->SetTitle("Efficiency (Data/MC)");
		PreSel_Correction_TOF_Graph[S]->SetTitle("R binning");
		PreSel_Correction_TOF_Graph[S]->GetYaxis()->SetRangeUser(0.9,1.2);

		PreSel_Correction_TOF_Graph[S]->GetXaxis()->SetTitleSize(0.045);
                PreSel_Correction_TOF_Graph[S]->GetYaxis()->SetTitleSize(0.045);

		PreSel_Correction_TOF_Graph[S]->Draw("AP");
		PreSel_CorrectionFit_TOF_Graph[S]->Draw("P4Csame");
		{
			TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
			leg->AddEntry(PreSel_Correction_TOF_Graph[S],"Efficiency correction", "ep");
			leg->AddEntry(PreSel_CorrectionFit_TOF_Graph[S],"Param.", "l");
			leg->SetLineWidth(2);
			leg->Draw("same");
		}


		c21[S]->cd(3);
		gPad->SetLogx();
		gPad->SetGridx();
		gPad->SetGridy();
		PreSel_Correction_NaF_Graph[S]=new TGraphErrors();
		PreSel_CorrectionFit_NaF_Graph[S]=new TGraphErrors();

		j=0;
		for(int i=0;i<nbinsNaF;i++) {
			if(PreSel_CorrectionFit_NaF[S] -> GetBinContent(i+1)>0){
				PreSel_Correction_NaF_Graph[S]->SetPoint(j,NaFPB.EkPerMassBinCent(i),PreSel_Correction_NaF[S] -> GetBinContent(i+1));
				PreSel_Correction_NaF_Graph[S]->SetPointError(j,0,PreSel_Correction_NaF[S] -> GetBinError(i+1));
				PreSel_CorrectionFit_NaF_Graph[S]->SetPoint(j,NaFPB.EkPerMassBinCent(i),PreSel_CorrectionFit_NaF[S] -> GetBinContent(i+1));
				PreSel_CorrectionFit_NaF_Graph[S]->SetPointError(j,0,PreSel_CorrectionFit_NaF[S] -> GetBinError(i+1));

				j++;
			}
		}
		PreSel_Correction_NaF_Graph[S]->SetLineColor(2);
		PreSel_CorrectionFit_NaF_Graph[S]->SetFillColor(2);
		PreSel_CorrectionFit_NaF_Graph[S]->SetFillStyle(3001);
		PreSel_CorrectionFit_NaF_Graph[S]->SetLineWidth(4);
		PreSel_Correction_NaF_Graph[S]->SetMarkerColor(2);
		PreSel_Correction_NaF_Graph[S]->SetMarkerStyle(8);
		PreSel_Correction_NaF_Graph[S]->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
		PreSel_Correction_NaF_Graph[S]->GetYaxis()->SetTitle("Efficiency (Data/MC)");
		PreSel_Correction_NaF_Graph[S]->SetTitle("R binning");
		PreSel_Correction_NaF_Graph[S]->GetYaxis()->SetRangeUser(0.9,1.2);

		PreSel_Correction_NaF_Graph[S]->GetXaxis()->SetTitleSize(0.045);
                PreSel_Correction_NaF_Graph[S]->GetYaxis()->SetTitleSize(0.045);

		PreSel_Correction_NaF_Graph[S]->Draw("AP");
		PreSel_CorrectionFit_NaF_Graph[S]->Draw("P4Csame");
		{
			TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
			leg->AddEntry(PreSel_Correction_NaF_Graph[S],"Efficiency correction", "ep");
			leg->AddEntry(PreSel_CorrectionFit_NaF_Graph[S],"Param.", "l");
			leg->SetLineWidth(2);
			leg->Draw("same");
		}



		c21[S]->cd(4);
		gPad->SetLogx();
		gPad->SetGridx();
		gPad->SetGridy();
		PreSel_Correction_Agl_Graph[S]=new TGraphErrors();
		PreSel_CorrectionFit_Agl_Graph[S]=new TGraphErrors();

		j=0;
		for(int i=0;i<nbinsAgl;i++) {
			if(PreSel_CorrectionFit_Agl[S] -> GetBinContent(i+1)>0){
				PreSel_Correction_Agl_Graph[S]->SetPoint(j,AglPB.EkPerMassBinCent(i),PreSel_Correction_Agl[S] -> GetBinContent(i+1));
				PreSel_Correction_Agl_Graph[S]->SetPointError(j,0,PreSel_Correction_Agl[S] -> GetBinError(i+1));
				PreSel_CorrectionFit_Agl_Graph[S]->SetPoint(j,AglPB.EkPerMassBinCent(i),PreSel_CorrectionFit_Agl[S] -> GetBinContent(i+1));
				PreSel_CorrectionFit_Agl_Graph[S]->SetPointError(j,0,PreSel_CorrectionFit_Agl[S] -> GetBinError(i+1));

				j++;
			}
		}
		PreSel_Correction_Agl_Graph[S]->SetLineColor(2);
		PreSel_CorrectionFit_Agl_Graph[S]->SetFillColor(2);
		PreSel_CorrectionFit_Agl_Graph[S]->SetFillStyle(3001);
		PreSel_CorrectionFit_Agl_Graph[S]->SetLineWidth(4);
		PreSel_Correction_Agl_Graph[S]->SetMarkerColor(2);
		PreSel_Correction_Agl_Graph[S]->SetMarkerStyle(8);
		PreSel_Correction_Agl_Graph[S]->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
		PreSel_Correction_Agl_Graph[S]->GetYaxis()->SetTitle("Efficiency (Data/MC)");
		PreSel_Correction_Agl_Graph[S]->SetTitle("R binning");
		PreSel_Correction_Agl_Graph[S]->GetYaxis()->SetRangeUser(0.9,1.2);

		PreSel_Correction_Agl_Graph[S]->GetXaxis()->SetTitleSize(0.045);
                PreSel_Correction_Agl_Graph[S]->GetYaxis()->SetTitleSize(0.045);

		PreSel_Correction_Agl_Graph[S]->Draw("AP");
		PreSel_CorrectionFit_Agl_Graph[S]->Draw("P4Csame");
		{
			TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
			leg->AddEntry(PreSel_Correction_Agl_Graph[S],"Efficiency correction", "ep");
			leg->AddEntry(PreSel_CorrectionFit_Agl_Graph[S],"Param.", "l");
			leg->SetLineWidth(2);
			leg->Draw("same");
		}














	}

	for(int S=0;S<3;S++) finalPlots.Add(c21[S]);
        finalPlots.writeObjsInFolder("DATA-driven Results/Data vs MC/Selections/Corrections");





	
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
