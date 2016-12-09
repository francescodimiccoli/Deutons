

void	DATApreSeleff_Plot(TH1 * LATpreSelDATA_R  ,
                           TH1 * preSelLATcorr    ,
                           TH1 * preSelLATcorr_fit
	){

	string tagli[3]= {"Matching TOF","Chi^2 R","1 Tr. Track"};
	string nome;
	string Legend[11]={"Lat. Zone 0","Lat. Zone 1","Lat. Zone 2","Lat. Zone 3","Lat. Zone 4","Lat. Zone 5","Lat. Zone 6","Lat. Zone 7","Lat. Zone 8","Lat. Zone 9","Lat. Zone 10" };

	TLegend* leig[3];
	TGraphErrors * Eff_preSelLAT[3][11];
	TCanvas *c14[4];
	for(int S=0; S<3; S++) {
		nome="Latitude Efficiency: "+tagli[S];
		c14[S]=new TCanvas(nome.c_str());
		c14[S]->Divide(2,1);
		c14[S]->cd(1);
		gPad->SetLogx();
		gPad->SetGridx();
		gPad->SetGridy();
		for(int l=0; l<11; l++) {
			Eff_preSelLAT[S][l]=new TGraphErrors();
			int point=0;
			for(int i=1; i<nbinsr; i++) {
					if(PRB.RigBinCent(i)>Rcut[l]){
					Eff_preSelLAT[S][l]->SetPoint(point,PRB.RigBinCent(i),LATpreSelDATA_R->GetBinContent(i+1,l+1,S+1));
					Eff_preSelLAT[S][l]->SetPointError(point,0,LATpreSelDATA_R->GetBinError(i+1,l+1,S+1));
					point++;
					}
			}
		}
		Eff_preSelLAT[S][10]->SetMarkerColor(1);
		Eff_preSelLAT[S][10]->SetLineColor(1);
		Eff_preSelLAT[S][10]->SetMarkerStyle(8);
		Eff_preSelLAT[S][10]->GetXaxis()->SetTitle("R [GV]");
		Eff_preSelLAT[S][10]->GetYaxis()->SetTitle("Efficiency");
		Eff_preSelLAT[S][10]->GetYaxis()->SetRangeUser(0.1,1);
		Eff_preSelLAT[S][10]->Draw("AP");
		for(int l=0; l<10; l++) {
			Eff_preSelLAT[S][l]->SetMarkerColor(55+2*l);
			Eff_preSelLAT[S][l]->SetLineColor(55+2*l);
			Eff_preSelLAT[S][l]->SetMarkerStyle(8);
			Eff_preSelLAT[S][l]->Draw("Psame");
		}
		{
			leig[S] =new TLegend(0.8, 0.1,0.98,0.95);
			for (int l=0; l<11; l++) leig[S]->AddEntry(Eff_preSelLAT[S][l],Legend[l].c_str(), "p");
			leig[S]->SetLineWidth(2);
			leig[S]->Draw("same");
		}

	}



	TGraphErrors *CorrLATpre[3];
	TGraphErrors *CorrLATpre_Spl[3];

	for(int S=0; S<3; S++) {
		c14[S]->cd(2);
		gPad->SetGridy();
		gPad->SetGridx();
		nome="Latitude Efficiency: "+tagli[S];
		CorrLATpre[S]=new TGraphErrors();
		CorrLATpre[S]->SetTitle("Latitude Efficiency Corr.");
		CorrLATpre[S]->GetXaxis()->SetTitle("Latitude");
		CorrLATpre[S]->GetYaxis()->SetTitle("Eff. Corr. Factor");
		CorrLATpre[S]->SetMarkerStyle(8);
		for(int i=1; i<11; i++) {
			CorrLATpre[S]->SetPoint(i-1,geomagC[i],preSelLATcorr->GetBinContent(i+1,S+1));
			CorrLATpre[S]->SetPointError(i-1,0,preSelLATcorr->GetBinError(i+1,S+1));
		}
		//CorrLATpre[S]->Fit(nome.c_str());
		CorrLATpre[S]->Draw("AP");

		nome="CorrLATpre_spl"+tagli[S];
		CorrLATpre_Spl[S]=new TGraphErrors(11);
		CorrLATpre_Spl[S]->SetName(tagli[S].c_str());
		int j=0;
		for(int i=1; i<11; i++) {
			CorrLATpre_Spl[S]->SetPoint(j,geomagC[i],preSelLATcorr_fit->GetBinContent(i+1,S+1));
			CorrLATpre_Spl[S]->SetPointError(j,0,preSelLATcorr_fit->GetBinError(i+1,S+1));
			j++;
		}
		CorrLATpre_Spl[S]->SetLineColor(2);
		CorrLATpre_Spl[S]->SetMarkerColor(2);
		CorrLATpre_Spl[S]->SetFillColor(2);
		CorrLATpre_Spl[S]->SetFillStyle(3001);
		CorrLATpre_Spl[S]->SetLineWidth(2);
		CorrLATpre_Spl[S]->Draw("PCsame");
	}


	for(int S=0; S<3; S++) finalPlots.Add(c14[S]);
	finalPlots.writeObjsInFolder("DATA-driven Results/Latitude effect/Clean-event Selections");
	for(int S=0; S<3; S++) finalPlots.Add(CorrLATpre_Spl[S]);
	finalPlots.writeObjsInFolder("Export");


}
