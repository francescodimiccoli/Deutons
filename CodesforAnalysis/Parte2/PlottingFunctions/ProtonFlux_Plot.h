


void	ProtonFlux_Plot( TH1 * ProtonsPrimaryFlux,
                         TH1 * ProtonsGeomagFlux ,
                         TH1 * P_pre_PrimaryFlux ,
                         TH1 * P_sel_PrimaryFlux ,
			TH1F * Tempi
	){


	TCanvas * c23 = new TCanvas("Protons Flux: Geo. Zones");
	TCanvas * c24 = new TCanvas("Primary Protons Flux");
	TCanvas * c25 = new TCanvas("Protons Flux: Eff. corr. Test");


	TGraphErrors * P_Fluxgeo[11];
	TGraphErrors * PFlux;
	TGraphErrors * PFluxNoQ;
	TGraphErrors * PFlux_pre;
	float potenza=2.7;

	c23->cd();
	gPad->SetLogx();
	gPad->SetLogy();
	gPad->SetGridx();
	gPad->SetGridy();
	string nome;
	for(int j=0; j<11; j++) {
		nome="Protons Flux: Geo. Zone "+to_string(j);
		P_Fluxgeo[j]=new TGraphErrors();
		P_Fluxgeo[j]->SetName(nome.c_str());
		for(int i=0; i<nbinsr; i++) {
			P_Fluxgeo[j]->SetPoint(i,PRB.EkPerMassBinCent(i),ProtonsGeomagFlux->GetBinContent(i+1,j+1)*pow(PRB.EkPerMassBinCent(i),potenza));
			P_Fluxgeo[j]->SetPointError(i,0,ProtonsGeomagFlux->GetBinError(i+1,j+1)*pow(PRB.EkPerMassBinCent(i),potenza));
		}
		P_Fluxgeo[j]->SetMarkerStyle(8);
		P_Fluxgeo[j]->SetMarkerColor(j-1);
		P_Fluxgeo[j]->SetLineColor(j-1);
		P_Fluxgeo[j]->SetLineWidth(2);
	}
	P_Fluxgeo[10]->SetTitle("Protons Flux: Geo. Zones");
	P_Fluxgeo[10]->GetXaxis()->SetTitle("Kin. En./nucl. [GeV/nucl.]");
	P_Fluxgeo[10]->GetYaxis()->SetTitle("Flux [(m^2 sr GeV/nucl.)^-1]");
	P_Fluxgeo[10]->GetXaxis()->SetTitleSize(0.045);
	P_Fluxgeo[10]->GetYaxis()->SetTitleSize(0.045);
	P_Fluxgeo[10]->GetYaxis()->SetRangeUser(1e-2,1e4);
	P_Fluxgeo[10]->Draw("AP");
	for(int j=0; j<11; j++) P_Fluxgeo[j]->Draw("Psame");



	c24->cd();
	gPad->SetLogx();
	gPad->SetLogy();
	gPad->SetGridx();
	gPad->SetGridy();
	PFlux=new TGraphErrors();
	PFluxNoQ = new TGraphErrors();
	for(int i=0; i<nbinsr; i++) {
		PFlux->SetPoint(i,PRB.EkPerMassBinCent(i),ProtonsPrimaryFlux->GetBinContent(i+1)*pow(PRB.EkPerMassBinCent(i),potenza));
		PFlux->SetPointError(i,0,ProtonsPrimaryFlux->GetBinError(i+1)*pow(PRB.EkPerMassBinCent(i),potenza));
		
		PFluxNoQ->SetPoint(i,PRB.EkPerMassBinCent(i),P_pre_PrimaryFlux->GetBinContent(i+1)*pow(PRB.EkPerMassBinCent(i),potenza));
                PFluxNoQ->SetPointError(i,0,P_pre_PrimaryFlux->GetBinError(i+1)*pow(PRB.EkPerMassBinCent(i),potenza));
	}
	PFlux->SetName("Protons Primary Flux");
	PFlux->SetMarkerStyle(8);
	PFlux->SetMarkerColor(2);
	PFluxNoQ->SetMarkerStyle(8);
        PFluxNoQ->SetMarkerColor(4);
	PFlux->SetTitle("Primary Protons Flux");
	PFlux->GetXaxis()->SetTitle("Kin. En./nucl. [GeV/nucl.]");
	PFlux->GetYaxis()->SetTitle("Flux [(m^2 sr GeV/nucl.)^-1]");
	PFlux->GetXaxis()->SetTitleSize(0.045);
	PFlux->GetYaxis()->SetTitleSize(0.045);
	PFlux->GetYaxis()->SetRangeUser(1e-2,1e4);
	PFlux->Draw("AP");
	PFluxNoQ->Draw("Psame");
	TGraph* galprop3P=new TGraph();
	TGraph* galprop3P2=new TGraph();
	float x,y=0;
	int j=0;
	{
		string nomefile="./Galprop/Tom/prot_100.dat";
		ifstream fp(nomefile.c_str());
		while (!fp.eof()) {
			fp>>x>>y;
			if(x/1e3>0.05&&x/1e3<=100)
				galprop3P->SetPoint(j,x/1e3,y*1e7*pow(x/1e3,potenza));
			j++;
		}
	}

	j=0;
	{
		string nomefile="./Galprop/Tom/prot_1500.dat";
		ifstream fp(nomefile.c_str());
		while (!fp.eof()) {
			fp>>x>>y;
			if(x/1e3>0.05&&x/1e3<=100)
				galprop3P2->SetPoint(j,x/1e3,y*1e7*pow(x/1e3,potenza));
			j++;
		}
	}

	galprop3P->Draw("sameC");
	galprop3P2->Draw("sameC");

	c25->cd();
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraphErrors * PFluxpre=new TGraphErrors();
	TGraphErrors * PFlux_=new TGraphErrors();
	PFlux_pre=new TGraphErrors();
	PFlux_pre->SetName("Protons Primary Flux (only pres.)");
	int p=0;
	for(int i=1; i<nbinsr; i++) {
		PFlux_->SetPoint(p,PRB.RigBinCent(i),1);
		PFlux_pre->SetPoint(p,PRB.EkPerMassBinCent(i),1);
		PFlux_->SetPointError(p,0,(P_sel_PrimaryFlux->GetBinError(i+1,2)+P_pre_PrimaryFlux->GetBinError(i+1,2))/P_pre_PrimaryFlux->GetBinContent(i+1,1));
		p++;
	}
	p=0;
	for(int i=1; i<nbinsr; i++) {
		PFluxpre->SetPoint(p,PRB.RigBinCent(i),P_sel_PrimaryFlux->GetBinContent(i+1)/P_pre_PrimaryFlux->GetBinContent(i+1,1));
		p++;
	}
	PFluxpre->SetMarkerStyle(8);
	PFluxpre->SetMarkerColor(2);
	PFlux_->SetMarkerSize(3);
        PFluxpre->SetMarkerSize(3);
	PFlux_->SetMarkerStyle(4);
	PFlux_->SetFillStyle(3002);
	PFlux_->SetFillColor(4);
	PFlux_->SetMarkerColor(2);
	PFluxpre->SetTitle("Primary Protons Flux");
	PFluxpre->GetXaxis()->SetTitle("R [GV]");
	PFluxpre->GetYaxis()->SetTitle("Fluxes ratio (Fullset / Clean Event)");
	PFluxpre->GetXaxis()->SetTitleSize(0.045);
	PFluxpre->GetYaxis()->SetTitleSize(0.045);
	PFluxpre->GetYaxis()->SetRangeUser(0.8,1.2);
	PFluxpre->SetTitle("Efficiency Corrections test");
	PFluxpre->Draw("AP");
	PFlux_->Draw("P4same");
	TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
        leg->AddEntry(PFluxpre,"Protons Flux (Clean Event + Likelihood + Distance)", "p");
	leg->AddEntry(PFlux_,"Protons Flux (Clean Event + Control Sample cuts + mass)", "p");
	leg->AddEntry(PFlux_,"Fluxes Uncertainty (Stat. + Syst.)", "f");
	leg->Draw("same");


	finalPlots.Add(c23);
	finalPlots.Add(c24);
	finalPlots.Add(c25);
	finalPlots.writeObjsInFolder("P Fluxes");

	finalPlots.Add(PFlux);
	finalPlots.Add(Tempi);
	finalPlots.writeObjsInFolder("Export/PFluxes");


	return;
}
