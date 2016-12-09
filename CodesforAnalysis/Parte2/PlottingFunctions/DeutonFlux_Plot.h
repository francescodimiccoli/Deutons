

void 	DeutonFlux_Plot(TH1 *DeutonsPrimaryFlux_TOF 	   ,
                        TH1 *DeutonsPrimaryFlux_NaF 	   ,
                        TH1 *DeutonsPrimaryFlux_Agl 	   ,
                        TH1 *DeutonsGeomagFlux_TOF  	   ,
                        TH1 *DeutonsGeomagFlux_NaF 	   ,
                        TH1 *DeutonsGeomagFlux_Agl  	   ,
                        TH1 *DeutonsPrimaryFlux_Dist_TOF,
                        TH1 *DeutonsPrimaryFlux_Dist_NaF,
                        TH1 *DeutonsPrimaryFlux_Dist_Agl,
                        TH1 *DeutonsGeomagFlux_Dist_TOF ,
                        TH1 *DeutonsGeomagFlux_Dist_NaF ,
                        TH1 *DeutonsGeomagFlux_Dist_Agl ,
                        TH1 *ProtonsPrimaryFlux_TOF 	   ,
                        TH1 *ProtonsPrimaryFlux_NaF 	   ,
                        TH1 *ProtonsPrimaryFlux_Agl 	   ,
                        TH1 *ProtonsPrimaryFlux_Dist_TOF,
                        TH1 *ProtonsPrimaryFlux_Dist_NaF,
                        TH1 *ProtonsPrimaryFlux_Dist_Agl,
			TH1 *DP_ratioTOF 		,
			TH1 *DP_ratioNaF 	,
                        TH1 *DP_ratioAgl 	,
                        TH1 * DP_ratioTOF_Dist,
                        TH1 *DP_ratioNaF_Dist,
                        TH1 *DP_ratioAgl_Dist,
			Flux * D_Flux,
			Flux * P_Flux,
			TH1 *ProtonsPrimaryFlux,
			TH2F * Corr_AcceptanceD_TOF,	
			TH2F * Corr_AcceptanceD_NaF,
			TH2F * Corr_AcceptanceD_Agl,
			
			TH1F * CountsTOF,
                        TH1F * CountsNaF,
                        TH1F * CountsAgl,
			
			TH1F *DStatTOF,
			TH1F *DStatNaF,
                        TH1F *DStatAgl,
                                
                        TH1F *DSystTOF,
                        TH1F *DSystNaF,
	                TH1F *DSystAgl


	){

	TCanvas * c32 = new TCanvas("Deutons Flux: Geo. Zones");
	TCanvas * c33 = new TCanvas("Exposure Time");
	TCanvas * c36 = new TCanvas("Deuterons Primaries counts");
	TCanvas * c34 = new TCanvas("Deuterons Flux: Primaries");
	TCanvas * c31 = new TCanvas("Protons Flux (analysis vs R bins)");
	TCanvas * c35 = new TCanvas("D/P ratio");

	float potenza=0;
	c33->Divide(1,2);
	gPad->SetGridx();
	gPad->SetGridy();

	TGraphErrors * espos_R=new TGraphErrors();
	TGraphErrors * esposd_TOF=new TGraphErrors();
	TGraphErrors * esposd_NaF=new TGraphErrors();
	TGraphErrors * esposd_Agl=new TGraphErrors();
	TGraphErrors * esposp_TOF=new TGraphErrors();
	TGraphErrors * esposp_NaF=new TGraphErrors();
	TGraphErrors * esposp_Agl=new TGraphErrors();
	for (int R=0;R<PRB.size();R++){
		espos_R->SetPoint(R,PRB.RigBinCent(R),P_Flux-> Exposure_R -> GetBinContent(R+1));
	}

	for(int m=0;m<ToFDB.size();m++){
		esposd_TOF->SetPoint(m,ToFDB.RigBinCent(m),D_Flux -> Exposure_TOF -> GetBinContent(m+1));
		esposp_TOF->SetPoint(m,ToFPB.RigBinCent(m),P_Flux -> Exposure_TOF -> GetBinContent(m+1));
		esposd_NaF->SetPoint(m,NaFDB.RigBinCent(m),D_Flux -> Exposure_NaF -> GetBinContent(m+1));
		esposp_NaF->SetPoint(m,NaFPB.RigBinCent(m),P_Flux -> Exposure_NaF -> GetBinContent(m+1));
		esposd_Agl->SetPoint(m,AglDB.RigBinCent(m),D_Flux -> Exposure_Agl -> GetBinContent(m+1));
		esposp_Agl->SetPoint(m,AglPB.RigBinCent(m),P_Flux -> Exposure_Agl -> GetBinContent(m+1));
	}
	c33->cd(1);
	gPad->SetLogy();
	gPad->SetGridx();
	gPad->SetGridy();
	TH2F * Frame = new TH2F("Frame","Frame",1000,0.1,20,1e3,1e4,1e7);	
	Frame->GetXaxis()->SetTitle("kin. En./nucl.");
	Frame->GetYaxis()->SetTitle("Exposure Time (sec)");
	Frame->Draw();
	esposp_TOF->SetMarkerStyle(8);
	esposp_TOF->SetMarkerColor(2);
	esposp_TOF->SetLineColor(2);
	esposd_TOF->SetMarkerStyle(8);
	esposd_TOF->SetMarkerColor(4);
	esposd_TOF->SetLineColor(4);
	esposd_TOF->GetXaxis()->SetTitle("kin. En./nucl.");
	esposd_TOF->GetYaxis()->SetTitle("Exposure Time");
	esposd_TOF->GetXaxis()->SetTitleSize(0.05);
	esposd_TOF->GetYaxis()->SetTitleSize(0.05);     
	esposd_TOF->SetTitle("TOF range");
	esposd_TOF->Draw("PCsame");
	//esposp_TOF->Draw("PCsame");
	
	esposp_NaF->SetMarkerStyle(8);
	esposp_NaF->SetMarkerColor(2);
	esposp_NaF->SetLineColor(2);
	esposd_NaF->SetMarkerStyle(4);
	esposd_NaF->SetMarkerColor(4);
	esposd_NaF->SetLineColor(4);
	esposd_NaF->GetXaxis()->SetTitle("kin. En./nucl.");
	esposd_NaF->GetYaxis()->SetTitle("Exposure Time");
	esposd_NaF->GetXaxis()->SetTitleSize(0.05);
	esposd_NaF->GetYaxis()->SetTitleSize(0.05);
	esposd_NaF->SetTitle("RICH NaF range");
	esposd_NaF->Draw("PCsame");
	//esposp_NaF->Draw("PCsame");
	
	esposp_Agl->SetMarkerStyle(8);
	esposp_Agl->SetMarkerColor(2);
	esposp_Agl->SetLineColor(2);
	esposd_Agl->SetMarkerStyle(3);
	esposd_Agl->SetMarkerColor(4);
	esposd_Agl->SetLineColor(4);
	esposd_Agl->GetXaxis()->SetTitle("kin. En./nucl.");
	esposd_Agl->GetYaxis()->SetTitle("Exposure Time");
	esposd_Agl->GetXaxis()->SetTitleSize(0.05);
	esposd_Agl->GetYaxis()->SetTitleSize(0.05);
	esposd_Agl->SetTitle("RICH Agl range");
	esposd_Agl->Draw("PCsame");
	//esposp_Agl->Draw("PCsame");

	c33->cd(2);
	gPad->SetLogy();
	gPad->SetGridx();
	gPad->SetGridy();
	gPad->SetLogx();
	espos_R->SetMarkerStyle(8);
	espos_R->SetMarkerColor(1);
	espos_R->SetLineColor(1);
	espos_R->Draw("APC");
	esposp_TOF->Draw("PCsame");
	esposp_NaF->Draw("PCsame");
	esposp_Agl->Draw("PCsame");
	esposd_TOF->Draw("PCsame");
	esposd_NaF->Draw("PCsame");
	esposd_Agl->Draw("PCsame");


	c32-> Divide(2,1);
	c32->cd(1);
	gPad->SetLogx();
	gPad->SetLogy();
	gPad->SetGridx();
	gPad->SetGridy();
	c32->cd(2);
	gPad->SetLogx();
	gPad->SetLogy();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraph* galprop3P=new TGraph();
	TGraph* galprop3P2=new TGraph();
	float x,y=0;
	int j=0;
	{
		string filename="./Galprop/Tom/deut_1500.dat";
		cout<<filename<<endl;
		ifstream fp(filename.c_str());
		while (!fp.eof()){
			fp>>x>>y;
			if(x/1e3>0.05&&x/1e3<=100)
				galprop3P->SetPoint(j,x/1e3,y*1e7*pow(x/1e3,potenza));
			j++;
		}
	}

	j=0;
	{
		string filename="./Galprop/Tom/deut_450.dat";
		cout<<filename<<endl;
		ifstream fp(filename.c_str());
		while (!fp.eof()){
			fp>>x>>y;
			if(x/1e3>0.05&&x/1e3<=100)
				galprop3P2->SetPoint(j,x/1e3,y*1e7*pow(x/1e3,potenza));
			j++;
		}
	}
	galprop3P->GetXaxis()->SetRangeUser(0.1,10);
	galprop3P->GetYaxis()->SetRangeUser(1e-3,1e3);

	c32->cd(1);
	galprop3P->SetTitle("Deutons Flux: Geo. Zones");
	galprop3P->GetXaxis()->SetTitle("Kin.En./nucl. [GeV/nucl.]");
	galprop3P ->GetYaxis()->SetTitle("Flux [(m^2 sr GeV/nucl.)^-1]");
	galprop3P ->GetXaxis()->SetTitleSize(0.045);
	galprop3P->GetYaxis()->SetTitleSize(0.045);
	galprop3P ->GetYaxis()->SetRangeUser(1e-2,1e4);

	galprop3P->Draw("AC");
	galprop3P2->Draw("sameC");
	c32->cd(2);
	gPad->SetTitle("Deutons Flux: Geo. Zone (Distance FIT)");
	galprop3P->Draw("AC");
	galprop3P2->Draw("sameC");


	TGraphErrors * D_FluxgeoTOF[11];
	TGraphErrors * D_FluxgeoNaF[11];
	TGraphErrors * D_FluxgeoAgl[11];	
	TGraphErrors * D_FluxgeoDistTOF[11];
	TGraphErrors * D_FluxgeoDistNaF[11];
	TGraphErrors * D_FluxgeoDistAgl[11];


	string nome;
	int p=0;
	for(int j=0;j<11;j++) {
		nome="Deutons Flux: Geo. Zone "+to_string(j);
		D_FluxgeoTOF[j]=new TGraphErrors();
		D_FluxgeoTOF[j]->SetName(nome.c_str());
		D_FluxgeoDistTOF[j]=new TGraphErrors();
		D_FluxgeoDistTOF[j]->SetName((nome + "Distance Fit").c_str());

		p=0;
		for(int m=0;m<nbinsToF;m++){
			D_FluxgeoTOF[j]->SetPoint(p,ToFDB.EkPerMassBinCent(m),DeutonsGeomagFlux_TOF->GetBinContent(m+1,j+1));
			D_FluxgeoTOF[j]->SetPointError(p,0,DeutonsGeomagFlux_TOF->GetBinError(m+1,j+1));
			D_FluxgeoDistTOF[j]->SetPoint(p,ToFDB.EkPerMassBinCent(m),DeutonsGeomagFlux_Dist_TOF->GetBinContent(m+1,j+1));
			D_FluxgeoDistTOF[j]->SetPointError(p,0,DeutonsGeomagFlux_Dist_TOF->GetBinError(m+1,j+1));
			p++;
		}
		D_FluxgeoTOF[j]->SetMarkerStyle(8);
		D_FluxgeoTOF[j]->SetMarkerSize(1.5);
		D_FluxgeoTOF[j]->SetMarkerColor(j-1);
		D_FluxgeoTOF[j]->SetLineColor(j-1);
		D_FluxgeoTOF[j]->SetLineWidth(2);

		D_FluxgeoDistTOF[j]->SetMarkerStyle(8);
		D_FluxgeoDistTOF[j]->SetMarkerSize(1.5);
		D_FluxgeoDistTOF[j]->SetMarkerColor(j-1);
		D_FluxgeoDistTOF[j]->SetLineColor(j-1);
		D_FluxgeoDistTOF[j]->SetLineWidth(2);

	}
	c32->cd(1);	
	D_FluxgeoTOF[10]->Draw("Psame");
	for(int j=0;j<11;j++) D_FluxgeoTOF[j]->Draw("Psame"); 

	c32->cd(2);
	D_FluxgeoDistTOF[10]->Draw("Psame");
	for(int j=0;j<11;j++) D_FluxgeoDistTOF[j]->Draw("Psame");

	p=0;
	for(int j=0;j<11;j++) {
		nome="Deutons Flux: Geo. Zone "+to_string(j);
		D_FluxgeoNaF[j]=new TGraphErrors();
		D_FluxgeoNaF[j]->SetName(nome.c_str());
		D_FluxgeoDistNaF[j]=new TGraphErrors();
		D_FluxgeoDistNaF[j]->SetName((nome + "Distance Fit").c_str());

		p=0;
		for(int m=0;m<nbinsToF;m++){
			D_FluxgeoNaF[j]->SetPoint(p,NaFDB.EkPerMassBinCent(m),DeutonsGeomagFlux_NaF->GetBinContent(m+1,j+1));
			D_FluxgeoNaF[j]->SetPointError(p,0,DeutonsGeomagFlux_NaF->GetBinError(m+1,j+1));
			D_FluxgeoDistNaF[j]->SetPoint(p,NaFDB.EkPerMassBinCent(m),DeutonsGeomagFlux_Dist_NaF->GetBinContent(m+1,j+1));
			D_FluxgeoDistNaF[j]->SetPointError(p,0,DeutonsGeomagFlux_Dist_NaF->GetBinError(m+1,j+1));
			p++;
		}
		D_FluxgeoNaF[j]->SetMarkerStyle(4);
		D_FluxgeoNaF[j]->SetMarkerSize(1.5);
		D_FluxgeoNaF[j]->SetMarkerColor(j-1);
		D_FluxgeoNaF[j]->SetLineColor(j-1);
		D_FluxgeoNaF[j]->SetLineWidth(2);

		D_FluxgeoDistNaF[j]->SetMarkerStyle(4);
		D_FluxgeoDistNaF[j]->SetMarkerSize(1.5);
		D_FluxgeoDistNaF[j]->SetMarkerColor(j-1);
		D_FluxgeoDistNaF[j]->SetLineColor(j-1);
		D_FluxgeoDistNaF[j]->SetLineWidth(2);

	}
	c32->cd(1);
	D_FluxgeoNaF[10]->Draw("Psame");
	for(int j=0;j<11;j++) D_FluxgeoNaF[j]->Draw("Psame");

	c32->cd(2);
	D_FluxgeoDistNaF[10]->Draw("Psame");
	for(int j=0;j<11;j++) D_FluxgeoDistNaF[j]->Draw("Psame");

	p=0;
	for(int j=0;j<11;j++) {
		nome="Deutons Flux: Geo. Zone "+to_string(j);
		D_FluxgeoAgl[j]=new TGraphErrors();
		D_FluxgeoAgl[j]->SetName(nome.c_str());
		D_FluxgeoDistAgl[j]=new TGraphErrors();
		D_FluxgeoDistAgl[j]->SetName((nome + "Distance Fit").c_str());

		p=0;
		for(int m=0;m<nbinsToF;m++){
			D_FluxgeoAgl[j]->SetPoint(p,AglDB.EkPerMassBinCent(m),DeutonsGeomagFlux_Agl->GetBinContent(m+1,j+1));
			D_FluxgeoAgl[j]->SetPointError(p,0,DeutonsGeomagFlux_Agl->GetBinError(m+1,j+1));
			D_FluxgeoDistAgl[j]->SetPoint(p,AglDB.EkPerMassBinCent(m),DeutonsGeomagFlux_Dist_Agl->GetBinContent(m+1,j+1));
			D_FluxgeoDistAgl[j]->SetPointError(p,0,DeutonsGeomagFlux_Dist_Agl->GetBinError(m+1,j+1));
			p++;
		}
		D_FluxgeoAgl[j]->SetMarkerStyle(3);
		D_FluxgeoAgl[j]->SetMarkerSize(1.5);
		D_FluxgeoAgl[j]->SetMarkerColor(j-1);
		D_FluxgeoAgl[j]->SetLineColor(j-1);
		D_FluxgeoAgl[j]->SetLineWidth(2);

		D_FluxgeoDistAgl[j]->SetMarkerStyle(3);
		D_FluxgeoDistAgl[j]->SetMarkerSize(1.5);
		D_FluxgeoDistAgl[j]->SetMarkerColor(j-1);
		D_FluxgeoDistAgl[j]->SetLineColor(j-1);
		D_FluxgeoDistAgl[j]->SetLineWidth(2);

	}
	c32->cd(1);
	D_FluxgeoAgl[10]->Draw("Psame");
	for(int j=0;j<11;j++) D_FluxgeoAgl[j]->Draw("Psame");

	c32->cd(2);
	D_FluxgeoDistAgl[10]->Draw("Psame");
	for(int j=0;j<11;j++) D_FluxgeoDistAgl[j]->Draw("Psame");

	TGraphErrors * D_FluxTOF;
	TGraphErrors * D_FluxNaF;
	TGraphErrors * D_FluxAgl;
	TGraphErrors * D_FluxDistTOF;
	TGraphErrors * D_FluxDistNaF;
	TGraphErrors * D_FluxDistAgl;

	c34-> Divide(2,1);
	c34->cd(1);
	gPad->SetLogx();
	gPad->SetLogy();
	gPad->SetGridx();
	gPad->SetGridy();
	c34->cd(2);
	gPad->SetLogx();
	gPad->SetLogy();
	gPad->SetGridx();
	gPad->SetGridy();

	c34->cd(1);
	galprop3P->Draw("AC");
	galprop3P2->Draw("sameC");
	c34->cd(2);
	galprop3P->Draw("AC");
	galprop3P2->Draw("sameC");

	nome="Deutons Flux: Primaries TOF" ;
	D_FluxTOF=new TGraphErrors();
	D_FluxTOF->SetName(nome.c_str());
	D_FluxDistTOF=new TGraphErrors();
	D_FluxDistTOF->SetName((nome + "Distance Fit").c_str());

	p=0;
	for(int m=0;m<nbinsToF;m++){
		D_FluxTOF->SetPoint(p,ToFDB.EkPerMassBinCent(m),DeutonsPrimaryFlux_TOF->GetBinContent(m+1));
		D_FluxTOF->SetPointError(p,0,DeutonsPrimaryFlux_TOF->GetBinError(m+1));
		D_FluxDistTOF->SetPoint(p,ToFDB.EkPerMassBinCent(m),DeutonsPrimaryFlux_Dist_TOF->GetBinContent(m+1));
		D_FluxDistTOF->SetPointError(p,0,DeutonsPrimaryFlux_Dist_TOF->GetBinError(m+1));
		p++;
	}
	D_FluxTOF->SetMarkerStyle(8);
	D_FluxTOF->SetMarkerSize(1.5);
	D_FluxTOF->SetMarkerColor(j-1);
	D_FluxTOF->SetLineColor(j-1);
	D_FluxTOF->SetLineWidth(2);

	D_FluxDistTOF->SetMarkerStyle(8);
	D_FluxDistTOF->SetMarkerSize(1.5);
	D_FluxDistTOF->SetMarkerColor(j-1);
	D_FluxDistTOF->SetLineColor(j-1);
	D_FluxDistTOF->SetLineWidth(2);

	c34->cd(1);
	D_FluxTOF->Draw("Psame");

	c34->cd(2);
	D_FluxDistTOF->Draw("Psame");

	nome="Deutons Flux: Primaries NaF" ;
	D_FluxNaF=new TGraphErrors();
	D_FluxNaF->SetName(nome.c_str());
	D_FluxDistNaF=new TGraphErrors();
	D_FluxDistNaF->SetName((nome + "Distance Fit").c_str());

	p=0;
	for(int m=0;m<nbinsNaF;m++){
		D_FluxNaF->SetPoint(p,NaFDB.EkPerMassBinCent(m),DeutonsPrimaryFlux_NaF->GetBinContent(m+1));
		D_FluxNaF->SetPointError(p,0,DeutonsPrimaryFlux_NaF->GetBinError(m+1));
		D_FluxDistNaF->SetPoint(p,NaFDB.EkPerMassBinCent(m),DeutonsPrimaryFlux_Dist_NaF->GetBinContent(m+1));
		D_FluxDistNaF->SetPointError(p,0,DeutonsPrimaryFlux_Dist_NaF->GetBinError(m+1));
		p++;
	}
	D_FluxNaF->SetMarkerStyle(4);
	D_FluxNaF->SetMarkerSize(1.5);
	D_FluxNaF->SetMarkerColor(j-1);
	D_FluxNaF->SetLineColor(j-1);
	D_FluxNaF->SetLineWidth(2);

	D_FluxDistNaF->SetMarkerStyle(4);
	D_FluxDistNaF->SetMarkerSize(1.5);
	D_FluxDistNaF->SetMarkerColor(j-1);
	D_FluxDistNaF->SetLineColor(j-1);
	D_FluxDistNaF->SetLineWidth(2);

	c34->cd(1);
	D_FluxNaF->Draw("Psame");

	c34->cd(2);
	D_FluxDistNaF->Draw("Psame");

	nome="Deutons Flux: Primaries Agl" ;
	D_FluxAgl=new TGraphErrors();
	D_FluxAgl->SetName(nome.c_str());
	D_FluxDistAgl=new TGraphErrors();
	D_FluxDistAgl->SetName((nome + "Distance Fit").c_str());

	p=0;
	for(int m=0;m<nbinsToF;m++){
		D_FluxAgl->SetPoint(p,AglDB.EkPerMassBinCent(m),DeutonsPrimaryFlux_Agl->GetBinContent(m+1));
		D_FluxAgl->SetPointError(p,0,DeutonsPrimaryFlux_Agl->GetBinError(m+1));
		D_FluxDistAgl->SetPoint(p,AglDB.EkPerMassBinCent(m),DeutonsPrimaryFlux_Dist_Agl->GetBinContent(m+1));
		D_FluxDistAgl->SetPointError(p,0,DeutonsPrimaryFlux_Dist_Agl->GetBinError(m+1));
		p++;
	}
	D_FluxAgl->SetMarkerStyle(3);
	D_FluxAgl->SetMarkerSize(1.5);
	D_FluxAgl->SetMarkerColor(j-1);
	D_FluxAgl->SetLineColor(j-1);
	D_FluxAgl->SetLineWidth(2);

	D_FluxDistAgl->SetMarkerStyle(3);
	D_FluxDistAgl->SetMarkerSize(1.5);
	D_FluxDistAgl->SetMarkerColor(j-1);
	D_FluxDistAgl->SetLineColor(j-1);
	D_FluxDistAgl->SetLineWidth(2);

	c34->cd(1);
	D_FluxAgl->Draw("Psame");

	c34->cd(2);
	D_FluxDistAgl->Draw("Psame");


	c31->cd();
	gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetGridx();
        gPad->SetGridy();

	TGraphErrors * PFlux=new TGraphErrors();
        TGraphErrors * PFluxTOF=new TGraphErrors();
        TGraphErrors * PFluxNaF=new TGraphErrors();
        TGraphErrors * PFluxAgl=new TGraphErrors();

	for(int i=0; i<nbinsr; i++) {
                PFlux->SetPoint(i,PRB.EkPerMassBinCent(i),ProtonsPrimaryFlux->GetBinContent(i+1)*pow(PRB.EkPerMassBinCent(i),potenza));
                PFlux->SetPointError(i,0,ProtonsPrimaryFlux->GetBinError(i+1)*pow(PRB.EkPerMassBinCent(i),potenza));
	}

	for(int i=0; i<nbinsToF; i++) {
		PFluxTOF->SetPoint(i,ToFPB.EkPerMassBinCent(i),ProtonsPrimaryFlux_TOF->GetBinContent(i+1)*pow(ToFPB.EkPerMassBinCent(i),potenza));
		PFluxTOF->SetPointError(i,0,ProtonsPrimaryFlux_TOF->GetBinError(i+1)*pow(ToFPB.EkPerMassBinCent(i),potenza));
	}

	for(int i=0; i<nbinsNaF; i++) {
		PFluxNaF->SetPoint(i,NaFPB.EkPerMassBinCent(i),ProtonsPrimaryFlux_NaF->GetBinContent(i+1)*pow(NaFPB.EkPerMassBinCent(i),potenza));
		PFluxNaF->SetPointError(i,0,ProtonsPrimaryFlux_NaF->GetBinError(i+1)*pow(NaFPB.EkPerMassBinCent(i),potenza));
	}
	for(int i=0; i<nbinsAgl; i++) {
		PFluxAgl->SetPoint(i,AglPB.EkPerMassBinCent(i),ProtonsPrimaryFlux_Agl->GetBinContent(i+1)*pow(AglPB.EkPerMassBinCent(i),potenza));
		PFluxAgl->SetPointError(i,0,ProtonsPrimaryFlux_Agl->GetBinError(i+1)*pow(AglPB.EkPerMassBinCent(i),potenza));
	}

	PFlux->SetName("Protons Primary Flux");
	PFlux->SetMarkerStyle(8);
	PFlux->SetMarkerColor(2);
	PFluxTOF->SetMarkerStyle(8);
	PFluxTOF->SetMarkerColor(1);
	PFluxNaF->SetMarkerStyle(4);
	PFluxNaF->SetMarkerColor(1);
	PFluxAgl->SetMarkerStyle(3);
	PFluxAgl->SetMarkerColor(1);
	
	PFlux->SetTitle("Primary Protons Flux");
        PFlux->GetXaxis()->SetTitle("Kin. En./nucl. [GeV/nucl.]");
        PFlux->GetYaxis()->SetTitle("Flux [(m^2 sr GeV/nucl.)^-1]");
        PFlux->GetXaxis()->SetTitleSize(0.045);
        PFlux->GetYaxis()->SetTitleSize(0.045);
        PFlux->GetYaxis()->SetRangeUser(1e-2,1e4);
        PFlux->Draw("AP");
	PFluxTOF->Draw("Psame");
	PFluxNaF->Draw("Psame");
	PFluxAgl->Draw("Psame");




	c35->Divide(2,1);
	c35->cd(1);
	gPad->SetLogx();
	gPad->SetLogy();
	gPad->SetGridx();
	gPad->SetGridy();
	c35->cd(2);
	gPad->SetLogx();
	gPad->SetLogy();
	gPad->SetGridx();
	gPad->SetGridy();

	TGraph* galpropratio1=new TGraph();
	TGraph* galpropratio2=new TGraph();
	j=0;
	{
		string filename="./Galprop/Trotta2011/PDratio/500.dat";
		cout<<filename<<endl;
		ifstream fp(filename.c_str());
		while (!fp.eof()){
			fp>>x>>y;
			if(x/1e3>0.05&&x/1e3<=100)
				galpropratio1->SetPoint(j,x/1e3,y);
			j++;
		}
	}

	j=0;
	{
		string filename="./Galprop/Trotta2011/PDratio/1000.dat";
		cout<<filename<<endl;
		ifstream fp(filename.c_str());
		while (!fp.eof()){
			fp>>x>>y;
			if(x/1e3>0.05&&x/1e3<=100)
				galpropratio2->SetPoint(j,x/1e3,y);
			j++;
		}
	}
	galpropratio1->GetXaxis()->SetRangeUser(0.1,10);
	galpropratio1->GetYaxis()->SetRangeUser(1e-3,1e-1);
	galpropratio1->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
	galpropratio1->GetYaxis()->SetTitle("Flux ratio");
	c35->cd(1);
	galpropratio1->Draw("AC");
	galpropratio2->Draw("sameC");
	c35->cd(2);
	galpropratio1->Draw("AC");
	galpropratio2->Draw("sameC");

	TGraphErrors * PD_ratioTOF=new TGraphErrors();
	TGraphErrors * PD_ratioTOF_Dist=new TGraphErrors();
	p=0;
	for(int m=0;m<nbinsToF;m++){
		PD_ratioTOF->SetPoint(p,ToFDB.EkPerMassBinCent(m),DP_ratioTOF->GetBinContent(m+1));
		PD_ratioTOF_Dist->SetPoint(p,ToFDB.EkPerMassBinCent(m),DP_ratioTOF_Dist->GetBinContent(m+1));
		PD_ratioTOF->SetPointError(p,0,DP_ratioTOF->GetBinError(m+1));
		PD_ratioTOF_Dist->SetPointError(p,0,DP_ratioTOF_Dist->GetBinError(m+1));
		p++;
	}
	PD_ratioTOF->SetName("PD_ratioTOF");
	PD_ratioTOF->SetMarkerStyle(8);
	PD_ratioTOF->SetMarkerSize(1.5);
	PD_ratioTOF->SetMarkerColor(2);
	PD_ratioTOF->SetLineColor(2);
	PD_ratioTOF->SetLineWidth(2);

	PD_ratioTOF_Dist->SetMarkerStyle(8);
	PD_ratioTOF_Dist->SetMarkerSize(1.5);
	PD_ratioTOF_Dist->SetMarkerColor(2);
	PD_ratioTOF_Dist->SetLineColor(2);
	PD_ratioTOF_Dist->SetLineWidth(2);

	c35->cd(1);
	PD_ratioTOF->Draw("Psame"); 

	c35->cd(2);
	PD_ratioTOF_Dist->Draw("Psame");


	TGraphErrors * PD_ratioNaF=new TGraphErrors();
	TGraphErrors * PD_ratioNaF_Dist=new TGraphErrors();
	p=0;
	for(int m=0;m<nbinsToF;m++){
		PD_ratioNaF->SetPoint(p,NaFDB.EkPerMassBinCent(m),DP_ratioNaF->GetBinContent(m+1));
		PD_ratioNaF_Dist->SetPoint(p,NaFDB.EkPerMassBinCent(m),DP_ratioNaF_Dist->GetBinContent(m+1));
		PD_ratioNaF->SetPointError(p,0,DP_ratioNaF->GetBinError(m+1));
		PD_ratioNaF_Dist->SetPointError(p,0,DP_ratioNaF_Dist->GetBinError(m+1));
		p++;
	}
	PD_ratioNaF->SetName("PD_ratioNaF");
	PD_ratioNaF->SetMarkerStyle(4);
	PD_ratioNaF->SetMarkerSize(1.5);
	PD_ratioNaF->SetMarkerColor(2);
	PD_ratioNaF->SetLineColor(2);
	PD_ratioNaF->SetLineWidth(2);

	PD_ratioNaF_Dist->SetMarkerStyle(4);
	PD_ratioNaF_Dist->SetMarkerSize(1.5);
	PD_ratioNaF_Dist->SetMarkerColor(2);
	PD_ratioNaF_Dist->SetLineColor(2);
	PD_ratioNaF_Dist->SetLineWidth(2);

	c35->cd(1);
	PD_ratioNaF->Draw("Psame");

	c35->cd(2);
	PD_ratioNaF_Dist->Draw("Psame");

	TGraphErrors * PD_ratioAgl=new TGraphErrors();
	TGraphErrors * PD_ratioAgl_Dist=new TGraphErrors();
	p=0;
	for(int m=0;m<nbinsToF;m++){
		PD_ratioAgl->SetPoint(p,AglDB.EkPerMassBinCent(m),DP_ratioAgl->GetBinContent(m+1));
		PD_ratioAgl_Dist->SetPoint(p,AglDB.EkPerMassBinCent(m),DP_ratioAgl_Dist->GetBinContent(m+1));
		PD_ratioAgl->SetPointError(p,0,DP_ratioAgl->GetBinError(m+1));
		PD_ratioAgl_Dist->SetPointError(p,0,DP_ratioAgl_Dist->GetBinError(m+1));
		p++;
	}
	PD_ratioAgl->SetName("PD_ratioAgl");
	PD_ratioAgl->SetMarkerStyle(3);
	PD_ratioAgl->SetMarkerSize(1.5);
	PD_ratioAgl->SetMarkerColor(2);
	PD_ratioAgl->SetLineColor(2);
	PD_ratioAgl->SetLineWidth(2);

	PD_ratioAgl_Dist->SetMarkerStyle(3);
	PD_ratioAgl_Dist->SetMarkerSize(1.5);
	PD_ratioAgl_Dist->SetMarkerColor(2);
	PD_ratioAgl_Dist->SetLineColor(2);
	PD_ratioAgl_Dist->SetLineWidth(2);

	c35->cd(1);
	PD_ratioAgl->Draw("Psame");

	c35->cd(2);
	PD_ratioAgl_Dist->Draw("Psame");





	TCanvas *e1 = new TCanvas("Flux Errors Breakdown");
	e1-> Divide(3,1);

	e1->cd(1);
	gPad->SetGridx();
	gPad->SetGridy();
	gPad->SetLogy();

	TH1F * ErrTOTD_TOF=(TH1F *)DStatTOF->Clone();
	for(int iR=0;iR<ErrTOTD_TOF->GetNbinsX();iR++) 
		ErrTOTD_TOF->SetBinContent(iR+1, DeutonsPrimaryFlux_TOF->GetBinError(iR+1)/DeutonsPrimaryFlux_TOF->GetBinContent(iR+1));

	TH1F * ErrAccD_TOF=(TH1F *)DStatTOF->Clone();
	for(int iR=0;iR<ErrAccD_TOF->GetNbinsX();iR++) 
		ErrAccD_TOF->SetBinContent(iR+1, Corr_AcceptanceD_TOF->GetBinError(iR+1,1)/Corr_AcceptanceD_TOF->GetBinContent(iR+1,1));

	ErrTOTD_TOF->GetXaxis()->SetTitle("Bin nr.");
        ErrTOTD_TOF->GetYaxis()->SetTitle("Relative Error");
        ErrTOTD_TOF->GetXaxis()->SetTitleSize(0.045);
        ErrTOTD_TOF->GetYaxis()->SetTitleSize(0.045);
        ErrTOTD_TOF->SetTitle("Flux TOTAL Error breakdown (TOF bins)");
        ErrTOTD_TOF->GetYaxis()->SetRangeUser(1e-4,1);
	ErrTOTD_TOF->SetLineColor(1);
	ErrTOTD_TOF->SetLineWidth(7);

	ErrAccD_TOF->SetLineColor(2);
        ErrAccD_TOF->SetLineWidth(4);
        DStatTOF->SetLineColor(3);
        DStatTOF->SetLineWidth(4);
        DSystTOF->SetLineColor(4);
        DSystTOF->SetLineWidth(4);

	DSystTOF->Smooth(2);
	ErrTOTD_TOF->Smooth(2);



	ErrTOTD_TOF->Draw("hist");
        ErrAccD_TOF->Draw("hist same");
        DStatTOF   ->Draw("hist same");
        DSystTOF   ->Draw("hist same");

	{       TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(  ErrTOTD_TOF  ,"TOTAL Error", "l");
                leg->AddEntry(  ErrAccD_TOF  ,"Acceptance Error", "l");
                leg->AddEntry(  DStatTOF     ,"T. Fit Stat. Error", "l");
                leg->AddEntry(  DSystTOF     ,"T. Fit Syst. Error", "l");
                
		leg->Draw("same");

        }


	e1->cd(2);
	gPad->SetGridx();
	gPad->SetGridy();
	gPad->SetLogy();

	TH1F * ErrTOTD_NaF=(TH1F *)DStatNaF->Clone();
	for(int iR=0;iR<ErrTOTD_NaF->GetNbinsX();iR++) 
		ErrTOTD_NaF->SetBinContent(iR+1, DeutonsPrimaryFlux_NaF->GetBinError(iR+1)/DeutonsPrimaryFlux_NaF->GetBinContent(iR+1));

	TH1F * ErrAccD_NaF=(TH1F *)DStatNaF->Clone();
	for(int iR=0;iR<ErrAccD_NaF->GetNbinsX();iR++) 
		ErrAccD_NaF->SetBinContent(iR+1, Corr_AcceptanceD_NaF->GetBinError(iR+1,1)/Corr_AcceptanceD_NaF->GetBinContent(iR+1,1));

	ErrTOTD_NaF->GetXaxis()->SetTitle("Bin nr.");
        ErrTOTD_NaF->GetYaxis()->SetTitle("Relative Error");
        ErrTOTD_NaF->GetXaxis()->SetTitleSize(0.045);
        ErrTOTD_NaF->GetYaxis()->SetTitleSize(0.045);
        ErrTOTD_NaF->SetTitle("Flux TOTAL Error breakdown (NaF bins)");
        ErrTOTD_NaF->GetYaxis()->SetRangeUser(1e-4,1);
	ErrTOTD_NaF->SetLineColor(1);
	ErrTOTD_NaF->SetLineWidth(7);

	ErrAccD_NaF->SetLineColor(2);
        ErrAccD_NaF->SetLineWidth(4);
        DStatNaF->SetLineColor(3);
        DStatNaF->SetLineWidth(4);
        DSystNaF->SetLineColor(4);
        DSystNaF->SetLineWidth(4);

	DSystNaF->Smooth(2);
	ErrTOTD_NaF->Smooth(2);

	ErrTOTD_NaF->Draw("hist");
        ErrAccD_NaF->Draw("hist same");
        DStatNaF   ->Draw("hist same");
        DSystNaF   ->Draw("hist same");

	{       TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(  ErrTOTD_NaF  ,"TOTAL Error", "l");
                leg->AddEntry(  ErrAccD_NaF  ,"Acceptance Error", "l");
                leg->AddEntry(  DStatNaF     ,"T. Fit Stat. Error", "l");
                leg->AddEntry(  DSystNaF     ,"T. Fit Syst. Error", "l");
                
		leg->Draw("same");

        }


	e1->cd(3);
	gPad->SetGridx();
	gPad->SetGridy();
	gPad->SetLogy();

	TH1F * ErrTOTD_Agl=(TH1F *)DStatAgl->Clone();
	for(int iR=0;iR<ErrTOTD_Agl->GetNbinsX();iR++) 
		ErrTOTD_Agl->SetBinContent(iR+1, DeutonsPrimaryFlux_Agl->GetBinError(iR+1)/DeutonsPrimaryFlux_Agl->GetBinContent(iR+1));

	TH1F * ErrAccD_Agl=(TH1F *)DStatAgl->Clone();
	for(int iR=0;iR<ErrAccD_Agl->GetNbinsX();iR++) 
		ErrAccD_Agl->SetBinContent(iR+1, Corr_AcceptanceD_Agl->GetBinError(iR+1,1)/Corr_AcceptanceD_Agl->GetBinContent(iR+1,1));

	ErrTOTD_Agl->GetXaxis()->SetTitle("Bin nr.");
        ErrTOTD_Agl->GetYaxis()->SetTitle("Relative Error");
        ErrTOTD_Agl->GetXaxis()->SetTitleSize(0.045);
        ErrTOTD_Agl->GetYaxis()->SetTitleSize(0.045);
        ErrTOTD_Agl->SetTitle("Flux TOTAL Error breakdown (Agl bins)");
        ErrTOTD_Agl->GetYaxis()->SetRangeUser(1e-4,1);
	ErrTOTD_Agl->SetLineColor(1);
	ErrTOTD_Agl->SetLineWidth(7);

	ErrAccD_Agl->SetLineColor(2);
        ErrAccD_Agl->SetLineWidth(4);
        DStatAgl->SetLineColor(3);
        DStatAgl->SetLineWidth(4);
        DSystAgl->SetLineColor(4);
        DSystAgl->SetLineWidth(4);
	
	DSystAgl->Smooth(1);
	ErrAccD_Agl->Smooth(4);
	ErrTOTD_Agl->Smooth(2);

	ErrTOTD_Agl->Draw("hist");
        ErrAccD_Agl->Draw("hist same");
        DStatAgl   ->Draw("hist same");
        DSystAgl   ->Draw("hist same");

	{       TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(  ErrTOTD_Agl  ,"TOTAL Error", "l");
                leg->AddEntry(  ErrAccD_Agl  ,"Acceptance Error", "l");
                leg->AddEntry(  DStatAgl     ,"T. Fit Stat. Error", "l");
                leg->AddEntry(  DSystAgl     ,"T. Fit Syst. Error", "l");
                
		leg->Draw("same");

        }


	c36->Divide(3,1);
	c36->cd(1);
	
	gPad->SetGridx();
	gPad->SetGridy();

	CountsTOF->SetLineColor(4);
	CountsTOF->SetLineWidth(4);
	CountsTOF->SetMarkerStyle(8);
        CountsTOF->SetMarkerColor(4);
	CountsTOF->SetMarkerSize(2);
	CountsTOF->GetXaxis()->SetTitle("Nr. bin");
	CountsTOF->GetYaxis()->SetTitle("Extracted Counts");
	CountsTOF->Draw();

	c36->cd(2);

        gPad->SetGridx();
        gPad->SetGridy();

        CountsNaF->SetLineColor(4);
        CountsNaF->SetLineWidth(4);
	CountsNaF->SetMarkerStyle(8);
        CountsNaF->SetMarkerColor(4);
        CountsNaF->SetMarkerSize(2);
	CountsNaF->GetXaxis()->SetTitle("Nr. bin");
        CountsNaF->GetYaxis()->SetTitle("Extracted Counts");

        CountsNaF->Draw();


	c36->cd(3);

        gPad->SetGridx();
        gPad->SetGridy();

        CountsAgl->SetLineColor(4);
        CountsAgl->SetLineWidth(4);
	CountsAgl->SetMarkerStyle(8);
        CountsAgl->SetMarkerColor(4);
        CountsAgl->SetMarkerSize(2);
	CountsAgl->GetXaxis()->SetTitle("Nr. bin");
        CountsAgl->GetYaxis()->SetTitle("Extracted Counts");


        CountsAgl->Draw();







	finalPlots.Add(c32);
        finalPlots.Add(c33);
        finalPlots.Add(c36);
	finalPlots.Add(c34);
	finalPlots.Add(e1);
        finalPlots.writeObjsInFolder("D Fluxes");


	finalPlots.Add(c31);
	finalPlots.Add(c35);
	finalPlots.writeObjsInFolder("D over P ratio");

	finalPlots.Add(D_FluxTOF);
	finalPlots.Add(D_FluxNaF);
	finalPlots.Add(D_FluxAgl);
	finalPlots.writeObjsInFolder("Export/DFluxes");
	
	finalPlots.Add(DStatTOF);
	finalPlots.Add(DStatNaF);
	finalPlots.Add(DStatAgl);
	finalPlots.Add(DSystTOF);
	finalPlots.Add(DSystNaF);
	finalPlots.Add(DSystAgl);
	finalPlots.writeObjsInFolder("Export/DErr");



	finalPlots.Add(PD_ratioTOF);
        finalPlots.Add(PD_ratioNaF);
        finalPlots.Add(PD_ratioAgl);
	finalPlots.writeObjsInFolder("Export/DP_ratios");

	return;
	}
