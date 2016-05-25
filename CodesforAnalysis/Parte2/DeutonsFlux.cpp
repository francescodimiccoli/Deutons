using namespace std;


void DeutonFlux() {
	string nomefile="../Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
	TFile * file1 = TFile::Open(nomefile.c_str(),"READ");

	TH2F * esposizionegeo_R    = (TH2F*)file1->Get(      "esposizionegeo"        );
	TH2F * esposizionedgeoTOF  = (TH2F*)file1->Get(      "esposizionedgeo"       );
	TH2F * esposizionedgeoNaF  = (TH2F*)file1->Get(      "esposizionedgeoNaF"    );
	TH2F * esposizionedgeoAgl  = (TH2F*)file1->Get(      "esposizionedgeoAgl"    );
	TH2F * esposizionepgeoTOF  = (TH2F*)file1->Get(      "esposizionepgeo"       );
	TH2F * esposizionepgeoNaF  = (TH2F*)file1->Get(      "esposizionepgeoNaF"    );
	TH2F * esposizionepgeoAgl  = (TH2F*)file1->Get(      "esposizionepgeoAgl"    );

	Flux * D_Flux         = new Flux(file1, "D_Flux"         ,"Results","Corr_AcceptanceD",1);
	Flux * D_Flux_geo     = new Flux(file1, "D_Flux_geo"     ,"Results","Geomag_AcceptanceD",11);

	Flux * D_Flux_Dist         = new Flux(file1, "D_Flux_Dist"         ,"Results","Corr_AcceptanceD",1);
	Flux * D_Flux_geo_Dist     = new Flux(file1, "D_Flux_geo_Dist"     ,"Results","Geomag_AcceptanceD",11);

	Flux * P_Flux         = new Flux(file1, "P_Flux"         ,"Results","Corr_AcceptanceP",1);
	Flux * P_Flux_Dist    = new Flux(file1, "P_Flux_Dist"    ,"Results","Corr_AcceptanceP",1);

	cout<<"*************** DEUTONS FLUXES CALCULATION *******************"<<endl;
	
	//Evaluation additional Fit Error
	
	TH1F * SystR  ;
        TH1F * SystTOF;
        TH1F * SystNaF;
        TH1F * SystAgl;
	
	if(D_Flux -> Counts_R  ) SystR   = (TH1F*) D_Flux -> Counts_R   ->Clone();
	if(D_Flux -> Counts_TOF) SystTOF = (TH1F*) D_Flux -> Counts_TOF ->Clone();	
	if(D_Flux -> Counts_NaF) SystNaF = (TH1F*) D_Flux -> Counts_NaF ->Clone();
	if(D_Flux -> Counts_Agl) SystAgl = (TH1F*) D_Flux -> Counts_Agl ->Clone();
	
	if(D_Flux -> Counts_R  )SystR   -> Add(	(TH1F*) D_Flux_Dist -> Counts_R  ->Clone()	,-1); 
        if(D_Flux -> Counts_TOF)SystTOF -> Add(	(TH1F*) D_Flux_Dist -> Counts_TOF->Clone()	,-1);
        if(D_Flux -> Counts_NaF)SystNaF -> Add(	(TH1F*) D_Flux_Dist -> Counts_NaF->Clone()	,-1);
        if(D_Flux -> Counts_Agl)SystAgl -> Add(	(TH1F*) D_Flux_Dist -> Counts_Agl->Clone()	,-1);

	D_Flux      -> Add_SystFitError(1,SystR,SystTOF,SystNaF,SystAgl);
	D_Flux_Dist -> Add_SystFitError(1,SystR,SystTOF,SystNaF,SystAgl);

	//

	D_Flux         -> Set_Exposure_Time (esposizionegeo_R,esposizionedgeoTOF,esposizionedgeoNaF,esposizionedgeoAgl);
	D_Flux_geo     -> Set_Exposure_Time (Tempi);

	D_Flux_Dist    -> Set_Exposure_Time (esposizionegeo_R,esposizionedgeoTOF,esposizionedgeoNaF,esposizionedgeoAgl);
	D_Flux_geo_Dist-> Set_Exposure_Time (Tempi);

	P_Flux         -> Set_Exposure_Time (esposizionegeo_R,esposizionepgeoTOF,esposizionepgeoNaF,esposizionepgeoAgl);
	P_Flux_Dist    -> Set_Exposure_Time (esposizionegeo_R,esposizionepgeoTOF,esposizionepgeoNaF,esposizionepgeoAgl);

	enum {protons, deutons};

	D_Flux         -> Eval_Flux(1 , deutons, 2 );
	D_Flux_geo     -> Eval_Flux(11, deutons, 2 );

	D_Flux_Dist    -> Eval_Flux(1 , deutons, 2 );
	D_Flux_geo_Dist-> Eval_Flux(11, deutons, 2 );

	P_Flux         -> Eval_Flux(1 , protons );
	P_Flux_Dist    -> Eval_Flux(1 , protons );


	//Fit on Mass
	TH1F * DeutonsPrimaryFlux_TOF 		= (TH1F *)D_Flux     -> Flux_TOF ;
	TH1F * DeutonsPrimaryFlux_NaF 		= (TH1F *)D_Flux     -> Flux_NaF ;
	TH1F * DeutonsPrimaryFlux_Agl 		= (TH1F *)D_Flux     -> Flux_Agl ;

	TH2F * DeutonsGeomagFlux_TOF  		= (TH2F *)D_Flux_geo -> Flux_TOF ;
	TH2F * DeutonsGeomagFlux_NaF 		= (TH2F *)D_Flux_geo -> Flux_NaF ;
	TH2F * DeutonsGeomagFlux_Agl  		= (TH2F *)D_Flux_geo -> Flux_Agl ;

	//Fit on Distance
	TH1F * DeutonsPrimaryFlux_Dist_TOF 	= (TH1F *)D_Flux_Dist-> Flux_TOF ;
	TH1F * DeutonsPrimaryFlux_Dist_NaF 	= (TH1F *)D_Flux_Dist-> Flux_NaF ;
	TH1F * DeutonsPrimaryFlux_Dist_Agl 	= (TH1F *)D_Flux_Dist-> Flux_Agl ;

	TH2F * DeutonsGeomagFlux_Dist_TOF  	= (TH2F *)D_Flux_geo_Dist-> Flux_TOF ;
	TH2F * DeutonsGeomagFlux_Dist_NaF 	= (TH2F *)D_Flux_geo_Dist-> Flux_NaF ;
	TH2F * DeutonsGeomagFlux_Dist_Agl  	= (TH2F *)D_Flux_geo_Dist-> Flux_Agl ;

	//for D/P ratio: P Flux
	//Fit on Mass
	TH1F * ProtonsPrimaryFlux_TOF 		= (TH1F *)P_Flux     -> Flux_TOF ;
	TH1F * ProtonsPrimaryFlux_NaF 		= (TH1F *)P_Flux     -> Flux_NaF ;
	TH1F * ProtonsPrimaryFlux_Agl 		= (TH1F *)P_Flux     -> Flux_Agl ;

	//Fit on Distance
	TH1F * ProtonsPrimaryFlux_Dist_TOF 	= (TH1F *)P_Flux_Dist-> Flux_TOF ;
	TH1F * ProtonsPrimaryFlux_Dist_NaF 	= (TH1F *)P_Flux_Dist-> Flux_NaF ;
	TH1F * ProtonsPrimaryFlux_Dist_Agl 	= (TH1F *)P_Flux_Dist-> Flux_Agl ;

	//D/P ratio
	//Fit on Mass	
	TH1F * DP_ratioTOF = (TH1F *)D_Flux     -> Flux_TOF -> Clone();
	TH1F * DP_ratioNaF = (TH1F *)D_Flux     -> Flux_NaF -> Clone();
	TH1F * DP_ratioAgl = (TH1F *)D_Flux     -> Flux_Agl -> Clone();

	DP_ratioTOF ->  Divide((TH1F *)P_Flux     -> Flux_TOF			);
	DP_ratioNaF ->  Divide((TH1F *)P_Flux     -> Flux_NaF			);
	DP_ratioAgl ->  Divide((TH1F *)P_Flux     -> Flux_Agl			);

	//Fit on Distance
	TH1F * DP_ratioTOF_Dist = (TH1F *)D_Flux_Dist-> Flux_TOF -> Clone();
	TH1F * DP_ratioNaF_Dist = (TH1F *)D_Flux_Dist-> Flux_NaF -> Clone();
	TH1F * DP_ratioAgl_Dist = (TH1F *)D_Flux_Dist-> Flux_Agl -> Clone();

	DP_ratioTOF_Dist ->  Divide((TH1F *)P_Flux_Dist-> Flux_TOF                   );
	DP_ratioNaF_Dist ->  Divide((TH1F *)P_Flux_Dist-> Flux_NaF                   );
	DP_ratioAgl_Dist ->  Divide((TH1F *)P_Flux_Dist-> Flux_Agl                   );	


	cout<<"*** Updating P1 file ****"<<endl;

	nomefile="../Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
	file1 = TFile::Open(nomefile.c_str(),"UPDATE");
	file1->cd("Results/Fluxes");	

	DeutonsPrimaryFlux_TOF 	    	->Write("DeutonsPrimaryFlux_TOF"); 
	DeutonsPrimaryFlux_NaF 	  	->Write("DeutonsPrimaryFlux_NaF");
	DeutonsPrimaryFlux_Agl 	  	->Write("DeutonsPrimaryFlux_Agl");

	DeutonsGeomagFlux_TOF  	  	->Write("DeutonsGeomagFlux_TOF");
	DeutonsGeomagFlux_NaF 	  	->Write("DeutonsGeomagFlux_NaF");
	DeutonsGeomagFlux_Agl  	  	->Write("DeutonsGeomagFlux_Agl");

	DeutonsPrimaryFlux_Dist_TOF	->Write("DeutonsPrimaryFlux_Dist_TOF");
	DeutonsPrimaryFlux_Dist_NaF	->Write("DeutonsPrimaryFlux_Dist_NaF");
	DeutonsPrimaryFlux_Dist_Agl	->Write("DeutonsPrimaryFlux_Dist_Agl");

	DeutonsGeomagFlux_Dist_TOF 	->Write("DeutonsGeomagFlux_Dist_TOF");
	DeutonsGeomagFlux_Dist_NaF 	->Write("DeutonsGeomagFlux_Dist_NaF");
	DeutonsGeomagFlux_Dist_Agl 	->Write("DeutonsGeomagFlux_Dist_Agl");

	ProtonsPrimaryFlux_TOF 	  	->Write("ProtonsPrimaryFlux_TOF");
	ProtonsPrimaryFlux_NaF 	  	->Write("ProtonsPrimaryFlux_NaF");
	ProtonsPrimaryFlux_Agl 	  	->Write("ProtonsPrimaryFlux_Agl");

	ProtonsPrimaryFlux_Dist_TOF	->Write("ProtonsPrimaryFlux_Dist_TOF");
	ProtonsPrimaryFlux_Dist_NaF	->Write("ProtonsPrimaryFlux_Dist_NaF");
	ProtonsPrimaryFlux_Dist_Agl	->Write("ProtonsPrimaryFlux_Dist_Agl");

	DP_ratioTOF 	 -> Write("DP_ratioTOF"		);
	DP_ratioNaF 	 -> Write("DP_ratioNaF"		);
	DP_ratioAgl 	 -> Write("DP_ratioAgl"		);	

	DP_ratioTOF_Dist -> Write("DP_ratioTOF_Dist"	);
	DP_ratioNaF_Dist -> Write("DP_ratioNaF_Dist"	);
	DP_ratioAgl_Dist -> Write("DP_ratioAgl_Dist"	);

	file1->Write();
	file1->Close();	

	TCanvas * c32 = new TCanvas("Deutons Flux: Geo. Zones");
	TCanvas * c33 = new TCanvas("Exposure Time");
	TCanvas * c34 = new TCanvas("Deutons Flux: Primaries");
	TCanvas * c35 = new TCanvas("D/P ratio");

	float potenza=0;
	c33->Divide(1,3);
	gPad->SetGridx();
	gPad->SetGridy();
	TGraphErrors * esposd_TOF=new TGraphErrors();
	TGraphErrors * esposd_NaF=new TGraphErrors();
	TGraphErrors * esposd_Agl=new TGraphErrors();
	TGraphErrors * esposp_TOF=new TGraphErrors();
	TGraphErrors * esposp_NaF=new TGraphErrors();
	TGraphErrors * esposp_Agl=new TGraphErrors();
	for(int m=0;m<nbinsbeta;m++){
		esposd_TOF->SetPoint(m,ToFPB.EkBinCent(m),D_Flux -> Exposure_TOF -> GetBinContent(m+1));
		esposp_TOF->SetPoint(m,ToFPB.EkBinCent(m),P_Flux -> Exposure_TOF -> GetBinContent(m+1));
		esposd_NaF->SetPoint(m,NaFPB.EkBinCent(m),D_Flux -> Exposure_NaF -> GetBinContent(m+1));
		esposp_NaF->SetPoint(m,NaFPB.EkBinCent(m),P_Flux -> Exposure_NaF -> GetBinContent(m+1));
		esposd_Agl->SetPoint(m,AglPB.EkBinCent(m),D_Flux -> Exposure_Agl -> GetBinContent(m+1));
		esposp_Agl->SetPoint(m,AglPB.EkBinCent(m),P_Flux -> Exposure_Agl -> GetBinContent(m+1));
	}
	c33->cd(1);
	gPad->SetLogy();
	gPad->SetGridx();
	gPad->SetGridy();
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
	esposd_TOF->Draw("APC");
	esposp_TOF->Draw("PCsame");
	c33->cd(2);
	gPad->SetLogy();
	gPad->SetGridx();
	gPad->SetGridy();
	esposp_NaF->SetMarkerStyle(8);
	esposp_NaF->SetMarkerColor(2);
	esposp_NaF->SetLineColor(2);
	esposd_NaF->SetMarkerStyle(8);
	esposd_NaF->SetMarkerColor(4);
	esposd_NaF->SetLineColor(4);
	esposd_NaF->GetXaxis()->SetTitle("kin. En./nucl.");
	esposd_NaF->GetYaxis()->SetTitle("Exposure Time");
	esposd_NaF->GetXaxis()->SetTitleSize(0.05);
	esposd_NaF->GetYaxis()->SetTitleSize(0.05);
	esposd_NaF->SetTitle("RICH NaF range");
	esposd_NaF->Draw("APC");
	esposp_NaF->Draw("PCsame");
	c33->cd(3);
	gPad->SetLogy();
	gPad->SetGridx();
	gPad->SetGridy();
	esposp_Agl->SetMarkerStyle(8);
	esposp_Agl->SetMarkerColor(2);
	esposp_Agl->SetLineColor(2);
	esposd_Agl->SetMarkerStyle(8);
	esposd_Agl->SetMarkerColor(4);
	esposd_Agl->SetLineColor(4);
	esposd_Agl->GetXaxis()->SetTitle("kin. En./nucl.");
	esposd_Agl->GetYaxis()->SetTitle("Exposure Time");
	esposd_Agl->GetXaxis()->SetTitleSize(0.05);
	esposd_Agl->GetYaxis()->SetTitleSize(0.05);
	esposd_Agl->SetTitle("RICH Agl range");
	esposd_Agl->Draw("APC");
	esposp_Agl->Draw("PCsame");

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
		string nomefile="./Galprop/Trotta2011/Def/new_D500.txt";
		cout<<nomefile<<endl;
		ifstream fp(nomefile.c_str());
		while (!fp.eof()){
			fp>>x>>y;
			if(x/1e3>0.05&&x/1e3<=100)
				galprop3P->SetPoint(j,x/1e3,y*1e7*pow(x/1e3,potenza));
			j++;
		}
	}

	j=0;
	{
		string nomefile="./Galprop/Trotta2011/Def/new_D1250.txt";
		cout<<nomefile<<endl;
		ifstream fp(nomefile.c_str());
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
		for(int m=1;m<nbinsToF;m++){
			D_FluxgeoTOF[j]->SetPoint(p,ToFPB.EkBinCent(m),DeutonsGeomagFlux_TOF->GetBinContent(m+1,j+1));
			D_FluxgeoTOF[j]->SetPointError(p,0,DeutonsGeomagFlux_TOF->GetBinError(m+1,j+1));
			D_FluxgeoDistTOF[j]->SetPoint(p,ToFPB.EkBinCent(m),DeutonsGeomagFlux_Dist_TOF->GetBinContent(m+1,j+1));
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
		for(int m=1;m<nbinsToF;m++){
			D_FluxgeoNaF[j]->SetPoint(p,NaFPB.EkBinCent(m),DeutonsGeomagFlux_NaF->GetBinContent(m+1,j+1));
			D_FluxgeoNaF[j]->SetPointError(p,0,DeutonsGeomagFlux_NaF->GetBinError(m+1,j+1));
			D_FluxgeoDistNaF[j]->SetPoint(p,NaFPB.EkBinCent(m),DeutonsGeomagFlux_Dist_NaF->GetBinContent(m+1,j+1));
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
		for(int m=1;m<nbinsToF;m++){
			D_FluxgeoAgl[j]->SetPoint(p,AglPB.EkBinCent(m),DeutonsGeomagFlux_Agl->GetBinContent(m+1,j+1));
			D_FluxgeoAgl[j]->SetPointError(p,0,DeutonsGeomagFlux_Agl->GetBinError(m+1,j+1));
			D_FluxgeoDistAgl[j]->SetPoint(p,AglPB.EkBinCent(m),DeutonsGeomagFlux_Dist_Agl->GetBinContent(m+1,j+1));
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

	nome="Deutons Flux: Primaries" ;
	D_FluxTOF=new TGraphErrors();
	D_FluxTOF->SetName(nome.c_str());
	D_FluxDistTOF=new TGraphErrors();
	D_FluxDistTOF->SetName((nome + "Distance Fit").c_str());

	p=0;
	for(int m=1;m<nbinsToF;m++){
		D_FluxTOF->SetPoint(p,ToFPB.EkBinCent(m),DeutonsPrimaryFlux_TOF->GetBinContent(m+1));
		D_FluxTOF->SetPointError(p,0,DeutonsPrimaryFlux_TOF->GetBinError(m+1));
		D_FluxDistTOF->SetPoint(p,ToFPB.EkBinCent(m),DeutonsPrimaryFlux_Dist_TOF->GetBinContent(m+1));
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

	nome="Deutons Flux: Primaries" ;
	D_FluxNaF=new TGraphErrors();
	D_FluxNaF->SetName(nome.c_str());
	D_FluxDistNaF=new TGraphErrors();
	D_FluxDistNaF->SetName((nome + "Distance Fit").c_str());

	p=0;
	for(int m=1;m<nbinsToF;m++){
		D_FluxNaF->SetPoint(p,NaFPB.EkBinCent(m),DeutonsPrimaryFlux_NaF->GetBinContent(m+1));
		D_FluxNaF->SetPointError(p,0,DeutonsPrimaryFlux_NaF->GetBinError(m+1));
		D_FluxDistNaF->SetPoint(p,NaFPB.EkBinCent(m),DeutonsPrimaryFlux_Dist_NaF->GetBinContent(m+1));
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

	nome="Deutons Flux: Primaries" ;
	D_FluxAgl=new TGraphErrors();
	D_FluxAgl->SetName(nome.c_str());
	D_FluxDistAgl=new TGraphErrors();
	D_FluxDistAgl->SetName((nome + "Distance Fit").c_str());

	p=0;
	for(int m=1;m<nbinsToF;m++){
		D_FluxAgl->SetPoint(p,AglPB.EkBinCent(m),DeutonsPrimaryFlux_Agl->GetBinContent(m+1));
		D_FluxAgl->SetPointError(p,0,DeutonsPrimaryFlux_Agl->GetBinError(m+1));
		D_FluxDistAgl->SetPoint(p,AglPB.EkBinCent(m),DeutonsPrimaryFlux_Dist_Agl->GetBinContent(m+1));
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
		string nomefile="./Galprop/Trotta2011/PDratio/500.dat";
		cout<<nomefile<<endl;
		ifstream fp(nomefile.c_str());
		while (!fp.eof()){
			fp>>x>>y;
			if(x/1e3>0.05&&x/1e3<=100)
				galpropratio1->SetPoint(j,x/1e3,y);
			j++;
		}
	}

	j=0;
	{
		string nomefile="./Galprop/Trotta2011/PDratio/1000.dat";
		cout<<nomefile<<endl;
		ifstream fp(nomefile.c_str());
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
	for(int m=1;m<nbinsToF;m++){
		PD_ratioTOF->SetPoint(p,ToFPB.EkBinCent(m),DP_ratioTOF->GetBinContent(m+1));
		PD_ratioTOF_Dist->SetPoint(p,ToFPB.EkBinCent(m),DP_ratioTOF_Dist->GetBinContent(m+1));
		PD_ratioTOF->SetPointError(p,0,DP_ratioTOF->GetBinError(m+1));
		PD_ratioTOF_Dist->SetPointError(p,0,DP_ratioTOF_Dist->GetBinError(m+1));
		p++;
	}
	PD_ratioTOF->SetMarkerStyle(8);
	PD_ratioTOF->SetMarkerSize(1.5);
	PD_ratioTOF->SetMarkerColor(4);
	PD_ratioTOF->SetLineColor(4);
	PD_ratioTOF->SetLineWidth(2);

	PD_ratioTOF_Dist->SetMarkerStyle(8);
	PD_ratioTOF_Dist->SetMarkerSize(1.5);
	PD_ratioTOF_Dist->SetMarkerColor(4);
	PD_ratioTOF_Dist->SetLineColor(4);
	PD_ratioTOF_Dist->SetLineWidth(2);

	c35->cd(1);
	PD_ratioTOF->Draw("Psame"); 

	c35->cd(2);
	PD_ratioTOF_Dist->Draw("Psame");


	TGraphErrors * PD_ratioNaF=new TGraphErrors();
	TGraphErrors * PD_ratioNaF_Dist=new TGraphErrors();
	p=0;
	for(int m=1;m<nbinsToF;m++){
		PD_ratioNaF->SetPoint(p,NaFPB.EkBinCent(m),DP_ratioNaF->GetBinContent(m+1));
		PD_ratioNaF_Dist->SetPoint(p,NaFPB.EkBinCent(m),DP_ratioNaF_Dist->GetBinContent(m+1));
		PD_ratioNaF->SetPointError(p,0,DP_ratioNaF->GetBinError(m+1));
		PD_ratioNaF_Dist->SetPointError(p,0,DP_ratioNaF_Dist->GetBinError(m+1));
		p++;
	}
	PD_ratioNaF->SetMarkerStyle(4);
	PD_ratioNaF->SetMarkerSize(1.5);
	PD_ratioNaF->SetMarkerColor(4);
	PD_ratioNaF->SetLineColor(4);
	PD_ratioNaF->SetLineWidth(2);

	PD_ratioNaF_Dist->SetMarkerStyle(4);
	PD_ratioNaF_Dist->SetMarkerSize(1.5);
	PD_ratioNaF_Dist->SetMarkerColor(4);
	PD_ratioNaF_Dist->SetLineColor(4);
	PD_ratioNaF_Dist->SetLineWidth(2);

	c35->cd(1);
	PD_ratioNaF->Draw("Psame");

	c35->cd(2);
	PD_ratioNaF_Dist->Draw("Psame");

	TGraphErrors * PD_ratioAgl=new TGraphErrors();
	TGraphErrors * PD_ratioAgl_Dist=new TGraphErrors();
	p=0;
	for(int m=1;m<nbinsToF;m++){
		PD_ratioAgl->SetPoint(p,AglPB.EkBinCent(m),DP_ratioAgl->GetBinContent(m+1));
		PD_ratioAgl_Dist->SetPoint(p,AglPB.EkBinCent(m),DP_ratioAgl_Dist->GetBinContent(m+1));
		PD_ratioAgl->SetPointError(p,0,DP_ratioAgl->GetBinError(m+1));
		PD_ratioAgl_Dist->SetPointError(p,0,DP_ratioAgl_Dist->GetBinError(m+1));
		p++;
	}
	PD_ratioAgl->SetMarkerStyle(3);
	PD_ratioAgl->SetMarkerSize(1.5);
	PD_ratioAgl->SetMarkerColor(4);
	PD_ratioAgl->SetLineColor(4);
	PD_ratioAgl->SetLineWidth(2);

	PD_ratioAgl_Dist->SetMarkerStyle(3);
	PD_ratioAgl_Dist->SetMarkerSize(1.5);
	PD_ratioAgl_Dist->SetMarkerColor(4);
	PD_ratioAgl_Dist->SetLineColor(4);
	PD_ratioAgl_Dist->SetLineWidth(2);

	c35->cd(1);
	PD_ratioAgl->Draw("Psame");

	c35->cd(2);
	PD_ratioAgl_Dist->Draw("Psame");
	
	cout<<"*** Updating Results file ***"<<endl;
	nomefile="./Final_plots/"+mese+".root";
	TFile *f_out=new TFile(nomefile.c_str(), "UPDATE");
	f_out->mkdir("D Fluxes");
	f_out->cd("D Fluxes");
	c33-> Write();
	c32->Write();
	c34-> Write();
	c35->Write();
	f_out->cd("Export");
	D_FluxTOF->Write("Deutons Primary Flux: TOF");
	D_FluxNaF->Write("Deutons Primary Flux: NaF");
	D_FluxAgl->Write("Deutons Primary Flux: Agl");
	f_out->Write();
	f_out->Close();



	return;
}
