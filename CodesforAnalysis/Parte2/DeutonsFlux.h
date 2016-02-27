using namespace std;

TCanvas * c32 = new TCanvas("Deutons Flux: Geo. Zones");
TCanvas * c33 = new TCanvas("Exposure Time");
TCanvas * c34 = new TCanvas("Deutons Flux: Primaries");
TCanvas * c35 = new TCanvas("D/P ratio");

TH3F * DFluxgeoTOF = new TH3F("DFluxsgeoTOF","DFluxsgeoTOF",18,0,18,11,0,11,2,0,2);
TH3F * DFluxgeoNaF = new TH3F("DFluxsgeoNaF","DFluxsgeoNaF",18,0,18,11,0,11,2,0,2);
TH3F * DFluxgeoAgl = new TH3F("DFluxsgeoAgl","DFluxsgeoAgl",18,0,18,11,0,11,2,0,2);

TH2F * DFluxTOF= new TH2F("DFluxTOF","DFluxTOF",18,0,18,2,0,2);
TH2F * DFluxNaF= new TH2F("DFluxNaF","DFluxNaF",18,0,18,2,0,2);
TH2F * DFluxAgl= new TH2F("DFluxAgl","DFluxAgl",18,0,18,2,0,2);

TH2F * PDratioTOF= new TH2F("PDratioTOF","PDratioTOF",18,0,18,2,0,2);
TH2F * PDratioNaF= new TH2F("PDratioNaF","PDratioNaF",18,0,18,2,0,2);
TH2F * PDratioAgl= new TH2F("PDratioAgl","PDratioAgl",18,0,18,2,0,2);

TGraphErrors * D_FluxgeoTOF[11];
TGraphErrors * D_FluxgeoNaF[11];
TGraphErrors * D_FluxgeoAgl[11];

TGraphErrors * D_FluxTOF;
TGraphErrors * D_FluxNaF;
TGraphErrors * D_FluxAgl;

TGraphErrors * PD_ratioTOF;
TGraphErrors * PD_ratioNaF;
TGraphErrors * PD_ratioAgl;

void DeutonFlux(TFile * file1){
	TH2F * esposizionepgeo = (TH2F*)file1->Get("esposizionepgeo");
	TH2F * esposizionepgeoNaF = (TH2F*)file1->Get("esposizionepgeoNaF");
	TH2F * esposizionepgeoAgl = (TH2F*)file1->Get("esposizionepgeoAgl");
	TH2F * esposizionedgeo = (TH2F*)file1->Get("esposizionedgeo");
	TH2F * esposizionedgeoNaF = (TH2F*)file1->Get("esposizionedgeoNaF");
	TH2F * esposizionedgeoAgl = (TH2F*)file1->Get("esposizionedgeoAgl");	

	cout<<"*************** Deutons Fluxes Calculation ******************"<<endl;
	float Esposd_TOF[18]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	float Esposd_NaF[18]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	float Esposd_Agl[18]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	float Esposp_TOF[18]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	float Esposp_NaF[18]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	float Esposp_Agl[18]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};	

	for(int m=0;m<18;m++)
		for(int i=1;i<11;i++) {
			Esposd_TOF[m]+=esposizionedgeo->GetBinContent(m+1,i);
			Esposp_TOF[m]+=esposizionepgeo->GetBinContent(m+1,i);
			Esposd_NaF[m]+=esposizionedgeoNaF->GetBinContent(m+1,i);
			Esposp_NaF[m]+=esposizionepgeoNaF->GetBinContent(m+1,i);
			Esposd_Agl[m]+=esposizionedgeoAgl->GetBinContent(m+1,i);
			Esposp_Agl[m]+=esposizionepgeoAgl->GetBinContent(m+1,i);
		}

	float AcceptLATpTOF[18]={0};
	float AcceptLATpNaF[18]={0};
	float AcceptLATpAgl[18]={0};
	float AcceptLATdTOF[18]={0};
	float AcceptLATdNaF[18]={0};
	float AcceptLATdAgl[18]={0};

	for(int m=0;m<18;m++){
		AcceptLATdTOF[m]=AcceptDzoneTOF->GetBinContent(m+1,1,2);
		AcceptLATdNaF[m]=AcceptDzoneNaF->GetBinContent(m+1,1,2);
		AcceptLATdAgl[m]=AcceptDzoneAgl->GetBinContent(m+1,1,2);
		AcceptLATpTOF[m]=AcceptPTOF->GetBinContent(m+1);
		AcceptLATpNaF[m]=AcceptPNaF->GetBinContent(m+1);
		AcceptLATpAgl[m]=AcceptPAgl->GetBinContent(m+1);
	}
	for(int m=0;m<18;m++){
		AcceptLATdTOF[m]/=CorrLATd_TOF_spl->GetBinContent(m+1,1);
		AcceptLATdNaF[m]/=CorrLATd_NaF_spl->GetBinContent(m+1,1);
		AcceptLATdAgl[m]/=CorrLATd_Agl_spl->GetBinContent(m+1,1);
		AcceptLATpTOF[m]/=CorrLATp_TOF_spl->GetBinContent(m+1,1);
		AcceptLATpNaF[m]/=CorrLATp_NaF_spl->GetBinContent(m+1,1);
		AcceptLATpAgl[m]/=CorrLATp_Agl_spl->GetBinContent(m+1,1);
	}
	float errore=0;	
	//Geo zone
	for(int l =0;l<11;l++)
		for(int m=0;m<18;m++)
			if(AcceptDzoneTOF->GetBinContent(m+1,l+1,2)>0&&Tempi->GetBinContent(l)>0){
				DFluxgeoTOF->SetBinContent(m+1,l+1,0,DCountsgeoTOF->GetBinContent(m+1,l+1,0)/(AcceptDzoneTOF->GetBinContent(m+1,l+1,2)*Tempi->GetBinContent(l)*deltaencinTOF[m]));
				//err stat
				errore=DCountsgeoTOF->GetBinContent(m+1,l+1,1)/DCountsgeoTOF->GetBinContent(m+1,l+1,0);
				//err lat corr
				//errore=errore+CorrLAT_tot_spl->GetBinContent(l+1,2);
				//err D vs MC
				//for(int S=0;S<3;S++) errore=errore+pow(PreDVSMC_P[S]->GetBinContent(i+1,2)/PreDVSMC_P[S]->GetBinContent(i+1,1),2);
				//errore=errore+pow(LikDVSMC_P_graph->GetBinContent(i+1,2)/LikDVSMC_P_graph->GetBinContent(i+1,1),2);

				//errore=pow(errore,0.5);
				if(errore>1) errore=0.5;
				if(errore>0) //cout<<m<<" "<<l<<" "<<errore<<endl;
				DFluxgeoTOF->SetBinContent(m+1,l+1,1,errore*DFluxgeoTOF->GetBinContent(m+1,l+1,0));
			}

	for(int l =0;l<11;l++)
		for(int m=0;m<18;m++)
			if(AcceptDzoneNaF->GetBinContent(m+1,l+1,2)>0&&Tempi->GetBinContent(l)>0){
				DFluxgeoNaF->SetBinContent(m+1,l+1,0,DCountsgeoNaF->GetBinContent(m+1,l+1,0)/(AcceptDzoneNaF->GetBinContent(m+1,l+1,2)*Tempi->GetBinContent(l)*deltaencinNaF[m]));
				//err stat
				errore=DCountsgeoNaF->GetBinContent(m+1,l+1,1)/DCountsgeoNaF->GetBinContent(m+1,l+1,0);
				if(errore>1) errore=0.5;
				if(errore>0) 
				DFluxgeoNaF->SetBinContent(m+1,l+1,1,errore*DFluxgeoNaF->GetBinContent(m+1,l+1,0));		
			}	
	for(int l =0;l<11;l++)
		for(int m=0;m<18;m++)
			if(AcceptDzoneAgl->GetBinContent(m+1,l+1,2)>0&&Tempi->GetBinContent(l)>0){
				DFluxgeoAgl->SetBinContent(m+1,l+1,0,DCountsgeoAgl->GetBinContent(m+1,l+1,0)/(AcceptDzoneAgl->GetBinContent(m+1,l+1,2)*Tempi->GetBinContent(l)*deltaencinAgl[m]));
				//err stat
				errore=DCountsgeoAgl->GetBinContent(m+1,l+1,1)/DCountsgeoAgl->GetBinContent(m+1,l+1,0);
				if(errore>1) errore=0.5;
				if(errore>0) 
				DFluxgeoAgl->SetBinContent(m+1,l+1,1,errore*DFluxgeoAgl->GetBinContent(m+1,l+1,0));		
			}


	//Primaries
	for(int m=0;m<18;m++)
		if( AcceptLATdTOF[m]>0&&Esposd_TOF[m]>0){	
			DFluxTOF->SetBinContent(m+1,0,DCountsgeoTOF->GetBinContent(m+1,12,0)/(AcceptLATdTOF[m]*Esposd_TOF[m]*deltaencinTOF[m]));
			//err stat
			errore=DCountsgeoTOF->GetBinContent(m+1,12,1)/DCountsgeoTOF->GetBinContent(m+1,12,0);
			if(errore>1) errore=0.5;
			if(errore>0) DFluxTOF->SetBinContent(m+1,1,errore*DFluxTOF->GetBinContent(m+1,0));
		}
	for(int m=0;m<18;m++)
		if( AcceptLATdNaF[m]>0&&Esposd_NaF[m]>0){       
			DFluxNaF->SetBinContent(m+1,0,DCountsgeoNaF->GetBinContent(m+1,12,0)/(AcceptLATdNaF[m]*Esposd_NaF[m]*deltaencinNaF[m]));
			//err stat
			errore=DCountsgeoNaF->GetBinContent(m+1,12,1)/DCountsgeoNaF->GetBinContent(m+1,12,0);
			if(errore>1) errore=0.5;
			if(errore>0) DFluxNaF->SetBinContent(m+1,1,errore*DFluxNaF->GetBinContent(m+1,0));		
		}
	for(int m=0;m<18;m++)
		if( AcceptLATdAgl[m]>0&&Esposd_Agl[m]>0){
			DFluxAgl->SetBinContent(m+1,0,DCountsgeoAgl->GetBinContent(m+1,12,0)/(AcceptLATdAgl[m]*Esposd_Agl[m]*deltaencinAgl[m]));
			//err stat
			errore=DCountsgeoAgl->GetBinContent(m+1,12,1)/DCountsgeoAgl->GetBinContent(m+1,12,0);
			if(errore>1) errore=0.5;
			if(errore>0) DFluxAgl->SetBinContent(m+1,1,errore*DFluxAgl->GetBinContent(m+1,0));	
		}

	//ratio
	for(int m=0;m<18;m++)
		if( AcceptLATdTOF[m]>0&&Esposd_TOF[m]>0&&AcceptLATpTOF[m]>0&&Esposp_TOF[m]>0&&PCountsTOF->GetBinContent(m+1,0)>0){
			PDratioTOF->SetBinContent(m+1,0,(DCountsgeoTOF->GetBinContent(m+1,12,0)/(AcceptLATdTOF[m]*Esposd_TOF[m]*deltaencinTOF[m]))/(PCountsTOF->GetBinContent(m+1,0)/(AcceptLATpTOF[m]*Esposp_TOF[m]*deltaencinTOF[m])));
			//err
			PDratioTOF->SetBinContent(m+1,1,pow(2,0.5)*DFluxTOF->GetBinContent(m+1,1)/DFluxTOF->GetBinContent(m+1,0)*PDratioTOF->GetBinContent(m+1,0));

		}
	for(int m=0;m<18;m++)
		if( AcceptLATdNaF[m]>0&&Esposd_NaF[m]>0&&AcceptLATpNaF[m]>0&&Esposp_NaF[m]>0){
			PDratioNaF->SetBinContent(m+1,0,(DCountsgeoNaF->GetBinContent(m+1,12,0)/(AcceptLATdNaF[m]*Esposd_NaF[m]*deltaencinNaF[m]))/(PCountsNaF->GetBinContent(m+1,0)/(AcceptLATpNaF[m]*Esposp_NaF[m]*deltaencinNaF[m])));
		//err
		PDratioNaF->SetBinContent(m+1,1,pow(2,0.5)*DFluxNaF->GetBinContent(m+1,1)/DFluxNaF->GetBinContent(m+1,0)*PDratioNaF->GetBinContent(m+1,0));
		}
	for(int m=0;m<18;m++)
		if( AcceptLATdAgl[m]>0&&Esposd_Agl[m]>0&&AcceptLATpAgl[m]>0&&Esposp_Agl[m]>0){
			PDratioAgl->SetBinContent(m+1,0,(DCountsgeoAgl->GetBinContent(m+1,12,0)/(AcceptLATdAgl[m]*Esposd_Agl[m]*deltaencinAgl[m]))/(PCountsAgl->GetBinContent(m+1,0)/(AcceptLATpAgl[m]*Esposp_Agl[m]*deltaencinAgl[m])));
		//err
		PDratioAgl->SetBinContent(m+1,1,pow(2,0.5)*DFluxAgl->GetBinContent(m+1,1)/DFluxAgl->GetBinContent(m+1,0)*PDratioAgl->GetBinContent(m+1,0));
		}



	potenza=0;
	c33->Divide(1,3);
	gPad->SetGridx();
	gPad->SetGridy();
	TGraphErrors * esposd_TOF=new TGraphErrors();
	TGraphErrors * esposd_NaF=new TGraphErrors();
	TGraphErrors * esposd_Agl=new TGraphErrors();
	TGraphErrors * esposp_TOF=new TGraphErrors();
	TGraphErrors * esposp_NaF=new TGraphErrors();
	TGraphErrors * esposp_Agl=new TGraphErrors();
	for(int m=0;m<18;m++){
		esposd_TOF->SetPoint(m,Ekincent[m],Esposd_TOF[m]);
		esposp_TOF->SetPoint(m,Ekincent[m],Esposp_TOF[m]);
		esposd_NaF->SetPoint(m,EkincentNaF[m],Esposd_NaF[m]);
		esposp_NaF->SetPoint(m,EkincentNaF[m],Esposp_NaF[m]);
		esposd_Agl->SetPoint(m,EkincentAgl[m],Esposd_Agl[m]);
		esposp_Agl->SetPoint(m,EkincentAgl[m],Esposp_Agl[m]);
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

	c32->cd();
	gPad->SetLogx();
	gPad->SetLogy();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraph* galprop3P=new TGraph();
	TGraph* galprop3P2=new TGraph();
	float x,y=0;
	int j=0;
	{
		string nomefile=percorso+"/CodesforAnalysis/Galprop/Trotta2011/Def/new_D500.txt";
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
		string nomefile=percorso+"/CodesforAnalysis/Galprop/Trotta2011/Def/new_D1250.txt";
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
	galprop3P->Draw("AC");
	galprop3P2->Draw("sameC");

	string nome;
	int p=0;
	for(int j=0;j<11;j++) {
		nome="Deutons Flux: Geo. Zone "+numero[j];
		D_FluxgeoTOF[j]=new TGraphErrors();
		D_FluxgeoTOF[j]->SetName(nome.c_str());
		p=0;
		for(int m=1;m<18;m++){
			D_FluxgeoTOF[j]->SetPoint(p,Ekincent[m],DFluxgeoTOF->GetBinContent(m+1,j+1,0));
			D_FluxgeoTOF[j]->SetPointError(p,0,DFluxgeoTOF->GetBinContent(m+1,j+1,1));
			p++;
		}
		D_FluxgeoTOF[j]->SetMarkerStyle(8);
		D_FluxgeoTOF[j]->SetMarkerSize(1.5);
		D_FluxgeoTOF[j]->SetMarkerColor(j-1);
		D_FluxgeoTOF[j]->SetLineColor(j-1);
		D_FluxgeoTOF[j]->SetLineWidth(2);
	}
	D_FluxgeoTOF[10]->SetTitle("Protons Flux: Geo. Zones");
	D_FluxgeoTOF[10]->GetXaxis()->SetTitle("R [GV]");
	D_FluxgeoTOF[10]->GetYaxis()->SetTitle("Flux [(m^2 sr GeV/nucl.)^-1]");
	D_FluxgeoTOF[10]->GetXaxis()->SetTitleSize(0.045);
	D_FluxgeoTOF[10]->GetYaxis()->SetTitleSize(0.045);
	D_FluxgeoTOF[10]->GetYaxis()->SetRangeUser(1e-2,1e4);
	D_FluxgeoTOF[10]->Draw("Psame");
	for(int j=0;j<11;j++) D_FluxgeoTOF[j]->Draw("Psame");	

	p=0;
	for(int j=0;j<11;j++) {
		nome="Deutons Flux: Geo. Zone "+numero[j];
		D_FluxgeoNaF[j]=new TGraphErrors();
		D_FluxgeoNaF[j]->SetName(nome.c_str());
		p=0;
		for(int m=1;m<18;m++){
			D_FluxgeoNaF[j]->SetPoint(p,EkincentNaF[m],DFluxgeoNaF->GetBinContent(m+1,j+1,0));
			D_FluxgeoNaF[j]->SetPointError(p,0,DFluxgeoNaF->GetBinContent(m+1,j+1,1));
			p++;
		}
		D_FluxgeoNaF[j]->SetMarkerStyle(4);
		D_FluxgeoNaF[j]->SetMarkerSize(1.5);
		D_FluxgeoNaF[j]->SetMarkerColor(j-1);
		D_FluxgeoNaF[j]->SetLineColor(j-1);
		D_FluxgeoNaF[j]->SetLineWidth(2);
	}
	D_FluxgeoNaF[10]->SetTitle("Protons Flux: Geo. Zones");
	D_FluxgeoNaF[10]->GetXaxis()->SetTitle("R [GV]");
	D_FluxgeoNaF[10]->GetYaxis()->SetTitle("Flux [(m^2 sr GeV/nucl.)^-1]");
	D_FluxgeoNaF[10]->GetXaxis()->SetTitleSize(0.045);
	D_FluxgeoNaF[10]->GetYaxis()->SetTitleSize(0.045);
	D_FluxgeoNaF[10]->GetYaxis()->SetRangeUser(1e-2,1e4);
	D_FluxgeoNaF[10]->Draw("Psame");
	for(int j=0;j<11;j++) D_FluxgeoNaF[j]->Draw("Psame");

	p=0;
	for(int j=0;j<11;j++) {
		nome="Deutons Flux: Geo. Zone "+numero[j];
		D_FluxgeoAgl[j]=new TGraphErrors();
		D_FluxgeoAgl[j]->SetName(nome.c_str());
		p=0;
		for(int m=1;m<18;m++){
			D_FluxgeoAgl[j]->SetPoint(p,EkincentAgl[m],DFluxgeoAgl->GetBinContent(m+1,j+1,0));
			D_FluxgeoAgl[j]->SetPointError(p,0,DFluxgeoAgl->GetBinContent(m+1,j+1,1));
			p++;
		}
		D_FluxgeoAgl[j]->SetMarkerStyle(3);
		D_FluxgeoAgl[j]->SetMarkerSize(1.5);
		D_FluxgeoAgl[j]->SetMarkerColor(j-1);
		D_FluxgeoAgl[j]->SetLineColor(j-1);
		D_FluxgeoAgl[j]->SetLineWidth(2);
	}
	D_FluxgeoAgl[10]->SetTitle("Protons Flux: Geo. Zones");
	D_FluxgeoAgl[10]->GetXaxis()->SetTitle("R [GV]");
	D_FluxgeoAgl[10]->GetYaxis()->SetTitle("Flux [(m^2 sr GeV/nucl.)^-1]");
	D_FluxgeoAgl[10]->GetXaxis()->SetTitleSize(0.045);
	D_FluxgeoAgl[10]->GetYaxis()->SetTitleSize(0.045);
	D_FluxgeoAgl[10]->GetYaxis()->SetRangeUser(1e-2,1e4);
	D_FluxgeoAgl[10]->Draw("Psame");
	for(int j=0;j<11;j++) D_FluxgeoAgl[j]->Draw("Psame");

	c34->cd();
	gPad->SetLogx();
	gPad->SetLogy();
	gPad->SetGridx();
	gPad->SetGridy();
	galprop3P->GetXaxis()->SetRangeUser(0.1,10);
	galprop3P->GetYaxis()->SetRangeUser(1e-3,1e3);
	galprop3P->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
        galprop3P->GetYaxis()->SetTitle("D Primary Flux [(m^2 sec sr GeV/nucl.)^-1]");
	galprop3P->Draw("AC");
	galprop3P2->Draw("sameC");
	nome="Deutons Primary Flux";
	D_FluxTOF=new TGraphErrors();
	D_FluxTOF->SetName(nome.c_str());
	p=0;
	for(int m=1;m<18;m++){
		D_FluxTOF->SetPoint(p,Ekincent[m],DFluxTOF->GetBinContent(m+1,0));
		cout<<DFluxTOF->GetBinContent(m+1,1)<<" "<<DFluxTOF->GetBinContent(m+1,0)<<endl;
		D_FluxTOF->SetPointError(p,0,DFluxTOF->GetBinContent(m+1,1));
		p++;
	}
	D_FluxTOF->SetMarkerStyle(8);
	D_FluxTOF->SetMarkerSize(1.5);
	D_FluxTOF->SetMarkerColor(4);
	D_FluxTOF->SetLineColor(4);
	D_FluxTOF->SetLineWidth(2);
	D_FluxTOF->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
	D_FluxTOF->GetYaxis()->SetTitle("Flux [(m^2 sr GeV/nucl.)^-1]");
	D_FluxTOF->GetXaxis()->SetTitleSize(0.045);
	D_FluxTOF->GetYaxis()->SetTitleSize(0.045);
	D_FluxTOF->GetYaxis()->SetRangeUser(1e-2,1e4);
	D_FluxTOF->Draw("Psame");

	D_FluxNaF=new TGraphErrors();
	D_FluxNaF->SetName(nome.c_str());
	p=0;
	for(int m=1;m<18;m++){
		D_FluxNaF->SetPoint(p,EkincentNaF[m],DFluxNaF->GetBinContent(m+1,0));
		D_FluxNaF->SetPointError(p,0,DFluxNaF->GetBinContent(m+1,1));
		p++;
	}
	D_FluxNaF->SetMarkerStyle(4);
	D_FluxNaF->SetMarkerSize(1.5);
	D_FluxNaF->SetMarkerColor(4);
	D_FluxNaF->SetLineColor(4);
	D_FluxNaF->SetLineWidth(2);
	D_FluxNaF->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
	D_FluxNaF->GetYaxis()->SetTitle("Flux [(m^2 sr GeV/nucl.)^-1]");
	D_FluxNaF->GetXaxis()->SetTitleSize(0.045);
	D_FluxNaF->GetYaxis()->SetTitleSize(0.045);
	D_FluxNaF->GetYaxis()->SetRangeUser(1e-2,1e4);
	D_FluxNaF->Draw("Psame");

	D_FluxAgl=new TGraphErrors();
	D_FluxAgl->SetName(nome.c_str());
	p=0;
	for(int m=1;m<18;m++){
		D_FluxAgl->SetPoint(p,EkincentAgl[m],DFluxAgl->GetBinContent(m+1,0));
		D_FluxAgl->SetPointError(p,0,DFluxAgl->GetBinContent(m+1,1));
		p++;
	}
	D_FluxAgl->SetMarkerStyle(3);
	D_FluxAgl->SetMarkerSize(1.5);
	D_FluxAgl->SetMarkerColor(4);
	D_FluxAgl->SetLineColor(4);
	D_FluxAgl->SetLineWidth(2);
	D_FluxAgl->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
	D_FluxAgl->GetYaxis()->SetTitle("Flux [(m^2 sr GeV/nucl.)^-1]");
	D_FluxAgl->GetXaxis()->SetTitleSize(0.045);
	D_FluxAgl->GetYaxis()->SetTitleSize(0.045);
	D_FluxAgl->GetYaxis()->SetRangeUser(1e-2,1e4);
	D_FluxAgl->Draw("Psame");

	c35->cd();
	gPad->SetLogx();
	gPad->SetLogy();
	gPad->SetGridx();
	gPad->SetGridy();	
	TGraph* galpropratio1=new TGraph();
	TGraph* galpropratio2=new TGraph();
	x,y=0;
	j=0;
	{
		string nomefile=percorso+"/CodesforAnalysis/Galprop/Trotta2011/PDratio/500.dat";
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
		string nomefile=percorso+"/CodesforAnalysis/Galprop/Trotta2011/PDratio/1000.dat";
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
	galpropratio1->Draw("AC");
	galpropratio2->Draw("sameC");
	float encinPamela[18]={0.120,0.132,0.144,0.158,0.173,0.190,0.208,0.228,0.250,0.274,0.300,0.329,0.361,0.395,0.433,0.475,0.520,0.570};
        float RatioPamela[17]={3.12,3.10,3.01,2.97,2.86,2.81,2.74,2.58,2.37,2.28,2.25,2.22,2.10,2.11,2.07,2.00,2.03};
        float ErrPamela[17]={0.19+0.20,0.15+0.18,0.13+0.16,0.12+0.14,0.11+0.13,0.10+0.12,0.09+0.11,0.08+0.10,0.07+0.09,0.06+0.08,0.06+0.08,0.06+0.08,0.06+0.07,0.05+0.07,0.05+0.07,0.05+0.07,0.06+0.07};
        TGraphAsymmErrors *pamelaratio=new TGraphAsymmErrors();

        for(int m=0;m<17;m++)
        {       pamelaratio->SetPoint(m,(encinPamela[m]+encinPamela[m+1])/2,RatioPamela[m]/100);
                pamelaratio->SetPointError(m,0,0,ErrPamela[m]/100,ErrPamela[m]/100);
        }
        pamelaratio->SetMarkerColor(1);
        pamelaratio->SetMarkerStyle(23);
	pamelaratio->SetMarkerSize(1.5);
        pamelaratio->SetLineColor(1);
	pamelaratio->Draw("Psame");
	PD_ratioTOF=new TGraphErrors();
	p=0;
	for(int m=1;m<18;m++){
		PD_ratioTOF->SetPoint(p,Ekincent[m],PDratioTOF->GetBinContent(m+1,0));
		PD_ratioTOF->SetPointError(p,0,PDratioTOF->GetBinContent(m+1,1));
		p++;
	}
	PD_ratioTOF->SetMarkerStyle(8);
	PD_ratioTOF->SetMarkerSize(1.5);
	PD_ratioTOF->SetMarkerColor(4);
	PD_ratioTOF->SetLineColor(4);
	PD_ratioTOF->SetLineWidth(2);
	PD_ratioTOF->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
	PD_ratioTOF->GetYaxis()->SetTitle("Flux ratio");
	PD_ratioTOF->GetXaxis()->SetTitleSize(0.045);
	PD_ratioTOF->GetYaxis()->SetTitleSize(0.045);
	PD_ratioTOF->GetYaxis()->SetRangeUser(1e-2,1e4);
	PD_ratioTOF->Draw("Psame");	

	PD_ratioNaF=new TGraphErrors();
	p=0;
	for(int m=1;m<18;m++){
		PD_ratioNaF->SetPoint(p,EkincentNaF[m],PDratioNaF->GetBinContent(m+1,0));
		PD_ratioNaF->SetPointError(p,0,PDratioNaF->GetBinContent(m+1,1));
		p++;
	}
	PD_ratioNaF->SetMarkerStyle(4);
	PD_ratioNaF->SetMarkerSize(1.5);
	PD_ratioNaF->SetMarkerColor(4);
	PD_ratioNaF->SetLineColor(4);
	PD_ratioNaF->SetLineWidth(2);
	PD_ratioNaF->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
	PD_ratioNaF->GetYaxis()->SetTitle("Flux ratio");
	PD_ratioNaF->GetXaxis()->SetTitleSize(0.045);
	PD_ratioNaF->GetYaxis()->SetTitleSize(0.045);
	PD_ratioNaF->GetYaxis()->SetRangeUser(1e-2,1e4);
	PD_ratioNaF->Draw("Psame");

	PD_ratioAgl=new TGraphErrors();
	p=0;
	for(int m=1;m<18;m++){
		PD_ratioAgl->SetPoint(p,EkincentAgl[m],PDratioAgl->GetBinContent(m+1,0));
		PD_ratioAgl->SetPointError(p,0,PDratioAgl->GetBinContent(m+1,1));
		p++;
	}
	PD_ratioAgl->SetMarkerStyle(3);
	PD_ratioAgl->SetMarkerSize(1.5);
	PD_ratioAgl->SetMarkerColor(4);
	PD_ratioAgl->SetLineColor(4);
	PD_ratioAgl->SetLineWidth(2);
	PD_ratioAgl->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
	PD_ratioAgl->GetYaxis()->SetTitle("Flux ratio");
	PD_ratioAgl->GetXaxis()->SetTitleSize(0.045);
	PD_ratioAgl->GetYaxis()->SetTitleSize(0.045);
	PD_ratioAgl->GetYaxis()->SetRangeUser(1e-2,1e4);
	PD_ratioAgl->Draw("Psame");
	P_Fluxgeo[10]->Draw("Psame");
}

