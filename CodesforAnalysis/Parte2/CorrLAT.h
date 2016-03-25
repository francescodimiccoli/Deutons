

void CorrLAT(TFile * file1){
	TH2F * esposizionepgeo 	   = (TH2F*)file1->Get(	"esposizionepgeo"	);
        TH2F * esposizionepgeoNaF  = (TH2F*)file1->Get(	"esposizionepgeoNaF"	);
        TH2F * esposizionepgeoAgl  = (TH2F*)file1->Get(	"esposizionepgeoAgl"	);
        TH2F * esposizionedgeo 	   = (TH2F*)file1->Get(	"esposizionedgeo"	);
        TH2F * esposizionedgeoNaF  = (TH2F*)file1->Get(	"esposizionedgeoNaF"	);
        TH2F * esposizionedgeoAgl  = (TH2F*)file1->Get(	"esposizionedgeoAgl"	);

	LATcorr * EffpreSelDATA     = new LATcorr(file1,"EffpreSelDATA"  ,"Results");	
	LATcorr * LikelihoodLATcorr = new LATcorr(file1,"EffLikDATA" 	 ,"Results");
        LATcorr * DistanceLATcorr   = new LATcorr(file1,"EffDistDATA"	 ,"Results");

	cout<<"******* TOTAL LAT. CORRECTION *************"<<endl;

	
	TH1F * PreLATCorr = ( (TH1F *)( (TH2F *)EffpreSelDATA ->  LATcorrR_fit ) -> ProjectionX("",1,1));

	PreLATCorr -> Multiply( ( (TH1F *)( (TH2F *)EffpreSelDATA ->  LATcorrR_fit ) -> ProjectionX("",2,2)));
	PreLATCorr -> Multiply( ( (TH1F *)( (TH2F *)EffpreSelDATA ->  LATcorrR_fit ) -> ProjectionX("",3,3)));		

	TH1F *  TOTLATCorr = (TH1F *) PreLATCorr -> Clone();
	
	TOTLATCorr  -> Multiply ( ( (TH1F *)LikelihoodLATcorr ->  LATcorrTOF_fit) );
	TOTLATCorr  -> Multiply ( ( (TH1F *)DistanceLATcorr   ->  LATcorrTOF_fit) );









/*	float CorrLAT_tot[11]={1,1,1,1,1,1,1,1,1,1,1};
	float CorrLAT_tot_err[11]={0};

	float CorrLAT_Pre[11]={1,1,1,1,1,1,1,1,1,1,1};
	float CorrLAT_Pre_err[11]={0};

	for(int m=0;m<11;m++) {
		for(int S=0;S<3;S++) CorrLAT_tot[m]*=CorrLAT_pre[S]->Eval(geomagC[m]);
		for(int S=0;S<3;S++) CorrLAT_Pre[m]*=CorrLAT_pre[S]->Eval(geomagC[m]);
		CorrLAT_tot[m]*=CorrLAT_Lik->Eval(geomagC[m]);
		CorrLAT_tot[m]*=CorrLAT_Dist->Eval(geomagC[m]);

		for(int S=0;S<3;S++) CorrLAT_tot_err[m]+=pow(CorrLATpre_spl[S]->GetBinContent(m+1,2)/CorrLATpre_spl[S]->GetBinContent(m+1,1),2);
		for(int S=0;S<3;S++) CorrLAT_Pre_err[m]+=pow(CorrLATpre_spl[S]->GetBinContent(m+1,2)/CorrLATpre_spl[S]->GetBinContent(m+1,1),2);
		CorrLAT_tot_err[m]+=pow(CorrLAT_Lik_spl->GetBinContent(m+1,2)/CorrLAT_Lik_spl->GetBinContent(m+1,1),2);
		CorrLAT_tot_err[m]+=pow(CorrLAT_Dist_spl->GetBinContent(m+1,2)/CorrLAT_Dist_spl->GetBinContent(m+1,1),2);
		CorrLAT_tot_err[m]=pow(CorrLAT_tot_err[m],0.5)*CorrLAT_tot[m];
		CorrLAT_Pre_err[m]=pow(CorrLAT_Pre_err[m],0.5)*CorrLAT_Pre[m];
	}		


	float CorrLATNaF_tot[11]={1,1,1,1,1,1,1,1,1,1,1};
	float CorrLATNaF_tot_err[11]={0};

	for(int m=0;m<11;m++) {
		for(int S=0;S<3;S++) CorrLATNaF_tot[m]*=CorrLAT_pre[S]->Eval(geomagC[m]);
		CorrLATNaF_tot[m]*=CorrLAT_LikNaF->Eval(geomagC[m]);
		CorrLATNaF_tot[m]*=CorrLAT_DistNaF->Eval(geomagC[m]);

		for(int S=0;S<3;S++) CorrLATNaF_tot_err[m]+=pow(CorrLATpre_spl[S]->GetBinContent(m+1,2)/CorrLATpre_spl[S]->GetBinContent(m+1,1),2);
		CorrLATNaF_tot_err[m]+=pow(CorrLAT_LikNaF_spl->GetBinContent(m+1,2)/CorrLAT_LikNaF_spl->GetBinContent(m+1,1),2);
		CorrLATNaF_tot_err[m]+=pow(CorrLAT_DistNaF_spl->GetBinContent(m+1,2)/CorrLAT_DistNaF_spl->GetBinContent(m+1,1),2);
		CorrLATNaF_tot_err[m]=pow(CorrLATNaF_tot_err[m],0.5)*CorrLATNaF_tot[m];
	}

	float CorrLATAgl_tot[11]={1,1,1,1,1,1,1,1,1,1,1};
	float CorrLATAgl_tot_err[11]={0};

	for(int m=0;m<11;m++) {
		for(int S=0;S<3;S++) CorrLATAgl_tot[m]*=CorrLAT_pre[S]->Eval(geomagC[m]);
		CorrLATAgl_tot[m]*=CorrLAT_LikAgl->Eval(geomagC[m]);
		CorrLATAgl_tot[m]*=CorrLAT_DistAgl->Eval(geomagC[m]);

		for(int S=0;S<3;S++) CorrLATAgl_tot_err[m]+=pow(CorrLATpre_spl[S]->GetBinContent(m+1,2)/CorrLATpre_spl[S]->GetBinContent(m+1,1),2);
		CorrLATAgl_tot_err[m]+=pow(CorrLAT_LikAgl_spl->GetBinContent(m+1,2)/CorrLAT_LikAgl_spl->GetBinContent(m+1,1),2);
		CorrLATAgl_tot_err[m]+=pow(CorrLAT_DistAgl_spl->GetBinContent(m+1,2)/CorrLAT_DistAgl_spl->GetBinContent(m+1,1),2);
		CorrLATAgl_tot_err[m]=pow(CorrLATAgl_tot_err[m],0.5)*CorrLATAgl_tot[m];
	}


	//METODO 2
	float Espos_R[43]={0};	
	float CorrLATtotM2[43]={0};
	float CorrLATpreM2[43]={0};
	float CorrLATtotM2_err[43]={0};
	float CorrLATpreM2_err[43]={0};
	for(int i=0;i<43;i++){
			for(int m=1;m<11;m++) Espos_R[i]+=esposizionegeo->GetBinContent(i+1,m);
			for(int m=0;m<11;m++) CorrLATtotM2[i]+=CorrLAT_tot[m]*esposizionegeo->GetBinContent(i+1,m)/Espos_R[i];
			for(int m=0;m<11;m++) CorrLATtotM2_err[i]+=pow(CorrLAT_tot_err[m]*esposizionegeo->GetBinContent(i+1,m)/Espos_R[i],2); 
			for(int m=0;m<11;m++) CorrLATpreM2[i]+=CorrLAT_Pre[m]*esposizionegeo->GetBinContent(i+1,m)/Espos_R[i];
			for(int m=0;m<11;m++) CorrLATpreM2_err[i]+=pow(CorrLAT_Pre_err[m]*esposizionegeo->GetBinContent(i+1,m)/Espos_R[i],2);
		}
	
	for(int i=0;i<43;i++) CorrLATtotM2_err[i]=pow(CorrLATtotM2_err[i],0.5);
	for(int i=0;i<43;i++) CorrLATpreM2_err[i]=pow(CorrLATpreM2_err[i],0.5);
	for(int i=0;i<43;i++) CorrLAT_totM2->SetBinContent(i+1,1,CorrLATtotM2[i]);
	for(int i=0;i<43;i++) CorrLAT_preM2->SetBinContent(i+1,1,CorrLATpreM2[i]);
	for(int i=0;i<43;i++) CorrLAT_totM2->SetBinContent(i+1,2,CorrLATtotM2_err[i]);
        for(int i=0;i<43;i++) CorrLAT_preM2->SetBinContent(i+1,2,CorrLATpreM2_err[i]);
	
	//bin beta
	//METODO 2	
	float Esposd_TOF[18]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        float Esposd_NaF[18]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        float Esposd_Agl[18]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        float Esposp_TOF[18]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        float Esposp_NaF[18]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        float Esposp_Agl[18]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	float CorrLATd_TOF[18]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        float CorrLATd_NaF[18]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        float CorrLATd_Agl[18]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        float CorrLATp_TOF[18]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        float CorrLATp_NaF[18]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        float CorrLATp_Agl[18]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};


        for(int m=0;m<18;m++)
                for(int i=1;i<11;i++) {
                        Esposd_TOF[m]+=esposizionedgeo->GetBinContent(m+1,i);
                        Esposp_TOF[m]+=esposizionepgeo->GetBinContent(m+1,i);
                        Esposd_NaF[m]+=esposizionedgeoNaF->GetBinContent(m+1,i);
                        Esposp_NaF[m]+=esposizionepgeoNaF->GetBinContent(m+1,i);
                        Esposd_Agl[m]+=esposizionedgeoAgl->GetBinContent(m+1,i);
                        Esposp_Agl[m]+=esposizionepgeoAgl->GetBinContent(m+1,i);
                        }

	for(int m=0;m<18;m++)	
		for(int i=1;i<11;i++) {
			CorrLATd_TOF[m]+=CorrLAT_tot[i]*esposizionedgeo->GetBinContent(m+1,i)/Esposd_TOF[m];
			CorrLATp_TOF[m]+=CorrLAT_tot[i]*esposizionepgeo->GetBinContent(m+1,i)/Esposp_TOF[m];
			CorrLATd_NaF[m]+=CorrLATNaF_tot[i]*esposizionedgeoNaF->GetBinContent(m+1,i)/Esposd_NaF[m];
                        CorrLATp_NaF[m]+=CorrLATNaF_tot[i]*esposizionepgeoNaF->GetBinContent(m+1,i)/Esposp_NaF[m];
			CorrLATd_Agl[m]+=CorrLATAgl_tot[i]*esposizionedgeoAgl->GetBinContent(m+1,i)/Esposd_Agl[m];
                        CorrLATp_Agl[m]+=CorrLATAgl_tot[i]*esposizionepgeoAgl->GetBinContent(m+1,i)/Esposp_Agl[m];				
		}
	
	for(int m=0;m<18;m++){
		CorrLATd_TOF_spl->SetBinContent(m+1,1,CorrLATd_TOF[m]);
		CorrLATp_TOF_spl->SetBinContent(m+1,1,CorrLATp_TOF[m]);
		CorrLATd_NaF_spl->SetBinContent(m+1,1,CorrLATd_NaF[m]);
                CorrLATp_NaF_spl->SetBinContent(m+1,1,CorrLATp_NaF[m]);
		CorrLATd_Agl_spl->SetBinContent(m+1,1,CorrLATd_Agl[m]);
                CorrLATp_Agl_spl->SetBinContent(m+1,1,CorrLATp_Agl[m]);
	}
*/

	TCanvas * c26 = new TCanvas("Latitude pile-up correction (R bins)");
	TCanvas * c26_bis = new TCanvas("Latitude pile-up correction (Beta bins)");


/*	
	c26->Divide(1,2);
	c26->cd(1);
	gPad->SetGridy();
        gPad->SetGridx();
	TGraphErrors *CorrLAT_tot_Spl=new TGraphErrors();
	for(int m=0;m<11;m++) {
		CorrLAT_tot_spl->SetBinContent(m+1,0,CorrLAT_tot[m]);
		CorrLAT_tot_spl->SetBinContent(m+1,1,CorrLAT_tot_err[m]);
		CorrLAT_tot_Spl->SetPoint(m,geomagC[m],CorrLAT_tot[m]);
		CorrLAT_tot_Spl->SetPointError(m,0,CorrLAT_tot_err[m]);
	}
	CorrLAT_tot_Spl->SetLineColor(1);
	CorrLAT_tot_Spl->SetMarkerColor(1);
	CorrLAT_tot_Spl->SetLineWidth(2);
	CorrLAT_tot_Spl->SetMarkerStyle(8);
	CorrLAT_tot_Spl->SetFillStyle(3002);
	CorrLAT_tot_Spl->GetXaxis()->SetTitle("Latitude");
        CorrLAT_tot_Spl->GetYaxis()->SetTitle("Eff. Corr. Factor");
	CorrLAT_tot_Spl->Draw("APC");
	
	c26->cd(2);
        gPad->SetGridy();
        gPad->SetGridx();
	gPad->SetLogx();
	TGraphErrors *CorrLAT_totM1_Spl=new TGraphErrors();
	for(int i=0;i<43;i++){
		CorrLAT_totM1_Spl->SetPoint(i,R_cent[i],CorrLAT_totM1->GetBinContent(i+1,1));
	}
	CorrLAT_totM1_Spl->SetLineColor(2);
        CorrLAT_totM1_Spl->SetMarkerColor(2);
        CorrLAT_totM1_Spl->SetLineWidth(2);
        CorrLAT_totM1_Spl->SetMarkerStyle(8);
        CorrLAT_totM1_Spl->SetFillStyle(3002);
	CorrLAT_totM1_Spl->GetXaxis()->SetTitle("R [GV]");
        CorrLAT_totM1_Spl->GetYaxis()->SetTitle("Eff. Corr. Factor");
	CorrLAT_totM1_Spl->Draw("APC");
	TGraphErrors *CorrLAT_totM2_Spl=new TGraphErrors();
        for(int i=0;i<43;i++){
                CorrLAT_totM2_Spl->SetPoint(i,R_cent[i],CorrLAT_totM2->GetBinContent(i+1,1));
        }
        CorrLAT_totM2_Spl->SetLineColor(4);
        CorrLAT_totM2_Spl->SetMarkerColor(4);
        CorrLAT_totM2_Spl->SetLineWidth(2);
        CorrLAT_totM2_Spl->SetMarkerStyle(8);
        CorrLAT_totM2_Spl->SetFillStyle(3002);
	CorrLAT_totM2_Spl->Draw("PCsame");

	c26_bis->Divide(1,3);
        c26_bis->cd(1);
        gPad->SetGridy();
        gPad->SetGridx();
	TGraphErrors * CorrLATp_TOF_Spl=new TGraphErrors();
	for(int m=0;m<18;m++){
			CorrLATp_TOF_Spl->SetPoint(m,Ekincent[m],CorrLATp_TOF_spl->GetBinContent(m+1,1));		
		}
	CorrLATp_TOF_Spl->SetLineColor(2);
        CorrLATp_TOF_Spl->SetMarkerColor(2);
        CorrLATp_TOF_Spl->SetLineWidth(2);
        CorrLATp_TOF_Spl->SetMarkerStyle(8);
        CorrLATp_TOF_Spl->SetFillStyle(3002);
        CorrLATp_TOF_Spl->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
        CorrLATp_TOF_Spl->GetYaxis()->SetTitle("Eff. Corr. Factor");
        CorrLATp_TOF_Spl->Draw("APC");
	TGraphErrors * CorrLATd_TOF_Spl=new TGraphErrors();
        for(int m=0;m<18;m++){
                        CorrLATd_TOF_Spl->SetPoint(m,Ekincent[m],CorrLATd_TOF_spl->GetBinContent(m+1,1));
                }
	CorrLATd_TOF_Spl->SetLineColor(4);
        CorrLATd_TOF_Spl->SetMarkerColor(4);
        CorrLATd_TOF_Spl->SetLineWidth(2);
        CorrLATd_TOF_Spl->SetMarkerStyle(8);
        CorrLATd_TOF_Spl->SetFillStyle(3002);		
	CorrLATd_TOF_Spl->Draw("PCsame");
	 
	c26_bis->cd(2);
        gPad->SetGridy();
        gPad->SetGridx();
        TGraphErrors * CorrLATp_NaF_Spl=new TGraphErrors();
        for(int m=0;m<18;m++){
                        CorrLATp_NaF_Spl->SetPoint(m,EkincentNaF[m],CorrLATp_NaF_spl->GetBinContent(m+1,1));
                }
        CorrLATp_NaF_Spl->SetLineColor(2);
        CorrLATp_NaF_Spl->SetMarkerColor(2);
        CorrLATp_NaF_Spl->SetLineWidth(2);
        CorrLATp_NaF_Spl->SetMarkerStyle(8);
        CorrLATp_NaF_Spl->SetFillStyle(3002);
        CorrLATp_NaF_Spl->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
        CorrLATp_NaF_Spl->GetYaxis()->SetTitle("Eff. Corr. Factor");
        CorrLATp_NaF_Spl->Draw("APC");
        TGraphErrors * CorrLATd_NaF_Spl=new TGraphErrors();
        for(int m=0;m<18;m++){
                        CorrLATd_NaF_Spl->SetPoint(m,EkincentNaF[m],CorrLATd_NaF_spl->GetBinContent(m+1,1));
                }
        CorrLATd_NaF_Spl->SetLineColor(4);
        CorrLATd_NaF_Spl->SetMarkerColor(4);
        CorrLATd_NaF_Spl->SetLineWidth(2);
        CorrLATd_NaF_Spl->SetMarkerStyle(8);
        CorrLATd_NaF_Spl->SetFillStyle(3002);
        CorrLATd_NaF_Spl->Draw("PCsame");

	c26_bis->cd(3);
        gPad->SetGridy();
        gPad->SetGridx();
        TGraphErrors * CorrLATp_Agl_Spl=new TGraphErrors();
        for(int m=0;m<18;m++){
                        CorrLATp_Agl_Spl->SetPoint(m,EkincentAgl[m],CorrLATp_Agl_spl->GetBinContent(m+1,1));
                }
        CorrLATp_Agl_Spl->SetLineColor(2);
        CorrLATp_Agl_Spl->SetMarkerColor(2);
        CorrLATp_Agl_Spl->SetLineWidth(2);
        CorrLATp_Agl_Spl->SetMarkerStyle(8);
        CorrLATp_Agl_Spl->SetFillStyle(3002);
        CorrLATp_Agl_Spl->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
        CorrLATp_Agl_Spl->GetYaxis()->SetTitle("Eff. Corr. Factor");
        CorrLATp_Agl_Spl->Draw("APC");
        TGraphErrors * CorrLATd_Agl_Spl=new TGraphErrors();
        for(int m=0;m<18;m++){
                        CorrLATd_Agl_Spl->SetPoint(m,EkincentAgl[m],CorrLATd_Agl_spl->GetBinContent(m+1,1));
                }
        CorrLATd_Agl_Spl->SetLineColor(4);
        CorrLATd_Agl_Spl->SetMarkerColor(4);
        CorrLATd_Agl_Spl->SetLineWidth(2);
        CorrLATd_Agl_Spl->SetMarkerStyle(8);
        CorrLATd_Agl_Spl->SetFillStyle(3002);
        CorrLATd_Agl_Spl->Draw("PCsame");
*/

	return;


}
