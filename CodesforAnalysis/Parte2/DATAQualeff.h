using namespace std;

TCanvas *c15=new TCanvas("Latitude Likelihood Efficiency");
TCanvas *c16=new TCanvas("Latitude Distance Efficiency");

TH2F * EffDistDATAP = new TH2F("EffDistDATAP","EffDistDATAP",43,0,43,11,0,11);
TH2F * EffLikDATAP = new TH2F("EffLikDATAP","EffLikDATAP",43,0,43,11,0,11);
TH2F * EffDistDATANaFP = new TH2F("EffDistDATANaFP","EffDistDATANaFP",43,0,43,11,0,11);
TH2F * EffLikDATANaFP = new TH2F("EffLikDATANaFP","EffLikDATANaFP",43,0,43,11,0,11);
TH2F * EffDistDATAAglP = new TH2F("EffDistDATAAglP","EffDistDATAAglP",43,0,43,11,0,11);
TH2F * EffLikDATAAglP = new TH2F("EffLikDATAAglP","EffLikDATAAglP",43,0,43,11,0,11);



TH2F * EffQualDATAP = new TH2F("EffQualDATAP","EffQualDATAP",43,0,43,11,0,11);
TH2F * EffQualDATANaFP = new TH2F("EffQualDATANaFP","EffQualDATANaFP",43,0,43,11,0,11);
TH2F * EffQualDATAAglP = new TH2F("EffQualDATAAglP","EffQualDATAAglP",43,0,43,11,0,11);

TF1 * CorrLAT_Lik;
TF1 * CorrLAT_Dist;

TH2F *CorrLAT_Lik_spl;
TGraphErrors *CorrLAT_Lik_Spl;

TH2F *CorrLAT_Dist_spl;
TGraphErrors *CorrLAT_Dist_Spl;


TF1 * CorrLAT_LikNaF;
TF1 * CorrLAT_DistNaF;

TH2F *CorrLAT_LikNaF_spl;
TGraphErrors *CorrLAT_LikNaF_Spl;

TH2F *CorrLAT_DistNaF_spl;
TGraphErrors *CorrLAT_DistNaF_Spl;


TF1 * CorrLAT_LikAgl;
TF1 * CorrLAT_DistAgl;

TH2F *CorrLAT_LikAgl_spl;
TGraphErrors *CorrLAT_LikAgl_Spl;

TH2F *CorrLAT_DistAgl_spl;
TGraphErrors *CorrLAT_DistAgl_Spl;


void DATAQualeff_Fill(TNtuple *ntupla, int l,int zona){

	int k = ntupla->GetEvent(l);
	if(Beta<=0||R<=0) return;
	if(!(EdepL1>0&&EdepL1<EdepL1beta->Eval(Beta)+0.1&&EdepL1>EdepL1beta->Eval(Beta)-0.1)) return;
	if(!(((int)Cutmask>>11)==0||((int)Cutmask>>11)==512)){
		for(int K=0;K<43;K++) if(R<bin[K+1]&&R>bin[K]&&R>Rcut[zona]) {EffQualDATAP->Fill(K,zona);}
		if(Dist5D_P<6) {
			for(int K=0;K<43;K++) if(R<bin[K+1]&&R>bin[K]&&R>Rcut[zona]) {EffDistDATAP->Fill(K,zona);}
		}
		if(Dist5D_P<6&Likcut){
			for(int K=0;K<43;K++) if(R<bin[K+1]&&R>bin[K]&&R>Rcut[zona]) {EffLikDATAP->Fill(K,zona);}
		}
	}
	if(((int)Cutmask>>11)==512){
		for(int K=0;K<43;K++) if(R<bin[K+1]&&R>bin[K]&&R>Rcut[zona]) {EffQualDATANaFP->Fill(K,zona);}
                if(Dist5D_P<6) {
                        for(int K=0;K<43;K++) if(R<bin[K+1]&&R>bin[K]&&R>Rcut[zona]) {EffDistDATANaFP->Fill(K,zona);}
                }
                if(Dist5D_P<6&Likcut){
                        for(int K=0;K<43;K++) if(R<bin[K+1]&&R>bin[K]&&R>Rcut[zona]) {EffLikDATANaFP->Fill(K,zona);}
                }
	}
	if(((int)Cutmask>>11)==0){
                for(int K=0;K<43;K++) if(R<bin[K+1]&&R>bin[K]&&R>Rcut[zona]) {EffQualDATAAglP->Fill(K,zona);}
                if(Dist5D_P<6) {
                        for(int K=0;K<43;K++) if(R<bin[K+1]&&R>bin[K]&&R>Rcut[zona]) {EffDistDATAAglP->Fill(K,zona);}
                }
                if(Dist5D_P<6&Likcut){
                        for(int K=0;K<43;K++) if(R<bin[K+1]&&R>bin[K]&&R>Rcut[zona]) {EffLikDATAAglP->Fill(K,zona);}
                }
        }

	return;
}


void DATAQualeff_Copy(TFile * file1){
	EffQualDATAP = (TH2F *)file1->Get("EffQualDATAP");
	EffDistDATAP = (TH2F *)file1->Get("EffDistDATAP");
	EffLikDATAP = (TH2F *)file1->Get("EffLikDATAP");
	EffQualDATANaFP = (TH2F *)file1->Get("EffQualDATANaFP");
        EffDistDATANaFP = (TH2F *)file1->Get("EffDistDATANaFP");
        EffLikDATANaFP = (TH2F *)file1->Get("EffLikDATANaFP");
	EffQualDATAAglP = (TH2F *)file1->Get("EffQualDATAAglP");
        EffDistDATAAglP = (TH2F *)file1->Get("EffDistDATAAglP");
        EffLikDATAAglP = (TH2F *)file1->Get("EffLikDATAAglP");
	cout<<EffQualDATAP<<" "<<EffDistDATAP<<" "<<EffLikDATAP<<" "<<EffQualDATANaFP<<" "<<EffDistDATANaFP<<" "<<EffLikDATANaFP<<" "<< EffQualDATAAglP<<" "<<EffDistDATAAglP<<" "<<EffLikDATAAglP<<endl;
	return;
}

void DATAQualeff_Write(){
        EffQualDATAP->Write();
        EffDistDATAP->Write();
        EffLikDATAP->Write();
        EffQualDATANaFP->Write();
        EffDistDATANaFP->Write();
        EffLikDATANaFP->Write();
        EffQualDATAAglP->Write();
        EffDistDATAAglP->Write();
        EffLikDATAAglP->Write();
        return;
}



void DATAQualeff(TFile * file1){
	cout<<"********************************************************** DATA QUALITY SEL. EFFICIENCIES ************************************************************"<<endl;		
	c15->Divide(2,3);
	c16->Divide(2,3);
	float EffLikDATA[43][11]={{0}};
	for(int l=0;l<11;l++)
		for(int i=1;i<43;i++) 
			if(EffQualDATAP->GetBinContent(i+1,l+1)>0)
				EffLikDATA[i][l]=EffLikDATAP->GetBinContent(i+1,l+1)/(float)EffDistDATAP->GetBinContent(i+1,l+1);

	float EffDistDATA[43][11]={0};
	for(int l=0;l<11;l++)
		for(int i=1;i<43;i++) 
			if(EffQualDATAP->GetBinContent(i+1,l+1)>0)
				EffDistDATA[i][l]=EffDistDATAP->GetBinContent(i+1,l+1)/(float)EffQualDATAP->GetBinContent(i+1,l+1);
	
	float EffLikDATANaF[43][11]={{0}};
        for(int l=0;l<11;l++)
                for(int i=1;i<43;i++)
                        if(EffQualDATANaFP->GetBinContent(i+1,l+1)>0)
                                EffLikDATANaF[i][l]=EffLikDATANaFP->GetBinContent(i+1,l+1)/(float)EffDistDATANaFP->GetBinContent(i+1,l+1);

        float EffDistDATANaF[43][11]={0};
        for(int l=0;l<11;l++)
                for(int i=1;i<43;i++)
                        if(EffQualDATANaFP->GetBinContent(i+1,l+1)>0)
                                EffDistDATANaF[i][l]=EffDistDATANaFP->GetBinContent(i+1,l+1)/(float)EffQualDATANaFP->GetBinContent(i+1,l+1);

	float EffLikDATAAgl[43][11]={{0}};
        for(int l=0;l<11;l++)
                for(int i=1;i<43;i++)
                        if(EffQualDATAAglP->GetBinContent(i+1,l+1)>0)
                                EffLikDATAAgl[i][l]=EffLikDATAAglP->GetBinContent(i+1,l+1)/(float)EffDistDATAAglP->GetBinContent(i+1,l+1);

        float EffDistDATAAgl[43][11]={0};
        for(int l=0;l<11;l++)
                for(int i=1;i<43;i++)
                        if(EffQualDATAAglP->GetBinContent(i+1,l+1)>0)
                                EffDistDATAAgl[i][l]=EffDistDATAAglP->GetBinContent(i+1,l+1)/(float)EffQualDATAAglP->GetBinContent(i+1,l+1);

	cout<<"************************************************************ LAT. Eff. CORRECTION ************************************************************"<<endl;

        double Tau_DLik[11]={0};
        double Tau_DLik_err[11]={0};
        float HELikeff1[11]={0};
        float HELikeff2[11]={0};
        float HELikeff[11]={0};
	double Tau_DDist[11]={0};
        double Tau_DDist_err[11]={0};
        float HEDisteff1[11]={0};
        float HEDisteff2[11]={0};
        float HEDisteff[11]={0};

                for(int i=1;i<11;i++)
                        for(int j=30;j<43;j++) {
                                HEDisteff1[i]+=EffDistDATAP->GetBinContent(j+1,i+1);
				HEDisteff2[i]+=EffQualDATAP->GetBinContent(j+1,i+1);
				HELikeff1[i]+=EffLikDATAP->GetBinContent(j+1,i+1);
                                HELikeff2[i]+=EffDistDATAP->GetBinContent(j+1,i+1);
				}
                for(int i=1;i<11;i++){
					HELikeff[i]=HELikeff1[i]/HELikeff2[i];
					HEDisteff[i]=HEDisteff1[i]/HEDisteff2[i];	
				}

                for(int i=1;i<11;i++) {
					Tau_DLik[i]=(HELikeff[1])/HELikeff[i];
                                        Tau_DLik_err[i]=pow(HELikeff1[i],0.5)/HELikeff1[i]*Tau_DLik[i];
					Tau_DDist[i]=(HEDisteff[1])/HEDisteff[i];
                                        Tau_DDist_err[i]=pow(HEDisteff1[i],0.5)/HEDisteff1[i]*Tau_DDist[i];
                                        }
	

	double Tau_DLikNaF[11]={0};
        double Tau_DLikNaF_err[11]={0};
        float HELikeffNaF1[11]={0};
        float HELikeffNaF2[11]={0};
        float HELikeffNaF[11]={0};
        double Tau_DDistNaF[11]={0};
        double Tau_DDistNaF_err[11]={0};
        float HEDisteffNaF1[11]={0};
        float HEDisteffNaF2[11]={0};
        float HEDisteffNaF[11]={0};

                for(int i=1;i<11;i++)
                        for(int j=30;j<43;j++) {
                                HEDisteffNaF1[i]+=EffDistDATANaFP->GetBinContent(j+1,i+1);
                                HEDisteffNaF2[i]+=EffQualDATANaFP->GetBinContent(j+1,i+1);
                                HELikeffNaF1[i]+=EffLikDATANaFP->GetBinContent(j+1,i+1);
                                HELikeffNaF2[i]+=EffDistDATANaFP->GetBinContent(j+1,i+1);
                                }
                for(int i=1;i<11;i++){
                                        HELikeffNaF[i]=HELikeffNaF1[i]/HELikeffNaF2[i];
                                        HEDisteffNaF[i]=HEDisteffNaF1[i]/HEDisteffNaF2[i];
                                }

                for(int i=1;i<11;i++) {
                                        Tau_DLikNaF[i]=(HELikeffNaF[1])/HELikeffNaF[i];
                                        Tau_DLikNaF_err[i]=pow(HELikeffNaF1[i],0.5)/HELikeffNaF1[i]*Tau_DLikNaF[i];
                                        Tau_DDistNaF[i]=(HEDisteffNaF[1])/HEDisteffNaF[i];
                                        Tau_DDistNaF_err[i]=pow(HEDisteffNaF1[i],0.5)/HEDisteffNaF1[i]*Tau_DDistNaF[i];
                                        }
	
	double Tau_DLikAgl[11]={0};
        double Tau_DLikAgl_err[11]={0};
        float HELikeffAgl1[11]={0};
        float HELikeffAgl2[11]={0};
        float HELikeffAgl[11]={0};
        double Tau_DDistAgl[11]={0};
        double Tau_DDistAgl_err[11]={0};
        float HEDisteffAgl1[11]={0};
        float HEDisteffAgl2[11]={0};
        float HEDisteffAgl[11]={0};

                for(int i=1;i<11;i++) 
                        for(int j=30;j<43;j++) {
                                HEDisteffAgl1[i]+=EffDistDATAAglP->GetBinContent(j+1,i+1);
                                HEDisteffAgl2[i]+=EffQualDATAAglP->GetBinContent(j+1,i+1);
                                HELikeffAgl1[i]+=EffLikDATAAglP->GetBinContent(j+1,i+1);
                                HELikeffAgl2[i]+=EffDistDATAAglP->GetBinContent(j+1,i+1);
                                }
                for(int i=1;i<11;i++){
                                        HELikeffAgl[i]=HELikeffAgl1[i]/HELikeffAgl2[i];
                                        HEDisteffAgl[i]=HEDisteffAgl1[i]/HEDisteffAgl2[i];
                                }

                for(int i=1;i<11;i++) {
                                        Tau_DLikAgl[i]=(HELikeffAgl[1])/HELikeffAgl[i];
                                        Tau_DLikAgl_err[i]=pow(HELikeffAgl1[i],0.5)/HELikeffAgl1[i]*Tau_DLikAgl[i];
                                        Tau_DDistAgl[i]=(HEDisteffAgl[1])/HEDisteffAgl[i];
                                        Tau_DDistAgl_err[i]=pow(HEDisteffAgl1[i],0.5)/HEDisteffAgl1[i]*Tau_DDistAgl[i];
                                        }

	cout<<"***** DRAWING *******"<<endl;
	c15->cd(1);
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraphErrors *EffDATALikP[11];
	for(int l=0;l<11;l++){
		EffDATALikP[l]=new TGraphErrors();
		int j=0;
		for(int i=1;i<43;i++) if(EffLikDATA[i][l]>0){
			EffDATALikP[l]->SetPoint(j,R_cent[i],EffLikDATA[i][l]);
			EffDATALikP[l]->SetPointError(j,0,pow(EffLikDATAP->GetBinContent(i+1,l+1),0.5)/EffLikDATAP->GetBinContent(i+1,l+1)*EffLikDATA[i][l]);
			j++;
		}	
	}
	cout<<endl;
	EffDATALikP[10]->SetMarkerColor(1);
	EffDATALikP[10]->SetMarkerStyle(8);
	EffDATALikP[10]->SetLineColor(1);
	EffDATALikP[10]->SetTitle("Likelihood Latitude Efficiency ");
	EffDATALikP[10]->GetXaxis()->SetTitle("R [GV]");
	EffDATALikP[10]->GetYaxis()->SetTitle("Efficiency");
	EffDATALikP[10]->GetXaxis()->SetTitleSize(0.045);
	EffDATALikP[10]->GetYaxis()->SetTitleSize(0.045);
	EffDATALikP[10]->Draw("AP");
	for(int l=0;l<10;l++){
		EffDATALikP[l]->SetMarkerColor(l);
		EffDATALikP[l]->SetMarkerStyle(8);	
		EffDATALikP[l]->SetLineColor(l);
		EffDATALikP[l]->Draw("Psame");
	}
	


	c15->cd(3);
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
        TGraphErrors *EffDATALikNaFP[11];
        for(int l=0;l<11;l++){
                EffDATALikNaFP[l]=new TGraphErrors();
                int j=0;
                for(int i=1;i<43;i++) if(EffLikDATANaF[i][l]>0){
                        EffDATALikNaFP[l]->SetPoint(j,R_cent[i],EffLikDATANaF[i][l]);
                        EffDATALikNaFP[l]->SetPointError(j,0,pow(EffLikDATANaFP->GetBinContent(i+1,l+1),0.5)/EffLikDATANaFP->GetBinContent(i+1,l+1)*EffLikDATANaF[i][l]);
                        j++;
                }
        }
        EffDATALikNaFP[10]->SetMarkerColor(1);
        EffDATALikNaFP[10]->SetMarkerStyle(8);
        EffDATALikNaFP[10]->SetLineColor(1);
        EffDATALikNaFP[10]->SetTitle("Likelihood Latitude Efficiency (NaF)");
        EffDATALikNaFP[10]->GetXaxis()->SetTitle("R [GV]");
        EffDATALikNaFP[10]->GetYaxis()->SetTitle("Efficiency");
        EffDATALikNaFP[10]->GetXaxis()->SetTitleSize(0.045);
        EffDATALikNaFP[10]->GetYaxis()->SetTitleSize(0.045);
        EffDATALikNaFP[10]->Draw("AP");
        for(int l=0;l<10;l++){
                EffDATALikNaFP[l]->SetMarkerColor(l);
                EffDATALikNaFP[l]->SetMarkerStyle(8);
                EffDATALikNaFP[l]->SetLineColor(l);
                EffDATALikNaFP[l]->Draw("Psame");
        }

        c15->cd(5);
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
        TGraphErrors *EffDATALikAglP[11];
        for(int l=0;l<11;l++){
                EffDATALikAglP[l]=new TGraphErrors();
                int j=0;
                for(int i=1;i<43;i++) if(EffLikDATAAgl[i][l]>0){
                        EffDATALikAglP[l]->SetPoint(j,R_cent[i],EffLikDATAAgl[i][l]);
                        EffDATALikAglP[l]->SetPointError(j,0,pow(EffLikDATAAglP->GetBinContent(i+1,l+1),0.5)/EffLikDATAAglP->GetBinContent(i+1,l+1)*EffLikDATAAgl[i][l]);
                        j++;
                }
        }
        EffDATALikAglP[10]->SetMarkerColor(1);
        EffDATALikAglP[10]->SetMarkerStyle(8);
        EffDATALikAglP[10]->SetLineColor(1);
        EffDATALikAglP[10]->SetTitle("Likelihood Latitude Efficiency (Agl)");
        EffDATALikAglP[10]->GetXaxis()->SetTitle("R [GV]");
        EffDATALikAglP[10]->GetYaxis()->SetTitle("Efficiency");
        EffDATALikAglP[10]->GetXaxis()->SetTitleSize(0.045);
        EffDATALikAglP[10]->GetYaxis()->SetTitleSize(0.045);
        EffDATALikAglP[10]->Draw("AP");
        for(int l=0;l<10;l++){
                EffDATALikAglP[l]->SetMarkerColor(l);
                EffDATALikAglP[l]->SetMarkerStyle(8);
                EffDATALikAglP[l]->SetLineColor(l);
                EffDATALikAglP[l]->Draw("Psame");
        }

	TGraphErrors *CorrLATLik;
        c15->cd(2);
        gPad->SetGridy();
        gPad->SetGridx();
        CorrLAT_Lik=new TF1("Likelihood Lat. Corr.","pol3");
        CorrLATLik=new TGraphErrors();
        CorrLATLik->SetTitle("Latitude Efficiency Corr.");
        CorrLATLik->GetXaxis()->SetTitle("Latitude");
        CorrLATLik->GetYaxis()->SetTitle("Eff. Corr. Factor");
        CorrLATLik->GetYaxis()->SetRangeUser(0.96,1.04);
	CorrLATLik->SetMarkerStyle(8);
        for(int i=1;i<11;i++) {
                        CorrLATLik->SetPoint(i-1,geomagC[i],Tau_DLik[i]);
                        CorrLATLik->SetPointError(i-1,0,Tau_DLik_err[i]);
                }
        CorrLATLik->Fit("Likelihood Lat. Corr.");
        CorrLATLik->Draw("AP");
	string nome="CorrLAT_Lik_spl";
        CorrLAT_Lik_spl=new TH2F(nome.c_str(),nome.c_str(),11,0,11,2,0,2);
        for(int i=0;i<11;i++) CorrLAT_Lik_spl->SetBinContent(i+1,1,CorrLAT_Lik->Eval(geomagC[i]));
        for(int i=0;i<11;i++) CorrLAT_Lik_spl->SetBinContent(i+1,2,Tau_DLik_err[i]);
	CorrLAT_Lik_Spl=new TGraphErrors("CorrLAT_Lik_Spl");
        CorrLAT_Lik_Spl->SetName("CorrLAT_Lik_Spl");
	for(int i=0;i<10;i++) {CorrLAT_Lik_Spl->SetPoint(i,geomagC[i+1],CorrLAT_Lik->Eval(geomagC[i+1]));
                                CorrLAT_Lik_Spl->SetPointError(i,0,Tau_DLik_err[i+1]);
                                }
        CorrLAT_Lik_Spl->SetLineColor(2);
        CorrLAT_Lik_Spl->SetMarkerColor(2);
        CorrLAT_Lik_Spl->SetFillColor(2);
        CorrLAT_Lik_Spl->SetFillStyle(3001);
        CorrLAT_Lik_Spl->SetLineWidth(2);
        CorrLAT_Lik_Spl->Draw("Lsame");


	TGraphErrors *CorrLATLikNaF;
        c15->cd(4);
        gPad->SetGridy();
        gPad->SetGridx();
        CorrLAT_LikNaF=new TF1("NaF Likelihood Lat. Corr.","pol3");
        CorrLATLikNaF=new TGraphErrors();
        CorrLATLikNaF->SetTitle("NaF Latitude Efficiency Corr.");
        CorrLATLikNaF->GetXaxis()->SetTitle("Latitude");
        CorrLATLikNaF->GetYaxis()->SetTitle("Eff. Corr. Factor");
        CorrLATLikNaF->GetYaxis()->SetRangeUser(0.96,1.04);
        CorrLATLikNaF->SetMarkerStyle(8);
        for(int i=1;i<11;i++) {
                        CorrLATLikNaF->SetPoint(i-1,geomagC[i],Tau_DLikNaF[i]);
                        CorrLATLikNaF->SetPointError(i-1,0,Tau_DLikNaF_err[i]);
                }
        CorrLATLikNaF->Fit("NaF Likelihood Lat. Corr.");
        CorrLATLikNaF->Draw("AP");
        nome="CorrLAT_LikNaF_spl";
        CorrLAT_LikNaF_spl=new TH2F(nome.c_str(),nome.c_str(),11,0,11,2,0,2);
        for(int i=0;i<11;i++) CorrLAT_LikNaF_spl->SetBinContent(i+1,1,CorrLAT_LikNaF->Eval(geomagC[i]));
        for(int i=0;i<11;i++) CorrLAT_LikNaF_spl->SetBinContent(i+1,2,Tau_DLikNaF_err[i]);
        CorrLAT_LikNaF_Spl=new TGraphErrors("CorrLAT_LikNaF_Spl");
        CorrLAT_LikNaF_Spl->SetName("CorrLAT_LikNaF_Spl");
        for(int i=0;i<10;i++) {CorrLAT_LikNaF_Spl->SetPoint(i,geomagC[i+1],CorrLAT_LikNaF->Eval(geomagC[i+1]));
                                CorrLAT_LikNaF_Spl->SetPointError(i,0,Tau_DLikNaF_err[i+1]);
                                }
        CorrLAT_LikNaF_Spl->SetLineColor(2);
        CorrLAT_LikNaF_Spl->SetMarkerColor(2);
        CorrLAT_LikNaF_Spl->SetFillColor(2);
        CorrLAT_LikNaF_Spl->SetFillStyle(3001);
        CorrLAT_LikNaF_Spl->SetLineWidth(2);
        CorrLAT_LikNaF_Spl->Draw("Lsame");

	TGraphErrors *CorrLATLikAgl;
        c15->cd(6);
        gPad->SetGridy();
        gPad->SetGridx();
        CorrLAT_LikAgl=new TF1("Agl Likelihood Lat. Corr.","pol3");
        CorrLATLikAgl=new TGraphErrors();
        CorrLATLikAgl->SetTitle("Agl Latitude Efficiency Corr.");
        CorrLATLikAgl->GetXaxis()->SetTitle("Latitude");
        CorrLATLikAgl->GetYaxis()->SetTitle("Eff. Corr. Factor");
        CorrLATLikAgl->GetYaxis()->SetRangeUser(0.96,1.04);
        CorrLATLikAgl->SetMarkerStyle(8);
        for(int i=1;i<11;i++) {
                        CorrLATLikAgl->SetPoint(i-1,geomagC[i],Tau_DLikAgl[i]);
                        CorrLATLikAgl->SetPointError(i-1,0,Tau_DLikAgl_err[i]);
                }
        CorrLATLikAgl->Fit("Agl Likelihood Lat. Corr.");
        CorrLATLikAgl->Draw("AP");
        nome="CorrLAT_LikAgl_spl";
        CorrLAT_LikAgl_spl=new TH2F(nome.c_str(),nome.c_str(),11,0,11,2,0,2);
        for(int i=0;i<11;i++) CorrLAT_LikAgl_spl->SetBinContent(i+1,1,CorrLAT_LikAgl->Eval(geomagC[i]));
        for(int i=0;i<11;i++) CorrLAT_LikAgl_spl->SetBinContent(i+1,2,Tau_DLikAgl_err[i]);
        CorrLAT_LikAgl_Spl=new TGraphErrors("CorrLAT_LikAgl_Spl");
        CorrLAT_LikAgl_Spl->SetName("CorrLAT_LikAgl_Spl");
        for(int i=0;i<10;i++) {CorrLAT_LikAgl_Spl->SetPoint(i,geomagC[i+1],CorrLAT_LikAgl->Eval(geomagC[i+1]));
                                CorrLAT_LikAgl_Spl->SetPointError(i,0,Tau_DLikAgl_err[i+1]);
                                }
        CorrLAT_LikAgl_Spl->SetLineColor(2);
        CorrLAT_LikAgl_Spl->SetMarkerColor(2);
        CorrLAT_LikAgl_Spl->SetFillColor(2);
        CorrLAT_LikAgl_Spl->SetFillStyle(3001);
        CorrLAT_LikAgl_Spl->SetLineWidth(2);
        CorrLAT_LikAgl_Spl->Draw("Lsame");


	TGraphErrors *CorrLATDistNaF;
        c16->cd(4);
        gPad->SetGridy();
        gPad->SetGridx();
        CorrLAT_DistNaF=new TF1("NaF Distance Lat. Corr.","pol3");
        CorrLATDistNaF=new TGraphErrors();
        CorrLATDistNaF->SetTitle("NaF Latitude Efficiency Corr.");
        CorrLATDistNaF->GetXaxis()->SetTitle("Latitude");
        CorrLATDistNaF->GetYaxis()->SetTitle("Eff. Corr. Factor");
        CorrLATDistNaF->GetYaxis()->SetRangeUser(0.96,1.04);
        CorrLATDistNaF->SetMarkerStyle(8);
        for(int i=1;i<11;i++) {
                        CorrLATDistNaF->SetPoint(i-1,geomagC[i],Tau_DDistNaF[i]);
                        CorrLATDistNaF->SetPointError(i-1,0,Tau_DDistNaF_err[i]);
                }
        CorrLATDistNaF->Fit("NaF Distance Lat. Corr.");
        CorrLATDistNaF->Draw("AP");
        nome="CorrLAT_DistNaF_spl";
        CorrLAT_DistNaF_spl=new TH2F(nome.c_str(),nome.c_str(),11,0,11,2,0,2);
        for(int i=0;i<11;i++) CorrLAT_DistNaF_spl->SetBinContent(i+1,1,CorrLAT_DistNaF->Eval(geomagC[i]));
        for(int i=0;i<11;i++) CorrLAT_DistNaF_spl->SetBinContent(i+1,2,Tau_DDistNaF_err[i]);
        CorrLAT_DistNaF_Spl=new TGraphErrors("CorrLAT_DistNaF_Spl");
        CorrLAT_DistNaF_Spl->SetName("CorrLAT_DistNaF_Spl");
        for(int i=0;i<10;i++) {CorrLAT_DistNaF_Spl->SetPoint(i,geomagC[i+1],CorrLAT_DistNaF->Eval(geomagC[i+1]));
                                CorrLAT_DistNaF_Spl->SetPointError(i,0,Tau_DDistNaF_err[i+1]);
                                }
        CorrLAT_DistNaF_Spl->SetLineColor(2);
        CorrLAT_DistNaF_Spl->SetMarkerColor(2);
        CorrLAT_DistNaF_Spl->SetFillColor(2);
        CorrLAT_DistNaF_Spl->SetFillStyle(3001);
        CorrLAT_DistNaF_Spl->SetLineWidth(2);
        CorrLAT_DistNaF_Spl->Draw("Lsame");

        TGraphErrors *CorrLATDistAgl;
        c16->cd(6);
        gPad->SetGridy();
        gPad->SetGridx();
        CorrLAT_DistAgl=new TF1("Agl Distance Lat. Corr.","pol3");
        CorrLATDistAgl=new TGraphErrors();
        CorrLATDistAgl->SetTitle("Agl Latitude Efficiency Corr.");
        CorrLATDistAgl->GetXaxis()->SetTitle("Latitude");
        CorrLATDistAgl->GetYaxis()->SetTitle("Eff. Corr. Factor");
        CorrLATDistAgl->GetYaxis()->SetRangeUser(0.96,1.04);
        CorrLATDistAgl->SetMarkerStyle(8);
        for(int i=1;i<11;i++) {
                        CorrLATDistAgl->SetPoint(i-1,geomagC[i],Tau_DDistAgl[i]);
                        CorrLATDistAgl->SetPointError(i-1,0,Tau_DDistAgl_err[i]);
                }
        CorrLATDistAgl->Fit("Agl Distance Lat. Corr.");
        CorrLATDistAgl->Draw("AP");
        nome="CorrLAT_DistAgl_spl";
        CorrLAT_DistAgl_spl=new TH2F(nome.c_str(),nome.c_str(),11,0,11,2,0,2);
        for(int i=0;i<11;i++) CorrLAT_DistAgl_spl->SetBinContent(i+1,1,CorrLAT_DistAgl->Eval(geomagC[i]));
        for(int i=0;i<11;i++) CorrLAT_DistAgl_spl->SetBinContent(i+1,2,Tau_DDistAgl_err[i]);
        CorrLAT_DistAgl_Spl=new TGraphErrors("CorrLAT_DistAgl_Spl");
        CorrLAT_DistAgl_Spl->SetName("CorrLAT_DistAgl_Spl");
        for(int i=0;i<10;i++) {CorrLAT_DistAgl_Spl->SetPoint(i,geomagC[i+1],CorrLAT_DistAgl->Eval(geomagC[i+1]));
                                CorrLAT_DistAgl_Spl->SetPointError(i,0,Tau_DDistAgl_err[i+1]);
                                }
        CorrLAT_DistAgl_Spl->SetLineColor(2);
        CorrLAT_DistAgl_Spl->SetMarkerColor(2);
        CorrLAT_DistAgl_Spl->SetFillColor(2);
        CorrLAT_DistAgl_Spl->SetFillStyle(3001);
        CorrLAT_DistAgl_Spl->SetLineWidth(2);
        CorrLAT_DistAgl_Spl->Draw("Lsame");


	c16->cd(1);
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraphErrors *EffDATADistP[11];
	for(int l=0;l<11;l++){
		EffDATADistP[l]=new TGraphErrors();
		int j=0;
		for(int i=1;i<43;i++) if(EffDistDATA[i][l]>0){
			EffDATADistP[l]->SetPoint(j,R_cent[i],EffDistDATA[i][l]);
			EffDATADistP[l]->SetPointError(j,0,pow(EffDistDATAP->GetBinContent(i+1,l+1),0.5)/EffDistDATAP->GetBinContent(i+1,l+1)*EffDistDATA[i][l]);
			j++;
		}
	}
	EffDATADistP[10]->SetMarkerColor(1);
	EffDATADistP[10]->SetMarkerStyle(8);
	EffDATADistP[10]->SetLineColor(1);
	EffDATADistP[10]->SetTitle("Distance Latitude Efficiency (TOF)");
	EffDATADistP[10]->GetXaxis()->SetTitle("R [GV]");
	EffDATADistP[10]->GetYaxis()->SetTitle("Efficiency");
	EffDATADistP[10]->GetXaxis()->SetTitleSize(0.045);
	EffDATADistP[10]->GetYaxis()->SetTitleSize(0.045);
	EffDATADistP[10]->Draw("AP");
	for(int l=0;l<10;l++){
		EffDATADistP[l]->SetMarkerColor(l); 
		EffDATADistP[l]->SetMarkerStyle(8);
		EffDATADistP[l]->SetLineColor(l);
		EffDATADistP[l]->Draw("Psame");
	}
	
	c16->cd(3);
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
        TGraphErrors *EffDATADistNaFP[11];
        for(int l=0;l<11;l++){
                EffDATADistNaFP[l]=new TGraphErrors();
                int j=0;
                for(int i=1;i<43;i++) if(EffDistDATANaF[i][l]>0){
                        EffDATADistNaFP[l]->SetPoint(j,R_cent[i],EffDistDATANaF[i][l]);
                        EffDATADistNaFP[l]->SetPointError(j,0,pow(EffDistDATANaFP->GetBinContent(i+1,l+1),0.5)/EffDistDATANaFP->GetBinContent(i+1,l+1)*EffDistDATANaF[i][l]);
                        j++;
                }
        }
        EffDATADistNaFP[10]->SetMarkerColor(1);
        EffDATADistNaFP[10]->SetMarkerStyle(8);
        EffDATADistNaFP[10]->SetLineColor(1);
        EffDATADistNaFP[10]->SetTitle("Distance Latitude Efficiency (NaF)");
        EffDATADistNaFP[10]->GetXaxis()->SetTitle("R [GV]");
        EffDATADistNaFP[10]->GetYaxis()->SetTitle("Efficiency");
        EffDATADistNaFP[10]->GetXaxis()->SetTitleSize(0.045);
        EffDATADistNaFP[10]->GetYaxis()->SetTitleSize(0.045);
        EffDATADistNaFP[10]->Draw("AP");
        for(int l=0;l<10;l++){
                EffDATADistNaFP[l]->SetMarkerColor(l);
                EffDATADistNaFP[l]->SetMarkerStyle(8);
                EffDATADistNaFP[l]->SetLineColor(l);
                EffDATADistNaFP[l]->Draw("Psame");
        }

	c16->cd(5);
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
        TGraphErrors *EffDATADistAglP[11];
        for(int l=0;l<11;l++){
                EffDATADistAglP[l]=new TGraphErrors();
                int j=0;
                for(int i=1;i<43;i++) if(EffDistDATAAgl[i][l]>0){
                        EffDATADistAglP[l]->SetPoint(j,R_cent[i],EffDistDATAAgl[i][l]);
                        EffDATADistAglP[l]->SetPointError(j,0,pow(EffDistDATAAglP->GetBinContent(i+1,l+1),0.5)/EffDistDATAAglP->GetBinContent(i+1,l+1)*EffDistDATAAgl[i][l]);
                        j++;
                }
        }
        EffDATADistAglP[10]->SetMarkerColor(1);
        EffDATADistAglP[10]->SetMarkerStyle(8);
        EffDATADistAglP[10]->SetLineColor(1);
        EffDATADistAglP[10]->SetTitle("Distance Latitude Efficiency (Agl)");
        EffDATADistAglP[10]->GetXaxis()->SetTitle("R [GV]");
        EffDATADistAglP[10]->GetYaxis()->SetTitle("Efficiency");
        EffDATADistAglP[10]->GetXaxis()->SetTitleSize(0.045);
        EffDATADistAglP[10]->GetYaxis()->SetTitleSize(0.045);
        EffDATADistAglP[10]->Draw("AP");
        for(int l=0;l<10;l++){
                EffDATADistAglP[l]->SetMarkerColor(l);
                EffDATADistAglP[l]->SetMarkerStyle(8);
                EffDATADistAglP[l]->SetLineColor(l);
                EffDATADistAglP[l]->Draw("Psame");
        }
             	
	TGraphErrors *CorrLATDist;
        c16->cd(2);
        gPad->SetGridy();
        gPad->SetGridx();
        CorrLAT_Dist=new TF1("Distance Lat. Corr.","pol4");
        CorrLATDist=new TGraphErrors();
        CorrLATDist->SetTitle("Latitude Efficiency Corr.");
        CorrLATDist->GetXaxis()->SetTitle("Latitude");
        CorrLATDist->GetYaxis()->SetTitle("Eff. Corr. Factor");
	CorrLATDist->GetYaxis()->SetRangeUser(0.96,1.04);
        CorrLATDist->SetMarkerStyle(8);
        for(int i=1;i<11;i++) {
                        CorrLATDist->SetPoint(i-1,geomagC[i],Tau_DDist[i]);
                        CorrLATDist->SetPointError(i-1,0,Tau_DDist_err[i]);
                }
        CorrLATDist->Fit("Distance Lat. Corr.");
        CorrLATDist->Draw("AP");
	nome="CorrLAT_Dist_spl";
        CorrLAT_Dist_spl=new TH2F(nome.c_str(),nome.c_str(),11,0,11,2,0,2);
        for(int i=0;i<11;i++) CorrLAT_Dist_spl->SetBinContent(i+1,1,CorrLAT_Dist->Eval(geomagC[i]));
	for(int i=0;i<11;i++) CorrLAT_Dist_spl->SetBinContent(i+1,2,Tau_DDist_err[i]);	
	CorrLAT_Dist_Spl=new TGraphErrors("CorrLAT_Dist_Spl");
        CorrLAT_Dist_Spl->SetName("CorrLAT_Dist_Spl");
	for(int i=0;i<10;i++) {CorrLAT_Dist_Spl->SetPoint(i,geomagC[i+1],CorrLAT_Dist->Eval(geomagC[i+1]));
                                CorrLAT_Dist_Spl->SetPointError(i,0,Tau_DDist_err[i+1]);
                                }
        CorrLAT_Dist_Spl->SetLineColor(2);
        CorrLAT_Dist_Spl->SetMarkerColor(2);
        CorrLAT_Dist_Spl->SetFillColor(2);
        CorrLAT_Dist_Spl->SetFillStyle(3001);
        CorrLAT_Dist_Spl->SetLineWidth(2);
        CorrLAT_Dist_Spl->Draw("Lsame");

}


