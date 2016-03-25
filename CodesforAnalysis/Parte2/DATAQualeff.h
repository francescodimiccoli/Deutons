using namespace std;


LATcorr * LikelihoodLATcorr = new LATcorr("EffLikDATA");
LATcorr * DistanceLATcorr   = new LATcorr("EffDistDATA");

void DATAQualeff_Fill(TNtuple *ntupla, int l,int zona){

	int k = ntupla->GetEvent(l);
	if(Beta<=0||R<=0) return;
	if(!(EdepL1>0&&EdepL1<EdepL1beta->Eval(Beta)+0.1&&EdepL1>EdepL1beta->Eval(Beta)-0.1)) return;
	//TOF
	if(!(((int)Cutmask>>11)==0||((int)Cutmask>>11)==512)){
		for(int K=0;K<43;K++) 
			if(R<bin[K+1]&&R>bin[K]&&R>Rcut[zona]) {
					DistanceLATcorr->beforeTOF->Fill(K,zona);				
					if(Dist5D_P<6) {
						DistanceLATcorr->afterTOF->Fill(K,zona);
						LikelihoodLATcorr->beforeTOF->Fill(K,zona);	
						}
					if(Dist5D_P<6&Likcut){
						LikelihoodLATcorr->afterTOF->Fill(K,zona);
						}
			}
		}
	//NaF	
	if(((int)Cutmask>>11)==512){
		for(int K=0;K<43;K++) 
			if(R<bin[K+1]&&R>bin[K]&&R>Rcut[zona]) {
					DistanceLATcorr->beforeNaF->Fill(K,zona);				
					if(Dist5D_P<6) {
						DistanceLATcorr->afterNaF->Fill(K,zona);
						LikelihoodLATcorr->beforeNaF->Fill(K,zona);	
						}
					if(Dist5D_P<6&Likcut){
						LikelihoodLATcorr->afterNaF->Fill(K,zona);
						}
			}
	}
	//Agl
	if(((int)Cutmask>>11)==0){
        		for(int K=0;K<43;K++) 
			if(R<bin[K+1]&&R>bin[K]&&R>Rcut[zona]) {
					DistanceLATcorr->beforeAgl->Fill(K,zona);				
					if(Dist5D_P<6) {
						DistanceLATcorr->afterAgl->Fill(K,zona);
						LikelihoodLATcorr->beforeAgl->Fill(K,zona);	
						}
					if(Dist5D_P<6&Likcut){
						LikelihoodLATcorr->afterAgl->Fill(K,zona);
						}
			}

        }

	return;
}

void DATAQualeff_Write(){
  LikelihoodLATcorr->Write();       
  DistanceLATcorr  ->Write();
  return;
}




void DATAQualeff(TFile * file1){

	LATcorr * LikelihoodLATcorr = new LATcorr(file1,"EffLikDATA");
	LATcorr * DistanceLATcorr   = new LATcorr(file1,"EffDistDATA");


	cout<<"****************************** DATA QUALITY SEL. EFFICIENCIES **************************************"<<endl;		
	
	LikelihoodLATcorr -> Eval_Efficiency(); 
        DistanceLATcorr   -> Eval_Efficiency();
	

	TH2F *EffDistDATATOF = (TH2F *) DistanceLATcorr     -> effTOF -> Clone();
        TH2F *EffLikDATATOF  = (TH2F *) LikelihoodLATcorr   -> effTOF -> Clone();
        TH2F *EffDistDATANaF = (TH2F *) DistanceLATcorr     -> effNaF -> Clone();
        TH2F *EffLikDATANaF  = (TH2F *) LikelihoodLATcorr   -> effNaF -> Clone();
        TH2F *EffDistDATAAgl = (TH2F *) DistanceLATcorr     -> effAgl -> Clone();
        TH2F *EffLikDATAAgl  = (TH2F *) LikelihoodLATcorr   -> effAgl -> Clone();



	cout<<"****************************** LAT. Eff. CORRECTION *************************************************"<<endl;

	LikelihoodLATcorr -> Eval_LATcorr(1); 
        DistanceLATcorr   -> Eval_LATcorr(1);

	TH2F *LikLATcorr_TOF  =	(TH2F *) LikelihoodLATcorr   -> LATcorrTOF -> Clone();
	TH2F *DistLATcorr_TOF =	(TH2F *) DistanceLATcorr     -> LATcorrTOF -> Clone();
	TH2F *LikLATcorr_NaF  =	(TH2F *) LikelihoodLATcorr   -> LATcorrNaF -> Clone();
	TH2F *DistLATcorr_NaF =	(TH2F *) DistanceLATcorr     -> LATcorrNaF -> Clone();
	TH2F *LikLATcorr_Agl  =	(TH2F *) LikelihoodLATcorr   -> LATcorrAgl -> Clone();
	TH2F *DistLATcorr_Agl =	(TH2F *) DistanceLATcorr     -> LATcorrAgl -> Clone();
	

        TH1F *LikLATcorr_TOF_fit  	= (TH1F *) LikelihoodLATcorr   -> LATcorrTOF_fit-> Clone();
        TH1F *DistLATcorr_TOF_fit	= (TH1F *) DistanceLATcorr     -> LATcorrTOF_fit-> Clone();
        TH1F *LikLATcorr_NaF_fit    	= (TH1F *) LikelihoodLATcorr   -> LATcorrNaF_fit-> Clone();
        TH1F *DistLATcorr_NaF_fit  	= (TH1F *) DistanceLATcorr     -> LATcorrNaF_fit-> Clone();
        TH1F *LikLATcorr_Agl_fit   	= (TH1F *) LikelihoodLATcorr   -> LATcorrAgl_fit-> Clone();
        TH1F *DistLATcorr_Agl_fit  	= (TH1F *) DistanceLATcorr     -> LATcorrAgl_fit-> Clone();
	
	
	cout<<"*** Updating P1 file ****"<<endl;
        string nomefile=percorso + "/Risultati/risultati/"+mese+"_"+frac+"_P1.root";
        file1 =TFile::Open(nomefile.c_str(),"UPDATE");
        if(!file1){
                nomefile=percorso + "/Risultati/"+mese+"/"+mese+"_"+frac+"_P1.root";
                file1 =TFile::Open(nomefile.c_str(),"UPDATE");
        }

        file1->cd("Results");

	EffDistDATATOF ->Write();
	EffLikDATATOF  ->Write();
	EffDistDATANaF ->Write();
	EffLikDATANaF  ->Write();
	EffDistDATAAgl ->Write();
	EffLikDATAAgl  ->Write();
        
	LikLATcorr_TOF  ->Write();
        DistLATcorr_TOF ->Write();
        LikLATcorr_NaF  ->Write();
        DistLATcorr_NaF ->Write();
        LikLATcorr_Agl  ->Write();
        DistLATcorr_Agl ->Write();

	LikLATcorr_TOF_fit   ->Write();
        DistLATcorr_TOF_fit  ->Write();
        LikLATcorr_NaF_fit   ->Write();
        DistLATcorr_NaF_fit  ->Write();
        LikLATcorr_Agl_fit   ->Write();
	DistLATcorr_Agl_fit  ->Write();
	
	file1->Write();
        file1->Close();
	


	TCanvas *c15=new TCanvas("Latitude Likelihood Efficiency");
	TCanvas *c16=new TCanvas("Latitude Distance Efficiency");

	c15->Divide(2,3);
	c16->Divide(2,3);

	c15->cd(1);
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraphErrors *EffDATALikP[11];
	for(int l=0;l<11;l++){
		EffDATALikP[l]=new TGraphErrors();
		int j=0;
		for(int i=1;i<43;i++) {
			EffDATALikP[l]->SetPoint(j,R_cent[i],EffLikDATATOF->GetBinContent(i+1,l+1));
			EffDATALikP[l]->SetPointError(j,0,EffLikDATATOF->GetBinError(i+1,l+1));
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
	EffDATALikP[10]->GetYaxis()->SetRangeUser(0.1,1.1);
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
                for(int i=1;i<43;i++) {
                        EffDATALikNaFP[l]->SetPoint(j,R_cent[i],EffLikDATANaF ->GetBinContent(i+1,l+1));
                        EffDATALikNaFP[l]->SetPointError(j,0,EffLikDATANaF ->GetBinError(i+1,l+1));
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
        EffDATALikNaFP[10]->GetYaxis()->SetRangeUser(0.1,1.1);
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
                for(int i=1;i<43;i++) {
                        EffDATALikAglP[l]->SetPoint(j,R_cent[i],EffLikDATAAgl ->GetBinContent(i+1,l+1));
                        EffDATALikAglP[l]->SetPointError(j,0,EffLikDATAAgl ->GetBinError(i+1,l+1));
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
        EffDATALikAglP[10]->GetYaxis()->SetRangeUser(0.1,1.1);
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
        CorrLATLik=new TGraphErrors();
        CorrLATLik->SetTitle("Latitude Efficiency Corr.");
        CorrLATLik->GetXaxis()->SetTitle("Latitude");
        CorrLATLik->GetYaxis()->SetTitle("Eff. Corr. Factor");
        CorrLATLik->GetYaxis()->SetRangeUser(0.96,1.04);
	CorrLATLik->SetMarkerStyle(8);
        for(int i=1;i<11;i++) {
                        CorrLATLik->SetPoint(i-1,geomagC[i],LikLATcorr_TOF->GetBinContent(i+1,1));
                        CorrLATLik->SetPointError(i-1,0,LikLATcorr_TOF->GetBinError(i+1,1));
                }
	CorrLATLik->Draw("AP");
	TGraphErrors *CorrLAT_Lik_Spl=new TGraphErrors("CorrLAT_Lik_Spl");
        CorrLAT_Lik_Spl->SetName("CorrLAT_Lik_Spl");
	for(int i=1;i<11;i++) { CorrLAT_Lik_Spl->SetPoint(i-1,geomagC[i],LikLATcorr_TOF_fit->GetBinContent(i+1));
                                CorrLAT_Lik_Spl->SetPointError(i-1,0,LikLATcorr_TOF_fit->GetBinError(i+1));
                                }
        CorrLAT_Lik_Spl->SetLineColor(2);
        CorrLAT_Lik_Spl->SetMarkerColor(2);
        CorrLAT_Lik_Spl->SetFillColor(2);
        CorrLAT_Lik_Spl->SetFillStyle(3001);
        CorrLAT_Lik_Spl->SetLineWidth(2);
        CorrLAT_Lik_Spl->Draw("Csame");


	TGraphErrors *CorrLATLikNaF;
        c15->cd(4);
        gPad->SetGridy();
        gPad->SetGridx();
        CorrLATLikNaF=new TGraphErrors();
        CorrLATLikNaF->SetTitle("Latitude Efficiency Corr.");
        CorrLATLikNaF->GetXaxis()->SetTitle("Latitude");
        CorrLATLikNaF->GetYaxis()->SetTitle("Eff. Corr. Factor");
        CorrLATLikNaF->GetYaxis()->SetRangeUser(0.96,1.04);
        CorrLATLikNaF->SetMarkerStyle(8);
        for(int i=1;i<11;i++) {
                        CorrLATLikNaF->SetPoint(i-1,geomagC[i],LikLATcorr_NaF->GetBinContent(i+1,1));
                        CorrLATLikNaF->SetPointError(i-1,0,LikLATcorr_NaF->GetBinError(i+1,1));
                }
        CorrLATLikNaF->Draw("AP");
        TGraphErrors *CorrLAT_LikNaF_Spl=new TGraphErrors("CorrLAT_LikNaF_Spl");
        CorrLAT_LikNaF_Spl->SetName("CorrLAT_LikNaF_Spl");
        for(int i=1;i<11;i++) { CorrLAT_LikNaF_Spl->SetPoint(i-1,geomagC[i],LikLATcorr_NaF_fit->GetBinContent(i+1));
                                CorrLAT_LikNaF_Spl->SetPointError(i-1,0,LikLATcorr_NaF_fit->GetBinError(i+1));
                                }
        CorrLAT_LikNaF_Spl->SetLineColor(2);
        CorrLAT_LikNaF_Spl->SetMarkerColor(2);
        CorrLAT_LikNaF_Spl->SetFillColor(2);
        CorrLAT_LikNaF_Spl->SetFillStyle(3001);
        CorrLAT_LikNaF_Spl->SetLineWidth(2);
        CorrLAT_LikNaF_Spl->Draw("Csame");

	TGraphErrors *CorrLATLikAgl;
        c15->cd(6);
        gPad->SetGridy();
        gPad->SetGridx();
        CorrLATLikAgl=new TGraphErrors();
        CorrLATLikAgl->SetTitle("Latitude Efficiency Corr.");
        CorrLATLikAgl->GetXaxis()->SetTitle("Latitude");
        CorrLATLikAgl->GetYaxis()->SetTitle("Eff. Corr. Factor");
        CorrLATLikAgl->GetYaxis()->SetRangeUser(0.96,1.04);
        CorrLATLikAgl->SetMarkerStyle(8);
        for(int i=1;i<11;i++) {
                        CorrLATLikAgl->SetPoint(i-1,geomagC[i],LikLATcorr_Agl->GetBinContent(i+1,1));
                        CorrLATLikAgl->SetPointError(i-1,0,LikLATcorr_Agl->GetBinError(i+1,1));
                }
        CorrLATLikAgl->Draw("AP");
        TGraphErrors *CorrLAT_LikAgl_Spl=new TGraphErrors("CorrLAT_LikAgl_Spl");
        CorrLAT_LikAgl_Spl->SetName("CorrLAT_LikAgl_Spl");
        for(int i=1;i<11;i++) { CorrLAT_LikAgl_Spl->SetPoint(i-1,geomagC[i],LikLATcorr_Agl_fit->GetBinContent(i+1));
                                CorrLAT_LikAgl_Spl->SetPointError(i-1,0,LikLATcorr_Agl_fit->GetBinError(i+1));
                                }
        CorrLAT_LikAgl_Spl->SetLineColor(2);
        CorrLAT_LikAgl_Spl->SetMarkerColor(2);
        CorrLAT_LikAgl_Spl->SetFillColor(2);
        CorrLAT_LikAgl_Spl->SetFillStyle(3001);
        CorrLAT_LikAgl_Spl->SetLineWidth(2);
        CorrLAT_LikAgl_Spl->Draw("Csame");
	


	c16->cd(1);
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraphErrors *EffDATADistP[11];
	for(int l=0;l<11;l++){
		EffDATADistP[l]=new TGraphErrors();
		int j=0;
		for(int i=1;i<43;i++) {
			EffDATADistP[l]->SetPoint(j,R_cent[i],EffDistDATATOF ->GetBinContent(i+1,l+1));
			EffDATADistP[l]->SetPointError(j,0, EffDistDATATOF ->GetBinError(i+1,l+1));
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
	EffDATADistP[10]->GetYaxis()->SetRangeUser(0.1,1.1);
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
                for(int i=1;i<43;i++) {
                        EffDATADistNaFP[l]->SetPoint(j,R_cent[i],EffDistDATANaF ->GetBinContent(i+1,l+1));
			EffDATADistNaFP[l]->SetPointError(j,0, EffDistDATANaF ->GetBinError(i+1,l+1));
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
        EffDATADistNaFP[10]->GetYaxis()->SetRangeUser(0.1,1.1);
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
                for(int i=1;i<43;i++) {
                        EffDATADistAglP[l]->SetPoint(j,R_cent[i],EffDistDATAAgl ->GetBinContent(i+1,l+1));
			EffDATADistAglP[l]->SetPointError(j,0, EffDistDATAAgl ->GetBinError(i+1,l+1));
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
        EffDATADistAglP[10]->GetYaxis()->SetRangeUser(0.1,1.1);
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
        CorrLATDist=new TGraphErrors();
        CorrLATDist->SetTitle("Latitude Efficiency Corr.");
        CorrLATDist->GetXaxis()->SetTitle("Latitude");
        CorrLATDist->GetYaxis()->SetTitle("Eff. Corr. Factor");
        CorrLATDist->GetYaxis()->SetRangeUser(0.96,1.04);
        CorrLATDist->SetMarkerStyle(8);
        for(int i=1;i<11;i++) {
                        CorrLATDist->SetPoint(i-1,geomagC[i],DistLATcorr_TOF->GetBinContent(i+1,1));
                        CorrLATDist->SetPointError(i-1,0,DistLATcorr_TOF->GetBinError(i+1,1));
                }
        CorrLATDist->Draw("AP");
        TGraphErrors *CorrLAT_Dist_Spl=new TGraphErrors("CorrLAT_Dist_Spl");
        CorrLAT_Dist_Spl->SetName("CorrLAT_Dist_Spl");
        for(int i=1;i<11;i++) { CorrLAT_Dist_Spl->SetPoint(i-1,geomagC[i],DistLATcorr_TOF_fit->GetBinContent(i+1,1));
                                CorrLAT_Dist_Spl->SetPointError(i-1,0,DistLATcorr_TOF_fit->GetBinError(i+1,1));
                                }
        CorrLAT_Dist_Spl->SetLineColor(2);
        CorrLAT_Dist_Spl->SetMarkerColor(2);
        CorrLAT_Dist_Spl->SetFillColor(2);
        CorrLAT_Dist_Spl->SetFillStyle(3001);
        CorrLAT_Dist_Spl->SetLineWidth(2);
        CorrLAT_Dist_Spl->Draw("Csame");


        TGraphErrors *CorrLATDistNaF;
        c16->cd(4);
        gPad->SetGridy();
        gPad->SetGridx();
        CorrLATDistNaF=new TGraphErrors();
        CorrLATDistNaF->SetTitle("Latitude Efficiency Corr.");
        CorrLATDistNaF->GetXaxis()->SetTitle("Latitude");
        CorrLATDistNaF->GetYaxis()->SetTitle("Eff. Corr. Factor");
        CorrLATDistNaF->GetYaxis()->SetRangeUser(0.96,1.04);
        CorrLATDistNaF->SetMarkerStyle(8);
        for(int i=1;i<11;i++) {
                        CorrLATDistNaF->SetPoint(i-1,geomagC[i],DistLATcorr_NaF->GetBinContent(i+1,1));
                        CorrLATDistNaF->SetPointError(i-1,0,DistLATcorr_NaF->GetBinError(i+1,1));
                }
        CorrLATDistNaF->Draw("AP");
        TGraphErrors *CorrLAT_DistNaF_Spl=new TGraphErrors("CorrLAT_DistNaF_Spl");
        CorrLAT_DistNaF_Spl->SetName("CorrLAT_DistNaF_Spl");
        for(int i=1;i<11;i++) { CorrLAT_DistNaF_Spl->SetPoint(i-1,geomagC[i],DistLATcorr_NaF_fit->GetBinContent(i+1));
                                CorrLAT_DistNaF_Spl->SetPointError(i-1,0,DistLATcorr_NaF_fit->GetBinError(i+1));
                                }
        CorrLAT_DistNaF_Spl->SetLineColor(2);
        CorrLAT_DistNaF_Spl->SetMarkerColor(2);
        CorrLAT_DistNaF_Spl->SetFillColor(2);
        CorrLAT_DistNaF_Spl->SetFillStyle(3001);
        CorrLAT_DistNaF_Spl->SetLineWidth(2);
        CorrLAT_DistNaF_Spl->Draw("Csame");

        TGraphErrors *CorrLATDistAgl;
        c16->cd(6);
        gPad->SetGridy();
        gPad->SetGridx();
        CorrLATDistAgl=new TGraphErrors();
        CorrLATDistAgl->SetTitle("Latitude Efficiency Corr.");
        CorrLATDistAgl->GetXaxis()->SetTitle("Latitude");
        CorrLATDistAgl->GetYaxis()->SetTitle("Eff. Corr. Factor");
        CorrLATDistAgl->GetYaxis()->SetRangeUser(0.96,1.04);
        CorrLATDistAgl->SetMarkerStyle(8);
        for(int i=1;i<11;i++) {
                        CorrLATDistAgl->SetPoint(i-1,geomagC[i],DistLATcorr_Agl->GetBinContent(i+1,1));
                        CorrLATDistAgl->SetPointError(i-1,0,DistLATcorr_Agl->GetBinError(i+1,1));
                }
        CorrLATDistAgl->Draw("AP");
        TGraphErrors *CorrLAT_DistAgl_Spl=new TGraphErrors("CorrLAT_DistAgl_Spl");
        CorrLAT_DistAgl_Spl->SetName("CorrLAT_DistAgl_Spl");
        for(int i=1;i<11;i++) { CorrLAT_DistAgl_Spl->SetPoint(i-1,geomagC[i],DistLATcorr_Agl_fit->GetBinContent(i+1));
                                CorrLAT_DistAgl_Spl->SetPointError(i-1,0,DistLATcorr_Agl_fit->GetBinError(i+1));
                                }
        CorrLAT_DistAgl_Spl->SetLineColor(2);
        CorrLAT_DistAgl_Spl->SetMarkerColor(2);
        CorrLAT_DistAgl_Spl->SetFillColor(2);
        CorrLAT_DistAgl_Spl->SetFillStyle(3001);
        CorrLAT_DistAgl_Spl->SetLineWidth(2);
        CorrLAT_DistAgl_Spl->Draw("Csame");

	cout<<"*** Updating Results file ***"<<endl;
        nomefile=percorso + "/CodesforAnalysis/Final_plots/"+mese+".root";
        TFile *f_out=new TFile(nomefile.c_str(), "UPDATE");
        f_out->mkdir("DATA-driven Results/Latitude effect/Quality");
	f_out->cd("DATA-driven Results/Latitude effect/Quality");
	c15->Write();
	c16->Write();
        f_out->Write();
        f_out->Close();

}


