using namespace std;

TCanvas * c31 = new TCanvas("Deutons Acceptance (R bins)");
TCanvas * c31_bis = new TCanvas("Deutons Acceptance (Beta bins)");
TCanvas * c31_tris = new TCanvas("Gen. Efficiency");

TH3F * AcceptDzoneTOF = new TH3F("AcceptDzoneTOF","AcceptDzoneTOF",18,0,18,11,0,11,6,0,6);
TH3F * AcceptDzoneNaF = new TH3F("AcceptDzoneNaF","AcceptDzoneNaF",18,0,18,11,0,11,6,0,6);
TH3F * AcceptDzoneAgl = new TH3F("AcceptDzoneAgl","AcceptDzoneAgl",18,0,18,11,0,11,6,0,6);

TH1F * AcceptPTOF = new TH1F("AcceptPTOF","AcceptPTOF",18,0,18);
TH1F * AcceptPNaF = new TH1F("AcceptPNaF","AcceptPNaF",18,0,18);
TH1F * AcceptPAgl = new TH1F("AcceptPAgl","AcceptPAgl",18,0,18);

TH2F * AcceptDTOF = new TH2F("AcceptDTOF","AcceptDTOF",18,0,18,6,0,6);
TH2F * AcceptDNaF = new TH2F("AcceptDNaF","AcceptDNaF",18,0,18,6,0,6);
TH2F * AcceptDAgl = new TH2F("AcceptDAgl","AcceptDAgl",18,0,18,6,0,6);

void AcceptanceD(TFile * file1){

	cout<<"****************** ACCEPTANCE CALCULATION ******************"<<endl;
	//Deutons
	float AcceptgeoD[43][6]={{0}};
	float AcceptSelMCD[43][6]={{0}};
	float AcceptpreMCD[43][6]={{0}};

	float AcceptgeoDTOF[18][6]={{0}};
        float AcceptSelMCDTOF[18][6]={{0}};
        float AcceptpreMCDTOF[18][6]={{0}};
	float AcceptSelMCDTOFzone[18][11][6]={{{0}}};
	
	float AcceptgeoDNaF[18][6]={{0}};
        float AcceptSelMCDNaF[18][6]={{0}};
        float AcceptpreMCDNaF[18][6]={{0}};
	float AcceptSelMCDNaFzone[18][11][6]={{{0}}};
	
	float AcceptgeoDAgl[18][6]={{0}};
        float AcceptSelMCDAgl[18][6]={{0}};
        float AcceptpreMCDAgl[18][6]={{0}};
	float AcceptSelMCDAglzone[18][11][6]={{{0}}};
	
	
	float eventiprova_D[6]={0};
	for(int h=0;h<6;h++) for(int i=1;i<43;i++) eventiprova_D[h]+=EffpreselMCD1_R->GetBinContent(i+1,h+1);
	float triggerbin_D[43][6];
	for(int h=0;h<6;h++) for(int i=1;i<43;i++)
		triggerbin_D[i][h]=(pow(0.0242236931,-1))*eventiprova_D[h]*(log(bin[i+1])-log(bin[i]))/(log(20)-log(0.5));
	
	float triggerbetabin_D[18][6]={{0}};
	for(int h=0;h<6;h++) for(int m=0;m<18;m++) 
	triggerbetabin_D[m][h]=(pow(0.0242236931,-1))*eventiprova_D[h]*(log(BetabinsR_D[m+1])-log(BetabinsR_D[m]))/(log(20)-log(0.5)); 

	float triggerbetabinNaF_D[18][6]={{0}};
        for(int h=0;h<6;h++) for(int m=0;m<18;m++)
        triggerbetabinNaF_D[m][h]=(pow(0.0242236931,-1))*eventiprova_D[h]*(log(BetabinsNaFR_D[m+1])-log(BetabinsNaFR_D[m]))/(log(20)-log(0.5));
	
	float triggerbetabinAgl_D[18][6]={{0}};
        for(int h=0;h<6;h++) for(int m=0;m<18;m++)
        triggerbetabinAgl_D[m][h]=(pow(0.0242236931,-1))*eventiprova_D[h]*(log(BetabinsAglR_D[m+1])-log(BetabinsAglR_D[m]))/(log(20)-log(0.5));

	for(int h=0;h<6;h++) for(int i=0;i<43;i++) 
	AcceptgeoD[i][h]=EffTriggerMCD_R_TH2F->GetBinContent(i+1,h+1)*EffTOF_MCD_R_TH2F->GetBinContent(i+1,h+1)*EffTriggMCD1_R->GetBinContent(i+1,h+1)/triggerbin_D[i][h]*47.78;
	for(int h=0;h<6;h++) for(int i=0;i<43;i++)
	AcceptpreMCD[i][h]=EffPreMCD_R_TH2F->GetBinContent(i+1,h+1)*EffpreselMCD1_R->GetBinContent(i+1,h+1)/triggerbin_D[i][h]*47.78;
	for(int h=0;h<6;h++) for(int i=0;i<43;i++)
        AcceptSelMCD[i][h]=AcceptpreMCD[i][h]*EffMCDistD_TH2F->GetBinContent(i+1,h+1);
	
	
	for(int h=0;h<6;h++) for(int m=0;m<18;m++)
        AcceptpreMCDTOF[m][h]=EffPreMCD_TH2F->GetBinContent(m+1,h+1)*EffpreselMCD1->GetBinContent(m+1,h+1)/triggerbetabin_D[m][h]*47.78;
	for(int h=0;h<6;h++) for(int m=0;m<18;m++)
	AcceptSelMCDTOF[m][h]=AcceptpreMCDTOF[m][h]*EffMCDistD_Beta_TH2F->GetBinContent(m+1,h+1);
	for(int h=0;h<6;h++) for(int i=0;i<11;i++) for(int m=0;m<18;m++){
		AcceptSelMCDTOFzone[m][i][h]=AcceptSelMCDTOF[m][h];
		//data-driven corrections	
		//latitude
		for(int S=0;S<3;S++) AcceptSelMCDTOFzone[m][i][h]/=CorrLAT_pre[S]->Eval(geomagC[i]);	
		AcceptSelMCDTOFzone[m][i][h]/=CorrLAT_Lik->Eval(geomagC[i])*CorrLAT_Dist->Eval(geomagC[i]);		
		//DvsMC
	}


	for(int h=0;h<6;h++) for(int m=0;m<18;m++)
        AcceptpreMCDNaF[m][h]=EffPreMCDNaF_TH2F->GetBinContent(m+1,h+1)*EffpreselMCD1NaF->GetBinContent(m+1,h+1)/triggerbetabinNaF_D[m][h]*47.78;
        for(int h=0;h<6;h++) for(int m=0;m<18;m++)
        AcceptSelMCDNaF[m][h]=AcceptpreMCDNaF[m][h]*EffMCDistD_BetaNaF_TH2F->GetBinContent(m+1,h+1);
	for(int h=0;h<6;h++) for(int i=0;i<11;i++) for(int m=0;m<18;m++){
		AcceptSelMCDNaFzone[m][i][h]=AcceptSelMCDNaF[m][h];
		//data-driven corrections	
		//latitude
		for(int S=0;S<3;S++) AcceptSelMCDNaFzone[m][i][h]/=CorrLAT_pre[S]->Eval(geomagC[i]);	
		AcceptSelMCDNaFzone[m][i][h]/=CorrLAT_LikNaF->Eval(geomagC[i])*CorrLAT_DistNaF->Eval(geomagC[i]);		
		//DvsMC
	}

	

	for(int h=0;h<6;h++) for(int m=0;m<18;m++)
        AcceptpreMCDAgl[m][h]=EffPreMCDAgl_TH2F->GetBinContent(m+1,h+1)*EffpreselMCD1Agl->GetBinContent(m+1,h+1)/triggerbetabinAgl_D[m][h]*47.78;
        for(int h=0;h<6;h++) for(int m=0;m<18;m++)
        AcceptSelMCDAgl[m][h]=AcceptpreMCDAgl[m][h]*EffMCDistD_BetaAgl_TH2F->GetBinContent(m+1,h+1);
	for(int h=0;h<6;h++) for(int i=0;i<11;i++) for(int m=0;m<18;m++){
		AcceptSelMCDAglzone[m][i][h]=AcceptSelMCDAgl[m][h];
		//data-driven corrections	
		//latitude
		for(int S=0;S<3;S++) AcceptSelMCDAglzone[m][i][h]/=CorrLAT_pre[S]->Eval(geomagC[i]);	
		AcceptSelMCDAglzone[m][i][h]/=CorrLAT_LikAgl->Eval(geomagC[i])*CorrLAT_DistAgl->Eval(geomagC[i]);		
		//DvsMC
	}

	//Protons
        float AcceptgeoPTOF[18]={0};
        float AcceptSelMCPTOF[18]={0};
        float AcceptpreMCPTOF[18]={0};

        float AcceptgeoPNaF[18]={0};
        float AcceptSelMCPNaF[18]={0};
        float AcceptpreMCPNaF[18]={0};

        float AcceptgeoPAgl[18]={0};
        float AcceptSelMCPAgl[18]={0};
        float AcceptpreMCPAgl[18]={0};
	
	float eventiprova_P=0;
        for(int i=0;i<43;i++) eventiprova_P+=EffpreselMCP1_R->GetBinContent(i+1);
        
	float triggerbetabin_P[18]={0};
	for(int m=0;m<18;m++)
        triggerbetabin_P[m]=(pow(0.0308232619,-1))*eventiprova_P*(log(BetabinsR_P[m+1])-log(BetabinsR_P[m]))/(log(100)-log(0.5));

	float triggerbetabinNaF_P[18]={0};
        for(int m=0;m<18;m++)
        triggerbetabinNaF_P[m]=(pow(0.0308232619,-1))*eventiprova_P*(log(BetabinsNaFR_P[m+1])-log(BetabinsNaFR_P[m]))/(log(100)-log(0.5));

	float triggerbetabinAgl_P[18]={0};
        for(int m=0;m<18;m++)
        triggerbetabinAgl_P[m]=(pow(0.0308232619,-1))*eventiprova_P*(log(BetabinsAglR_P[m+1])-log(BetabinsAglR_P[m]))/(log(100)-log(0.5));

	for(int m=0;m<18;m++)
        AcceptpreMCPTOF[m]=EffPreMCP_TH1F->GetBinContent(m+1)*EffpreselMCP1->GetBinContent(m+1)/triggerbetabin_P[m]*47.78;
	for(int m=0;m<18;m++)
	AcceptSelMCPTOF[m]=AcceptpreMCPTOF[m]*EffMCDistP_Beta_TH1F->GetBinContent(m+1);
	
	for(int m=0;m<18;m++)
        AcceptpreMCPNaF[m]=EffPreMCPNaF_TH1F->GetBinContent(m+1)*EffpreselMCP1NaF->GetBinContent(m+1)/triggerbetabinNaF_P[m]*47.78;
        for(int m=0;m<18;m++)
        AcceptSelMCPNaF[m]=AcceptpreMCPNaF[m]*EffMCDistP_BetaNaF_TH1F->GetBinContent(m+1);
	
	for(int m=0;m<18;m++)
        AcceptpreMCPAgl[m]=EffPreMCPAgl_TH1F->GetBinContent(m+1)*EffpreselMCP1Agl->GetBinContent(m+1)/triggerbetabinAgl_P[m]*47.78;
        for(int m=0;m<18;m++)
        AcceptSelMCPAgl[m]=AcceptpreMCPAgl[m]*EffMCDistP_BetaAgl_TH1F->GetBinContent(m+1);	

	
	string MCLegend[7]={"protons.B800","d.pl1.0_520_GG_Blic","d.pl1.0_520_GG_BlicDPMJet","d.pl1.0_520_GG_QMD","d.pl1.0_520_Shen_Blic","d.pl1.0_520_Shen_BlicDPMJet","d.pl1.0_520_Shen_QMD"};
	c31->cd();
	gPad->SetLogx();
        gPad->SetLogy();
	gPad->SetGridx();
        gPad->SetGridy();
	TGraphErrors * AccgeoD[6];
	TGraphErrors * AccPreMCD[6];
	TGraphErrors * AccSelMCD[6];
	for(int h=0;h<6;h++){
		AccgeoD[h]= new TGraphErrors();
		AccPreMCD[h]= new TGraphErrors();
		AccSelMCD[h]= new TGraphErrors();
		int p=0;
		for(int i=0;i<43;i++) if(AcceptgeoD[i][h]>0) {AccgeoD[h]->SetPoint(p,R_cent[i],AcceptgeoD[i][h]);p++;}	
		p=0;
		for(int i=0;i<43;i++) if(AcceptpreMCD[i][h]>0) {AccPreMCD[h]->SetPoint(p,R_cent[i],AcceptpreMCD[i][h]);p++;}
		p=0;
		for(int i=0;i<43;i++) if(AcceptSelMCD[i][h]>0) {AccSelMCD[h]->SetPoint(p,R_cent[i],AcceptSelMCD[i][h]);p++;}
		for(int m=0;m<18;m++) AcceptDTOF->SetBinContent(m+1,h+1,AcceptSelMCDTOF[m][h]);
		AccgeoD[h]->SetMarkerStyle(8);
		AccgeoD[h]->SetMarkerColor(4);
		AccgeoD[h]->SetMarkerSize(2);
		AccgeoD[h]->SetLineColor(4);
		AccgeoD[h]->SetLineWidth(1);
		AccgeoD[h]->SetMarkerStyle(h+3);
		AccPreMCD[h]->SetMarkerStyle(8);
		AccPreMCD[h]->SetMarkerColor(4);
		AccPreMCD[h]->SetMarkerSize(2);
		AccPreMCD[h]->SetLineColor(4);
		AccPreMCD[h]->SetLineWidth(1);
		AccPreMCD[h]->SetMarkerStyle(h+3);
		AccSelMCD[h]->SetMarkerStyle(8);
                AccSelMCD[h]->SetMarkerColor(4);
         	AccSelMCD[h]->SetMarkerSize(2);
                AccSelMCD[h]->SetLineColor(4);
                AccSelMCD[h]->SetLineWidth(1);
                AccSelMCD[h]->SetMarkerStyle(h+3);
		AccgeoD[h]->SetTitle("Final effective Acceptance");
		AccgeoD[h]->GetXaxis()->SetTitle("R [GV]");
		AccgeoD[h]->GetYaxis()->SetTitle("Acceptance [m^2 sr]");
		AccgeoD[h]->GetXaxis()->SetTitleSize(0.045);
		AccgeoD[h]->GetYaxis()->SetTitleSize(0.045);
		AccgeoD[h]->GetYaxis()->SetRangeUser(1e-2,1.3);
	}
	{		TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(AccgeoD[0],MCLegend[1].c_str(), "ep");
		AccgeoD[0]->Draw("ACP");
		AccSelMCD[0]->Draw("Csame");
		AccPreMCD[0]->Draw("Csame");
	for(int h=1;h<6;h++){
		AccgeoD[h]->Draw("CPsame");
		AccPreMCD[h]->Draw("CPsame");
		AccSelMCD[h]->Draw("CPsame");		
	}
	leg->Draw("same");	
	}
	
	c31_bis->cd();
        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetGridx();
        gPad->SetGridy();
        TGraphErrors * AccgeoDbeta[6];
        TGraphErrors * AccPreMCDbeta[6];
        TGraphErrors * AccSelMCDbeta[6];
        TGraphErrors * AccSelMCPbeta=new TGraphErrors();
	for(int h=0;h<6;h++){
                AccgeoDbeta[h]= new TGraphErrors();
                AccPreMCDbeta[h]= new TGraphErrors();
                AccSelMCDbeta[h]= new TGraphErrors();
                int p=0;
                /*for(int i=0;i<43;i++) if(AcceptgeoD[i][h]>0) {AccgeoDbeta[h]->SetPoint(p,R_cent[i],AcceptgeoD[i][h]);p++;}
                p=0;
                for(int i=0;i<43;i++) if(AcceptpreMCD[i][h]>0) {AccPreMCDbeta[h]->SetPoint(p,R_cent[i],AcceptpreMCD[i][h]);p++;}
                p=0;*/
                for(int m=0;m<18;m++) if(AcceptSelMCDTOF[m][h]>0) {AccSelMCDbeta[h]->SetPoint(p,Ekincent[m],AcceptSelMCDTOFzone[m][1][h]);p++;}
		for(int m=0;m<18;m++) for(int i=0;i<11;i++)AcceptDzoneTOF->SetBinContent(m+1,i+1,h,AcceptSelMCDTOFzone[m][i][h]);
		if(h==0) 
		{
			p=0;
			for(int m=0;m<18;m++) if(AcceptSelMCPTOF[m]>0) {AccSelMCPbeta->SetPoint(p,Ekincent[m],AcceptSelMCPTOF[m]);p++;}
			for(int m=0;m<18;m++) AcceptPTOF->SetBinContent(m+1,AcceptSelMCPTOF[m]);
			AccSelMCDbeta[0]->SetPoint(p,50,0.001);
			AccSelMCPbeta->SetMarkerStyle(8);
                	AccSelMCPbeta->SetMarkerColor(2);
                	AccSelMCPbeta->SetMarkerSize(2);
                	AccSelMCPbeta->SetLineColor(2);
                	AccSelMCPbeta->SetLineWidth(1);
		}
		AccSelMCDbeta[h]->SetMarkerStyle(8);
                AccSelMCDbeta[h]->SetMarkerColor(4);
                AccSelMCDbeta[h]->SetMarkerSize(2);
                AccSelMCDbeta[h]->SetLineColor(4);
		AccSelMCDbeta[h]->SetLineStyle(2);
                AccSelMCDbeta[h]->SetLineWidth(1);
                AccSelMCDbeta[h]->SetMarkerStyle(h+3);
                AccSelMCDbeta[0]->SetTitle("Final effective Acceptance");
                AccSelMCDbeta[0]->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
                AccSelMCDbeta[0]->GetYaxis()->SetTitle("Acceptance [m^2 sr]");
                AccSelMCDbeta[0]->GetXaxis()->SetTitleSize(0.045);
                AccSelMCDbeta[0]->GetYaxis()->SetTitleSize(0.045);
                AccSelMCDbeta[0]->GetYaxis()->SetRangeUser(1e-2,1.3);
        }
	{	TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(AccSelMCDbeta[0],MCLegend[1].c_str(), "ep");
		AccSelMCDbeta[0]->Draw("AP");
		AccSelMCPbeta->SetLineStyle(2);
		AccSelMCPbeta->Draw("CPsame");
        for(int h=1;h<6;h++){
                leg->AddEntry(AccSelMCDbeta[h],MCLegend[h+1].c_str(), "ep");
		AccSelMCDbeta[h]->Draw("CPsame");
        }
        	leg->AddEntry(AccSelMCPbeta,"Protons B800","ep");		
		leg->Draw("same");
	}

        TGraphErrors * AccgeoDbetaNaF[6];
        TGraphErrors * AccPreMCDbetaNaF[6];
        TGraphErrors * AccSelMCDbetaNaF[6];
        TGraphErrors * AccSelMCPbetaNaF=new TGraphErrors();
	for(int h=0;h<6;h++){
                AccgeoDbetaNaF[h]= new TGraphErrors();
                AccPreMCDbetaNaF[h]= new TGraphErrors();
                AccSelMCDbetaNaF[h]= new TGraphErrors();
                int p=0;
                /*for(int i=0;i<43;i++) if(AcceptgeoD[i][h]>0) {AccgeoDbeta[h]->SetPoint(p,R_cent[i],AcceptgeoD[i][h]);p++;}
                p=0;
                for(int i=0;i<43;i++) if(AcceptpreMCD[i][h]>0) {AccPreMCDbeta[h]->SetPoint(p,R_cent[i],AcceptpreMCD[i][h]);p++;}
                p=0;*/
                for(int m=0;m<18;m++) if(AcceptSelMCDNaF[m][h]>0) {AccSelMCDbetaNaF[h]->SetPoint(p,EkincentNaF[m],AcceptSelMCDNaFzone[m][1][h]);p++;}
		for(int m=0;m<18;m++) for(int i=0;i<11;i++)AcceptDzoneNaF->SetBinContent(m+1,i+1,h,AcceptSelMCDNaFzone[m][i][h]);
                for(int m=0;m<18;m++) AcceptDNaF->SetBinContent(m+1,h+1,AcceptSelMCDNaF[m][h]);
		if(h==0)
                {
                        p=0;
			for(int m=0;m<18;m++) if(AcceptSelMCPNaF[m]>0) {AccSelMCPbetaNaF->SetPoint(p,EkincentNaF[m],AcceptSelMCPNaF[m]);p++;}
                        for(int m=0;m<18;m++) AcceptPNaF->SetBinContent(m+1,AcceptSelMCPNaF[m]);
                        AccSelMCPbetaNaF->SetMarkerStyle(8);
                        AccSelMCPbetaNaF->SetMarkerColor(2);
                        AccSelMCPbetaNaF->SetMarkerSize(2);
                        AccSelMCPbetaNaF->SetLineColor(2);
                        AccSelMCPbetaNaF->SetLineWidth(1);
                }
		/*AccgeoDbeta[h]->SetMarkerStyle(8);
                AccgeoDbeta[h]->SetMarkerColor(4);
                AccgeoDbeta[h]->SetMarkerSize(2);
                AccgeoDbeta[h]->SetLineColor(4);
                AccgeoDbeta[h]->SetLineWidth(1);
                AccgeoDbeta[h]->SetMarkerStyle(h+3);
                AccPreMCDbeta[h]->SetMarkerStyle(8);
                AccPreMCDbeta[h]->SetMarkerColor(4);
                AccPreMCDbeta[h]->SetMarkerSize(2);
                AccPreMCDbeta[h]->SetLineColor(4);
                AccPreMCDbeta[h]->SetLineWidth(1);
                AccPreMCDbeta[h]->SetMarkerStyle(h+3);*/
                AccSelMCDbetaNaF[h]->SetMarkerStyle(8);
                AccSelMCDbetaNaF[h]->SetMarkerColor(4);
                AccSelMCDbetaNaF[h]->SetMarkerSize(2);
                AccSelMCDbetaNaF[h]->SetLineColor(4);
                AccSelMCDbetaNaF[h]->SetLineStyle(2);
		AccSelMCDbetaNaF[h]->SetLineWidth(1);
                AccSelMCDbetaNaF[h]->SetMarkerStyle(h+3);
                AccSelMCDbetaNaF[h]->SetTitle("Deutons Acceptance");
                AccSelMCDbetaNaF[h]->GetXaxis()->SetTitle("R [GV]");
                AccSelMCDbetaNaF[h]->GetYaxis()->SetTitle("Acceptance [m^2 sr]");
                AccSelMCDbetaNaF[h]->GetXaxis()->SetTitleSize(0.045);
                AccSelMCDbetaNaF[h]->GetYaxis()->SetTitleSize(0.045);
                AccSelMCDbetaNaF[h]->GetYaxis()->SetRangeUser(1e-2,1.3);
        }
                AccSelMCDbetaNaF[0]->Draw("CPsame");
		AccSelMCPbetaNaF->SetLineStyle(2);
		AccSelMCPbetaNaF->Draw("CPsame");
        for(int h=1;h<6;h++){
                AccSelMCDbetaNaF[h]->Draw("CPsame");
        }
	
	TGraphErrors * AccgeoDbetaAgl[6];
        TGraphErrors * AccPreMCDbetaAgl[6];
        TGraphErrors * AccSelMCDbetaAgl[6];
	TGraphErrors * AccSelMCPbetaAgl=new TGraphErrors();
	for(int h=0;h<6;h++){
                AccgeoDbetaAgl[h]= new TGraphErrors();
                AccPreMCDbetaAgl[h]= new TGraphErrors();
                AccSelMCDbetaAgl[h]= new TGraphErrors();
                int p=0;
                /*for(int i=0;i<43;i++) if(AcceptgeoD[i][h]>0) {AccgeoDbeta[h]->SetPoint(p,R_cent[i],AcceptgeoD[i][h]);p++;}
                p=0;
                for(int i=0;i<43;i++) if(AcceptpreMCD[i][h]>0) {AccPreMCDbeta[h]->SetPoint(p,R_cent[i],AcceptpreMCD[i][h]);p++;}
                p=0;*/
                for(int m=0;m<18;m++) if(AcceptSelMCDAgl[m][h]>0) {AccSelMCDbetaAgl[h]->SetPoint(p,EkincentAgl[m],AcceptSelMCDAglzone[m][1][h]);p++;}
		for(int m=0;m<18;m++) for(int i=0;i<11;i++)AcceptDzoneAgl->SetBinContent(m+1,i+1,h,AcceptSelMCDAglzone[m][i][h]);
                for(int m=0;m<18;m++) AcceptDAgl->SetBinContent(m+1,h+1,AcceptSelMCDAgl[m][h]); 
		if(h==0)
                {
                        p=0;
			for(int m=0;m<18;m++) if(AcceptSelMCPAgl[m]>0) {AccSelMCPbetaAgl->SetPoint(p,EkincentAgl[m],AcceptSelMCPAgl[m]);p++;}
                        for(int m=0;m<18;m++) AcceptPAgl->SetBinContent(m+1,AcceptSelMCPAgl[m]);
                        AccSelMCPbetaAgl->SetMarkerStyle(8);
                        AccSelMCPbetaAgl->SetMarkerColor(2);
                        AccSelMCPbetaAgl->SetMarkerSize(2);
                        AccSelMCPbetaAgl->SetLineColor(2);
                        AccSelMCPbetaAgl->SetLineWidth(1);
                }

		/*AccgeoDbeta[h]->SetMarkerStyle(8);
                AccgeoDbeta[h]->SetMarkerColor(4);
                AccgeoDbeta[h]->SetMarkerSize(2);
                AccgeoDbeta[h]->SetLineColor(4);
                AccgeoDbeta[h]->SetLineWidth(1);
                AccgeoDbeta[h]->SetMarkerStyle(h+3);
                AccPreMCDbeta[h]->SetMarkerStyle(8);
                AccPreMCDbeta[h]->SetMarkerColor(4);
                AccPreMCDbeta[h]->SetMarkerSize(2);
                AccPreMCDbeta[h]->SetLineColor(4);
                AccPreMCDbeta[h]->SetLineWidth(1);
                AccPreMCDbeta[h]->SetMarkerStyle(h+3);*/
                AccSelMCDbetaAgl[h]->SetMarkerStyle(8);
                AccSelMCDbetaAgl[h]->SetMarkerColor(4);
                AccSelMCDbetaAgl[h]->SetMarkerSize(2);
                AccSelMCDbetaAgl[h]->SetLineColor(4);
                AccSelMCDbetaAgl[h]->SetLineWidth(1);
                AccSelMCDbetaAgl[h]->SetLineStyle(2);
		AccSelMCDbetaAgl[h]->SetMarkerStyle(h+3);
                AccSelMCDbetaAgl[h]->SetTitle("Deutons Acceptance");
                AccSelMCDbetaAgl[h]->GetXaxis()->SetTitle("R [GV]");
                AccSelMCDbetaAgl[h]->GetYaxis()->SetTitle("Acceptance [m^2 sr]");
                AccSelMCDbetaAgl[h]->GetXaxis()->SetTitleSize(0.045);
                AccSelMCDbetaAgl[h]->GetYaxis()->SetTitleSize(0.045);
                AccSelMCDbetaAgl[h]->GetYaxis()->SetRangeUser(1e-2,1.3);
        }
                AccSelMCDbetaAgl[0]->Draw("CPsame");
		AccSelMCPbetaAgl->SetLineStyle(2);
        	AccSelMCPbetaAgl->Draw("CPsame");
	for(int h=1;h<6;h++){
                AccSelMCDbetaAgl[h]->Draw("CPsame");
        }
	c31_tris->Divide(2,1);
	c31_tris->cd(1);
	gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetGridx();
        gPad->SetGridy();
        TGraphErrors * EffgenbetaD[6];
	int p =0;
	for(int h=0;h<6;h++){
		EffgenbetaD[h]=new TGraphErrors();
		p=0;
		for(int i=0;i<43;i++) {EffgenbetaD[h]->SetPoint(p,encindeut[i],EffpreselMCD1_R->GetBinContent(i+1,h+1)/triggerbin_D[i][h]);p++;}			
	
	EffgenbetaD[h]->SetMarkerStyle(8);
        EffgenbetaD[h]->SetMarkerColor(4);
        EffgenbetaD[h]->SetMarkerSize(1);
        EffgenbetaD[h]->SetLineColor(4);
        EffgenbetaD[h]->SetLineWidth(1);
        EffgenbetaD[h]->SetMarkerStyle(h+3);
        EffgenbetaD[h]->SetTitle("");
        EffgenbetaD[h]->GetXaxis()->SetTitle("R [GV]");
        EffgenbetaD[h]->GetYaxis()->SetTitle("Gen. Eff.");
        EffgenbetaD[h]->GetXaxis()->SetTitleSize(0.045);
        EffgenbetaD[h]->GetYaxis()->SetTitleSize(0.045);
	}
	EffgenbetaD[0]->Draw("ACP");
        for(int h=1;h<6;h++){
                EffgenbetaD[h]->Draw("CPsame");
        }
	TGraphErrors * EffgenbetaDTOF[6];
        p=0;
        for(int h=0;h<6;h++){
                EffgenbetaDTOF[h]=new TGraphErrors();
                p=0;
                for(int i=0;i<18;i++) {EffgenbetaDTOF[h]->SetPoint(p,Ekincent[i],EffpreselMCD1->GetBinContent(i+1,h+1)/triggerbetabin_D[i][h]);p++;}
        EffgenbetaDTOF[h]->SetMarkerStyle(8);
        EffgenbetaDTOF[h]->SetMarkerColor(4);
        EffgenbetaDTOF[h]->SetMarkerSize(2);
        EffgenbetaDTOF[h]->SetLineColor(4);
        EffgenbetaDTOF[h]->SetLineWidth(1);
        EffgenbetaDTOF[h]->SetMarkerStyle(h+3);
        EffgenbetaDTOF[h]->SetTitle("");
        EffgenbetaDTOF[h]->GetXaxis()->SetTitle("R [GV]");
        EffgenbetaDTOF[h]->GetYaxis()->SetTitle("Gen. Eff.");
        EffgenbetaDTOF[h]->GetXaxis()->SetTitleSize(0.045);
        EffgenbetaDTOF[h]->GetYaxis()->SetTitleSize(0.045);
        }
        EffgenbetaDTOF[0]->Draw("CPsame");
        for(int h=1;h<6;h++){
                EffgenbetaDTOF[h]->Draw("CPsame");
        }
	
	TGraphErrors * EffgenbetaDNaF[6];
        p=0;
        for(int h=0;h<6;h++){
                EffgenbetaDNaF[h]=new TGraphErrors();
                p=0;
                for(int i=0;i<18;i++) {EffgenbetaDNaF[h]->SetPoint(p,EkincentNaF[i],EffpreselMCD1NaF->GetBinContent(i+1,h+1)/triggerbetabinNaF_D[i][h]);p++;}

        EffgenbetaDNaF[h]->SetMarkerStyle(8);
        EffgenbetaDNaF[h]->SetMarkerColor(4);
        EffgenbetaDNaF[h]->SetMarkerSize(2);
        EffgenbetaDNaF[h]->SetLineColor(4);
        EffgenbetaDNaF[h]->SetLineWidth(1);
        EffgenbetaDNaF[h]->SetMarkerStyle(h+3);
        EffgenbetaDNaF[h]->SetTitle("");
        EffgenbetaDNaF[h]->GetXaxis()->SetTitle("R [GV]");
        EffgenbetaDNaF[h]->GetYaxis()->SetTitle("Gen. Eff.");
        EffgenbetaDNaF[h]->GetXaxis()->SetTitleSize(0.045);
        EffgenbetaDNaF[h]->GetYaxis()->SetTitleSize(0.045);
        }
        EffgenbetaDNaF[0]->Draw("CPsame");
        for(int h=1;h<6;h++){
                EffgenbetaDNaF[h]->Draw("CPsame");
        }
	TGraphErrors * EffgenbetaDAgl[6];
        p=0;
        for(int h=0;h<6;h++){
                EffgenbetaDAgl[h]=new TGraphErrors();
                p=0;
                for(int i=0;i<18;i++) {EffgenbetaDAgl[h]->SetPoint(p,EkincentAgl[i],EffpreselMCD1Agl->GetBinContent(i+1,h+1)/triggerbetabinAgl_D[i][h]);p++;}

        EffgenbetaDAgl[h]->SetMarkerStyle(8);
        EffgenbetaDAgl[h]->SetMarkerColor(4);
        EffgenbetaDAgl[h]->SetMarkerSize(2);
        EffgenbetaDAgl[h]->SetLineColor(4);
        EffgenbetaDAgl[h]->SetLineWidth(1);
        EffgenbetaDAgl[h]->SetMarkerStyle(h+3);
        EffgenbetaDAgl[h]->SetTitle("");
        EffgenbetaDAgl[h]->GetXaxis()->SetTitle("R [GV]");
        EffgenbetaDAgl[h]->GetYaxis()->SetTitle("Gen. Eff.");
        EffgenbetaDAgl[h]->GetXaxis()->SetTitleSize(0.045);
        EffgenbetaDAgl[h]->GetYaxis()->SetTitleSize(0.045);
        }
        EffgenbetaDAgl[0]->Draw("CPsame");
        for(int h=1;h<6;h++){
                EffgenbetaDAgl[h]->Draw("CPsame");
        }
	c31_tris->cd(2);
        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetGridx();
        gPad->SetGridy();
        float eventiprova=0;
        for(int i=0;i<43;i++) eventiprova+=EffpreselMCP1_R->GetBinContent(i+1);
        float triggerbin=(pow(0.0308232619,-1))*eventiprova/43;
	//float triggerbin[43];
	//for(int i=0;i<43;i++) triggerbin[i]=(pow(0.0308232619,-1))*eventiprova_P*(log(bin[i+1])-log(bin[i]))/(log(100)-log(0.5));
	TGraphErrors * EffgenbetaP;
	p =0;
	EffgenbetaP=new TGraphErrors();
	p=0;
	for(int i=0;i<43;i++) {EffgenbetaP->SetPoint(p,encinprot[i],EffpreselMCP1_R->GetBinContent(i+1)/triggerbin);p++;}

	EffgenbetaP->SetMarkerStyle(8);
	EffgenbetaP->SetMarkerColor(2);
	EffgenbetaP->SetMarkerSize(1);
	EffgenbetaP->SetLineColor(2);
	EffgenbetaP->SetLineWidth(1);
	EffgenbetaP->SetTitle("");
	EffgenbetaP->GetXaxis()->SetTitle("R [GV]");
	EffgenbetaP->GetYaxis()->SetTitle("Gen. Eff.");
	EffgenbetaP->GetXaxis()->SetTitleSize(0.045);
	EffgenbetaP->GetYaxis()->SetTitleSize(0.045);
	EffgenbetaP->Draw("ACP");
	
	TGraphErrors * EffgenbetaPTOF;
        p=0;
	EffgenbetaPTOF=new TGraphErrors();
	p=0;
	for(int i=0;i<18;i++) {EffgenbetaPTOF->SetPoint(p,Ekincent[i],EffpreselMCP1->GetBinContent(i+1)/triggerbetabin_P[i]);p++;}
	EffgenbetaPTOF->SetMarkerStyle(8);
	EffgenbetaPTOF->SetMarkerColor(2);
	EffgenbetaPTOF->SetMarkerSize(2);
	EffgenbetaPTOF->SetLineColor(2);
	EffgenbetaPTOF->SetLineWidth(1);
	EffgenbetaPTOF->SetTitle("");
	EffgenbetaPTOF->GetXaxis()->SetTitle("R [GV]");
	EffgenbetaPTOF->GetYaxis()->SetTitle("Gen. Eff.");
	EffgenbetaPTOF->GetXaxis()->SetTitleSize(0.045);
	EffgenbetaPTOF->GetYaxis()->SetTitleSize(0.045);
	EffgenbetaPTOF->Draw("CPsame");
		
	TGraphErrors * EffgenbetaPNaF;
        p=0;
        EffgenbetaPNaF=new TGraphErrors();
        p=0;
        for(int i=0;i<18;i++) {EffgenbetaPNaF->SetPoint(p,EkincentNaF[i],EffpreselMCP1NaF->GetBinContent(i+1)/triggerbetabinNaF_P[i]);p++;}
        EffgenbetaPNaF->SetMarkerStyle(8);
        EffgenbetaPNaF->SetMarkerColor(2);
        EffgenbetaPNaF->SetMarkerSize(2);
        EffgenbetaPNaF->SetLineColor(2);
        EffgenbetaPNaF->SetLineWidth(1);
        EffgenbetaPNaF->SetTitle("");
        EffgenbetaPNaF->GetXaxis()->SetTitle("R [GV]");
        EffgenbetaPNaF->GetYaxis()->SetTitle("Gen. Eff.");
        EffgenbetaPNaF->GetXaxis()->SetTitleSize(0.045);
        EffgenbetaPNaF->GetYaxis()->SetTitleSize(0.045);
        EffgenbetaPNaF->Draw("CPsame");
	
	TGraphErrors * EffgenbetaPAgl;
        p=0;
        EffgenbetaPAgl=new TGraphErrors();
        p=0;
        for(int i=0;i<18;i++) {EffgenbetaPAgl->SetPoint(p,EkincentAgl[i],EffpreselMCP1Agl->GetBinContent(i+1)/triggerbetabinAgl_P[i]);p++;}
        EffgenbetaPAgl->SetMarkerStyle(8);
        EffgenbetaPAgl->SetMarkerColor(2);
        EffgenbetaPAgl->SetMarkerSize(2);
        EffgenbetaPAgl->SetLineColor(2);
        EffgenbetaPAgl->SetLineWidth(1);
        EffgenbetaPAgl->SetTitle("");
        EffgenbetaPAgl->GetXaxis()->SetTitle("R [GV]");
        EffgenbetaPAgl->GetYaxis()->SetTitle("Gen. Eff.");
        EffgenbetaPAgl->GetXaxis()->SetTitleSize(0.045);
        EffgenbetaPAgl->GetYaxis()->SetTitleSize(0.045);
        EffgenbetaPAgl->Draw("CPsame");

}

