using namespace std;

TCanvas *c4=new TCanvas("Preselections Efficiency (R bins)");
TCanvas *c4_bis=new TCanvas("Preselections Efficiency (Beta bins)");

TH1F * EffpreselMCP1=new TH1F("EffpreselMCP1","EffpreselMCP1",18,0,18);
TH1F * EffpreselMCP2=new TH1F("EffpreselMCP2","EffpreselMCP2",18,0,18);
TH1F * EffpreselMCP1NaF=new TH1F("EffpreselMCP1NaF","EffpreselMCP1NaF",18,0,18);
TH1F * EffpreselMCP2NaF=new TH1F("EffpreselMCP2NaF","EffpreselMCP2NaF",18,0,18);
TH1F * EffpreselMCP1Agl=new TH1F("EffpreselMCP1Agl","EffpreselMCP1Agl",18,0,18);
TH1F * EffpreselMCP2Agl=new TH1F("EffpreselMCP2Agl","EffpreselMCP2Agl",18,0,18);

TH1F * EffpreselMCP1_R=new TH1F("EffpreselMCP1_R","EffpreselMCP1_R",43,0,43);
TH1F * EffpreselMCP2_R=new TH1F("EffpreselMCP2_R","EffpreselMCP2_R",43,0,43);
TH2F * EffpreselMCD1=new TH2F("EffpreselMCD1","EffpreselMCD1",18,0,18,6,0,6);
TH2F * EffpreselMCD2=new TH2F("EffpreselMCD2","EffpreselMCD2",18,0,18,6,0,6);
TH2F * EffpreselMCD1NaF=new TH2F("EffpreselMCD1NaF","EffpreselMCD1NaF",18,0,18,6,0,6);
TH2F * EffpreselMCD2NaF=new TH2F("EffpreselMCD2NaF","EffpreselMCD2NaF",18,0,18,6,0,6);
TH2F * EffpreselMCD1Agl=new TH2F("EffpreselMCD1Agl","EffpreselMCD1Agl",18,0,18,6,0,6);
TH2F * EffpreselMCD2Agl=new TH2F("EffpreselMCD2Agl","EffpreselMCD2Agl",18,0,18,6,0,6);
TH2F * EffpreselMCD1_R=new TH2F("EffpreselMCD1_R","EffpreselMCD1_R",43,0,43,6,0,6);
TH2F * EffpreselMCD2_R=new TH2F("EffpreselMCD2_R","EffpreselMCD2_R",43,0,43,6,0,6);

TH1F *EffPreMCP_R_TH1F = new TH1F("EffPreMCP_R_TH1F","EffPreMCP_R_TH1F",43,0,43);
TH1F *EffPreMCP_TH1F = new TH1F("EffPreMCP_TH1F","EffPreMCP_TH1F",18,0,18);
TH1F *EffPreMCPNaF_TH1F = new TH1F("EffPreMCPNaF_TH1F","EffPreMCPNaF_TH1F",18,0,18);
TH1F *EffPreMCPAgl_TH1F = new TH1F("EffPreMCPAgl_TH1F","EffPreMCPAgl_TH1F",18,0,18);
TH2F *EffPreMCD_R_TH2F = new TH2F("EffPreMCD_R_TH2F","EffPreMCD_R_TH2F",43,0,43,6,0,6);
TH2F *EffPreMCD_TH2F = new TH2F("EffPreMCD_TH2F","EffPreMCD_TH2F",18,0,18,6,0,6);
TH2F *EffPreMCDNaF_TH2F = new TH2F("EffPreMCDNaF_TH2F","EffPreMCDNaF_TH2F",18,0,18,6,0,6);
TH2F *EffPreMCDAgl_TH2F = new TH2F("EffPreMCDAgl_TH2F","EffPreMCDAgl_TH2F",18,0,18,6,0,6);



void MCpreseff_Fill(TNtuple *ntupla, int l){
	int k = ntupla->GetEvent(l);
	
	if(Massa_gen<1&&Massa_gen>0.5) {
		for(int M=0;M<43;M++) if(fabs(Momento_gen)<bin[M+1]&&fabs(Momento_gen)>bin[M]) 
			EffpreselMCP1_R->Fill(M);
		if(Unbias==0&&((int)Cutmask&187)==187&&Beta_pre>0&&R_pre>0) 
			for(int M=0;M<43;M++) if(fabs(R_pre)<bin[M+1]&&fabs(R_pre)>bin[M]) 	
				EffpreselMCP2_R->Fill(M);

		for(int m=0;m<18;m++)  if(Var3>BetaP[m]&&Var3<=BetaP[m+1]) 
			EffpreselMCP1->Fill(m);
		if(Unbias==0&&((int)Cutmask&187)==187&&Beta_pre>0&&R_pre>0) 
			for(int m=0;m<18;m++)  if(Var>BetaP[m]&&Var<=BetaP[m+1])		
				EffpreselMCP2->Fill(m);	

		for(int m=0;m<18;m++)  if(Var3>BetaNaFP[m]&&Var3<=BetaNaFP[m+1])
			EffpreselMCP1NaF->Fill(m);
		if(Unbias==0&&((int)Cutmask&187)==187&&Beta_pre>0&&R_pre>0)
			for(int m=0;m<18;m++)  if((((int)Cutmask)>>11)==512&&Var2>BetaNaFP[m]&&Var2<=BetaNaFP[m+1])
				EffpreselMCP2NaF->Fill(m);	

		for(int m=0;m<18;m++)  if(Var3>BetaAglP[m]&&Var3<=BetaAglP[m+1])
			EffpreselMCP1Agl->Fill(m);
		if(Unbias==0&&((int)Cutmask&187)==187&&Beta_pre>0&&R_pre>0)
			for(int m=0;m<18;m++)  if((((int)Cutmask)>>11)==0&&Var2>BetaAglP[m]&&Var2<=BetaAglP[m+1])
				EffpreselMCP2Agl->Fill(m);

	}				 

	if(Massa_gen>1&&Massa_gen<2) {
		for(int M=0;M<43;M++) if(fabs(Momento_gen)<bin[M+1]&&fabs(Momento_gen)>bin[M]) 
			EffpreselMCD1_R->Fill(M,(int)(10000*Massa_gen-18570));
		if(((int)Cutmask&187)==187&&Beta_pre>0&&Unbias==0&&R_pre>0) 
			for(int M=0;M<43;M++) if(fabs(R_pre)<bin[M+1]&&fabs(R_pre)>bin[M])	
				EffpreselMCD2_R->Fill(M,(int)(10000*Massa_gen-18570));

		for(int m=0;m<18;m++) if(Var3>BetaD[m]&&Var3<=BetaD[m+1])
			EffpreselMCD1->Fill(m,(int)(10000*Massa_gen-18570));
		for(int m=0;m<18;m++)  if(Var>BetaD[m]&&Var<=BetaD[m+1])	
			if(((int)Cutmask&187)==187&&Beta_pre>0&&Unbias==0&&R_pre>0) 
				EffpreselMCD2->Fill(m,(int)(10000*Massa_gen-18570));

		for(int m=0;m<18;m++) if(Var3>BetaNaFD[m]&&Var3<=BetaNaFD[m+1])
			EffpreselMCD1NaF->Fill(m,(int)(10000*Massa_gen-18570));
		for(int m=0;m<18;m++)  if((((int)Cutmask)>>11)==512&&Var2>BetaNaFD[m]&&Var2<=BetaNaFD[m+1])
			if(((int)Cutmask&187)==187&&Beta_pre>0&&Unbias==0&&R_pre>0)
				EffpreselMCD2NaF->Fill(m,(int)(10000*Massa_gen-18570));

		for(int m=0;m<18;m++) if(Var3>BetaAglD[m]&&Var3<=BetaAglD[m+1])
			EffpreselMCD1Agl->Fill(m,(int)(10000*Massa_gen-18570));
		for(int m=0;m<18;m++)  if((((int)Cutmask)>>11)==0&&Var2>BetaAglD[m]&&Var2<=BetaAglD[m+1])
			if(((int)Cutmask&187)==187&&Beta_pre>0&&Unbias==0&&R_pre>0)
				EffpreselMCD2Agl->Fill(m,(int)(10000*Massa_gen-18570));

		
	}

	return;
}

void MCpreeff_Copy(TFile * file){
	EffpreselMCP1= (TH1F*) file->Get("EffpreselMCP1");
	EffpreselMCP1NaF= (TH1F*) file->Get("EffpreselMCP1NaF");
	EffpreselMCP1Agl= (TH1F*) file->Get("EffpreselMCP1Agl");
	EffpreselMCD1= (TH2F*) file->Get("EffpreselMCD1");
	EffpreselMCD1NaF= (TH2F*) file->Get("EffpreselMCD1NaF");
	EffpreselMCD1Agl= (TH2F*) file->Get("EffpreselMCD1Agl");
	EffpreselMCP2= (TH1F*) file->Get("EffpreselMCP2");
	EffpreselMCP2NaF= (TH1F*) file->Get("EffpreselMCP2NaF");
	EffpreselMCP2Agl= (TH1F*) file->Get("EffpreselMCP2Agl");
	EffpreselMCD2= (TH2F*) file->Get("EffpreselMCD2");
	EffpreselMCD2NaF= (TH2F*) file->Get("EffpreselMCD2NaF");
	EffpreselMCD2Agl= (TH2F*) file->Get("EffpreselMCD2Agl");
	EffpreselMCP1_R =(TH1F*) file->Get("EffpreselMCP1_R");
	EffpreselMCD1_R =(TH2F*) file->Get("EffpreselMCD1_R");
	EffpreselMCP2_R =(TH1F*) file->Get("EffpreselMCP2_R");
	EffpreselMCD2_R =(TH2F*) file->Get("EffpreselMCD2_R");

	return;	
}

void MCpreeff_Write(){
        EffpreselMCP1->Write();
        EffpreselMCP1NaF->Write();
        EffpreselMCP1Agl->Write();
        EffpreselMCD1->Write();
        EffpreselMCD1NaF->Write();
        EffpreselMCD1Agl->Write();
        EffpreselMCP2->Write();
        EffpreselMCP2NaF->Write();
        EffpreselMCP2Agl->Write();
        EffpreselMCD2->Write();
        EffpreselMCD2NaF->Write();
        EffpreselMCD2Agl->Write();
        EffpreselMCP1_R->Write();
        EffpreselMCD1_R->Write();
        EffpreselMCP2_R->Write();
        EffpreselMCD2_R->Write();
        return;
}


void MCpreeff(TFile * file1){

	string numero[18]={"0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","18"};
	string tagli[10]={"Trigger","3of4 TOF","TRD Segments","Rigidity exists","Chi^2 R","Matching TOF","Matching TRD","In TRD Accept.","1 Particle","1 Tr. Track"};
	string nome;

	cout<<"**** MC PRESELECTIONS EFFICIENCY (FULL SET) ****"<<endl;
	c4_bis->Divide(3,1);
	float EffpreselMCP[18]={0};
	for(int i=0;i<18;i++) if(EffpreselMCP1->GetBinContent(i+1)>0) if(EffpreselMCP2->GetBinContent(i+1)<EffpreselMCP1->GetBinContent(i+1))
		EffpreselMCP[i]=EffpreselMCP2->GetBinContent(i+1)/(float)EffpreselMCP1->GetBinContent(i+1);
	
	float EffpreselMCPNaF[18]={0};
        for(int i=0;i<18;i++) if(EffpreselMCP1NaF->GetBinContent(i+1)>0) if(EffpreselMCP2NaF->GetBinContent(i+1)<EffpreselMCP1NaF->GetBinContent(i+1))
                EffpreselMCPNaF[i]=EffpreselMCP2NaF->GetBinContent(i+1)/(float)EffpreselMCP1NaF->GetBinContent(i+1);
	
	float EffpreselMCPAgl[18]={0};
        for(int i=0;i<18;i++) if(EffpreselMCP1Agl->GetBinContent(i+1)>0) if(EffpreselMCP2Agl->GetBinContent(i+1)<EffpreselMCP1Agl->GetBinContent(i+1))
                EffpreselMCPAgl[i]=EffpreselMCP2Agl->GetBinContent(i+1)/(float)EffpreselMCP1Agl->GetBinContent(i+1);

	float EffpreselMCD[18][6]={{0}};
	for(int i=0;i<18;i++) for(int h=0;h<6;h++) if(EffpreselMCD2->GetBinContent(i+1,h+1)<EffpreselMCD1->GetBinContent(i+1,h+1))
		EffpreselMCD[i][h]=EffpreselMCD2->GetBinContent(i+1,h+1)/(float)EffpreselMCD1->GetBinContent(i+1,h+1);
	
	float EffpreselMCDNaF[18][6]={{0}};
        for(int i=0;i<18;i++) for(int h=0;h<6;h++) if(EffpreselMCD2NaF->GetBinContent(i+1,h+1)<EffpreselMCD1NaF->GetBinContent(i+1,h+1))
                EffpreselMCDNaF[i][h]=EffpreselMCD2NaF->GetBinContent(i+1,h+1)/(float)EffpreselMCD1NaF->GetBinContent(i+1,h+1);
	
	float EffpreselMCDAgl[18][6]={{0}};
        for(int i=0;i<18;i++) for(int h=0;h<6;h++) if(EffpreselMCD2Agl->GetBinContent(i+1,h+1)<EffpreselMCD1Agl->GetBinContent(i+1,h+1))
                EffpreselMCDAgl[i][h]=EffpreselMCD2Agl->GetBinContent(i+1,h+1)/(float)EffpreselMCD1Agl->GetBinContent(i+1,h+1);
	
	float EffpreselMCP_R[43]={0};
	for(int i=1;i<43;i++) EffpreselMCP_R[i]=EffpreselMCP2_R->GetBinContent(i+1)/(float)EffpreselMCP1_R->GetBinContent(i+1);
	float EffpreselMCD_R[43][6]={{0}};
	for(int i=4;i<43;i++) for(int h=0;h<6;h++) if(EffpreselMCD1_R->GetBinContent(i+1,h+1)>EffpreselMCD2_R->GetBinContent(i+1,h+1))
		EffpreselMCD_R[i][h]=EffpreselMCD2_R->GetBinContent(i+1,h+1)/(float)EffpreselMCD1_R->GetBinContent(i+1,h+1);

	c4->cd();
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	string MCLegend[7]={"protons.B800","d.pl1.0_520_GG_Blic","d.pl1.0_520_GG_BlicDPMJet","d.pl1.0_520_GG_QMD","d.pl1.0_520_Shen_Blic","d.pl1.0_520_Shen_BlicDPMJet","d.pl1.0_520_Shen_QMD"};
	TGraph * EffPreMCP_R = new TGraph();
	EffPreMCP_R->SetTitle(MCLegend[0].c_str());
	for(int i=0;i<43;i++) EffPreMCP_R->SetPoint(i,R_cent[i],EffpreselMCP_R[i]);
	for(int i=0;i<43;i++) EffPreMCP_R_TH1F->SetBinContent(i+1,EffpreselMCP_R[i]);
	TGraph * EffPreMCD_R[6];
	EffPreMCP_R->SetMarkerColor(2);
	EffPreMCP_R->SetMarkerStyle(8);
	EffPreMCP_R->SetLineColor(2);
	EffPreMCP_R->SetLineWidth(2);
	EffPreMCP_R->SetTitle("Preselections Efficiency MC (R bins)");
	EffPreMCP_R->GetXaxis()->SetTitle("R [GV]");
	EffPreMCP_R->GetYaxis()->SetTitle("Pres. Efficiency");
	EffPreMCP_R->GetXaxis()->SetTitleSize(0.045);
	EffPreMCP_R->GetYaxis()->SetTitleSize(0.045);
	{
		EffPreMCP_R->Draw("ACP");
		TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
		leg->AddEntry(EffPreMCP_R,MCLegend[0].c_str(), "ep");

		for(int h=0;h<6;h++){
			EffPreMCD_R[h]= new TGraph();
			EffPreMCD_R[h]->SetTitle(MCLegend[h+1].c_str());
			for(int i=1;i<43;i++) EffPreMCD_R[h]->SetPoint(i,R_cent[i],EffpreselMCD_R[i][h]);
			for(int i=1;i<43;i++) EffPreMCD_R_TH2F->SetBinContent(i+1,h+1,EffpreselMCD_R[i][h]);
			leg->AddEntry(EffPreMCD_R[h],MCLegend[h+1].c_str(), "ep");
			EffPreMCD_R[h]->SetMarkerColor(4);
			EffPreMCD_R[h]->SetMarkerStyle(h+3);
			EffPreMCD_R[h]->SetMarkerSize(2);
			EffPreMCD_R[h]->SetLineColor(4);
			EffPreMCD_R[h]->SetLineWidth(2);
			EffPreMCD_R[h]->Draw("Psame");
			leg->Draw();
		}
	}

	c4_bis->cd(1);
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraph * EffPreMCP = new TGraph();
	for(int i=0;i<18;i++) EffPreMCP->SetPoint(i,Ekincent[i],EffpreselMCP[i]);
	for(int i=0;i<18;i++) EffPreMCP_TH1F->SetBinContent(i+1,EffpreselMCP[i]);
	TGraph * EffPreMCD[6];
	EffPreMCP->SetMarkerColor(2);
	EffPreMCP->SetMarkerStyle(8);
	EffPreMCP->SetLineColor(2);
	EffPreMCP->SetLineWidth(2);
	EffPreMCP->SetTitle("Preselections Efficiency MC (Beta bins)");
	EffPreMCP->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
	EffPreMCP->GetYaxis()->SetTitle("Pres. Efficiency");
	EffPreMCP->GetXaxis()->SetTitleSize(0.045);
	EffPreMCP->GetYaxis()->SetTitleSize(0.045);
	{
		EffPreMCP->Draw("ACP");
		TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
		leg->AddEntry(EffPreMCP,MCLegend[0].c_str(), "ep");

		for(int h=0;h<6;h++){
			EffPreMCD[h]= new TGraph();
			for(int i=0;i<18;i++) EffPreMCD[h]->SetPoint(i,Ekincent[i],EffpreselMCD[i][h]);
			for(int i=0;i<18;i++) EffPreMCD_TH2F->SetBinContent(i+1,h+1,EffpreselMCD[i][h]);
			EffPreMCD[h]->SetMarkerColor(4);
			EffPreMCD[h]->SetMarkerStyle(h+3);
			leg->AddEntry(EffPreMCD[h],MCLegend[h+1].c_str(), "ep");
			EffPreMCD[h]->SetMarkerSize(2);
			EffPreMCD[h]->SetLineColor(4);
			EffPreMCD[h]->SetLineWidth(2);
			EffPreMCD[h]->Draw("Psame");
			leg->Draw();
		}
	}
	
	c4_bis->cd(2);
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
        TGraph * EffPreMCPNaF = new TGraph();
        for(int i=0;i<18;i++) EffPreMCPNaF->SetPoint(i,EkincentNaF[i],EffpreselMCPNaF[i]);
        for(int i=0;i<18;i++) EffPreMCPNaF_TH1F->SetBinContent(i+1,EffpreselMCPNaF[i]);
        TGraph * EffPreMCDNaF[6];
        EffPreMCPNaF->SetMarkerColor(2);
        EffPreMCPNaF->SetMarkerStyle(8);
        EffPreMCPNaF->SetLineColor(2);
        EffPreMCPNaF->SetLineWidth(2);
        EffPreMCPNaF->SetTitle("Preselections Efficiency MC (Beta bins NaF)");
        EffPreMCPNaF->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
        EffPreMCPNaF->GetYaxis()->SetTitle("Pres. Efficiency");
        EffPreMCPNaF->GetXaxis()->SetTitleSize(0.045);
        EffPreMCPNaF->GetYaxis()->SetTitleSize(0.045);
        {
                EffPreMCPNaF->Draw("ACP");
                TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(EffPreMCPNaF,MCLegend[0].c_str(), "ep");

                for(int h=0;h<6;h++){
                        EffPreMCDNaF[h]= new TGraph();
                        for(int i=0;i<18;i++) EffPreMCDNaF[h]->SetPoint(i,EkincentNaF[i],EffpreselMCDNaF[i][h]);
                        for(int i=0;i<18;i++) EffPreMCDNaF_TH2F->SetBinContent(i+1,h+1,EffpreselMCDNaF[i][h]);
                        EffPreMCDNaF[h]->SetMarkerColor(4);
                        EffPreMCDNaF[h]->SetMarkerStyle(h+3);
                        leg->AddEntry(EffPreMCDNaF[h],MCLegend[h+1].c_str(), "ep");
                        EffPreMCDNaF[h]->SetMarkerSize(2);
                        EffPreMCDNaF[h]->SetLineColor(4);
                        EffPreMCDNaF[h]->SetLineWidth(2);
                        EffPreMCDNaF[h]->Draw("Psame");
                        leg->Draw();
                }
        }

	c4_bis->cd(3);
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
        TGraph * EffPreMCPAgl = new TGraph();
        for(int i=0;i<18;i++) EffPreMCPAgl->SetPoint(i,EkincentAgl[i],EffpreselMCPAgl[i]);
        for(int i=0;i<18;i++) EffPreMCPAgl_TH1F->SetBinContent(i+1,EffpreselMCPAgl[i]);
        TGraph * EffPreMCDAgl[6];
        EffPreMCPAgl->SetMarkerColor(2);
        EffPreMCPAgl->SetMarkerStyle(8);
        EffPreMCPAgl->SetLineColor(2);
        EffPreMCPAgl->SetLineWidth(2);
        EffPreMCPAgl->SetTitle("Preselections Efficiency MC (Beta bins Agl)");
        EffPreMCPAgl->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
        EffPreMCPAgl->GetYaxis()->SetTitle("Pres. Efficiency");
        EffPreMCPAgl->GetXaxis()->SetTitleSize(0.045);
        EffPreMCPAgl->GetYaxis()->SetTitleSize(0.045);
        {
                EffPreMCPAgl->Draw("ACP");
                TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(EffPreMCPAgl,MCLegend[0].c_str(), "ep");

                for(int h=0;h<6;h++){
                        EffPreMCDAgl[h]= new TGraph();
                        for(int i=0;i<18;i++) EffPreMCDAgl[h]->SetPoint(i,EkincentAgl[i],EffpreselMCDAgl[i][h]);
                        for(int i=0;i<18;i++) EffPreMCDAgl_TH2F->SetBinContent(i+1,h+1,EffpreselMCDAgl[i][h]);
                        EffPreMCDAgl[h]->SetMarkerColor(4);
                        EffPreMCDAgl[h]->SetMarkerStyle(h+3);
                        leg->AddEntry(EffPreMCDAgl[h],MCLegend[h+1].c_str(), "ep");
                        EffPreMCDAgl[h]->SetMarkerSize(2);
                        EffPreMCDAgl[h]->SetLineColor(4);
                        EffPreMCDAgl[h]->SetLineWidth(2);
                        EffPreMCDAgl[h]->Draw("Psame");
                        leg->Draw();
                }
        }
	
}
