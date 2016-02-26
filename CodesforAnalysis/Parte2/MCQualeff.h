using namespace std;

TCanvas *c5=new TCanvas("Likelihood Efficiency (R bins)");
TCanvas *c6=new TCanvas("Distance Efficiency (R bins)");
TCanvas *c5_bis=new TCanvas("Likelihood Efficiency (Beta bins)");
TCanvas *c6_bis=new TCanvas("Distance Efficiency (Beta bins)");	

TH1F * EffDistMCP = new TH1F("EffDistMCP","EffDistMCP",43,0,43);
TH1F * EffLikMCP = new TH1F("EffLikMCP","EffLikMCP",43,0,43);
TH1F * EffDistMCP_Beta = new TH1F("EffDistMCP_Beta","EffDistMCP_Beta",18,0,18);
TH1F * EffLikMCP_Beta = new TH1F("EffLikMCP_Beta","EffLikMCP_Beta",18,0,18);
TH1F * EffDistMCP_BetaNaF = new TH1F("EffDistMCP_BetaNaF","EffDistMCP_BetaNaF",18,0,18);
TH1F * EffLikMCP_BetaNaF = new TH1F("EffLikMCP_BetaNaF","EffLikMCP_BetaNaF",18,0,18);
TH1F * EffDistMCP_BetaAgl = new TH1F("EffDistMCP_BetaAgl","EffDistMCP_BetaAgl",18,0,18);
TH1F * EffLikMCP_BetaAgl = new TH1F("EffLikMCP_BetaAgl","EffLikMCP_BetaAgl",18,0,18);

TH2F * EffDistMCD = new TH2F("EffDistMCD","EffDistMCD",43,0,43,6,0,6);
TH2F * EffLikMCD = new TH2F("EffLikMCD","EffLikMCD",43,0,43,6,0,6);
TH2F * EffDistMCD_Beta = new TH2F("EffDistMCD_Beta","EffDistMCD_Beta",18,0,18,6,0,6);
TH2F * EffLikMCD_Beta = new TH2F("EffLikMCD_Beta","EffLikMCD_Beta",18,0,18,6,0,6);
TH2F * EffDistMCD_BetaNaF = new TH2F("EffDistMCD_BetaNaF","EffDistMCD_BetaNaF",18,0,18,6,0,6);
TH2F * EffLikMCD_BetaNaF = new TH2F("EffLikMCD_BetaNaF","EffLikMCD_BetaNaF",18,0,18,6,0,6);
TH2F * EffDistMCD_BetaAgl = new TH2F("EffDistMCD_BetaAgl","EffDistMCD_BetaAgl",18,0,18,6,0,6);
TH2F * EffLikMCD_BetaAgl = new TH2F("EffLikMCD_BetaAgl","EffLikMCD_BetaAgl",18,0,18,6,0,6);

TH1F * EffQualMCP = new TH1F("EffQualMCP","EffQualMCP",43,0,43);
TH1F * EffQualMCP_Beta = new TH1F("EffQualMCP_Beta","EffQualMCP_Beta",18,0,18);
TH1F * EffQualMCP_BetaNaF = new TH1F("EffQualMCP_BetaNaF","EffQualMCP_BetaNaF",18,0,18);
TH1F * EffQualMCP_BetaAgl = new TH1F("EffQualMCP_BetaAgl","EffQualMCP_BetaAgl",18,0,18);
TH2F * EffQualMCD = new TH2F("EffQualMCD","EffQualMCD",43,0,43,6,0,6);
TH2F * EffQualMCD_Beta = new TH2F("EffQualMCD_Beta","EffQualMCD_Beta",18,0,18,6,0,6);
TH2F * EffQualMCD_BetaNaF = new TH2F("EffQualMCD_BetaNaF","EffQualMCD_BetaNaF",18,0,18,6,0,6);
TH2F * EffQualMCD_BetaAgl = new TH2F("EffQualMCD_BetaAgl","EffQualMCD_BetaAgl",18,0,18,6,0,6);

TH1F * EffMCLikP_TH1F = new TH1F("EffMCLikP_TH1F","EffMCLikP_TH1F",43,0,43);
TH2F * EffMCLikD_TH2F = new TH2F("EffMCLikD_TH2F","EffMCLikD_TH2F",43,0,43,6,0,6);
TH1F * EffMCLikP_Beta_TH1F = new TH1F("EffMCLikP_Beta_TH1F","EffMCLikP_Beta_TH1F",18,0,18);
TH2F * EffMCLikD_Beta_TH2F = new TH2F("EffMCLikD_Beta_TH2F","EffMCLikD_Beta_TH2F",18,0,18,6,0,6);
TH1F * EffMCLikP_BetaNaF_TH1F = new TH1F("EffMCLikP_BetaNaF_TH1F","EffMCLikP_BetaNaF_TH1F",18,0,18);
TH2F * EffMCLikD_BetaNaF_TH2F = new TH2F("EffMCLikD_BetaNaF_TH2F","EffMCLikD_BetaNaF_TH2F",18,0,18,6,0,6);
TH1F * EffMCLikP_BetaAgl_TH1F = new TH1F("EffMCLikP_BetaAgl_TH1F","EffMCLikP_BetaAgl_TH1F",18,0,18);
TH2F * EffMCLikD_BetaAgl_TH2F = new TH2F("EffMCLikD_BetaAgl_TH2F","EffMCLikD_BetaAgl_TH2F",18,0,18,6,0,6);

TH1F * EffMCDistP_TH1F = new TH1F("EffMCDistP_TH1F","EffMCDistP_TH1F",43,0,43);
TH2F * EffMCDistD_TH2F = new TH2F("EffMCDistD_TH2F","EffMCDistD_TH2F",43,0,43,6,0,6);
TH1F * EffMCDistP_Beta_TH1F = new TH1F("EffMCDistP_Beta_TH1F","EffMCDistP_Beta_TH1F",18,0,18);
TH2F * EffMCDistD_Beta_TH2F = new TH2F("EffMCDistD_Beta_TH2F","EffMCDistD_Beta_TH2F",18,0,18,6,0,6);
TH1F * EffMCDistP_BetaNaF_TH1F = new TH1F("EffMCDistP_BetaNaF_TH1F","EffMCDistP_BetaNaF_TH1F",18,0,18);
TH2F * EffMCDistD_BetaNaF_TH2F = new TH2F("EffMCDistD_BetaNaF_TH2F","EffMCDistD_BetaNaF_TH2F",18,0,18,6,0,6);
TH1F * EffMCDistP_BetaAgl_TH1F = new TH1F("EffMCDistP_BetaAgl_TH1F","EffMCDistP_BetaAgl_TH1F",18,0,18);
TH2F * EffMCDistD_BetaAgl_TH2F = new TH2F("EffMCDistD_BetaAgl_TH2F","EffMCDistD_BetaAgl_TH2F",18,0,18,6,0,6);

void MCQualeff_Fill(TNtuple *ntupla, int l){
	
	int k = ntupla->GetEvent(l);
	if(Beta<=0||R<=0||Unbias!=0) return;
	if(Massa_gen<1){
		for(int K=0;K<43;K++) if(R<bin[K+1]&&R>bin[K]) {EffQualMCP->Fill(K);}
                for(int m=0;m<18;m++) if(Var<BetaP[m+1]&&Var>BetaP[m]) EffQualMCP_Beta->Fill(m);
		for(int m=0;m<18;m++)
                                if((((int)Cutmask)>>11)==512&&Var2>BetaNaFP[m]&&Var2<=BetaNaFP[m+1]) EffQualMCP_BetaNaF->Fill(m);
		for(int m=0;m<18;m++)
                                if((((int)Cutmask)>>11)==0&&Var2>BetaAglP[m]&&Var2<=BetaAglP[m+1]) EffQualMCP_BetaAgl->Fill(m);		

		if(Dist5D_P<6&&Likcut) for(int K=0;K<43;K++) if(R<bin[K+1]&&R>bin[K]) EffDistMCP->Fill(K);
		if(Distcut&&Likcut) {
			for(int m=0;m<18;m++) if(Var<BetaP[m+1]&&Var>BetaP[m]) EffDistMCP_Beta->Fill(m);
			for(int m=0;m<18;m++) if((((int)Cutmask)>>11)==512&&Var2>BetaNaFP[m]&&Var2<=BetaNaFP[m+1]) EffDistMCP_BetaNaF->Fill(m);
                        for(int m=0;m<18;m++) if((((int)Cutmask)>>11)==0&&Var2>BetaAglP[m]&&Var2<=BetaAglP[m+1]) EffDistMCP_BetaAgl->Fill(m);			
		}
		if(Likcut){
			for(int K=0;K<43;K++) if(R<bin[K+1]&&R>bin[K]) {EffLikMCP->Fill(K);}
			for(int m=0;m<18;m++) if(Var<BetaP[m+1]&&Var>BetaP[m]) EffLikMCP_Beta->Fill(m);
			for(int m=0;m<18;m++) if((((int)Cutmask)>>11)==512&&Var2>BetaNaFP[m]&&Var2<=BetaNaFP[m+1]) EffLikMCP_BetaNaF->Fill(m);
			for(int m=0;m<18;m++) if((((int)Cutmask)>>11)==0&&Var2>BetaAglP[m]&&Var2<=BetaAglP[m+1]) EffLikMCP_BetaAgl->Fill(m);
		}
	}
	if(Massa_gen<2&&Massa_gen>1){
                for(int K=0;K<43;K++) if(R<bin[K+1]&&R>bin[K]) {EffQualMCD->Fill(K,(int)(10000*Massa_gen-18570));}
                for(int m=0;m<18;m++) if(Var<BetaD[m+1]&&Var>BetaD[m]) EffQualMCD_Beta->Fill(m,(int)(10000*Massa_gen-18570));
		for(int m=0;m<18;m++)
                                if((((int)Cutmask)>>11)==512&&Var2>BetaNaFD[m]&&Var2<=BetaNaFD[m+1]) EffQualMCD_BetaNaF->Fill(m,(int)(10000*Massa_gen-18570));
		for(int m=0;m<18;m++)
                                if((((int)Cutmask)>>11)==0&&Var2>BetaAglD[m]&&Var2<=BetaAglD[m+1]) EffQualMCD_BetaAgl->Fill(m,(int)(10000*Massa_gen-18570));

		if(Dist5D<6&&Likcut) for(int K=0;K<43;K++) if(R<bin[K+1]&&R>bin[K]) EffDistMCD->Fill(K,(int)(10000*Massa_gen-18570));
		if(Distcut&&Likcut){       
			for(int m=0;m<18;m++) if(Var<BetaD[m+1]&&Var>BetaD[m]) EffDistMCD_Beta->Fill(m,(int)(10000*Massa_gen-18570));
			for(int m=0;m<18;m++) 
                                if((((int)Cutmask)>>11)==512&&Var2>BetaNaFD[m]&&Var2<=BetaNaFD[m+1]) EffDistMCD_BetaNaF->Fill(m,(int)(10000*Massa_gen-18570));
                        for(int m=0;m<18;m++) 
                                if((((int)Cutmask)>>11)==0&&Var2>BetaAglD[m]&&Var2<=BetaAglD[m+1]) EffDistMCD_BetaAgl->Fill(m,(int)(10000*Massa_gen-18570));
		}
                if(Likcut){
                        for(int K=0;K<43;K++) if(R<bin[K+1]&&R>bin[K]) {EffLikMCD->Fill(K,(int)(10000*Massa_gen-18570));}
                        for(int m=0;m<18;m++) if(Var<BetaD[m+1]&&Var>BetaD[m]) EffLikMCD_Beta->Fill(m,(int)(10000*Massa_gen-18570));
                	for(int m=0;m<18;m++) 
				if((((int)Cutmask)>>11)==512&&Var2>BetaNaFD[m]&&Var2<=BetaNaFD[m+1]) EffLikMCD_BetaNaF->Fill(m,(int)(10000*Massa_gen-18570)); 
                        for(int m=0;m<18;m++) 
				if((((int)Cutmask)>>11)==0&&Var2>BetaAglD[m]&&Var2<=BetaAglD[m+1]) EffLikMCD_BetaAgl->Fill(m,(int)(10000*Massa_gen-18570));
	
		}
        }
	
	return;
}


void MCQualeff_Copy(TFile * file1){
	EffQualMCP = (TH1F *)file1->Get("EffQualMCP");
	EffQualMCP_Beta = (TH1F *)file1->Get("EffQualMCP_Beta");
	EffQualMCP_BetaNaF = (TH1F *)file1->Get("EffQualMCP_BetaNaF");
	EffQualMCP_BetaAgl = (TH1F *)file1->Get("EffQualMCP_BetaAgl");	
	EffQualMCD = (TH2F *)file1->Get("EffQualMCD");
        EffQualMCD_Beta = (TH2F *)file1->Get("EffQualMCD_Beta");
	EffQualMCD_BetaNaF = (TH2F *)file1->Get("EffQualMCD_BetaNaF");
	EffQualMCD_BetaAgl = (TH2F *)file1->Get("EffQualMCD_BetaAgl");
	EffDistMCP = (TH1F *)file1->Get("EffDistMCP");
	EffLikMCP = (TH1F *)file1->Get("EffLikMCP");
	EffDistMCP_Beta = (TH1F *)file1->Get("EffDistMCP_Beta");
	EffLikMCP_Beta = (TH1F *)file1->Get("EffLikMCP_Beta");
	EffDistMCP_BetaNaF = (TH1F *)file1->Get("EffDistMCP_BetaNaF");
        EffLikMCP_BetaNaF = (TH1F *)file1->Get("EffLikMCP_BetaNaF");
	EffDistMCP_BetaAgl = (TH1F *)file1->Get("EffDistMCP_BetaAgl");
        EffLikMCP_BetaAgl = (TH1F *)file1->Get("EffLikMCP_BetaAgl");
	EffDistMCD = (TH2F *)file1->Get("EffDistMCD");
        EffLikMCD = (TH2F *)file1->Get("EffLikMCD");
        EffDistMCD_Beta = (TH2F *)file1->Get("EffDistMCD_Beta");
        EffLikMCD_Beta = (TH2F *)file1->Get("EffLikMCD_Beta");
	EffDistMCD_BetaNaF = (TH2F *)file1->Get("EffDistMCD_BetaNaF");
        EffLikMCD_BetaNaF = (TH2F *)file1->Get("EffLikMCD_BetaNaF");
	EffDistMCD_BetaAgl = (TH2F *)file1->Get("EffDistMCD_BetaAgl");
        EffLikMCD_BetaAgl = (TH2F *)file1->Get("EffLikMCD_BetaAgl");

	return;
}

void MCQualeff_Write(){
        EffQualMCP->Write(); 
        EffQualMCP_Beta ->Write();
        EffQualMCP_BetaNaF ->Write();
        EffQualMCP_BetaAgl ->Write();
        EffQualMCD ->Write();
        EffQualMCD_Beta->Write();
        EffQualMCD_BetaNaF->Write();
        EffQualMCD_BetaAgl->Write();
        EffDistMCP ->Write();
        EffLikMCP->Write();
        EffDistMCP_Beta ->Write();
        EffLikMCP_Beta ->Write();
        EffDistMCP_BetaNaF ->Write();
        EffLikMCP_BetaNaF->Write();
        EffDistMCP_BetaAgl ->Write();
        EffLikMCP_BetaAgl ->Write();
        EffDistMCD->Write();
        EffLikMCD ->Write();
        EffDistMCD_Beta ->Write();
        EffLikMCD_Beta ->Write();
        EffDistMCD_BetaNaF ->Write();
        EffLikMCD_BetaNaF ->Write();
        EffDistMCD_BetaAgl->Write();
        EffLikMCD_BetaAgl ->Write();

        return;
}


void MCQualeff(TFile * file1){
	cout<<"**** MC QUALITY SEL. EFFICIENCIES ****"<<endl;		

	float EffLikP[43]={0};
	for(int i=1;i<43;i++) EffLikP[i]=EffLikMCP->GetBinContent(i+1)/(float)EffQualMCP->GetBinContent(i+1);
	float EffLikD[43][6]={{0}};
	for(int i=1;i<43;i++) for(int h=0;h<6;h++) if(EffQualMCD->GetBinContent(i+1,h+1)>0&&EffLikMCD->GetBinContent(i+1,h+1)<=EffQualMCD->GetBinContent(i+1,h+1))
		EffLikD[i][h]=EffLikMCD->GetBinContent(i+1,h+1)/(float)EffQualMCD->GetBinContent(i+1,h+1);
	float EffDistP[43]={0};
	for(int i=1;i<43;i++) EffDistP[i]=EffDistMCP->GetBinContent(i+1)/(float)EffQualMCP->GetBinContent(i+1);
	float EffDistD[43][6]={{0}};
	for(int i=1;i<43;i++) for(int h=0;h<6;h++) if(EffQualMCD->GetBinContent(i+1,h+1)>0&&EffDistMCD->GetBinContent(i+1,h+1)<=EffQualMCD->GetBinContent(i+1,h+1)) 
		EffDistD[i][h]=EffDistMCD->GetBinContent(i+1,h+1)/(float)EffQualMCD->GetBinContent(i+1,h+1);

	float EffLikP_Beta[18]={0};
	for(int i=0;i<18;i++) if(EffQualMCP_Beta->GetBinContent(i+1)>0) if(EffLikMCP_Beta->GetBinContent(i+1)<=EffQualMCP_Beta->GetBinContent(i+1))
		EffLikP_Beta[i]=EffLikMCP_Beta->GetBinContent(i+1)/(float)EffQualMCP_Beta->GetBinContent(i+1);
	float EffLikD_Beta[18][6]={{0}};
	for(int i=0;i<18;i++) for(int h=0;h<6;h++) 
		if(EffQualMCD_Beta->GetBinContent(i+1,h+1)>0) if(EffLikMCD_Beta->GetBinContent(i+1,h+1)<=EffQualMCD_Beta->GetBinContent(i+1,h+1))
			EffLikD_Beta[i][h]=EffLikMCD_Beta->GetBinContent(i+1,h+1)/(float)EffQualMCD_Beta->GetBinContent(i+1,h+1);	
	float EffLikP_BetaNaF[18]={0};
        for(int i=0;i<18;i++) if(EffQualMCP_BetaNaF->GetBinContent(i+1)>0) if(EffLikMCP_BetaNaF->GetBinContent(i+1)<=EffQualMCP_BetaNaF->GetBinContent(i+1))
                EffLikP_BetaNaF[i]=EffLikMCP_BetaNaF->GetBinContent(i+1)/(float)EffQualMCP_BetaNaF->GetBinContent(i+1);
        float EffLikD_BetaNaF[18][6]={{0}};
        for(int i=0;i<18;i++) for(int h=0;h<6;h++)
                if(EffQualMCD_BetaNaF->GetBinContent(i+1,h+1)>0) if(EffLikMCD_BetaNaF->GetBinContent(i+1,h+1)<=EffQualMCD_BetaNaF->GetBinContent(i+1,h+1))
                        EffLikD_BetaNaF[i][h]=EffLikMCD_BetaNaF->GetBinContent(i+1,h+1)/(float)EffQualMCD_BetaNaF->GetBinContent(i+1,h+1);
	float EffLikP_BetaAgl[18]={0};
        for(int i=0;i<18;i++) if(EffQualMCP_BetaAgl->GetBinContent(i+1)>0) if(EffLikMCP_BetaAgl->GetBinContent(i+1)<=EffQualMCP_BetaAgl->GetBinContent(i+1))
                EffLikP_BetaAgl[i]=EffLikMCP_BetaAgl->GetBinContent(i+1)/(float)EffQualMCP_BetaAgl->GetBinContent(i+1);
        float EffLikD_BetaAgl[18][6]={{0}};
        for(int i=0;i<18;i++) for(int h=0;h<6;h++)
                if(EffQualMCD_BetaAgl->GetBinContent(i+1,h+1)>0) if(EffLikMCD_BetaAgl->GetBinContent(i+1,h+1)<=EffQualMCD_BetaAgl->GetBinContent(i+1,h+1))
                        EffLikD_BetaAgl[i][h]=EffLikMCD_BetaAgl->GetBinContent(i+1,h+1)/(float)EffQualMCD_BetaAgl->GetBinContent(i+1,h+1);

	
	float EffDistP_Beta[18]={0};
        for(int i=0;i<18;i++) if(EffQualMCP_Beta->GetBinContent(i+1)>0) if(EffDistMCP_Beta->GetBinContent(i+1)<=EffQualMCP_Beta->GetBinContent(i+1))
                EffDistP_Beta[i]=EffDistMCP_Beta->GetBinContent(i+1)/(float)EffQualMCP_Beta->GetBinContent(i+1);
	float EffDistD_Beta[18][6]={{0}};
        for(int i=0;i<18;i++) for(int h=0;h<6;h++)
                if(EffQualMCD_Beta->GetBinContent(i+1,h+1)>0) if(EffDistMCD_Beta->GetBinContent(i+1,h+1)<=EffQualMCD_Beta->GetBinContent(i+1,h+1))
                        EffDistD_Beta[i][h]=EffDistMCD_Beta->GetBinContent(i+1,h+1)/(float)EffQualMCD_Beta->GetBinContent(i+1,h+1);
        float EffDistP_BetaNaF[18]={0};
        for(int i=0;i<18;i++) if(EffQualMCP_BetaNaF->GetBinContent(i+1)>0) if(EffDistMCP_BetaNaF->GetBinContent(i+1)<=EffQualMCP_BetaNaF->GetBinContent(i+1))
                EffDistP_BetaNaF[i]=EffDistMCP_BetaNaF->GetBinContent(i+1)/(float)EffQualMCP_BetaNaF->GetBinContent(i+1);
        float EffDistD_BetaNaF[18][6]={{0}};
        for(int i=0;i<18;i++) for(int h=0;h<6;h++)
                if(EffQualMCD_BetaNaF->GetBinContent(i+1,h+1)>0) if(EffDistMCD_BetaNaF->GetBinContent(i+1,h+1)<=EffQualMCD_BetaNaF->GetBinContent(i+1,h+1))
                        EffDistD_BetaNaF[i][h]=EffDistMCD_BetaNaF->GetBinContent(i+1,h+1)/(float)EffQualMCD_BetaNaF->GetBinContent(i+1,h+1);
        float EffDistP_BetaAgl[18]={0};
        for(int i=0;i<18;i++) if(EffQualMCP_BetaAgl->GetBinContent(i+1)>0) if(EffDistMCP_BetaAgl->GetBinContent(i+1)<=EffQualMCP_BetaAgl->GetBinContent(i+1))
                EffDistP_BetaAgl[i]=EffDistMCP_BetaAgl->GetBinContent(i+1)/(float)EffQualMCP_BetaAgl->GetBinContent(i+1);
        float EffDistD_BetaAgl[18][6]={{0}};
        for(int i=0;i<18;i++) for(int h=0;h<6;h++)
                if(EffQualMCD_BetaAgl->GetBinContent(i+1,h+1)>0) if(EffDistMCD_BetaAgl->GetBinContent(i+1,h+1)<=EffQualMCD_BetaAgl->GetBinContent(i+1,h+1))
                        EffDistD_BetaAgl[i][h]=EffDistMCD_BetaAgl->GetBinContent(i+1,h+1)/(float)EffQualMCD_BetaAgl->GetBinContent(i+1,h+1);
	
		

	c5->cd();
	string MCLegend[7]={"protons.B800","d.pl1.0_520_GG_Blic","d.pl1.0_520_GG_BlicDPMJet","d.pl1.0_520_GG_QMD","d.pl1.0_520_Shen_Blic","d.pl1.0_520_Shen_BlicDPMJet","d.pl1.0_520_Shen_QMD"};
	gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
	TGraph *EffMCLikP= new TGraph();
	TGraph *EffMCLikD[6];
	int j=0;
	for(int i=1;i<43;i++) {EffMCLikP->SetPoint(j,R_cent[i],EffLikP[i]);j++;}	
	for(int i=1;i<43;i++) EffMCLikP_TH1F->SetBinContent(i+1,EffLikP[i]);
	EffMCLikP->SetMarkerColor(2);
        EffMCLikP->SetMarkerStyle(8);
        EffMCLikP->SetLineColor(2);
        EffMCLikP->SetLineWidth(2);
        EffMCLikP->SetTitle("Likelihood Efficiency MC on top of Pres. (R bins)");
        EffMCLikP->GetXaxis()->SetTitle("R [GV]");
        EffMCLikP->GetYaxis()->SetTitle("Efficiency");
        EffMCLikP->GetXaxis()->SetTitleSize(0.045);
        EffMCLikP->GetYaxis()->SetTitleSize(0.045);
	{
                EffMCLikP->GetYaxis()->SetRangeUser(0,1);
		EffMCLikP->Draw("ACP");
                TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(EffMCLikP,MCLegend[0].c_str(), "ep");

                for(int h=0;h<6;h++){
                        EffMCLikD[h]= new TGraph();
                        EffMCLikD[h]->SetTitle(MCLegend[h+1].c_str());
                        for(int i=1;i<43;i++) EffMCLikD[h]->SetPoint(i,R_cent[i],EffLikD[i][h]);
			for(int i=1;i<43;i++) EffMCLikD_TH2F->SetBinContent(i+1,h+1,EffLikD[i][h]);
			//leg->AddEntry(EffMCLikD[h],MCLegend[h+1].c_str(), "ep");
                        EffMCLikD[h]->SetMarkerColor(4);
                        EffMCLikD[h]->SetMarkerStyle(h+3);
                        EffMCLikD[h]->SetMarkerSize(2);
                        EffMCLikD[h]->SetLineColor(4);
                        EffMCLikD[h]->SetLineWidth(1);
                        EffMCLikD[h]->Draw("Psame");
                        leg->Draw();
                } 
        }

	c5_bis->Divide(3,1);
	c5_bis->cd(1);
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
        TGraph *EffMCLikP_Beta= new TGraph();
	TGraph *EffMCLikD_Beta[6];
        j=0;
        for(int i=0;i<18;i++) {EffMCLikP_Beta->SetPoint(j,Ekincent[i],EffLikP_Beta[i]);j++;}
	for(int i=0;i<18;i++) EffMCLikP_Beta_TH1F->SetBinContent(i+1,EffLikP_Beta[i]);
        EffMCLikP_Beta->SetMarkerColor(2);
        EffMCLikP_Beta->SetMarkerStyle(8);
        EffMCLikP_Beta->SetLineColor(2);
        EffMCLikP_Beta->SetLineWidth(2);
        EffMCLikP_Beta->SetTitle("Likelihood Efficiency MC on top of Pres. (Beta bins)");
        EffMCLikP_Beta->GetXaxis()->SetTitle("Kin. En.  [GeV/nucl.]");
        EffMCLikP_Beta->GetYaxis()->SetTitle("Efficiency");
        EffMCLikP_Beta->GetXaxis()->SetTitleSize(0.045);
        EffMCLikP_Beta->GetYaxis()->SetTitleSize(0.045);
	{
                EffMCLikP_Beta->GetYaxis()->SetRangeUser(0,1);
		EffMCLikP_Beta->Draw("ACP");
                TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(EffMCLikP_Beta,MCLegend[0].c_str(), "ep");

                for(int h=0;h<6;h++){
                        EffMCLikD_Beta[h]= new TGraph();
                        EffMCLikD_Beta[h]->SetTitle(MCLegend[h+1].c_str());
                        for(int i=1;i<18;i++) EffMCLikD_Beta[h]->SetPoint(i,Ekincent[i],EffLikD_Beta[i][h]);
			for(int i=1;i<18;i++) EffMCLikD_Beta_TH2F->SetBinContent(i+1,h+1,EffLikD_Beta[i][h]);
                        leg->AddEntry(EffMCLikD_Beta[h],MCLegend[h+1].c_str(), "ep");
                        EffMCLikD_Beta[h]->SetMarkerColor(4);
                        EffMCLikD_Beta[h]->SetMarkerStyle(h+3);
                        EffMCLikD_Beta[h]->SetMarkerSize(2);
                        EffMCLikD_Beta[h]->SetLineColor(4);
                        EffMCLikD_Beta[h]->SetLineWidth(1);
                        EffMCLikD_Beta[h]->Draw("Psame");
                        leg->Draw();
                }
        }
	c5_bis->cd(2);
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
        TGraph *EffMCLikP_BetaNaF= new TGraph();
        TGraph *EffMCLikD_BetaNaF[6];
        j=0;
        for(int i=0;i<18;i++) {EffMCLikP_BetaNaF->SetPoint(j,EkincentNaF[i],EffLikP_BetaNaF[i]);j++;}
        for(int i=0;i<18;i++) EffMCLikP_BetaNaF_TH1F->SetBinContent(i+1,EffLikP_BetaNaF[i]);
        EffMCLikP_BetaNaF->SetMarkerColor(2);
        EffMCLikP_BetaNaF->SetMarkerStyle(8);
        EffMCLikP_BetaNaF->SetLineColor(2);
        EffMCLikP_BetaNaF->SetLineWidth(2);
        EffMCLikP_BetaNaF->SetTitle("Likelihood Efficiency MC on top of Pres. (Beta bins)");
        EffMCLikP_BetaNaF->GetXaxis()->SetTitle("Kin. En.  [GeV/nucl.]");
        EffMCLikP_BetaNaF->GetYaxis()->SetTitle("Efficiency");
        EffMCLikP_BetaNaF->GetXaxis()->SetTitleSize(0.045);
        EffMCLikP_BetaNaF->GetYaxis()->SetTitleSize(0.045);
        {
                EffMCLikP_BetaNaF->GetYaxis()->SetRangeUser(0,1);
		EffMCLikP_BetaNaF->Draw("ACP");
                TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(EffMCLikP_BetaNaF,MCLegend[0].c_str(), "ep");

                for(int h=0;h<6;h++){
                        EffMCLikD_BetaNaF[h]= new TGraph();
                        EffMCLikD_BetaNaF[h]->SetTitle(MCLegend[h+1].c_str());
                        for(int i=1;i<18;i++) EffMCLikD_BetaNaF[h]->SetPoint(i,EkincentNaF[i],EffLikD_BetaNaF[i][h]);
                        for(int i=1;i<18;i++) EffMCLikD_BetaNaF_TH2F->SetBinContent(i+1,h+1,EffLikD_BetaNaF[i][h]);
                        leg->AddEntry(EffMCLikD_BetaNaF[h],MCLegend[h+1].c_str(), "ep");
                        EffMCLikD_BetaNaF[h]->SetMarkerColor(4);
                        EffMCLikD_BetaNaF[h]->SetMarkerStyle(h+3);
                        EffMCLikD_BetaNaF[h]->SetMarkerSize(2);
                        EffMCLikD_BetaNaF[h]->SetLineColor(4);
                        EffMCLikD_BetaNaF[h]->SetLineWidth(1);
                        EffMCLikD_BetaNaF[h]->Draw("Psame");
                        leg->Draw();
                }
        }
	c5_bis->cd(3);
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
        TGraph *EffMCLikP_BetaAgl= new TGraph();
        TGraph *EffMCLikD_BetaAgl[6];
        j=0;
        for(int i=0;i<18;i++) {EffMCLikP_BetaAgl->SetPoint(j,EkincentAgl[i],EffLikP_BetaAgl[i]);j++;}
        for(int i=0;i<18;i++) EffMCLikP_BetaAgl_TH1F->SetBinContent(i+1,EffLikP_BetaAgl[i]);
        EffMCLikP_BetaAgl->SetMarkerColor(2);
        EffMCLikP_BetaAgl->SetMarkerStyle(8);
        EffMCLikP_BetaAgl->SetLineColor(2);
        EffMCLikP_BetaAgl->SetLineWidth(2);
        EffMCLikP_BetaAgl->SetTitle("Likelihood Efficiency MC on top of Pres. (Beta bins)");
        EffMCLikP_BetaAgl->GetXaxis()->SetTitle("Kin. En.  [GeV/nucl.]");
        EffMCLikP_BetaAgl->GetYaxis()->SetTitle("Efficiency");
        EffMCLikP_BetaAgl->GetXaxis()->SetTitleSize(0.045);
        EffMCLikP_BetaAgl->GetYaxis()->SetTitleSize(0.045);
        {
                EffMCLikP_BetaAgl->GetYaxis()->SetRangeUser(0,1);
		EffMCLikP_BetaAgl->Draw("ACP");
                TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(EffMCLikP_BetaAgl,MCLegend[0].c_str(), "ep");

                for(int h=0;h<6;h++){
                        EffMCLikD_BetaAgl[h]= new TGraph();
                        EffMCLikD_BetaAgl[h]->SetTitle(MCLegend[h+1].c_str());
                        for(int i=1;i<18;i++) EffMCLikD_BetaAgl[h]->SetPoint(i,EkincentAgl[i],EffLikD_BetaAgl[i][h]);
                        for(int i=1;i<18;i++) EffMCLikD_BetaAgl_TH2F->SetBinContent(i+1,h+1,EffLikD_BetaAgl[i][h]);
                        leg->AddEntry(EffMCLikD_BetaAgl[h],MCLegend[h+1].c_str(), "ep");
                        EffMCLikD_BetaAgl[h]->SetMarkerColor(4);
                        EffMCLikD_BetaAgl[h]->SetMarkerStyle(h+3);
                        EffMCLikD_BetaAgl[h]->SetMarkerSize(2);
                        EffMCLikD_BetaAgl[h]->SetLineColor(4);
                        EffMCLikD_BetaAgl[h]->SetLineWidth(1);
                        EffMCLikD_BetaAgl[h]->Draw("Psame");
                        leg->Draw();
                }
        }	
        c6->cd();
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
        TGraph *EffMCDistP= new TGraph();
	TGraph *EffMCDistD[6];
        j=0;
        for(int i=1;i<43;i++) {EffMCDistP->SetPoint(j,R_cent[i],EffDistP[i]);j++;}
	for(int i=1;i<43;i++) EffMCDistP_TH1F->SetBinContent(i+1,EffDistP[i]);
        EffMCDistP->SetMarkerColor(2);
        EffMCDistP->SetMarkerStyle(8);
        EffMCDistP->SetLineColor(2);
        EffMCDistP->SetLineWidth(2);
        EffMCDistP->SetTitle("Distance Efficiency MC on top of Pres. (R bins)");
        EffMCDistP->GetXaxis()->SetTitle("R [GV]");
        EffMCDistP->GetYaxis()->SetTitle("Efficiency");
        EffMCDistP->GetXaxis()->SetTitleSize(0.045);
        EffMCDistP->GetYaxis()->SetTitleSize(0.045);
	{
                EffMCDistP->GetYaxis()->SetRangeUser(0,1);
		EffMCDistP->Draw("ACP");
                TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(EffMCDistP,MCLegend[0].c_str(), "ep");

                for(int h=0;h<6;h++){
                        EffMCDistD[h]= new TGraph();
                        EffMCDistD[h]->SetTitle(MCLegend[h+1].c_str());
                        for(int i=1;i<43;i++) EffMCDistD[h]->SetPoint(i,R_cent[i],EffDistD[i][h]);
                        for(int i=1;i<43;i++) EffMCDistD_TH2F->SetBinContent(i+1,h+1,EffDistD[i][h]);
			leg->AddEntry(EffMCDistD[h],MCLegend[h+1].c_str(), "ep");
                        EffMCDistD[h]->SetMarkerColor(4);
                        EffMCDistD[h]->SetMarkerStyle(h+3);
                        EffMCDistD[h]->SetMarkerSize(2);
                        EffMCDistD[h]->SetLineColor(4);
                        EffMCDistD[h]->SetLineWidth(1);
                        EffMCDistD[h]->Draw("Psame");
                        leg->Draw();
                }
        }

        c6_bis->Divide(3,1);
	c6_bis->cd(1);
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
        TGraph *EffMCDistP_Beta= new TGraph();
        TGraph *EffMCDistD_Beta[6];
	j=0;
        for(int i=0;i<18;i++) {EffMCDistP_Beta->SetPoint(j,Ekincent[i],EffDistP_Beta[i]);j++;}
        for(int i=0;i<18;i++) EffMCDistP_Beta_TH1F->SetBinContent(i+1,EffDistP_Beta[i]);
	EffMCDistP_Beta->SetMarkerColor(2);
        EffMCDistP_Beta->SetMarkerStyle(8);
        EffMCDistP_Beta->SetLineColor(2);
        EffMCDistP_Beta->SetLineWidth(2);
        EffMCDistP_Beta->SetTitle("Distance Efficiency MC on top of Pres. (Beta bins)");
        EffMCDistP_Beta->GetXaxis()->SetTitle("Kin. En.  [GeV/nucl.]");
        EffMCDistP_Beta->GetYaxis()->SetTitle("Efficiency");
        EffMCDistP_Beta->GetXaxis()->SetTitleSize(0.045);
        EffMCDistP_Beta->GetYaxis()->SetTitleSize(0.045);
	{
                EffMCDistP_Beta->GetYaxis()->SetRangeUser(0,1);
		EffMCDistP_Beta->Draw("ACP");
                TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(EffMCDistP_Beta,MCLegend[0].c_str(), "ep");

                for(int h=0;h<6;h++){
                        EffMCDistD_Beta[h]= new TGraph();
                        EffMCDistD_Beta[h]->SetTitle(MCLegend[h+1].c_str());
                        for(int i=1;i<18;i++) EffMCDistD_Beta[h]->SetPoint(i,Ekincent[i],EffDistD_Beta[i][h]);
			for(int i=1;i<18;i++) EffMCDistD_Beta_TH2F->SetBinContent(i+1,h+1,EffDistD_Beta[i][h]);
                        leg->AddEntry(EffMCDistD_Beta[h],MCLegend[h+1].c_str(), "ep");
                        EffMCDistD_Beta[h]->SetMarkerColor(4);
                        EffMCDistD_Beta[h]->SetMarkerStyle(h+3);
                        EffMCDistD_Beta[h]->SetMarkerSize(2);
                        EffMCDistD_Beta[h]->SetLineColor(4);
                        EffMCDistD_Beta[h]->SetLineWidth(1);
                        EffMCDistD_Beta[h]->Draw("Psame");
                        leg->Draw();
                }
        }

	c6_bis->cd(2);
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
        TGraph *EffMCDistP_BetaNaF= new TGraph();
        TGraph *EffMCDistD_BetaNaF[6];
        j=0;
        for(int i=0;i<18;i++) {EffMCDistP_BetaNaF->SetPoint(j,EkincentNaF[i],EffDistP_BetaNaF[i]);j++;}
        for(int i=0;i<18;i++) EffMCDistP_BetaNaF_TH1F->SetBinContent(i+1,EffDistP_BetaNaF[i]);
        EffMCDistP_BetaNaF->SetMarkerColor(2);
        EffMCDistP_BetaNaF->SetMarkerStyle(8);
        EffMCDistP_BetaNaF->SetLineColor(2);
        EffMCDistP_BetaNaF->SetLineWidth(2);
        EffMCDistP_BetaNaF->SetTitle("Distelihood Efficiency MC on top of Pres. (Beta bins)");
        EffMCDistP_BetaNaF->GetXaxis()->SetTitle("Kin. En.  [GeV/nucl.]");
        EffMCDistP_BetaNaF->GetYaxis()->SetTitle("Efficiency");
        EffMCDistP_BetaNaF->GetXaxis()->SetTitleSize(0.045);
        EffMCDistP_BetaNaF->GetYaxis()->SetTitleSize(0.045);
        {
                EffMCDistP_BetaNaF->GetYaxis()->SetRangeUser(0,1);
		EffMCDistP_BetaNaF->Draw("ACP");
                TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(EffMCDistP_BetaNaF,MCLegend[0].c_str(), "ep");

                for(int h=0;h<6;h++){
                        EffMCDistD_BetaNaF[h]= new TGraph();
                        EffMCDistD_BetaNaF[h]->SetTitle(MCLegend[h+1].c_str());
                        for(int i=1;i<18;i++) EffMCDistD_BetaNaF[h]->SetPoint(i,EkincentNaF[i],EffDistD_BetaNaF[i][h]);
                        for(int i=1;i<18;i++) EffMCDistD_BetaNaF_TH2F->SetBinContent(i+1,h+1,EffDistD_BetaNaF[i][h]);
                        leg->AddEntry(EffMCDistD_BetaNaF[h],MCLegend[h+1].c_str(), "ep");
                        EffMCDistD_BetaNaF[h]->SetMarkerColor(4);
                        EffMCDistD_BetaNaF[h]->SetMarkerStyle(h+3);
                        EffMCDistD_BetaNaF[h]->SetMarkerSize(2);
                        EffMCDistD_BetaNaF[h]->SetLineColor(4);
                        EffMCDistD_BetaNaF[h]->SetLineWidth(1);
                        EffMCDistD_BetaNaF[h]->Draw("Psame");
                        leg->Draw();
                }
        }
        c6_bis->cd(3);
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
        TGraph *EffMCDistP_BetaAgl= new TGraph();
        TGraph *EffMCDistD_BetaAgl[6];
        j=0;
        for(int i=0;i<18;i++) {EffMCDistP_BetaAgl->SetPoint(j,EkincentAgl[i],EffDistP_BetaAgl[i]);j++;}
        for(int i=0;i<18;i++) EffMCDistP_BetaAgl_TH1F->SetBinContent(i+1,EffDistP_BetaAgl[i]);
        EffMCDistP_BetaAgl->SetMarkerColor(2);
        EffMCDistP_BetaAgl->SetMarkerStyle(8);
        EffMCDistP_BetaAgl->SetLineColor(2);
        EffMCDistP_BetaAgl->SetLineWidth(2);
        EffMCDistP_BetaAgl->SetTitle("Distelihood Efficiency MC on top of Pres. (Beta bins)");
        EffMCDistP_BetaAgl->GetXaxis()->SetTitle("Kin. En.  [GeV/nucl.]");
        EffMCDistP_BetaAgl->GetYaxis()->SetTitle("Efficiency");
        EffMCDistP_BetaAgl->GetXaxis()->SetTitleSize(0.045);
        EffMCDistP_BetaAgl->GetYaxis()->SetTitleSize(0.045);
        {
                EffMCDistP_BetaAgl->GetYaxis()->SetRangeUser(0,1);
		EffMCDistP_BetaAgl->Draw("ACP");
                TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(EffMCDistP_BetaAgl,MCLegend[0].c_str(), "ep");

                for(int h=0;h<6;h++){
                        EffMCDistD_BetaAgl[h]= new TGraph();
                        EffMCDistD_BetaAgl[h]->SetTitle(MCLegend[h+1].c_str());
                        for(int i=1;i<18;i++) EffMCDistD_BetaAgl[h]->SetPoint(i,EkincentAgl[i],EffDistD_BetaAgl[i][h]);
                        for(int i=1;i<18;i++) EffMCDistD_BetaAgl_TH2F->SetBinContent(i+1,h+1,EffDistD_BetaAgl[i][h]);
                        leg->AddEntry(EffMCDistD_BetaAgl[h],MCLegend[h+1].c_str(), "ep");
                        EffMCDistD_BetaAgl[h]->SetMarkerColor(4);
                        EffMCDistD_BetaAgl[h]->SetMarkerStyle(h+3);
                        EffMCDistD_BetaAgl[h]->SetMarkerSize(2);
                        EffMCDistD_BetaAgl[h]->SetLineColor(4);
                        EffMCDistD_BetaAgl[h]->SetLineWidth(1);
                        EffMCDistD_BetaAgl[h]->Draw("Psame");
                        leg->Draw();
                }
        }

}


