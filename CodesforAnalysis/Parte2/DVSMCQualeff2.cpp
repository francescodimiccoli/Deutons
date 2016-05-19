using namespace std;

TH1F * EffDistMCvsDP1 = new TH1F("EffDistMCvsDP1","EffDistMCvsDP1",nbinsr,0,nbinsr);
TH1F * EffDistMCvsDP2 = new TH1F("EffDistMCvsDP2","EffDistMCvsDP2",nbinsr,0,nbinsr);
TH1F * EffLik2MCvsDP1 = new TH1F("EffLik2MCvsDP1","EffLik2MCvsDP1",nbinsr,0,nbinsr);
TH1F * EffLik2MCvsDP2 = new TH1F("EffLik2MCvsDP2","EffLik2MCvsDP2",nbinsr,0,nbinsr);

TH2F * EffDistMCvsDP_D_1 = new TH2F("EffDistMCvsDP_D_1","EffDistMCvsDP_D_1",nbinsr,0,nbinsr,11,0,11);
TH2F * EffDistMCvsDP_D_2 = new TH2F("EffDistMCvsDP_D_2","EffDistMCvsDP_D_2",nbinsr,0,nbinsr,11,0,11);
TH2F * EffLik2MCvsDP_D_1 = new TH2F("EffLik2MCvsDP_D_1","EffLik2MCvsDP_D_1",nbinsr,0,nbinsr,11,0,11);
TH2F * EffLik2MCvsDP_D_2 = new TH2F("EffLik2MCvsDP_D_2","EffLik2MCvsDP_D_2",nbinsr,0,nbinsr,11,0,11);


void DVSMCQualeff2_D_Fill(TNtuple *ntupla, int l,int zona){

	 ntupla->GetEvent(l);
	if(Beta<=0||R<=0||R<1.2*Rcutoff||Beta>protons->Eval(R)+0.1||Beta<protons->Eval(R)-0.1) return;
	if((R>Rcut[zona]&&zona<10)||(zona==10)) {
		if(Herejcut){
			//if(Likcut){
			for(int K=0;K<nbinsr;K++) if(R<bin[K+1]&&R>bin[K]) {EffDistMCvsDP_D_2->Fill(K,zona);}
			        if(Dist5D_P<6) for(int K=0;K<nbinsr;K++) if(R<bin[K+1]&&R>bin[K]) {EffDistMCvsDP_D_1->Fill(K,zona);}
			//}
			if(Dist5D_P<6) {
				for(int K=0;K<nbinsr;K++) if(R<bin[K+1]&&R>bin[K]) {EffLik2MCvsDP_D_1->Fill(K,zona);}
				if(Likcut)
					for(int K=0;K<nbinsr;K++) if(R<bin[K+1]&&R>bin[K]) {EffLik2MCvsDP_D_2->Fill(K,zona);}
			}

		}
	}
}

void DVSMCQualeff2_Fill(TNtuple *ntupla, int l){
                                
	 ntupla->GetEvent(l);    
	if(Beta<=0||R<=0||Beta>protons->Eval(R)+0.1||Beta<protons->Eval(R)-0.1) return;
	if(Massa_gen<1){
	if(Herejcut){
		//if(Likcut){
			for(int K=0;K<nbinsr;K++) if(R<bin[K+1]&&R>bin[K]) {EffDistMCvsDP2->Fill(K);}
			if(Dist5D_P<6) for(int K=0;K<nbinsr;K++) if(R<bin[K+1]&&R>bin[K]) {EffDistMCvsDP1->Fill(K);}
		//}
		if(Dist5D_P<6) {
			for(int K=0;K<nbinsr;K++) if(R<bin[K+1]&&R>bin[K]) {EffLik2MCvsDP1->Fill(K);}
			if(Likcut)
				for(int K=0;K<nbinsr;K++) if(R<bin[K+1]&&R>bin[K]) {EffLik2MCvsDP2->Fill(K);}
		}

	}
	}
}

void DVSMCQualeff2_Write(){
        EffDistMCvsDP_D_1->Write(); 
        EffDistMCvsDP_D_2->Write();
        EffLik2MCvsDP_D_1->Write();
        EffLik2MCvsDP_D_2->Write();
 	EffDistMCvsDP1->Write();
        EffDistMCvsDP2->Write();
        EffLik2MCvsDP1->Write();
        EffLik2MCvsDP2->Write();
        return;
}





TCanvas *c20=new TCanvas("Data vs MC: Likelihood");
TCanvas *c21=new TCanvas("Data vs MC: Distance");
TGraphErrors *LikDVSMC_P_Graph;
TGraphErrors *DistDVSMC_P_Graph;
TH2F * LikDVSMC_P_graph=new TH2F("LikDVSMC_P_graph","LikDVSMC_P_graph",nbinsr,0,nbinsr,2,0,2);
TH2F * DistDVSMC_P_graph=new TH2F("DistDVSMC_P_graph","DistDVSMC_P_graph",nbinsr,0,nbinsr,2,0,2);


void DVSMCQualeff2(TFile * file1){
	TH2F * EffDistMCvsDP_D_1 = (TH2F *)file1->Get("EffDistMCvsDP_D_1");
	TH2F * EffDistMCvsDP_D_2 = (TH2F *)file1->Get("EffDistMCvsDP_D_2");
	TH2F * EffLik2MCvsDP_D_1 = (TH2F *)file1->Get("EffLik2MCvsDP_D_1");
	TH2F * EffLik2MCvsDP_D_2 = (TH2F *)file1->Get("EffLik2MCvsDP_D_2");
	TH1F * EffDistMCvsDP1 = (TH1F *)file1->Get("EffDistMCvsDP1");
	TH1F * EffDistMCvsDP2 = (TH1F *)file1->Get("EffDistMCvsDP2");
	TH1F * EffLik2MCvsDP1 = (TH1F *)file1->Get("EffLik2MCvsDP1");
	TH1F * EffLik2MCvsDP2 = (TH1F *)file1->Get("EffLik2MCvsDP2");

	cout<<"************************************************************ MC QUALITY SEL. EFFICIENCIES ************************************************************"<<endl;		
	float EffLik2MCvsDP[nbinsr][11]={{0}};
	for(int l=0;l<11;l++)
	for(int i=1;i<nbinsr;i++) if(EffLik2MCvsDP_D_2->GetBinContent(i+1,l+1)>100&&EffLik2MCvsDP_D_1->GetBinContent(i+1,l+1)>100)
		EffLik2MCvsDP[i][l]=EffLik2MCvsDP_D_2->GetBinContent(i+1,l+1)/(float)EffLik2MCvsDP_D_1->GetBinContent(i+1,l+1)*CorrLAT_Lik->Eval(geomagC[l]);

	float EffDistMCvsDP[nbinsr][11]={{0}};
	for(int l=0;l<11;l++)
	for(int i=1;i<nbinsr;i++) if(EffDistMCvsDP_D_2->GetBinContent(i+1,l+1)>100&&EffDistMCvsDP_D_1->GetBinContent(i+1,l+1)>100)
		EffDistMCvsDP[i][l]=EffDistMCvsDP_D_1->GetBinContent(i+1,l+1)/(float)EffDistMCvsDP_D_2->GetBinContent(i+1,l+1)*CorrLAT_Dist->Eval(geomagC[l]);
	
	float EffLik2MCvsDP_MC[nbinsr]={0};
        for(int i=1;i<nbinsr;i++) if(EffLik2MCvsDP2->GetBinContent(i+1)>100&&EffLik2MCvsDP1->GetBinContent(i+1)>100)
                EffLik2MCvsDP_MC[i]=EffLik2MCvsDP2->GetBinContent(i+1)/(float)EffLik2MCvsDP1->GetBinContent(i+1);

        float EffDistMCvsDP_MC[nbinsr]={0};
        for(int i=1;i<nbinsr;i++) if(EffDistMCvsDP2->GetBinContent(i+1)>100&&EffDistMCvsDP1->GetBinContent(i+1)>100)
                EffDistMCvsDP_MC[i]=EffDistMCvsDP1->GetBinContent(i+1)/(float)EffDistMCvsDP2->GetBinContent(i+1);

	//WEIGHTED MEAN over LAT
	float EffDistMCvsDP_mean[nbinsr]={0};
	float EffLik2MCvsDP_mean[nbinsr]={0};
	float EffDistMCvsDP_meanerr[nbinsr]={0};
        float EffLik2MCvsDP_meanerr[nbinsr]={0};
	TGraphErrors * EffDistMCvsDP_Mean=new TGraphErrors();
	TGraphErrors * EffLik2MCvsDP_Mean=new TGraphErrors();
	float denom=0;
	int p=0;
	for(int i=1;i<nbinsr;i++) {
		for(int l=0;l<11;l++) {
			if(EffLik2MCvsDP_D_2->GetBinContent(i+1,l+1)>100&&EffLik2MCvsDP_D_1->GetBinContent(i+1,l+1)>100){
				denom=pow(EffLik2MCvsDP_D_1->GetBinContent(i+1,l+1),0.5)/EffLik2MCvsDP_D_1->GetBinContent(i+1,l+1)*EffLik2MCvsDP[i][l]/EffLik2MCvsDP_MC[i];
				EffLik2MCvsDP_mean[i]+=(EffLik2MCvsDP[i][l]/EffLik2MCvsDP_MC[i])/pow(denom,2);
				EffLik2MCvsDP_meanerr[i]+=1/pow(denom,2);
			}
		}
		EffLik2MCvsDP_mean[i]=EffLik2MCvsDP_mean[i]/EffLik2MCvsDP_meanerr[i];
		EffLik2MCvsDP_meanerr[i]=pow(1/EffLik2MCvsDP_meanerr[i],0.5);
		if(EffLik2MCvsDP_mean[i]>0&&EffLik2MCvsDP_meanerr[i]>0){
		EffLik2MCvsDP_Mean->SetPoint(p,R_cent[i],EffLik2MCvsDP_mean[i]);
		EffLik2MCvsDP_Mean->SetPointError(p,0,EffLik2MCvsDP_meanerr[i]);
		p++;
		}
	}
	
	denom=0;
	p=0;
        for(int i=1;i<nbinsr;i++) {
                for(int l=0;l<11;l++) {
                        if(EffDistMCvsDP_D_2->GetBinContent(i+1,l+1)>100&&EffDistMCvsDP_D_1->GetBinContent(i+1,l+1)>100){
                                denom=pow(EffDistMCvsDP_D_1->GetBinContent(i+1,l+1),0.5)/EffDistMCvsDP_D_1->GetBinContent(i+1,l+1)*EffDistMCvsDP[i][l]/EffDistMCvsDP_MC[i];
                                EffDistMCvsDP_mean[i]+=(EffDistMCvsDP[i][l]/EffDistMCvsDP_MC[i])/pow(denom,2);
                                EffDistMCvsDP_meanerr[i]+=1/pow(denom,2);
                        }
                }
                EffDistMCvsDP_mean[i]=EffDistMCvsDP_mean[i]/EffDistMCvsDP_meanerr[i];
                EffDistMCvsDP_meanerr[i]=pow(1/EffDistMCvsDP_meanerr[i],0.5);
                if(EffDistMCvsDP_mean[i]>0&&EffDistMCvsDP_meanerr[i]>0){
		EffDistMCvsDP_Mean->SetPoint(p,R_cent[i],EffDistMCvsDP_mean[i]);
                EffDistMCvsDP_Mean->SetPointError(p,0,EffDistMCvsDP_meanerr[i]);
		p++;
		}
        }
	//SMOOTHING & FIT
	float ratioL_smooth[nbinsr]={0};
	float ratioD_smooth[nbinsr]={0};

	for(int j=1;j<nbinsr;j++){
		ratioL_smooth[j]=EffLik2MCvsDP_mean[j];
		ratioD_smooth[j]=EffDistMCvsDP_mean[j];
	}
	ratioL_smooth[0]=ratioL_smooth[1];
	ratioD_smooth[0]=ratioD_smooth[1];
	
        TGraphErrors *EffMCvsDLikP_MC=new TGraphErrors();
        int j=0;
        for(int i=1;i<nbinsr;i++) if(EffLik2MCvsDP_MC[i]>0)
                {EffMCvsDLikP_MC->SetPoint(j,R_cent[i],EffLik2MCvsDP_MC[i]/EffLik2MCvsDP_MC[i]);j++;}
	
        TGraphErrors *EffMCvsDDistP_MC= new TGraphErrors();
        j=0;
        for(int i=1;i<nbinsr;i++) if(EffDistMCvsDP_MC[i]>0)
                {EffMCvsDDistP_MC->SetPoint(j,R_cent[i],EffDistMCvsDP_MC[i]/EffDistMCvsDP_MC[i]);j++;}

	TF1 *polL=new TF1("polL","pol3");
	TF1 *polD=new TF1("polD","pol3");
	
	EffLik2MCvsDP_Mean->Fit("polL","","",8,70);
	EffDistMCvsDP_Mean->Fit("polD","","",8,70);

	for(int j=3;j<nbinsr;j++){
		if(R_cent[j]<8){
		ratioL_smooth[j]=(ratioL_smooth[j]+ratioL_smooth[j-1]+ratioL_smooth[j+1])/3;
		ratioD_smooth[j]=(ratioD_smooth[j]+ratioD_smooth[j-1]+ratioD_smooth[j+1])/3;
		}
		else{
		ratioL_smooth[j]=polL->Eval(R_cent[j]);
		ratioD_smooth[j]=polD->Eval(R_cent[j]);
		}
	}
	///errore fit
        float errorefitL=FitError(ratioL_smooth,EffLik2MCvsDP_mean,EffLik2MCvsDP_meanerr,nbinsr,3);
        float errorefitD=FitError(ratioD_smooth,EffDistMCvsDP_mean,EffDistMCvsDP_meanerr,nbinsr,3);
	///////////
	

	c20->cd();
	gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
	LikDVSMC_P_Graph=new TGraphErrors();
	LikDVSMC_P_Graph->SetName("LikDVSMC_P_Graph");
	j=0;
	
	for(int i=1;i<nbinsr;i++) {
		if(EffLik2MCvsDP_MC[i]>0){
		LikDVSMC_P_Graph->SetPoint(j,R_cent[i],ratioL_smooth[i]);
		LikDVSMC_P_graph->SetBinContent(i+1,1,ratioL_smooth[i]);
		LikDVSMC_P_Graph->SetPointError(j,0,errorefitL);
		LikDVSMC_P_graph->SetBinContent(i+1,2,errorefitL);
	j++;
	}
	}
	LikDVSMC_P_Graph->SetLineColor(4);
	LikDVSMC_P_Graph->SetFillColor(4);
	LikDVSMC_P_Graph->SetFillStyle(3001);
        LikDVSMC_P_Graph->SetLineWidth(4);
	EffLik2MCvsDP_Mean->SetMarkerColor(2);
        EffLik2MCvsDP_Mean->SetMarkerStyle(4);
        EffLik2MCvsDP_Mean->SetLineColor(2);
        EffLik2MCvsDP_Mean->SetLineWidth(2);
        EffMCvsDLikP_MC->SetMarkerColor(2);
        EffMCvsDLikP_MC->SetMarkerStyle(8);
        EffMCvsDLikP_MC->SetLineColor(2);
        EffMCvsDLikP_MC->SetLineWidth(2);
	EffLik2MCvsDP_Mean->SetTitle("Likelihood Efficiency DATA/MC");
        EffLik2MCvsDP_Mean->GetXaxis()->SetTitle("R [GV]");
        EffLik2MCvsDP_Mean->GetYaxis()->SetTitle("Efficiency");
        EffLik2MCvsDP_Mean->GetXaxis()->SetTitleSize(0.045);
        EffLik2MCvsDP_Mean->GetYaxis()->SetTitleSize(0.045);
	EffLik2MCvsDP_Mean->GetYaxis()->SetRangeUser(0.7,1.4);
	{
                EffLik2MCvsDP_Mean->Draw("AP");
		EffMCvsDLikP_MC->Draw("CPsame");
                LikDVSMC_P_Graph->Draw("L4same");
		TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                //leg->AddEntry(EffMCLikP,MCLegend[0].c_str(), "ep");

        }

	c21->cd();
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
	DistDVSMC_P_Graph=new TGraphErrors();
        DistDVSMC_P_Graph->SetName("DistDVSMC_P_Graph");
	j=0;
        for(int i=1;i<nbinsr;i++) {
                if(EffLik2MCvsDP_MC[i]>0){
                DistDVSMC_P_Graph->SetPoint(j,R_cent[i],ratioD_smooth[i]);
               	DistDVSMC_P_graph->SetBinContent(i+1,1,ratioD_smooth[i]);
		DistDVSMC_P_Graph->SetPointError(j,0,errorefitD);
        	DistDVSMC_P_graph->SetBinContent(i+1,2,errorefitD);
	j++;
        }
        }
        DistDVSMC_P_Graph->SetLineColor(4);
        DistDVSMC_P_Graph->SetFillColor(4);
        DistDVSMC_P_Graph->SetFillStyle(3001);
        DistDVSMC_P_Graph->SetLineWidth(4);
	EffDistMCvsDP_Mean->SetMarkerColor(2);
        EffDistMCvsDP_Mean->SetMarkerStyle(4);
        EffDistMCvsDP_Mean->SetLineColor(2);
        EffDistMCvsDP_Mean->SetLineWidth(2);
	EffMCvsDDistP_MC->SetMarkerColor(2);
        EffMCvsDDistP_MC->SetMarkerStyle(8);
        EffMCvsDDistP_MC->SetLineColor(2);
        EffMCvsDDistP_MC->SetLineWidth(2);
        EffDistMCvsDP_Mean->SetTitle("Distance Efficiency DATA/MC");
        EffDistMCvsDP_Mean->GetXaxis()->SetTitle("R [GV]");
        EffDistMCvsDP_Mean->GetYaxis()->SetTitle("Efficiency");
        EffDistMCvsDP_Mean->GetXaxis()->SetTitleSize(0.045);
        EffDistMCvsDP_Mean->GetYaxis()->SetTitleSize(0.045);
	EffDistMCvsDP_Mean->GetYaxis()->SetRangeUser(0.7,1.4);
        {
                EffDistMCvsDP_Mean->Draw("AP");
		EffMCvsDDistP_MC->Draw("CPsame");
                DistDVSMC_P_Graph->Draw("L4same");
		TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                //leg->AddEntry(EffMCDistP,MCLegend[0].c_str(), "ep");

        }	


}

