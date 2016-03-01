TCanvas * c36=new TCanvas("Sigma E. dep. Track vs TOF");
TCanvas * c37=new TCanvas("Eff. He cut");

TH2F * HecutMC_P=new TH2F("HecutMC_P","HecutMC_P",1000,0,40,1000,0,40);
TH1F * HecutMC_P1=new TH1F("HecutMC_P1","HecutMC_P1",43,0,43);
TH1F * HecutMC_P2=new TH1F("HecutMC_P2","HecutMC_P2",43,0,43);


TH2F * Hecut_D=new TH2F("Hecut_D","Hecut_D",1000,0,40,1000,0,40);

void HecutMC_Fill(TNtuple *ntupla,int l){
	 int k = ntupla->GetEvent(l);
	if(Beta<=0||R<=0) return;
	float EdepTOFud=(EdepTOFU+EdepTOFD)/2;
	if(Massa_gen<1){
		HecutMC_P->Fill(fabs(EdepTOFbeta->Eval(Beta)-EdepTOFud)/(pow(EdepTOFbeta->Eval(Beta),2)*etofu->Eval(Beta)),fabs(EdepTrackbeta->Eval(Beta)-EdepTrack)/(pow(EdepTrackbeta->Eval(Beta),2)*etrack->Eval(Beta)));
		for(int K=0;K<43;K++) if(R<bin[K+1]&&R>bin[K]) HecutMC_P1->Fill(K);
		if(Herejcut)
					for(int K=0;K<43;K++) if(R<bin[K+1]&&R>bin[K]) HecutMC_P2->Fill(K);	
	}

}

void HecutD_Fill(TNtuple *ntupla,int l){
         int k = ntupla->GetEvent(l);
        if(Beta<=0||R<=0) return;
        float EdepTOFud=(EdepTOFU+EdepTOFD)/2;
                Hecut_D->Fill(fabs(EdepTOFbeta->Eval(Beta)-EdepTOFud)/(pow(EdepTOFbeta->Eval(Beta),2)*etofu->Eval(Beta)),fabs(EdepTrackbeta->Eval(Beta)-EdepTrack)/(pow(EdepTrackbeta->Eval(Beta),2)*etrack->Eval(Beta)));

}


void HecutMC_Copy(TFile * file1){
	HecutMC_P=(TH2F*)file1->Get("HecutMC_P");
	Hecut_D=(TH2F*)file1->Get("Hecut_D");	
	HecutMC_P1=(TH1F*)file1->Get("HecutMC_P1");
	HecutMC_P2=(TH1F*)file1->Get("HecutMC_P2");
}

void HecutMC_Write(){
	HecutMC_P->Write();
	Hecut_D->Write();
	HecutMC_P1->Write();
	HecutMC_P2->Write();
}



void Hecut(TFile * file1){
	float EffHecut[43]={0};
	for(int K=0;K<43;K++) EffHecut[K]=1-HecutMC_P2->GetBinContent(K+1)/HecutMC_P1->GetBinContent(K+1);
	c36->cd();
	gPad->SetLogz();
	gPad->SetGridx();
        gPad->SetGridy();
	HecutMC_P->SetMarkerColor(2);
	Hecut_D->GetXaxis()->SetTitle("E.dep TOF (|meas -teo|). (# of sigmas)");
	Hecut_D->GetYaxis()->SetTitle("E.dep Track (|meas -teo|). (# of sigmas)");
	Hecut_D->Draw("col");
	//HecutMC_P->Draw("same");	
	c37->cd();
	gPad->SetGridx();
	gPad->SetGridy();
	gPad->SetLogx();
	TGraphErrors *effHecut=new TGraphErrors();
	for(int K=0;K<43;K++) {effHecut->SetPoint(K,R_cent[K],EffHecut[K]);
			       effHecut->SetPointError(K,0,(pow(HecutMC_P2->GetBinContent(K+1),0.5)/HecutMC_P2->GetBinContent(K+1))*(1-EffHecut[K]));}
	effHecut->SetMarkerColor(2);
	effHecut->SetLineColor(2);
	effHecut->SetMarkerStyle(8);
	effHecut->GetXaxis()->SetTitle("");
	effHecut->Draw("AP");			

}

