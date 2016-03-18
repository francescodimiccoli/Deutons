TCanvas * c36=new TCanvas("Sigma E. dep. Track vs TOF");
TCanvas * c36_bis=new TCanvas("Sigma E. dep. Track vs TOF (MC)");

TCanvas * c37=new TCanvas("Eff. He cut");

TH2F * HecutMC_P=new TH2F("HecutMC_P","HecutMC_P",1000,0,40,1000,0,40);
TH1F * HecutMC_P1=new TH1F("HecutMC_P1","HecutMC_P1",43,0,43);
TH1F * HecutMC_P2=new TH1F("HecutMC_P2","HecutMC_P2",43,0,43);
TH2F * HecutMC_He=new TH2F("HecutMC_He","HecutMC_He",1000,0,40,1000,0,40);
TH1F * HecutMC_He1=new TH1F("HecutMC_He1","HecutMC_He1",43,0,43);
TH1F * HecutMC_He2=new TH1F("HecutMC_He2","HecutMC_He2",43,0,43);


TH2F * Hecut_D=new TH2F("Hecut_D","Hecut_D",1000,0,40,1000,0,40);

void HecutMC_Fill(TNtuple *ntupla,int l){
	 int k = ntupla->GetEvent(l);
	if(Beta<=0||R<=0) return;
	float EdepTOFud=(EdepTOFU+EdepTOFD)/2;
	if(Massa_gen<1){
		HecutMC_P->Fill(fabs(EdepTOFbeta->Eval(Beta)-EdepTOFud)/(pow(EdepTOFbeta->Eval(Beta),2)*etofu->Eval(Beta)),fabs(EdepTrackbeta->Eval(Beta)-EdepTrack)/(pow(EdepTrackbeta->Eval(Beta),2)*etrack->Eval(Beta)));
		for(int K=0;K<43;K++) if(Momento_gen<bin[K+1]&&Momento_gen>bin[K]) HecutMC_P1->Fill(K);
		if(Herejcut)
					for(int K=0;K<43;K++) if(Momento_gen<bin[K+1]&&Momento_gen>bin[K]) HecutMC_P2->Fill(K);	
	}
	if(Massa_gen>2){
		HecutMC_He->Fill(fabs(EdepTOFbeta->Eval(Beta)-EdepTOFud)/(pow(EdepTOFbeta->Eval(Beta),2)*etofu->Eval(Beta)),fabs(EdepTrackbeta->Eval(Beta)-EdepTrack)/(pow(EdepTrackbeta->Eval(Beta),2)*etrack->Eval(Beta)));	
		for(int K=0;K<43;K++) if(Momento_gen<bin[K+1]&&Momento_gen>bin[K]) HecutMC_He1->Fill(K);
                if(Herejcut)
                                        for(int K=0;K<43;K++) if(Momento_gen<bin[K+1]&&Momento_gen>bin[K]) HecutMC_He2->Fill(K);
	}

}

void HecutD_Fill(TNtuple *ntupla,int l){
         int k = ntupla->GetEvent(l);
        if(Beta<=0||R<=0) return;
        float EdepTOFud=(EdepTOFU+EdepTOFD)/2;
                Hecut_D->Fill(fabs(EdepTOFbeta->Eval(Beta)-EdepTOFud)/(pow(EdepTOFbeta->Eval(Beta),2)*etofu->Eval(Beta)),fabs(EdepTrackbeta->Eval(Beta)-EdepTrack)/(pow(EdepTrackbeta->Eval(Beta),2)*etrack->Eval(Beta)));

}


void HecutMC_Write(){
	HecutMC_P->Write();
	HecutMC_He->Write();
	Hecut_D->Write();
	HecutMC_P1->Write();
	HecutMC_P2->Write();
	HecutMC_He1->Write();
        HecutMC_He2->Write();
}



void Hecut(TFile * file1){
	TH2F* HecutMC_P=(TH2F*)file1->Get("HecutMC_P");
	TH2F* HecutMC_He=(TH2F*)file1->Get("HecutMC_He");
	TH2F* Hecut_D=(TH2F*)file1->Get("Hecut_D");
	TH1F*  HecutMC_P1=(TH1F*)file1->Get("HecutMC_P1");
	TH1F* HecutMC_P2=(TH1F*)file1->Get("HecutMC_P2");
	TH1F* HecutMC_He1=(TH1F*)file1->Get("HecutMC_He1");
	TH1F* HecutMC_He2=(TH1F*)file1->Get("HecutMC_He2");

	cout<<"*************** He control sample cut *******************"<<endl;
	float EffHecut[43]={0};
	for(int K=0;K<43;K++) if(HecutMC_P2->GetBinContent(K+1)<HecutMC_P1->GetBinContent(K+1)) 
					EffHecut[K]=HecutMC_P2->GetBinContent(K+1)/HecutMC_P1->GetBinContent(K+1);
	float EffHecutHe[43]={0};
        for(int K=0;K<43;K++) if(HecutMC_He2->GetBinContent(K+1)<HecutMC_He1->GetBinContent(K+1))
				EffHecutHe[K]=HecutMC_He2->GetBinContent(K+1)/HecutMC_He1->GetBinContent(K+1);
	
	c36->cd();
	gPad->SetLogz();
	gPad->SetGridx();
        gPad->SetGridy();
	HecutMC_P->SetMarkerColor(2);
	Hecut_D->GetXaxis()->SetTitle("E.dep TOF (|meas -teo|). (# of sigmas)");
	Hecut_D->GetYaxis()->SetTitle("E.dep Track (|meas -teo|). (# of sigmas)");
	Hecut_D->Draw("col");
	
	c36_bis->cd();
	gPad->SetLogz();
        gPad->SetGridx();
        gPad->SetGridy();
	HecutMC_P->SetMarkerColor(2);
	HecutMC_He->SetMarkerColor(3);
	HecutMC_He->GetXaxis()->SetTitle("E.dep TOF (|meas -teo|). (# of sigmas)");
        HecutMC_He->GetYaxis()->SetTitle("E.dep Track (|meas -teo|). (# of sigmas)");
	HecutMC_He->Draw("col");
	//HecutMC_P->Draw("same");	
	c37->cd();
	gPad->SetGridx();
	gPad->SetGridy();
	gPad->SetLogx();
	TGraphErrors *effHecut=new TGraphErrors();
	for(int K=0;K<43;K++) {effHecut->SetPoint(K,R_cent[K],EffHecut[K]);
			       effHecut->SetPointError(K,0,(pow(HecutMC_P2->GetBinContent(K+1),0.5)/HecutMC_P2->GetBinContent(K+1))*(EffHecut[K]));}
	TGraphErrors *effHecutHe=new TGraphErrors();
        for(int K=0;K<43;K++) {effHecutHe->SetPoint(K,R_cent[K],EffHecutHe[K]);
                               effHecutHe->SetPointError(K,0,(pow(HecutMC_He2->GetBinContent(K+1),0.5)/HecutMC_He2->GetBinContent(K+1))*(EffHecutHe[K]));}
	effHecut->SetMarkerColor(2);
	effHecut->SetLineColor(2);
	effHecut->SetMarkerStyle(8);
	effHecutHe->SetMarkerColor(2);
        effHecutHe->SetLineColor(2);
        effHecutHe->SetMarkerStyle(8);
	effHecut->GetXaxis()->SetTitle("R [GV]");
	effHecut->GetYaxis()->SetRangeUser(0.88,1.05);
	effHecut->Draw("AP");			
	effHecutHe->Draw("Psame");

}

