
TH2F * Hecut_D=new TH2F("Hecut_D","Hecut_D",1000,0,40,1000,0,40);
TH2F * HecutMC_P=new TH2F("HecutMC_P","HecutMC_P",1000,0,40,1000,0,40);
TH2F * HecutMC_He=new TH2F("HecutMC_He","HecutMC_He",1000,0,40,1000,0,40);

Efficiency * HecutMCP = new Efficiency("HecutMCP");
Efficiency * HecutMCHe = new Efficiency("HecutMCHe");

void HecutMC_Fill(TNtuple *ntupla,int l){
	 int k = ntupla->GetEvent(l);
	if(Beta<=0||R<=0) return;
	float EdepTOFud=(EdepTOFU+EdepTOFD)/2;
	if(Massa_gen<1){
		HecutMC_P->Fill(fabs(EdepTOFbeta->Eval(Beta)-EdepTOFud)/(pow(EdepTOFbeta->Eval(Beta),2)*etofu->Eval(Beta)),fabs(EdepTrackbeta->Eval(Beta)-EdepTrack)/(pow(EdepTrackbeta->Eval(Beta),2)*etrack->Eval(Beta)));
		for(int K=0;K<43;K++) 
			if(Momento_gen<bin[K+1]&&Momento_gen>bin[K]){
					 HecutMCP->beforeR->Fill(K);
					 if(Herejcut)
						 HecutMCP->afterR->Fill(K);	
				}
	}
	if(Massa_gen>2){
		HecutMC_He->Fill(fabs(EdepTOFbeta->Eval(Beta)-EdepTOFud)/(pow(EdepTOFbeta->Eval(Beta),2)*etofu->Eval(Beta)),fabs(EdepTrackbeta->Eval(Beta)-EdepTrack)/(pow(EdepTrackbeta->Eval(Beta),2)*etrack->Eval(Beta)));	
		for(int K=0;K<43;K++) 
			if(Momento_gen<bin[K+1]&&Momento_gen>bin[K]) {
					HecutMCHe->beforeR->Fill(K);
                			if(Herejcut)
                                        	HecutMCHe->afterR->Fill(K);
			}
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
	HecutMCP->Write();
	HecutMCHe->Write();
}



void Hecut(TFile * file1){
	TH2F* HecutMC_P =(TH2F*)file1->Get("HecutMC_P");
	TH2F* HecutMC_He=(TH2F*)file1->Get("HecutMC_He");
	TH2F* Hecut_D   =(TH2F*)file1->Get("Hecut_D");
	
	Efficiency * HecutMCP = new Efficiency(file1,"HecutMCP");
	Efficiency * HecutMCHe = new Efficiency(file1,"HecutMCHe");

	cout<<"*************** He control sample cut Efficiency on P*******************"<<endl;
	
	HecutMCP->UpdateErrorbars();
	HecutMCHe->UpdateErrorbars();
	
	TH1F * HecutMCP_TH1F = 	(TH1F *)HecutMCP ->afterR->Clone();
        TH1F * HecutMCHe_TH1F=  (TH1F *)HecutMCHe->afterR->Clone();
				
	HecutMCP_TH1F -> Divide( HecutMCP ->beforeR );
	HecutMCHe_TH1F-> Divide( HecutMCHe->beforeR );
	

	cout<<"*** Updating P1 file ****"<<endl;
        string nomefile=percorso + "/Risultati/risultati/"+mese+"_"+frac+"_P1.root";
        file1 =TFile::Open(nomefile.c_str(),"UPDATE");
        if(!file1){
                nomefile=percorso + "/Risultati/"+mese+"/"+mese+"_"+frac+"_P1.root";
                file1 =TFile::Open(nomefile.c_str(),"UPDATE");
        }
	file1->cd("Results");
	HecutMC_P 	->Write();
	HecutMC_He	->Write();
	Hecut_D   	->Write();
	HecutMCP_TH1F 	->Write();
	HecutMCHe_TH1F	->Write();
	
	TCanvas * c36	 =new TCanvas("Sigma E. dep. Track vs TOF");
	TCanvas * c36_bis=new TCanvas("Sigma E. dep. Track vs TOF (MC)");
	TCanvas * c37	=new TCanvas("Eff. He cut");

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
	for(int K=0;K<43;K++) {effHecut->SetPoint(K,R_cent[K], HecutMCP_TH1F->GetBinContent(K+1));
			       effHecut->SetPointError(K,0,    HecutMCP_TH1F->GetBinError(K+1));}
	TGraphErrors *effHecutHe=new TGraphErrors();
        for(int K=0;K<43;K++) {effHecutHe->SetPoint(K,R_cent[K], HecutMCHe_TH1F->GetBinContent(K+1));
                               effHecutHe->SetPointError(K,0,    HecutMCHe_TH1F->GetBinError(K+1));}
	effHecut->SetMarkerColor(2);
	effHecut->SetLineColor(2);
	effHecut->SetMarkerStyle(8);
	effHecutHe->SetMarkerColor(2);
        effHecutHe->SetLineColor(2);
        effHecutHe->SetMarkerStyle(8);
	effHecut->GetXaxis()->SetTitle("R [GV]");
	effHecut->GetYaxis()->SetRangeUser(0.88,1.05);
	effHecut->Draw("AP");			
	//effHecutHe->Draw("Psame");

	cout<<"*** Updating Results file ***"<<endl;
        nomefile=percorso + "/CodesforAnalysis/Final_plots/"+mese+".root";
        TFile *f_out=new TFile(nomefile.c_str(), "UPDATE");
        f_out->mkdir("MC Results/He control sample cut");
        f_out->cd("MC Results/He control sample cut");		
	c36	->Write();	 
 	c36_bis ->Write();	
 	c37	->Write();
	f_out->Write();
	f_out->Close();
}	
