
TH2F * Hecut_D=new TH2F("Hecut_D","Hecut_D",1000,0,40,1000,0,40);
TH2F * HecutMC_P=new TH2F("HecutMC_P","HecutMC_P",1000,0,40,1000,0,40);
TH2F * HecutMC_He=new TH2F("HecutMC_He","HecutMC_He",1000,0,40,1000,0,40);

Efficiency * HecutMCP = new Efficiency("HecutMCP");
Efficiency * HecutMCHe = new Efficiency("HecutMCHe");

TemplateFIT * HeliumContaminationTOF = new TemplateFIT("HeliumContaminationTOF",1,0,60);
TemplateFIT * HeliumContaminationNaF = new TemplateFIT("HeliumContaminationNaF",1,0,60);
TemplateFIT * HeliumContaminationAgl = new TemplateFIT("HeliumContaminationAgl",1,0,60);

void HecutMC_Fill(TNtuple *ntupla,int l) {
	 ntupla->GetEvent(l);
	if(Tup.Beta<=0||Tup.R<=0) return;
	float EdepTOFud=(Tup.EdepTOFU+Tup.EdepTOFD)/2;
	int Kbin=RB.GetRBin(Tup.Momento_gen);

	if(Massa_gen<1) {
		HecutMC_P->Fill( fabs(EdepTOFbeta  ->Eval(Tup.Beta)-EdepTOFud) / (pow(EdepTOFbeta  ->Eval(Tup.Beta),2)*etofu ->Eval(Tup.Beta)) ,
				fabs(EdepTrackbeta->Eval(Tup.Beta)-Tup.EdepTrack) / (pow(EdepTrackbeta->Eval(Tup.Beta),2)*etrack->Eval(Tup.Beta)) );
		HecutMCP->beforeR->Fill(Kbin);
		if(Herejcut) HecutMCP->afterR->Fill(Kbin);

		if(Betastrongcut&&Likcut){		
			if(cmask.isOnlyFromToF()&&Tup.R<3)    HeliumContaminationTOF -> TemplateP -> Fill(Tup.Dist5D_P,0);
			if(cmask.isFromNaF()&&Tup.R<6)		  HeliumContaminationNaF -> TemplateP -> Fill(Tup.Dist5D_P,0);	
			if(cmask.isFromAgl()&&Tup.R<14)  	  HeliumContaminationAgl -> TemplateP -> Fill(Tup.Dist5D_P,0);
		}	
	}
	if(Massa_gen>1&&Massa_gen<2){
		if(Betastrongcut&&Likcut){
			if(cmask.isOnlyFromToF()&&Tup.R<3)    HeliumContaminationTOF -> TemplateD -> Fill(Tup.Dist5D_P,0);
			if(cmask.isFromNaF()&&Tup.R<6)		  HeliumContaminationNaF -> TemplateD -> Fill(Tup.Dist5D_P,0);
			if(cmask.isFromAgl()&&Tup.R<14) 		  HeliumContaminationAgl -> TemplateD -> Fill(Tup.Dist5D_P,0);
		}
		
	}
	if(Massa_gen>2) {
		HecutMC_He->Fill(fabs(EdepTOFbeta  ->Eval(Tup.Beta)-EdepTOFud)/(pow(EdepTOFbeta  ->Eval(Tup.Beta),2)*etofu ->Eval(Tup.Beta)),
				fabs(EdepTrackbeta->Eval(Tup.Beta)-Tup.EdepTrack)/(pow(EdepTrackbeta->Eval(Tup.Beta),2)*etrack->Eval(Tup.Beta)));
		HecutMCHe->beforeR->Fill(Kbin);
		if(Herejcut) HecutMCHe->afterR->Fill(Kbin);

		if(Betastrongcut&&Likcut){
			if(cmask.isOnlyFromToF()&&Tup.R<3)    HeliumContaminationTOF -> TemplateHe -> Fill(Tup.Dist5D_P,0);
			if(cmask.isFromNaF()&&Tup.R<6)		  HeliumContaminationNaF -> TemplateHe -> Fill(Tup.Dist5D_P,0);
			if(cmask.isFromAgl()&&Tup.R<14) 		  HeliumContaminationAgl -> TemplateHe -> Fill(Tup.Dist5D_P,0);
		}
	}
}

void HecutD_Fill(TNtuple *ntupla,int l) {
	 ntupla->GetEvent(l);
	if(Tup.Beta<=0||Tup.R<=0) return;
	if(!Tup.R>1.2*Tup.Rcutoff) return;
	float EdepTOFud=(Tup.EdepTOFU+Tup.EdepTOFD)/2;
	Hecut_D->Fill(fabs(EdepTOFbeta->Eval(Tup.Beta)-EdepTOFud)/(pow(EdepTOFbeta->Eval(Tup.Beta),2)*etofu->Eval(Tup.Beta)),fabs(EdepTrackbeta->Eval(Tup.Beta)-Tup.EdepTrack)/(pow(EdepTrackbeta->Eval(Tup.Beta),2)*etrack->Eval(Tup.Beta)));

	if(Betastrongcut&&Likcut){
		if(cmask.isOnlyFromToF()&&Tup.R<3)  HeliumContaminationTOF -> DATA -> Fill(Tup.Dist5D_P,0);
		if(cmask.isFromNaF()&&Tup.R<6)		HeliumContaminationNaF -> DATA -> Fill(Tup.Dist5D_P,0);
		if(cmask.isFromAgl()&&Tup.R<14)  	HeliumContaminationAgl -> DATA -> Fill(Tup.Dist5D_P,0);	
	}
}


void HecutMC_Write() {
	HecutMC_P->Write();
	HecutMC_He->Write();
	Hecut_D->Write();
	HecutMCP->Write();
	HecutMCHe->Write();
	HeliumContaminationTOF -> Write();
	HeliumContaminationNaF -> Write();	
	HeliumContaminationAgl -> Write();
}



void Hecut(TFile * file1) {
	TH2F* HecutMC_P =(TH2F*)file1->Get("HecutMC_P");
	TH2F* HecutMC_He=(TH2F*)file1->Get("HecutMC_He");
	TH2F* Hecut_D   =(TH2F*)file1->Get("Hecut_D");

	Efficiency * HecutMCP = new Efficiency(file1,"HecutMCP");
	Efficiency * HecutMCHe = new Efficiency(file1,"HecutMCHe");

	TemplateFIT * HeliumContaminationTOF = new TemplateFIT(file1,"HeliumContaminationTOF","HeliumContaminationTOF");
	TemplateFIT * HeliumContaminationNaF = new TemplateFIT(file1,"HeliumContaminationNaF","HeliumContaminationNaF");
	TemplateFIT * HeliumContaminationAgl = new TemplateFIT(file1,"HeliumContaminationAgl","HeliumContaminationAgl");

	cout<<"*************** He control sample cut Efficiency on P*******************"<<endl;

	HecutMCP->Eval_Efficiency();
	HecutMCHe->Eval_Efficiency();

	TH1F * HecutMCP_TH1F = 	(TH1F *)HecutMCP ->effR->Clone();
	TH1F * HecutMCHe_TH1F=  (TH1F *)HecutMCHe->effR->Clone();

	cout<<"*************** He Contamination ******************"<<endl;

	 HeliumContaminationTOF -> SetFitConstraints(0.0,1,0.0,1,0.00,1);
         HeliumContaminationNaF -> SetFitConstraints(0.0,1,0.0,1,0.00,1);
         HeliumContaminationAgl -> SetFitConstraints(0.0,1,0.0,1,0.00,1);
	
	HeliumContaminationTOF-> TemplateFits();
	HeliumContaminationNaF-> TemplateFits();
	HeliumContaminationAgl-> TemplateFits();

	float HeCountsTOF= HeliumContaminationTOF->GetResult_He(0)->Integral(0,20);//HeliumContaminationTOF->GetResult_He(0)->FindBin(4));
	float PCountsTOF = HeliumContaminationTOF->GetResult_P(0)->Integral(0,20);//HeliumContaminationTOF->GetResult_P(0)->FindBin(4)); 
	float DCountsTOF = HeliumContaminationTOF->GetResult_D(0)->Integral(0,20);
	float HeCont_TOF = HeCountsTOF/(PCountsTOF+DCountsTOF);

	float HeCountsNaF= HeliumContaminationNaF->GetResult_He(0)->Integral(0,20);//HeliumContaminationNaF->GetResult_He(0)->FindBin(4));
        float PCountsNaF = HeliumContaminationNaF->GetResult_P(0)->Integral(0,20);//HeliumContaminationNaF->GetResult_P(0)->FindBin(4));
        float DCountsNaF = HeliumContaminationNaF->GetResult_D(0)->Integral(0,20);
	float HeCont_NaF = HeCountsNaF/(PCountsNaF+DCountsNaF);
	
	float HeCountsAgl= HeliumContaminationAgl->GetResult_He(0)->Integral(0,20);//HeliumContaminationAgl->GetResult_He(0)->FindBin(4));
        float PCountsAgl = HeliumContaminationAgl->GetResult_P(0)->Integral(0,20);//HeliumContaminationAgl->GetResult_P(0)->FindBin(4));
        float DCountsAgl = HeliumContaminationNaF->GetResult_D(0)->Integral(0,20);
	float HeCont_Agl = HeCountsAgl/(PCountsAgl+DCountsAgl);	

	cout<<"*** Updating P1 file ****"<<endl;
	string nomefile="../Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
	file1 =TFile::Open(nomefile.c_str(),"UPDATE");

	file1->cd("Results");
	HecutMC_P 	->Write();
	HecutMC_He	->Write();
	Hecut_D   	->Write();
	HecutMCP_TH1F 	->Write();
	HecutMCHe_TH1F	->Write();
	file1->Write();
	file1->Close();

	TCanvas * c36	=new TCanvas("Sigma E. dep. Track vs TOF");
	TCanvas * c36_bis    =new TCanvas("Sigma E. dep. Track vs TOF (MC)");
	TCanvas * c37	=new TCanvas("Eff. He cut");
	TCanvas * c38 	=new TCanvas("He fragm. contamination");

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
	for(int K=0; K<nbinsr; K++) {
		effHecut->SetPoint(K,RB.RigBinCent(K), HecutMCP_TH1F->GetBinContent(K+1));
		effHecut->SetPointError(K,0,    HecutMCP_TH1F->GetBinError(K+1));
	}
	TGraphErrors *effHecutHe=new TGraphErrors();
	for(int K=0; K<nbinsr; K++) {
		effHecutHe->SetPoint(K,RB.RigBinCent(K), HecutMCHe_TH1F->GetBinContent(K+1));
		effHecutHe->SetPointError(K,0,    HecutMCHe_TH1F->GetBinError(K+1));
	}
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

	c38 -> Divide(3,1);
	c38 -> cd(1);
	gPad -> SetGridx();
        gPad -> SetGridy();
        gPad -> SetLogy();
	TH1F * P_TOF = 	HeliumContaminationTOF   -> GetResult_P(0);  
	TH1F * D_TOF =  HeliumContaminationTOF   -> GetResult_D(0);  
	TH1F * He_TOF=  HeliumContaminationTOF   -> GetResult_He(0); 
	
	HeliumContaminationTOF   -> GetResult_Data(0) -> GetXaxis() -> SetTitle("Distance from D");
	HeliumContaminationTOF   -> GetResult_Data(0) -> SetTitle(("He fragm. cont. : " + to_string(HeCont_TOF)).c_str());	
	P_TOF -> SetFillColor(2);
	D_TOF -> SetFillColor(4);
	He_TOF-> SetFillColor(3);
	HeliumContaminationTOF   -> GetResult_Data(0) -> SetMarkerStyle(8);
	
	P_TOF  -> SetFillStyle(3001);
        D_TOF  -> SetFillStyle(3001);
        He_TOF -> SetFillStyle(3001);

	HeliumContaminationTOF   -> GetResult_Data(0) -> Draw("P");
	P_TOF -> Draw("same");
        D_TOF  -> Draw("same");
        He_TOF -> Draw("same");

	c38 -> cd(2);
	gPad -> SetGridx();
        gPad -> SetGridy();
        gPad -> SetLogy();
	TH1F * P_NaF = 	HeliumContaminationNaF   -> GetResult_P(0); 
        TH1F * D_NaF =  HeliumContaminationNaF   -> GetResult_D(0); 
	TH1F * He_NaF=  HeliumContaminationNaF   -> GetResult_He(0);
	
	HeliumContaminationNaF   -> GetResult_Data(0) -> GetXaxis() -> SetTitle("Distance from D");
	HeliumContaminationNaF   -> GetResult_Data(0) -> SetTitle(("He fragm. cont. : " + to_string(HeCont_NaF)).c_str());
	P_NaF  -> SetFillColor(2);
        D_NaF  -> SetFillColor(4);
        He_NaF -> SetFillColor(3);
        HeliumContaminationNaF   -> GetResult_Data(0) -> SetMarkerStyle(8);

        P_NaF  -> SetFillStyle(3001);
        D_NaF  -> SetFillStyle(3001);
        He_NaF -> SetFillStyle(3001);

        HeliumContaminationNaF   -> GetResult_Data(0) -> Draw("P");
        P_NaF -> Draw("same");
        D_NaF  -> Draw("same");
        He_NaF -> Draw("same");

	c38 -> cd(3);
	gPad -> SetGridx();
	gPad -> SetGridy();
	gPad -> SetLogy();
	TH1F * P_Agl =  HeliumContaminationAgl   -> GetResult_P(0);
        TH1F * D_Agl =  HeliumContaminationAgl   -> GetResult_D(0);
        TH1F * He_Agl=  HeliumContaminationAgl   -> GetResult_He(0);

	HeliumContaminationAgl   -> GetResult_Data(0) -> GetXaxis() -> SetTitle("Distance from D");
	HeliumContaminationAgl   -> GetResult_Data(0) -> SetTitle(("He fragm. cont. : " + to_string(HeCont_Agl)).c_str());
	P_Agl  -> SetFillColor(2);
        D_Agl  -> SetFillColor(4);
        He_Agl -> SetFillColor(3);
        HeliumContaminationAgl   -> GetResult_Data(0) -> SetMarkerStyle(8);

        P_Agl  -> SetFillStyle(3001);
        D_Agl  -> SetFillStyle(3001);
        He_Agl -> SetFillStyle(3001);

        HeliumContaminationAgl   -> GetResult_Data(0) -> Draw("P");
        P_Agl  -> Draw("same");
        D_Agl   -> Draw("same");
        He_Agl  -> Draw("same");





	cout<<"*** Updating Results file ***"<<endl;
	nomefile="./Final_plots/"+mese+".root";
	TFile *f_out=new TFile(nomefile.c_str(), "RECREATE");
	f_out->mkdir("MC Results/He related cuts");
	f_out->cd("MC Results/He related cuts");
	c36	->Write();
	c36_bis ->Write();
	c37	->Write();
	c38	->Write();

	f_out->Write();
	f_out->Close();
}


