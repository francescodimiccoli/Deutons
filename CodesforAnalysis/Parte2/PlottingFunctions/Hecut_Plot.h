void Hecut_Plot(

	TH1 * HecutMC_He	, 
        TH1 * Hecut_D	 ,
        TH1 * HecutMCP_TH1F, 
        TH1 * HecutMCHe_TH1F,
	TH1 * HecutMC_P     ,
	TemplateFIT * HeliumContaminationTOF,
	TemplateFIT * HeliumContaminationNaF,
        TemplateFIT * HeliumContaminationAgl,
	float   HeCont_TOF,float HeCont_NaF,float HeCont_Agl

){


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
                effHecut->SetPoint(K,PRB.RigBinCent(K), HecutMCP_TH1F->GetBinContent(K+1));
                effHecut->SetPointError(K,0,    HecutMCP_TH1F->GetBinError(K+1));
        }
        TGraphErrors *effHecutHe=new TGraphErrors();
        for(int K=0; K<nbinsr; K++) {
                effHecutHe->SetPoint(K,PRB.RigBinCent(K), HecutMCHe_TH1F->GetBinContent(K+1));
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


	bool recreate = true;

	finalPlots.Add(c36	 );
	finalPlots.Add(c36_bis);
	finalPlots.Add(c37	 );
	finalPlots.Add(c38	 );

        finalPlots.writeObjsInFolder("MC Results/He related cuts",recreate);

}
