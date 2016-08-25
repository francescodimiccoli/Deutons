void Hecut_Plot(

	TH1 * HecutMC_He	, 
        TH1 * Hecut_D	 ,
        TH1 * HecutMCP_TH1F, 
        TH1 * HecutMCHe_TH1F,
	TH1 * HecutMC_P,
	TH1 * fragmeffTOF,
	TH1 * fragmeffNaF,
	TH1 * fragmeffAgl, 
	TH1 * ContaminationTOF,
	TH1 * ContaminationNaF,
	TH1 * ContaminationAgl    
){


	TCanvas * c36	=new TCanvas("Sigma E. dep. Track vs TOF");
	TCanvas * c36_bis    =new TCanvas("Sigma E. dep. Track vs TOF (MC)");
	TCanvas * c37	=new TCanvas("Eff. He Control sample cut");
	TCanvas * c38   =new TCanvas("Helium fragmentation");
	TCanvas * c39   =new TCanvas("Helium expected contamination");
	
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

	c38->cd();
        gPad->SetGridx();
        gPad->SetGridy();
        gPad->SetLogx();
	TH2F * Frame = new TH2F("Helium Fragmentation (He->D,P,T)","Helium Fragmentation (He->D,P,T)",1000,0.5,30,1000,0,1);
        TGraphErrors *efffragmTOF=new TGraphErrors();
        for(int K=0; K<nbinsToF; K++) {
                efffragmTOF->SetPoint(K,ToFDB.RigBinCent(K),fragmeffTOF->GetBinContent(K+1));
                efffragmTOF->SetPointError(K,0,   fragmeffTOF->GetBinError(K+1));
        }
	TGraphErrors *efffragmNaF=new TGraphErrors();
	for(int K=0; K<nbinsNaF; K++) {
                efffragmNaF->SetPoint(K,NaFDB.RigBinCent(K),fragmeffNaF->GetBinContent(K+1));
                efffragmNaF->SetPointError(K,0,   fragmeffNaF->GetBinError(K+1));
        }
	TGraphErrors *efffragmAgl=new TGraphErrors();
	for(int K=0; K<nbinsAgl; K++) {
                efffragmAgl->SetPoint(K,AglDB.RigBinCent(K),fragmeffAgl->GetBinContent(K+1));
                efffragmAgl->SetPointError(K,0,   fragmeffAgl->GetBinError(K+1));
        }

	efffragmTOF->SetMarkerStyle(8);
	efffragmTOF->SetMarkerColor(3);
	efffragmTOF->SetLineColor(3);
	efffragmNaF->SetMarkerStyle(3);
	efffragmNaF->SetMarkerColor(3);
	efffragmNaF->SetLineColor(3);
	efffragmAgl->SetMarkerStyle(4);
	efffragmAgl->SetMarkerColor(3);
	efffragmAgl->SetLineColor(3);

	Frame -> GetYaxis() -> SetRangeUser(0,0.3);	
	Frame -> GetYaxis() -> SetTitle("He -> D,P,T / He -> He");
	Frame -> GetXaxis() -> SetTitle("R [GV]"); 
	Frame->Draw();

	efffragmTOF->Draw("cpsame");
	efffragmNaF->Draw("cpsame");
	efffragmAgl->Draw("cpsame");

	c39->cd();
        gPad->SetGridx();
        gPad->SetGridy();
        gPad->SetLogx();
        TH2F * Frame2 = new TH2F("He -> (D,P,T) Contamination estimation","He -> (D,P,T) Contamination estimation",1000,0.5,30,1000,0,1);
        TGraphErrors *contTOF=new TGraphErrors();
        for(int K=0; K<nbinsToF; K++) {
                contTOF->SetPoint(K,ToFDB.RigBinCent(K),ContaminationTOF->GetBinContent(K+1));
                contTOF->SetPointError(K,0,   ContaminationTOF->GetBinError(K+1));
        }
        TGraphErrors *contNaF=new TGraphErrors();
        for(int K=0; K<nbinsNaF; K++) {
                contNaF->SetPoint(K,NaFDB.RigBinCent(K),ContaminationNaF->GetBinContent(K+1));
                contNaF->SetPointError(K,0,   ContaminationNaF->GetBinError(K+1));
        }
        TGraphErrors *contAgl=new TGraphErrors();
        for(int K=0; K<nbinsAgl; K++) {
                contAgl->SetPoint(K,AglDB.RigBinCent(K),ContaminationAgl->GetBinContent(K+1));
                contAgl->SetPointError(K,0,   ContaminationAgl->GetBinError(K+1));
        }

        contTOF->SetMarkerStyle(8);
        contTOF->SetMarkerColor(3);
        contTOF->SetLineColor(3);
        contNaF->SetMarkerStyle(3);
        contNaF->SetMarkerColor(3);
        contNaF->SetLineColor(3);
        contAgl->SetMarkerStyle(4);
        contAgl->SetMarkerColor(3);
        contAgl->SetLineColor(3);

        Frame2 -> GetYaxis() -> SetRangeUser(0,0.3);
        Frame2 -> GetYaxis() -> SetTitle("Expected He -> (Q=1) / Q=1");
        Frame2 -> GetXaxis() -> SetTitle("R [GV]");
        Frame2->Draw();

        contTOF->Draw("cpsame");
        contNaF->Draw("cpsame");
        contAgl->Draw("cpsame");




	bool recreate = true;

	finalPlots.Add(c36	 );
	finalPlots.Add(c36_bis);
	finalPlots.Add(c37	 );

        finalPlots.writeObjsInFolder("MC Results/He related cuts/Control Sample cuts",recreate);

	finalPlots.Add(c38       );
	finalPlots.Add(c39       );
	
	finalPlots.writeObjsInFolder("MC Results/He related cuts/He fragm.");

}
