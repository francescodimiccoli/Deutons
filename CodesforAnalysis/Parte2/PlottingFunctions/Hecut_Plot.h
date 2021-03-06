void Hecut_Plot(

	TH1 * L1TOF_DATA 	, 
        TH1 * L1TOF_DATAcutoff, 
        TH1 * L1NaF_DATA 	 ,
        TH1 * L1NaF_DATAcutoff ,
        TH1 * L1Agl_DATA 	 ,
        TH1 * L1Agl_DATAcutoff ,

	TH1 * L1TOFs_DATA 	, 
        TH1 * L1TOFs_DATAcutoff, 
        TH1 * L1NaFs_DATA 	 ,
        TH1 * L1NaFs_DATAcutoff ,
        TH1 * L1Agls_DATA 	 ,
        TH1 * L1Agls_DATAcutoff ,

	TH1 * HecutMC_He	, 
        TH1 * Hecut_D	 ,
        TH1 * HecutMCP_TH1F, 
        TH1 * HecutMCHe_TH1F,
	TH1 * HecutMC_P,

	TH1 *fragmTRDTOF_MC,
        TH1 *fragmTRDNaF_MC,
	TH1 *fragmTRDAgl_MC,
	TH1 *fragmTRDTOF_D,	
	TH1 *fragmTRDNaF_D,
        TH1 *fragmTRDAgl_D,

	TH1 * fragmeffTOF,
	TH1 * fragmeffNaF,
	TH1 * fragmeffAgl, 
	TH1 * ContaminationTOF,
	TH1 * ContaminationNaF,
	TH1 * ContaminationAgl    
){


	
	TCanvas * c35 	=new TCanvas("L1 E. dep.sigma");
	TCanvas * c35_1=new TCanvas("L1 E. dep. TOF");
	TCanvas * c35_2=new TCanvas("L1 E. dep. NaF");
	TCanvas * c35_3=new TCanvas("L1 E. dep. Agl");
	TCanvas * c36	=new TCanvas("Sigma E. dep. Track vs TOF");
	TCanvas * c36_bis    =new TCanvas("Sigma E. dep. Track vs TOF (MC)");
	TCanvas * c37	=new TCanvas("Eff. He Control sample cut");
	TCanvas * c38   =new TCanvas("Helium fragmentation above L1");
	TCanvas * c38_bis    =new TCanvas("Helium fragmentation in TRD");
	TCanvas * c39   =new TCanvas("Helium expected contamination");
	

	c35->Divide(3,1);
	c35->cd(1);
	gPad->SetLogy();
        gPad->SetGridx();
        gPad->SetGridy();

	L1TOFs_DATA -> SetLineColor(4);
	L1TOFs_DATAcutoff ->SetLineColor(1);
	L1TOFs_DATA -> SetLineWidth(3);
        L1TOFs_DATAcutoff ->SetLineWidth(3);
	L1TOFs_DATA -> SetTitle("TOF range");
	L1TOFs_DATA -> GetXaxis() -> SetTitle("L1 E. dep. - Q=1(teo.)");

	L1TOFs_DATA -> Draw();
	L1TOFs_DATAcutoff -> Draw("same");

	c35->cd(2);
        gPad->SetLogy();
        gPad->SetGridx();
        gPad->SetGridy();

	L1NaFs_DATA -> SetLineColor(4);
        L1NaFs_DATAcutoff ->SetLineColor(1);
        L1NaFs_DATA -> SetLineWidth(3);
        L1NaFs_DATAcutoff ->SetLineWidth(3);
        L1NaFs_DATA -> SetTitle("TOF range");
        L1NaFs_DATA -> GetXaxis() -> SetTitle("L1 E. dep. - Q=1(teo.)");

	
        L1NaFs_DATA -> Draw();
        L1NaFs_DATAcutoff -> Draw("same");

	c35->cd(3);
        gPad->SetLogy();
        gPad->SetGridx();
        gPad->SetGridy();

	L1Agls_DATA -> SetLineColor(4);
        L1Agls_DATAcutoff ->SetLineColor(1);
        L1Agls_DATA -> SetLineWidth(3);
        L1Agls_DATAcutoff ->SetLineWidth(3);
        L1Agls_DATA -> SetTitle("TOF range");
        L1Agls_DATA -> GetXaxis() -> SetTitle("L1 E. dep. - Q=1(teo.)");

        L1Agls_DATA -> Draw();
        L1Agls_DATAcutoff -> Draw("same");

	c35_1->Divide(6,3);
	int binsTOF = L1TOF_DATA->GetNbinsY();
	TH1F * SlicesTOF[binsTOF];
	TH1F * SlicesTOFcutoff[binsTOF];
	for(int i=0;i<binsTOF;i++){ 
		c35_1->cd(i+1);
		gPad->SetLogy();
        	gPad->SetGridx();
        	gPad->SetGridy();
		SlicesTOF[i]       = ProjectionXtoTH1F( (TH2F *)L1TOF_DATA,("SliceTOF"  + to_string(i)).c_str(),i+1,i+1);
		SlicesTOFcutoff[i] = ProjectionXtoTH1F( (TH2F *)L1TOF_DATAcutoff,("SliceTOFc"  + to_string(i)).c_str(),i+1,i+1);
		SlicesTOF[i] -> SetLineColor(4);
		SlicesTOFcutoff[i] -> SetLineWidth(2);
		SlicesTOFcutoff[i] -> SetLineColor(1);
                SlicesTOF[i] -> SetLineWidth(2);
		SlicesTOF[i] -> SetTitle(("Edep L1 TOF: bin" + to_string(i)).c_str());
		SlicesTOF[i] ->  GetXaxis() -> SetTitle("L1 E. dep. [keV]");
		SlicesTOF[i] -> Draw();	
		SlicesTOFcutoff[i] -> Draw("same");	
	}
	c35_2->Divide(6,3);
        int binsNaF = L1NaF_DATA->GetNbinsY();
        TH1F * SlicesNaF[binsNaF];
	TH1F * SlicesNaFcutoff[binsNaF];
	for(int i=0;i<binsNaF;i++){
                c35_2->cd(i+1);
                gPad->SetLogy();
                gPad->SetGridx();
                gPad->SetGridy();
                SlicesNaF[i] = ProjectionXtoTH1F( (TH2F *)L1NaF_DATA,("SliceNaF"  + to_string(i)).c_str(),i+1,i+1);
                SlicesNaFcutoff[i] = ProjectionXtoTH1F( (TH2F *)L1NaF_DATAcutoff,("SliceNaFc"  + to_string(i)).c_str(),i+1,i+1);
		SlicesNaF[i] -> SetLineColor(4);
                SlicesNaF[i] -> SetLineWidth(2);
		SlicesNaFcutoff[i] -> SetLineColor(1);
                SlicesNaFcutoff[i] -> SetLineWidth(2);
                SlicesNaF[i] -> SetTitle(("Edep L1 NaF: bin" + to_string(i)).c_str());
                SlicesNaF[i] ->  GetXaxis() -> SetTitle("L1 E. dep. [keV]");
                SlicesNaF[i] -> Draw();
		SlicesNaFcutoff[i] -> Draw("same");
        }
	c35_3->Divide(6,3);
        int binsAgl = L1Agl_DATA->GetNbinsY();
        TH1F * SlicesAgl[binsAgl];
	TH1F * SlicesAglcutoff[binsAgl];
	for(int i=0;i<binsAgl;i++){
                c35_3->cd(i+1);
                gPad->SetLogy();
                gPad->SetGridx();
                gPad->SetGridy();
                SlicesAgl[i]   = ProjectionXtoTH1F( (TH2F *)L1Agl_DATA,("SliceAgl"  + to_string(i)).c_str(),i+1,i+1);
                SlicesAglcutoff[i]   = ProjectionXtoTH1F( (TH2F *)L1Agl_DATAcutoff,("SliceAglc"  + to_string(i)).c_str(),i+1,i+1);
		SlicesAgl[i]  -> SetLineColor(4);
                SlicesAgl[i]  -> SetLineWidth(2);
                SlicesAglcutoff[i]  -> SetLineColor(1);
                SlicesAglcutoff[i]  -> SetLineWidth(2);
		SlicesAgl[i]  -> SetTitle(("Edep L1 Agl: bin" + to_string(i)).c_str());
                SlicesAgl[i]  ->  GetXaxis() -> SetTitle("L1 E. dep. [keV]");
                SlicesAgl[i]  -> Draw();
		SlicesAglcutoff[i]  -> Draw("same");
        }



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
	TH2F * Frame = new TH2F("Total Helium Fragmentation (MC) (He->D,P,T)","Total Helium Fragmentation (MC) (He->D,P,T)",10000,0.5,30,10000,0,0.3);
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


	c38_bis->cd();
        gPad->SetGridx();
        gPad->SetGridy();
        gPad->SetLogx();
	TH2F * FrameTRD = new TH2F("Helium Fragmentation in TRD (He->D,P,T)","Helium Fragmentation in TRD (He->D,P,T)",10000,0.5,30,10000,0,0.3);
        TGraphErrors *TRDfragmTOF=new TGraphErrors();
        TGraphErrors *TRDfragmTOFMC=new TGraphErrors();
	for(int K=0; K<nbinsToF; K++) {
                TRDfragmTOF->SetPoint(K,ToFDB.RigBinCent(K),(fragmTRDTOF_D->GetBinContent(K)+fragmTRDTOF_D->GetBinContent(K+1))/2);
                TRDfragmTOF->SetPointError(K,0,   fragmTRDTOF_D->GetBinError(K+1));
		TRDfragmTOFMC->SetPoint(K,ToFDB.RigBinCent(K),(fragmTRDTOF_MC->GetBinContent(K+1) + fragmTRDTOF_MC->GetBinContent(K+1))/2);
                TRDfragmTOFMC->SetPointError(K,0,   fragmTRDTOF_MC->GetBinError(K+1));
        }
	TGraphErrors *TRDfragmNaF=new TGraphErrors();
	TGraphErrors *TRDfragmNaFMC=new TGraphErrors();
	for(int K=0; K<nbinsNaF; K++) {
                TRDfragmNaF->SetPoint(K,NaFDB.RigBinCent(K),(fragmTRDNaF_D->GetBinContent(K) + fragmTRDNaF_D->GetBinContent(K+1))/2);
                TRDfragmNaF->SetPointError(K,0,   fragmTRDNaF_D->GetBinError(K+1));
		TRDfragmNaFMC->SetPoint(K,NaFDB.RigBinCent(K),(fragmTRDNaF_MC->GetBinContent(K) + fragmTRDNaF_MC->GetBinContent(K+1))/2 );
                TRDfragmNaFMC->SetPointError(K,0,   fragmTRDNaF_MC->GetBinError(K+1));
        }
	TGraphErrors *TRDfragmAgl=new TGraphErrors();
	TGraphErrors *TRDfragmAglMC=new TGraphErrors();
	for(int K=0; K<nbinsAgl; K++) {
                TRDfragmAgl->SetPoint(K,AglDB.RigBinCent(K),(fragmTRDAgl_D->GetBinContent(K)+fragmTRDAgl_D->GetBinContent(K+1))/2);
                TRDfragmAgl->SetPointError(K,0,   fragmTRDAgl_D->GetBinError(K+1));
        	TRDfragmAglMC->SetPoint(K,AglDB.RigBinCent(K),(fragmTRDAgl_MC->GetBinContent(K)+fragmTRDAgl_MC->GetBinContent(K+1))/2);
                TRDfragmAglMC->SetPointError(K,0,   fragmTRDAgl_MC->GetBinError(K+1));
	}

	TRDfragmTOF->SetMarkerStyle(8);
	TRDfragmTOF->SetMarkerColor(3);
	TRDfragmTOF->SetLineColor(3);
	TRDfragmNaF->SetMarkerStyle(3);
	TRDfragmNaF->SetMarkerColor(3);
	TRDfragmNaF->SetLineColor(3);
	TRDfragmAgl->SetMarkerStyle(4);
	TRDfragmAgl->SetMarkerColor(3);
	TRDfragmAgl->SetLineColor(3);

	TRDfragmTOFMC->SetMarkerStyle(8);
	TRDfragmTOFMC->SetMarkerColor(1);
	TRDfragmTOFMC->SetLineColor(1);
	TRDfragmNaFMC->SetMarkerStyle(3);
	TRDfragmNaFMC->SetMarkerColor(1);
	TRDfragmNaFMC->SetLineColor(1);
	TRDfragmAglMC->SetMarkerStyle(4);
	TRDfragmAglMC->SetMarkerColor(1);
	TRDfragmAglMC->SetLineColor(1);

	FrameTRD -> GetYaxis() -> SetRangeUser(0,0.3);	
	FrameTRD -> GetYaxis() -> SetTitle("He -> D,P,T / He -> He");
	FrameTRD -> GetXaxis() -> SetTitle("R [GV]"); 
	FrameTRD->Draw();

	TRDfragmTOF->Draw("cpsame");
	TRDfragmNaF->Draw("cpsame");
	TRDfragmAgl->Draw("cpsame");

	TRDfragmTOFMC->Draw("cpsame");
        TRDfragmNaFMC->Draw("cpsame");
        TRDfragmAglMC->Draw("cpsame");



	c39->cd();
        gPad->SetGridx();
        gPad->SetGridy();
        gPad->SetLogx();
        TH2F * Frame2 = new TH2F("He -> (D,P,T) Contamination estimation","He -> (D,P,T) Contamination estimation",10000,0.5,30,10000,0,0.3);
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

	finalPlots.Add(c35       );
	finalPlots.Add(c35_1     );
	finalPlots.Add(c35_2     );
	finalPlots.Add(c35_3     );
	finalPlots.writeObjsInFolder("He Fragmentation/Layer1 E.dep",recreate);

	finalPlots.Add(c36	 );
	finalPlots.Add(c36_bis);
	finalPlots.Add(c37	 );

        finalPlots.writeObjsInFolder("He Fragmentation/Control Sample cuts");

	finalPlots.Add(c38       );
	finalPlots.Add(c38_bis   );
	finalPlots.Add(c39       );
	
	finalPlots.writeObjsInFolder("He Fragmentation/He fragm.");

}
