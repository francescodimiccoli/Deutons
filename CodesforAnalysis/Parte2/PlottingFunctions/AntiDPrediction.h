#include "AntiDFluxes.h"



void AntiDpredictions_Plot(   TH1 * Acceptance_AntiDTOF,
                              TH1 * Acceptance_AntiDNaF,
                              TH1 * Acceptance_AntiDAgl,
				TH1F * RP_TOF ,
			        TH1F * RP_NaF ,
		                TH1F * RP_Agl, 
				float TOF_Threshold,	
                                float NaF_Threshold,
                                float Agl_Threshold,
				float TOF_Pbar_Threshold,
				float NaF_Pbar_Threshold,
				float Agl_Pbar_Threshold
){
	
	TCanvas * c1 = new TCanvas("Anti Deuterons Acceptance");

	c1->cd();
	gPad->SetGridx();
	gPad->SetLogx();
	gPad->SetGridy();
        gPad->SetLogy();

	TGraphErrors * AccAntiDTOF;
        int p=0;
        AccAntiDTOF=new TGraphErrors();
        p=0;
        for(int i=0;i<nbinsToF;i++) {AccAntiDTOF->SetPoint(p,ToFDB.EkPerMassBinCent(i),Acceptance_AntiDTOF  ->GetBinContent(i+1));p++;}
        AccAntiDTOF->SetPoint(p,50,0.000001);
	AccAntiDTOF->SetMarkerStyle(8);
        AccAntiDTOF->SetMarkerColor(4);
        AccAntiDTOF->SetMarkerSize(1.4);
        AccAntiDTOF->SetLineColor(4);
        AccAntiDTOF->SetLineWidth(1);
        AccAntiDTOF->SetTitle("");
        AccAntiDTOF->GetXaxis()->SetTitle("R [GV]");
        AccAntiDTOF->GetYaxis()->SetTitle("Gen. Eff.");
        AccAntiDTOF->GetXaxis()->SetTitleSize(0.045);
        AccAntiDTOF->GetYaxis()->SetTitleSize(0.045);
        AccAntiDTOF->Draw("AP");

	TGraphErrors * AccAntiDNaF;
        AccAntiDNaF=new TGraphErrors();
        p=0;
        for(int i=0;i<nbinsNaF;i++) {AccAntiDNaF->SetPoint(p,NaFDB.EkPerMassBinCent(i),Acceptance_AntiDNaF  ->GetBinContent(i+1));p++;}
	AccAntiDNaF->SetMarkerStyle(8);
        AccAntiDNaF->SetMarkerColor(4);
        AccAntiDNaF->SetMarkerSize(1.4);
        AccAntiDNaF->SetLineColor(4);
        AccAntiDNaF->SetLineWidth(1);
        AccAntiDNaF->SetTitle("");
        AccAntiDNaF->GetXaxis()->SetTitle("R [GV]");
        AccAntiDNaF->GetYaxis()->SetTitle("Gen. Eff.");
        AccAntiDNaF->GetXaxis()->SetTitleSize(0.045);
        AccAntiDNaF->GetYaxis()->SetTitleSize(0.045);
        AccAntiDNaF->Draw("Psame");


  	TGraphErrors * AccAntiDAgl;
        AccAntiDAgl=new TGraphErrors();
        p=0;
        for(int i=0;i<nbinsAgl;i++) {AccAntiDAgl->SetPoint(p,AglDB.EkPerMassBinCent(i),Acceptance_AntiDAgl  ->GetBinContent(i+1));p++;}
	AccAntiDAgl->SetMarkerStyle(8);
        AccAntiDAgl->SetMarkerColor(4);
        AccAntiDAgl->SetMarkerSize(1.4);
        AccAntiDAgl->SetLineColor(4);
        AccAntiDAgl->SetLineWidth(1);
        AccAntiDAgl->SetTitle("");
        AccAntiDAgl->GetXaxis()->SetTitle("R [GV]");
        AccAntiDAgl->GetYaxis()->SetTitle("Gen. Eff.");
        AccAntiDAgl->GetXaxis()->SetTitleSize(0.045);
        AccAntiDAgl->GetYaxis()->SetTitleSize(0.045);
        AccAntiDAgl->Draw("Psame");
	
	        TCanvas * c2 = new TCanvas("Anti Protons Rejection Power");
	float RP;
        c2->cd();
        gPad->SetGridx();
        gPad->SetLogx();
        gPad->SetGridy();
        gPad->SetLogy();

        TGraphErrors * RPTOF;
        RPTOF=new TGraphErrors();
        p=0;
	for(int i=0;i<nbinsToF;i++) {
					RP=Acceptance_AntiDTOF  ->GetBinContent(i+1)/RP_TOF  ->GetBinContent(i+1);
					if(RP>0) {RPTOF->SetPoint(p,ToFDB.EkPerMassBinCent(i),RP);p++;}}
        RPTOF->SetPoint(p,50,0.000001);
        RPTOF->SetMarkerStyle(8);
        RPTOF->SetMarkerColor(2);
        RPTOF->SetMarkerSize(1.4);
        RPTOF->SetLineColor(2);
        RPTOF->SetLineWidth(1);
        RPTOF->SetTitle("");
        RPTOF->GetXaxis()->SetTitle("R [GV]");
        RPTOF->GetYaxis()->SetTitle("Gen. Eff.");
        RPTOF->GetXaxis()->SetTitleSize(0.045);
        RPTOF->GetYaxis()->SetTitleSize(0.045);
        RPTOF->Draw("AP");

        TGraphErrors * RPNaF;
        RPNaF=new TGraphErrors();
        p=0;
        for(int i=0;i<nbinsNaF;i++) {	RP=Acceptance_AntiDNaF  ->GetBinContent(i+1)/RP_NaF  ->GetBinContent(i+1);
					if(RP>0){ RPNaF->SetPoint(p,NaFDB.EkPerMassBinCent(i),RP);p++;}}
        RPNaF->SetMarkerStyle(8);
        RPNaF->SetMarkerColor(2);
        RPNaF->SetMarkerSize(1.4);
        RPNaF->SetLineColor(2);
        RPNaF->SetLineWidth(1);
        RPNaF->SetTitle("");
        RPNaF->GetXaxis()->SetTitle("R [GV]");
        RPNaF->GetYaxis()->SetTitle("Gen. Eff.");
        RPNaF->GetXaxis()->SetTitleSize(0.045);
        RPNaF->GetYaxis()->SetTitleSize(0.045);
        RPNaF->Draw("Psame");


        TGraphErrors * RPAgl;
        RPAgl=new TGraphErrors();
        p=0;
        for(int i=0;i<nbinsAgl;i++) {	RP=Acceptance_AntiDAgl  ->GetBinContent(i+1)/RP_Agl  ->GetBinContent(i+1);
					if(RP>0){ RPAgl->SetPoint(p,AglDB.EkPerMassBinCent(i),RP);p++;}}
        RPAgl->SetMarkerStyle(8);
        RPAgl->SetMarkerColor(2);
        RPAgl->SetMarkerSize(1.4);
	RPAgl->SetLineColor(2);
        RPAgl->SetLineWidth(1);
        RPAgl->SetTitle("");
        RPAgl->GetXaxis()->SetTitle("R [GV]");
        RPAgl->GetYaxis()->SetTitle("Gen. Eff.");
        RPAgl->GetXaxis()->SetTitleSize(0.045);
        RPAgl->GetYaxis()->SetTitleSize(0.045);
	RPAgl->Draw("Psame");


	TCanvas * c3 = new TCanvas("Anti Deuterons Thresholds");

	c3->cd();
	gPad->SetGridx();
	gPad->SetLogx();
	gPad->SetGridy();
        gPad->SetLogy();

	TH2F * Frame = new TH2F("Antideuterons AMS-02 Sensitivity","Antideuterons AMS-02 Sensitivity",10000,0.1,50,10000,1e-8,1e-3);
	Frame ->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
	Frame ->GetYaxis()->SetTitle("Flux [(m2 s sr GeV/nucl)^-1]");
	Frame ->GetXaxis()->SetTitleSize(0.05);
        Frame ->GetYaxis()->SetTitleSize(0.05);
	
	TSpline3 * AntiD_secondaries = new TSpline3("Secondaries",AntiD_SecondaryX,AntiD_SecondaryY,28,0,100);
	TSpline3 * AntiD_LSP  = new TSpline3("LSP(SUSY) (X mass = 30 GeV)",AntiD_LSPX,AntiD_LSPY,23,0,20);
	TSpline3 * AntiD_grav = new TSpline3("gravitino (decay) (m = 50 GeV)",AntiD_gravX,AntiD_gravY,26,0,20);	
	TSpline3 * AntiD_HDM = new TSpline3("heavy DM (10 TeV)",AntiD_HDMX,AntiD_HDMY,25,0,20);

	TGraphErrors * ThresTOF;
	TGraphErrors * ThresNaF;
	TGraphErrors * ThresAgl;

        p=0;
        ThresTOF=new TGraphErrors();
        ThresNaF=new TGraphErrors();
        ThresAgl=new TGraphErrors();
        
	p=0;
        for(int i=0;i<nbinsToF;i++) {ThresTOF->SetPoint(p,ToFDB.EkPerMassBinCent(i),TOF_Threshold);p++;}
        p=0;
        for(int i=0;i<nbinsNaF;i++) {ThresNaF->SetPoint(p,NaFDB.EkPerMassBinCent(i),NaF_Threshold);p++;}
        p=0;
        for(int i=0;i<nbinsAgl;i++) {ThresAgl->SetPoint(p,AglDB.EkPerMassBinCent(i),Agl_Threshold);p++;}

	ThresTOF->SetLineWidth(2);
        ThresNaF->SetLineWidth(2);
	ThresAgl->SetLineWidth(2);

	ThresTOF->SetLineColor(4);
        ThresNaF->SetLineColor(4);
	ThresAgl->SetLineColor(4);


	TGraphErrors * BThresTOF;
	TGraphErrors * BThresNaF;
	TGraphErrors * BThresAgl;

        p=0;
        BThresTOF=new TGraphErrors();
        BThresNaF=new TGraphErrors();
        BThresAgl=new TGraphErrors();
        
	p=0;
        for(int i=0;i<nbinsToF;i++) {BThresTOF->SetPoint(p,ToFDB.EkPerMassBinCent(i),TOF_Pbar_Threshold);p++;}
        p=0;
        for(int i=0;i<nbinsNaF;i++) {BThresNaF->SetPoint(p,NaFDB.EkPerMassBinCent(i),NaF_Pbar_Threshold);p++;}
        p=0;
        for(int i=0;i<nbinsAgl;i++) {BThresAgl->SetPoint(p,AglDB.EkPerMassBinCent(i),Agl_Pbar_Threshold);p++;}

	BThresTOF->SetLineWidth(2);
        BThresNaF->SetLineWidth(2);
	BThresAgl->SetLineWidth(2);

	BThresTOF->SetLineColor(2);
        BThresNaF->SetLineColor(2);
	BThresAgl->SetLineColor(2);


	AntiD_secondaries ->SetLineWidth(3);
        AntiD_LSP         ->SetLineWidth(5);
        AntiD_grav        ->SetLineWidth(5);
        AntiD_HDM         ->SetLineWidth(5);


	AntiD_secondaries ->SetLineColor(1);
        AntiD_LSP         ->SetLineColor(4);
        AntiD_grav        ->SetLineColor(8);
        AntiD_HDM         ->SetLineColor(6);


	AntiD_secondaries ->SetLineStyle(1);
        AntiD_LSP         ->SetLineStyle(2);
	AntiD_grav        ->SetLineStyle(2);
        AntiD_HDM         ->SetLineStyle(2);
	
	Frame->Draw();
	AntiD_secondaries -> Draw("same");
	AntiD_LSP         ->Draw("same"); 
        AntiD_grav        ->Draw("same"); 
	AntiD_HDM         ->Draw("same");	


	ThresTOF->Draw("Lsame");
        ThresNaF->Draw("Lsame");
	ThresAgl->Draw("Lsame");

	BThresTOF->Draw("Lsame");
        BThresNaF->Draw("Lsame");
	BThresAgl->Draw("Lsame");


	        TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(		    AntiD_secondaries,AntiD_secondaries->GetTitle(), "l");
                leg->AddEntry(              AntiD_LSP        ,AntiD_LSP        ->GetTitle(), "l");
                leg->AddEntry(              AntiD_grav       ,AntiD_grav       ->GetTitle(), "l");
                leg->AddEntry(              AntiD_HDM        ,AntiD_HDM        ->GetTitle(), "l");
		leg->AddEntry(ThresTOF,"Sensitivity exclusion");
		leg->AddEntry(BThresTOF,"Signal/noise(antiP) < 1");
		leg->SetFillColor(0);
		leg->Draw("same");
        









	


	finalPlots.Add(c1);
        finalPlots.Add(c2);
	finalPlots.Add(c3);
	finalPlots.writeObjsInFolder("Anti_D Predictions/Predictions ");










	return;
}
