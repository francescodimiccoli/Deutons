
void DistanceCut_Plot(	TH1 * Sum_TOF,
			TH1 *Sum_NaF,
			TH1 *Sum_Agl,
			TGraph *P_TOF_Efficiency,
			TGraph *D_TOF_Efficiency,
			TGraph *P_NaF_Efficiency,
			TGraph *D_NaF_Efficiency, 	 
			TGraph *P_Agl_Efficiency,
			TGraph *D_Agl_Efficiency, 	 
			TGraph *Herej_TOF,
			TGraph *Herej_NaF, 
			TGraph *Herej_Agl, 
			TH1 *DistvsLikTOF_P ,
			TH1 *DistvsLikTOF_D ,
			TH1 *DistvsLikNaF_P ,
			TH1 *DistvsLikNaF_D ,
			TH1 *DistvsLikAgl_P ,
			TH1 *DistvsLikAgl_D ,
			TGraph *BadPrej_TOF,
			TGraph *BadPrej_NaF,
			TGraph *BadPrej_Agl,   	
			TGraph *BadPrejLik_TOF,
			TGraph *BadPrejLik_NaF,
			TGraph *BadPrejLik_Agl
	){
	TCanvas * c1 = new TCanvas("Distance Distributions TOF");
	c1->cd();
	gPad->SetLogy();
	gPad->SetLogx();
	gPad->SetGridy();
        gPad->SetGridx();
	Sum_TOF -> GetXaxis() -> SetTitle("Distance from D");
	Sum_TOF -> GetYaxis() -> SetTitle("Distance from P");
	Sum_TOF -> SetTitle("Distance Distribution");
	Sum_TOF -> Draw("col");

	TCanvas * c2 = new TCanvas("Distance Distributions NaF");
        c2->cd();
        gPad->SetLogy();
        gPad->SetLogx();
        gPad->SetGridy();
        gPad->SetGridx();
        Sum_NaF -> GetXaxis() -> SetTitle("Distance from D");
        Sum_NaF -> GetYaxis() -> SetTitle("Distance from P");
        Sum_NaF -> SetTitle("Distance Distribution");
        Sum_NaF -> Draw("col");

	TCanvas * c3 = new TCanvas("Distance Distributions Agl");
        c3->cd();
        gPad->SetLogy();
        gPad->SetLogx();
        gPad->SetGridy();
        gPad->SetGridx();
        Sum_Agl -> GetXaxis() -> SetTitle("Distance from D");
        Sum_Agl -> GetYaxis() -> SetTitle("Distance from P");
        Sum_Agl -> SetTitle("Distance Distribution");
        Sum_Agl -> Draw("col");

	TCanvas * c4 = new TCanvas("Distance cut Eff.");
	c4->Divide(3,1);
	c4->cd(1);
	gPad->SetLogx();
	gPad->SetGridy();
        gPad->SetGridx();
	P_TOF_Efficiency -> SetLineColor(2);
	P_TOF_Efficiency -> SetMarkerColor(2);	
	P_TOF_Efficiency -> SetLineWidth(4);
	P_TOF_Efficiency ->SetTitle("TOF");
	P_TOF_Efficiency ->GetXaxis()->SetTitle("cut value");
	P_TOF_Efficiency ->GetYaxis()->SetTitle("Efficiency");
	P_TOF_Efficiency ->Draw("APC");
	D_TOF_Efficiency -> SetLineColor(4);
        D_TOF_Efficiency -> SetMarkerColor(4);
        D_TOF_Efficiency -> SetLineWidth(4);
        D_TOF_Efficiency ->Draw("PCsame");

	c4->cd(2);
        gPad->SetLogx();
        gPad->SetGridy();
        gPad->SetGridx();
        P_NaF_Efficiency -> SetLineColor(2);
        P_NaF_Efficiency -> SetMarkerColor(2);
        P_NaF_Efficiency -> SetLineWidth(4);
        P_NaF_Efficiency ->SetTitle("NaF");
        P_NaF_Efficiency ->GetXaxis()->SetTitle("cut value");
        P_NaF_Efficiency ->GetYaxis()->SetTitle("Efficiency");
        P_NaF_Efficiency ->Draw("APC");
        D_NaF_Efficiency -> SetLineColor(4);
        D_NaF_Efficiency -> SetMarkerColor(4);
        D_NaF_Efficiency -> SetLineWidth(4);
        D_TOF_Efficiency ->Draw("PCsame");

	c4->cd(3);
        gPad->SetLogx();
        gPad->SetGridy();
        gPad->SetGridx();
        P_Agl_Efficiency -> SetLineColor(2);
        P_Agl_Efficiency -> SetMarkerColor(2);
        P_Agl_Efficiency -> SetLineWidth(4);
        P_Agl_Efficiency ->SetTitle("Agl");
        P_Agl_Efficiency ->GetXaxis()->SetTitle("cut value");
        P_Agl_Efficiency ->GetYaxis()->SetTitle("Efficiency");
        P_Agl_Efficiency ->Draw("APC");
        D_Agl_Efficiency -> SetLineColor(4);
        D_Agl_Efficiency -> SetMarkerColor(4);
        D_Agl_Efficiency -> SetLineWidth(4);
        D_Agl_Efficiency ->Draw("PCsame");
	
	TCanvas * c5 = new TCanvas("Distance cut He rej.");
        c5->Divide(3,1);
        c5->cd(1);
        gPad->SetLogx();
        gPad->SetGridy();
        gPad->SetGridx();
        Herej_TOF -> SetLineColor(3);
        Herej_TOF -> SetMarkerColor(3);
        Herej_TOF -> SetLineWidth(4);
        Herej_TOF ->SetTitle("TOF");
        Herej_TOF ->GetXaxis()->SetTitle("cut value");
        Herej_TOF ->GetYaxis()->SetTitle("He rejection");
        Herej_TOF ->Draw("APC");

	c5->cd(2);
        gPad->SetLogx();
        gPad->SetGridy();
        gPad->SetGridx();
        Herej_NaF -> SetLineColor(3);
        Herej_NaF -> SetMarkerColor(3);
        Herej_NaF -> SetLineWidth(4);
        Herej_NaF ->SetTitle("TOF");
        Herej_NaF ->GetXaxis()->SetTitle("cut value");
        Herej_NaF ->GetYaxis()->SetTitle("He rejection");
        Herej_NaF ->Draw("APC");

	c5->cd(3);
        gPad->SetLogx();
        gPad->SetGridy();
        gPad->SetGridx();
        Herej_Agl -> SetLineColor(3);
        Herej_Agl -> SetMarkerColor(3);
        Herej_Agl -> SetLineWidth(4);
        Herej_Agl ->SetTitle("TOF");
        Herej_Agl ->GetXaxis()->SetTitle("cut value");
        Herej_Agl ->GetYaxis()->SetTitle("He rejection");
        Herej_Agl ->Draw("APC");

	
	TCanvas * c6 = new TCanvas("Distance vs Likelihood");
        c6->Divide(3,1);
        c6->cd(1);
        gPad->SetLogy();
        gPad->SetGridy();
        gPad->SetGridx();
	DistvsLikTOF_P->SetMarkerColor(2);
	DistvsLikTOF_D->SetMarkerColor(4);
	DistvsLikTOF_D->SetMarkerStyle(8);
	DistvsLikTOF_P->SetMarkerStyle(8);
	DistvsLikTOF_D->SetMarkerSize(0.2);
        DistvsLikTOF_P->SetMarkerSize(0.3);
	DistvsLikTOF_D->SetTitle("Distance vs Likelihood: TOF");
	DistvsLikTOF_D->GetXaxis()->SetRangeUser(0,2.6);
	DistvsLikTOF_D->GetXaxis()->SetTitle("-log(1-Tup.LDiscriminant)");
	DistvsLikTOF_D->GetYaxis()->SetTitle("Distance from D");
	DistvsLikTOF_D->Draw();
	DistvsLikTOF_P->Draw("same");

        c6->cd(2);
        gPad->SetLogy();
        gPad->SetGridy();
        gPad->SetGridx();
	DistvsLikNaF_P->SetMarkerColor(2);
        DistvsLikNaF_D->SetMarkerColor(4);
	DistvsLikNaF_D->SetMarkerStyle(8);
        DistvsLikNaF_P->SetMarkerStyle(8);
        DistvsLikNaF_D->SetMarkerSize(0.3);
        DistvsLikNaF_P->SetMarkerSize(0.3);
	DistvsLikNaF_D->SetTitle("Distance vs Likelihood: NaF");
	DistvsLikNaF_D->GetXaxis()->SetTitle("-log(1-Tup.LDiscriminant)");
        DistvsLikNaF_D->GetYaxis()->SetTitle("Distance from D");
	DistvsLikNaF_D->Draw();
        DistvsLikNaF_P->Draw("same");


        c6->cd(3);
        gPad->SetLogy();
        gPad->SetGridy();
        gPad->SetGridx();
	DistvsLikAgl_P->SetMarkerColor(2);
        DistvsLikAgl_D->SetMarkerColor(4);
        DistvsLikAgl_D->SetMarkerStyle(8);
        DistvsLikAgl_P->SetMarkerStyle(8);
        DistvsLikAgl_D->SetMarkerSize(0.3);
        DistvsLikAgl_P->SetMarkerSize(0.3);
	DistvsLikAgl_D->SetTitle("Distance vs Likelihood: Agl");
	DistvsLikAgl_D->GetXaxis()->SetTitle("-log(1-Tup.LDiscriminant)");
        DistvsLikAgl_D->GetYaxis()->SetTitle("Distance from D");
	DistvsLikAgl_D->Draw();
        DistvsLikAgl_P->Draw("same");

	TCanvas * c7 = new TCanvas("Distance Bad P optimization");
        c7->Divide(3,1);
        c7->cd(1);
        gPad->SetLogx();
        gPad->SetGridy();
        gPad->SetGridx();
        BadPrej_TOF -> SetLineColor(3);
        BadPrej_TOF -> SetMarkerColor(3);
        BadPrej_TOF -> SetLineWidth(4);
        BadPrej_TOF ->SetTitle("TOF");
        BadPrej_TOF ->GetXaxis()->SetTitle("cut value");
        BadPrej_TOF ->GetYaxis()->SetTitle("Bad P rejection");
        BadPrej_TOF ->Draw("APC");

        c7->cd(2);
        gPad->SetLogx();
        gPad->SetGridy();
        gPad->SetGridx();
        BadPrej_NaF  -> SetLineColor(3);
        BadPrej_NaF  -> SetMarkerColor(3);
        BadPrej_NaF  -> SetLineWidth(4);
        BadPrej_NaF  ->SetTitle("NaF");
        BadPrej_NaF  ->GetXaxis()->SetTitle("cut value");
        BadPrej_NaF  ->GetYaxis()->SetTitle("Bad P rejection");
        BadPrej_NaF  ->Draw("APC");

        c7->cd(3);
        gPad->SetLogx();
        gPad->SetGridy();
        gPad->SetGridx();
        BadPrej_Agl -> SetLineColor(3);
        BadPrej_Agl -> SetMarkerColor(3);
        BadPrej_Agl -> SetLineWidth(4);
        BadPrej_Agl ->SetTitle("Agl");
        BadPrej_Agl ->GetXaxis()->SetTitle("cut value");
        BadPrej_Agl ->GetYaxis()->SetTitle("Bad P rejection");
        BadPrej_Agl ->Draw("APC");

	TCanvas * c8 = new TCanvas("Likelihood Bad P optimization");
        c8->Divide(3,1);
        c8->cd(1);
        gPad->SetLogx();
        gPad->SetGridy();
        gPad->SetGridx();
        BadPrejLik_TOF -> SetLineColor(3);
        BadPrejLik_TOF -> SetMarkerColor(3);
        BadPrejLik_TOF -> SetLineWidth(4);
        BadPrejLik_TOF ->SetTitle("TOF");
        BadPrejLik_TOF ->GetXaxis()->SetTitle("cut value");
        BadPrejLik_TOF ->GetYaxis()->SetTitle("Bad P rejection");
        BadPrejLik_TOF ->Draw("APC");

        c8->cd(2);
        gPad->SetLogx();
        gPad->SetGridy();
        gPad->SetGridx();
        BadPrejLik_NaF  -> SetLineColor(3);
        BadPrejLik_NaF  -> SetMarkerColor(3);
        BadPrejLik_NaF  -> SetLineWidth(4);
        BadPrejLik_NaF  ->SetTitle("NaF");
        BadPrejLik_NaF  ->GetXaxis()->SetTitle("cut value");
        BadPrejLik_NaF  ->GetYaxis()->SetTitle("Bad P rejection");
        BadPrejLik_NaF  ->Draw("APC");

        c8->cd(3);
        gPad->SetLogx();
        gPad->SetGridy();
        gPad->SetGridx();
        BadPrejLik_Agl -> SetLineColor(3);
        BadPrejLik_Agl -> SetMarkerColor(3);
        BadPrejLik_Agl -> SetLineWidth(4);
        BadPrejLik_Agl ->SetTitle("Agl");
        BadPrejLik_Agl ->GetXaxis()->SetTitle("cut value");
        BadPrejLik_Agl ->GetYaxis()->SetTitle("Bad P rejection");
        BadPrejLik_Agl ->Draw("APC");

	finalPlots.Add(c1);
	finalPlots.Add(c2);
	finalPlots.Add(c3);
	finalPlots.Add(c4);
	finalPlots.Add(c5);
	finalPlots.Add(c6);
	finalPlots.Add(c7);
	finalPlots.Add(c8);
	finalPlots.writeObjsInFolder("MC Results/Distance Cut");

	return;
}
