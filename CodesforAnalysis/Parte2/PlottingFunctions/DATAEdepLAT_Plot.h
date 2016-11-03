void    DATAEdepLAT_Plot(
			
			TH1F *EdepUTOF_coll[],
                        TH1F *EdepLTOF_coll[],
			TH1F *EdepTrack_coll[]

			){


	TCanvas *c = new TCanvas("E. dep. LAT distributions");
	c->Divide(3,1);
	
	c->cd(1);
	gPad->SetLogy();
	gPad->SetGridx();
	gPad->SetGridy();
	
	EdepUTOF_coll[1]->SetLineColor(1);
	EdepUTOF_coll[1]->SetLineWidth(2);
	EdepUTOF_coll[1]->SetTitle("Energy deposition Upper TOF");
	EdepUTOF_coll[1]->GetXaxis()->SetTitle("E. dep (|meas - teo|) [# of sigmas]");
	EdepUTOF_coll[1]->Draw("L");
	for(int i=2;i<11;i++) {
			EdepUTOF_coll[i]->SetLineColor(i);
			EdepUTOF_coll[i]->SetLineWidth(2);
			EdepUTOF_coll[i]->Draw("Lsame");
	}

	c->cd(2);
	gPad->SetLogy();
	gPad->SetGridx();
	gPad->SetGridy();
	
	EdepLTOF_coll[1]->SetLineColor(1);
	EdepLTOF_coll[1]->SetLineWidth(2);
	EdepLTOF_coll[1]->SetTitle("Energy deposition Lower TOF");
	EdepLTOF_coll[1]->GetXaxis()->SetTitle("E. dep (|meas - teo|) [# of sigmas]");
	EdepLTOF_coll[1]->Draw("L");
	TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
	leg->AddEntry(EdepLTOF_coll[1],("Geo. Zone" + to_string(1)).c_str(), "l");
	for(int i=2;i<11;i++) {
			EdepLTOF_coll[i]->SetLineColor(i);
			EdepLTOF_coll[i]->SetLineWidth(2);
			leg->AddEntry(EdepLTOF_coll[i],("Geo. Zone" + to_string(i)).c_str(), "l");
			EdepLTOF_coll[i]->Draw("Lsame");
	}
	leg->Draw("same");


	c->cd(3);
	gPad->SetLogy();
	gPad->SetGridx();
	gPad->SetGridy();
	
	EdepTrack_coll[1]->SetLineColor(1);
	EdepTrack_coll[1]->SetLineWidth(2);
	EdepTrack_coll[1]->SetTitle("Energy deposition Inner Tracker");
	EdepTrack_coll[1]->GetXaxis()->SetTitle("E. dep (|meas - teo|) [# of sigmas]");
	EdepTrack_coll[1]->Draw("L");
	for(int i=2;i<11;i++) {
			EdepTrack_coll[i]->SetLineColor(i);
			EdepTrack_coll[i]->SetLineWidth(2);
			EdepTrack_coll[i]->Draw("Lsame");
	}


	finalPlots.Add(c);
        finalPlots.writeObjsInFolder("DATA-driven Results/Latitude effect/Quality");

}
