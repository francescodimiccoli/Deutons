
void 	DVSMCTrackeff_Plot(TH1 *TrackerEfficiencyData,
	                   TH1 *TrackerEfficiencyMC,
			   TH1 *TrackerGlobalFactor ,
			   TH1 * ECALvsR_D ,					
        		   TH1 * ECALvsR_MC
	){
	
	
	TCanvas *c28= new TCanvas("R vs ECAL E.dep.");
	c28->Divide(1,2);
	c28->cd(1);
	gPad->SetLogx();
	gPad->SetLogy();
	gPad->SetLogz();
	ECALvsR_MC->SetTitle("Protons MC");
	ECALvsR_MC->GetXaxis()->SetTitle("R [GV]");
	ECALvsR_MC->GetYaxis()->SetTitle("ECAL E.dep.");
	ECALvsR_MC->Draw("col");
	c28->cd(2);
        gPad->SetLogx();
        gPad->SetLogy();
	gPad->SetLogz();
        ECALvsR_D->SetTitle("DATA");
        ECALvsR_D->GetXaxis()->SetTitle("R [GV]");
        ECALvsR_D->GetYaxis()->SetTitle("ECAL E.dep.");
        ECALvsR_D->Draw("col");

	TCanvas *c29= new TCanvas("Global Tracker Efficiency");
	c29 -> cd();
	TrackerEfficiencyMC    -> SetFillColor(2);
	TrackerEfficiencyData  -> SetFillColor(1);
	TrackerEfficiencyMC  ->SetBarWidth(0.5);
        TrackerEfficiencyData->SetBarWidth(0.5);
	
	TrackerEfficiencyMC    -> Draw("B");
        TrackerEfficiencyData  -> Draw("B,same");

	finalPlots.Add(c28);
        finalPlots.Add(c29);
	finalPlots.writeObjsInFolder("DATA-driven Results/Data vs MC/Tracker Efficiency");
	finalPlots.Add(TrackerEfficiencyData);
	finalPlots.Add(TrackerGlobalFactor);
	finalPlots.writeObjsInFolder("Export");
	
}
