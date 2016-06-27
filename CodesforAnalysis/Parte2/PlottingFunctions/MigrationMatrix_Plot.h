
void MigrationMatrix_Plot( TH1 * MigrMatrix){
	
	TCanvas * c27 = new TCanvas("Rigidity Migration matrix");
	c27->cd();
	gPad->SetLogz();
	MigrMatrix->GetXaxis()->SetTitle("n.bin (R meas)");
	MigrMatrix->GetYaxis()->SetTitle("n.bin (R gen)");
	MigrMatrix->Draw("col");


	finalPlots.Add(c27);
	finalPlots.writeObjsInFolder("MC Results");

}
