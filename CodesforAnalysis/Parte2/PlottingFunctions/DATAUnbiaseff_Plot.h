

void 	DATAUnbiaseff_Plot(TH1 * EffUnbDATA_R_TH1F , 
			   TH1 * TriggerGlobalFactor     

	){

	TCanvas *c12=new TCanvas ("DATA: Unb. Trigger Efficiency");

	c12->cd ();
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	string MCLegend[2]= {"protons","deutons"};
	TGraph * EffUnbDATA_R = new TGraph();
	EffUnbDATA_R->SetTitle (MCLegend[0].c_str() );
	for (int i=0; i<nbinsr; i++) EffUnbDATA_R->SetPoint (i,PRB.RigBinCent (i),EffUnbDATA_R_TH1F->GetBinContent (i+1) );
	EffUnbDATA_R->SetMarkerColor (2);
	EffUnbDATA_R->SetMarkerStyle (8);
	EffUnbDATA_R->SetLineColor (2);
	EffUnbDATA_R->SetLineWidth (2);
	EffUnbDATA_R->SetTitle ("Physical Trigg. Efficiency  (R bins)");
	EffUnbDATA_R->GetYaxis()->SetRangeUser(0,1);
	EffUnbDATA_R->GetXaxis()->SetTitle ("R [GV]");
	EffUnbDATA_R->GetYaxis()->SetTitle ("Efficiency");
	EffUnbDATA_R->GetXaxis()->SetTitleSize (0.045);
	EffUnbDATA_R->GetYaxis()->SetTitleSize (0.045);
	{
		EffUnbDATA_R->Draw ("ACP");
		TLegend* leg =new TLegend (0.4, 0.7,0.95,0.95);
		leg->AddEntry (EffUnbDATA_R,MCLegend[0].c_str(), "ep");

	}

	finalPlots.Add(c12);
	finalPlots.writeObjsInFolder("DATA-driven Results");
	finalPlots.Add(TriggerGlobalFactor);
	finalPlots.writeObjsInFolder("Export");	
}
