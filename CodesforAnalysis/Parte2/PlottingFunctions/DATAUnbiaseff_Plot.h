

void 	DATAUnbiaseff_Plot(TH1 * EffUnbDATA_R_TH1F , 
			   TH1 * EffUnbDATAQ_R_TH1F ,
			   TH1 * TriggerGlobalFactor,     
			   TH1 * TriggerGlobalFactorQ
	){

	TCanvas *c12=new TCanvas ("DATA: Unb. Trigger Efficiency");

	c12->cd ();
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	string MCLegend[2]= {"protons","deutons"};
	TGraphErrors * EffUnbDATA_R = new TGraphErrors();
	TGraphErrors * EffUnbDATAQ_R = new TGraphErrors();
	TGraphErrors * FactorQ_R = new TGraphErrors();
	EffUnbDATA_R->SetTitle (MCLegend[0].c_str() );
	for (int i=0; i<nbinsr; i++) {
			EffUnbDATA_R->SetPoint (i,PRB.RigBinCent (i),EffUnbDATA_R_TH1F->GetBinContent (i+1) );
			EffUnbDATAQ_R->SetPoint (i,PRB.RigBinCent (i),EffUnbDATAQ_R_TH1F->GetBinContent (i+1) );
			FactorQ_R->SetPoint (i,PRB.RigBinCent (i),TriggerGlobalFactorQ->GetBinContent(1));
			EffUnbDATA_R->SetPointError (i,0,EffUnbDATA_R_TH1F->GetBinError (i+1) );
                        EffUnbDATAQ_R->SetPointError (i,0,EffUnbDATAQ_R_TH1F->GetBinError (i+1) );
			FactorQ_R->SetPointError (i,0,TriggerGlobalFactorQ->GetBinError(1));
			}

	EffUnbDATA_R->SetMarkerColor (2);
	EffUnbDATA_R->SetMarkerStyle (8);
	EffUnbDATA_R->SetLineColor (2);
	EffUnbDATA_R->SetLineWidth (2);
	EffUnbDATAQ_R->SetMarkerColor (1);
        EffUnbDATAQ_R->SetMarkerStyle (8);
        EffUnbDATAQ_R->SetLineColor (1);
        EffUnbDATAQ_R->SetLineWidth (1);
	
	FactorQ_R->SetLineColor (2);
	FactorQ_R->SetLineWidth (2);
	FactorQ_R->SetFillColor (2);	
	FactorQ_R->SetFillStyle (3001);

	EffUnbDATA_R->SetTitle ("Physical Trigg. Efficiency  (R bins)");
	EffUnbDATA_R->GetYaxis()->SetRangeUser(0,1);
	EffUnbDATA_R->GetXaxis()->SetTitle ("R [GV]");
	EffUnbDATA_R->GetYaxis()->SetTitle ("Efficiency");
	EffUnbDATA_R->GetXaxis()->SetTitleSize (0.045);
	EffUnbDATA_R->GetYaxis()->SetTitleSize (0.045);
	{
		EffUnbDATA_R->Draw ("AP");
		EffUnbDATAQ_R->Draw ("Psame");
		FactorQ_R->Draw ("L4same");
		TLegend* leg =new TLegend (0.4, 0.7,0.95,0.95);
		leg->AddEntry (EffUnbDATA_R,"Trigger Efficiency (Preselections only)", "ep");
		leg->AddEntry (EffUnbDATAQ_R,"Trigger Efficiency (Full-set Selections)", "ep");
		leg->AddEntry (FactorQ_R,"Param.", "l");
		leg->Draw("same");
	}

	finalPlots.Add(c12);
	finalPlots.writeObjsInFolder("DATA-driven Results");
	finalPlots.Add(TriggerGlobalFactor);
	finalPlots.Add(TriggerGlobalFactorQ);
	finalPlots.writeObjsInFolder("Export");	
}
