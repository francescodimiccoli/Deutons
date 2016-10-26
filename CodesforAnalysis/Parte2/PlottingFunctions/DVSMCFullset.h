using namespace std;

void DVSMCFullset_Plot(
                TH2F *SystPlot_R,
                TH2F *SystPlot_TOF,
                TH2F * SystPlot_NaF,
                TH2F *SystPlot_Agl    ){

	TCanvas *c = new TCanvas("Systematics bands");
	c->Divide(2,2);
	
	c->cd(1);
	SystPlot_R->GetXaxis()->SetTitle("Bin nr.");
	SystPlot_R->GetXaxis()->SetTitle("Efficiency correction (Data/MC)");
	SystPlot_R->SetTitle("0<R<100 GV range");
	SystPlot_R->GetXaxis()->SetTitleSize(0.045);
	SystPlot_R->GetYaxis()->SetTitleSize(0.045);
	SystPlot_R->Draw("col");
	

	c->cd(2);
	SystPlot_TOF->GetXaxis()->SetTitle("Bin nr.");
	SystPlot_TOF->GetXaxis()->SetTitle("Efficiency correction (Data/MC)");
	SystPlot_TOF->SetTitle("TOF range");
	SystPlot_TOF->GetXaxis()->SetTitleSize(0.045);
	SystPlot_TOF->GetYaxis()->SetTitleSize(0.045);
	SystPlot_TOF->Draw("col");
	
	c->cd(3);
	SystPlot_NaF->GetXaxis()->SetTitle("Bin nr.");
	SystPlot_NaF->GetXaxis()->SetTitle("Efficiency correction (Data/MC)");
	SystPlot_NaF->SetTitle("NaF range");
	SystPlot_NaF->GetXaxis()->SetTitleSize(0.045);
	SystPlot_NaF->GetYaxis()->SetTitleSize(0.045);
	SystPlot_NaF->Draw("col");


	c->cd(4);
        SystPlot_Agl->GetXaxis()->SetTitle("Bin nr.");
        SystPlot_Agl->GetXaxis()->SetTitle("Efficiency correction (Data/MC)");
        SystPlot_Agl->SetTitle("Agl range");
        SystPlot_Agl->GetXaxis()->SetTitleSize(0.045);
        SystPlot_Agl->GetYaxis()->SetTitleSize(0.045);
        SystPlot_Agl->Draw("col");	


//	finalPlots.Add(c);
//        finalPlots.writeObjsInFolder("DATA-driven Results/Data vs MC/Protons");



	return;

}

