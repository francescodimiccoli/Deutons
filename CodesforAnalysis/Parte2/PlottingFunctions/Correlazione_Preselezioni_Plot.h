
void Correlazione_Preselezioni_Plot(TH2F * CorrelazionePreselezioni){

        TCanvas *c13=new TCanvas("Correlazione Selezioni");
        c13->cd();
        CorrelazionePreselezioni->Draw("col");

        finalPlots.Add(c13);
        finalPlots.writeObjsInFolder("MC Results");

        return;
}

