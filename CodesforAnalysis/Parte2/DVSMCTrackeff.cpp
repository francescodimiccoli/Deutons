using namespace std;

TH2F * ECALvsR_D=new TH2F("ECALvsR_D","ECALvsR_D",1000,0,100,1000,0,100);
TH2F * ECALvsR_MC=new TH2F("ECALvsR_MC","ECALvsR_MC",1000,0,100,1000,0,100);

void DVSMCTrackeff_D_Fill(TNtuple *ntupla, int l){
         ntupla->GetEvent(l);
        if(Tup.Unbias!=0||Tup.Beta_pre<=0) return;
		if(Tup.R_pre>1.2*Tup.Rcutoff){
			if(Tup.EdepTOFU<EdepTOFbeta->Eval(Tup.Beta_pre)+1&&Tup.EdepTOFU>EdepTOFbeta->Eval(Tup.Beta_pre)-1&&((int)Tup.Cutmask&187)==187&&Tup.EdepECAL>1)
			ECALvsR_D->Fill(Tup.R_pre,Tup.EdepECAL);
		}
        return;
}


void DVSMCTrackeff_Fill(TNtuple *ntupla, int l){
	 ntupla->GetEvent(l);
        if(Tup.Unbias!=0||Tup.Beta_pre<=0) return;
                if(Massa_gen<1&&Massa_gen>0.5){
                        if(Tup.EdepTOFU<EdepTOFbeta->Eval(Tup.Beta_pre)+1&&Tup.EdepTOFU>EdepTOFbeta->Eval(Tup.Beta_pre)-1&&((int)Tup.Cutmask&187)==187&&Tup.EdepECAL>1)
                        ECALvsR_MC->Fill(Tup.R_pre,Tup.EdepECAL);
                }
        return;

}


void DVSMCTrackeff_Write(){
        ECALvsR_D->Write(); 
        ECALvsR_MC->Write();
        return;
}

TCanvas *c28= new TCanvas("R vs ECAL E.dep.");
void DVSMCTrackeff(TFile * file){
	TH2F * ECALvsR_D =(TH2F*) file->Get("ECALvsR_D");
        TH2F * ECALvsR_MC =(TH2F*) file->Get("ECALvsR_MC");
	
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
}
