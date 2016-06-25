using namespace std;

TH2F * ECALvsR_D=new TH2F("ECALvsR_D","ECALvsR_D",1000,0,100,1000,0,100);
TH2F * ECALvsR_MC=new TH2F("ECALvsR_MC","ECALvsR_MC",1000,0,100,1000,0,100);

Efficiency * TrackerEfficiencyMCP = new Efficiency("TrackerEfficiencyMCP");
Efficiency * TrackerEfficiencyD   = new Efficiency("TrackerEfficiencyC"  );

void DVSMCTrackeff_D_Fill(){
	
	//cuts
	if(Tup.Unbias!=0) return;
	if(!(Tup.EdepTOFU<EdepTOFbeta->Eval(Tup.Beta_pre)+1&&Tup.EdepTOFU>EdepTOFbeta->Eval(Tup.Beta_pre)-1)) return;
	if(!(Tup.R_pre>1.2*Tup.Rcutoff))  return;	

	if(cmask.isPreselected()&&Tup.EdepECAL>1)
		ECALvsR_D->Fill(Tup.R_pre,Tup.EdepECAL);
        //R bins
	int Kbin=PRB.GetRBin (20) ;
	if(((int) Tup.Cutmask&3 ) == 3   && Tup.Beta_pre>0)            TrackerEfficiencyD -> beforeR -> Fill(Kbin);
	if(((int) Tup.Cutmask&11) == 11  && Tup.Beta_pre>0) 	       TrackerEfficiencyD -> afterR  -> Fill(Kbin); 	

	return;
}


void DVSMCTrackeff_Fill(){
        //cuts
	if(Tup.Unbias!=0||Tup.Beta_pre<=0) return;
	if(!(Tup.EdepTOFU<EdepTOFbeta->Eval(Tup.Beta_pre)+1&&Tup.EdepTOFU>EdepTOFbeta->Eval(Tup.Beta_pre)-1)) return;
	if(!(Massa_gen<1&&Massa_gen>0.5)) return;
 	if(!(Tup.R_pre>0)) return;	
	if(cmask.isPreselected()&&Tup.EdepECAL>1)
		ECALvsR_MC->Fill(Tup.R_pre,Tup.EdepECAL);	
        //R bins
        int Kbin=PRB.GetRBin (20) ;
	if(((int) Tup.Cutmask&3)  == 3   && Tup.Beta_pre>0) 	      TrackerEfficiencyMCP -> beforeR -> Fill(Kbin);
        if(((int) Tup.Cutmask&11) == 11  && Tup.Beta_pre>0)	      TrackerEfficiencyMCP -> afterR  -> Fill(Kbin);
        
	return;

}


void DVSMCTrackeff_Write(){
        ECALvsR_D->Write(); 
        ECALvsR_MC->Write();
        TrackerEfficiencyMCP->Write();
	TrackerEfficiencyD  ->Write();
	return;
}


void Set_GlobalDatavsMCCorr(TH1F * Correction, TH1F * DataEff, TH1F * MCEff){
	for(int i =0; i< Correction -> GetNbinsX(); i++) {
		Correction -> SetBinContent(i,DataEff -> GetBinContent(PRB.GetRBin(20)+1)/(float)MCEff -> GetBinContent(PRB.GetRBin(20)+1) );
		Correction -> SetBinError(i,DataEff -> GetBinError(PRB.GetRBin(20)+1));
	}
	return;
}

void DVSMCTrackeff(TFile * file){
	file->ReOpen("READ");

	Efficiency * TrackerEfficiencyMCP = new Efficiency(file,"TrackerEfficiencyMCP");
	Efficiency * TrackerEfficiencyD   = new Efficiency(file,"TrackerEfficiencyC"  );

	TH2F * ECALvsR_D =(TH2F*) file->Get("ECALvsR_D");
        TH2F * ECALvsR_MC =(TH2F*) file->Get("ECALvsR_MC");

	cout<<"*************** Tracker Eff: Data vs MC ***************"<<endl;	
 	TrackerEfficiencyMCP  -> Eval_Efficiency();
   	TrackerEfficiencyD    -> Eval_Efficiency();

	TH1F * TrackerEfficiencyMC    = (TH1F *)TrackerEfficiencyMCP  -> effR    ->    Clone();
	TH1F * TrackerEfficiencyData  = (TH1F *)TrackerEfficiencyD    -> effR    ->    Clone();		

	TH1F * TrackerGlobalFactor = new TH1F ("TrackerGlobalFactor","TrackerGlobalFactor",1,0,1);
	TrackerGlobalFactor -> SetBinContent(1,TrackerEfficiencyData -> GetBinContent(PRB.GetRBin(20)+1)/(float)TrackerEfficiencyMC -> GetBinContent(PRB.GetRBin(20)+1));
	TrackerGlobalFactor -> SetBinError(1,TrackerEfficiencyData -> GetBinError(PRB.GetRBin(20)+1));		
	
	cout<<"*** Updating P1 file ****"<<endl;
   	inputHistoFile->ReOpen("UPDATE");

	//inputHistoFile->mkdir("Results");
	inputHistoFile->cd("Results");
 
        TrackerEfficiencyData -> Write("TrackerEfficiencyData");
	TrackerGlobalFactor -> Write();
	
	
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
