using namespace std;

TH2F * ECALvsR_D=new TH2F("ECALvsR_D","ECALvsR_D",1000,0,100,1000,0,100);
TH2F * ECALvsR_MC=new TH2F("ECALvsR_MC","ECALvsR_MC",1000,0,100,1000,0,100);

Efficiency * TrakerEfficiencyMCP = new Efficiency("TrakerEfficiencyMCP");
Efficiency * TrakerEfficiencyD   = new Efficiency("TrakerEfficiencyC"  );

void DVSMCTrackeff_D_Fill(){
	
	//cuts
	if(Tup.Unbias!=0) return;
	if(!(Tup.EdepTOFU<EdepTOFbeta->Eval(Tup.Beta_pre)+1&&Tup.EdepTOFU>EdepTOFbeta->Eval(Tup.Beta_pre)-1)) return;
	if(!(Tup.R_pre>1.2*Tup.Rcutoff))  return;	

	if(cmask.isPreselected()&&Tup.EdepECAL>1)
		ECALvsR_D->Fill(Tup.R_pre,Tup.EdepECAL);
        //R bins
	int Kbin=PRB.GetRBin (Tup.R_pre) ;
	if(((int) Tup.Cutmask&3 ) == 3   && Tup.Beta_pre>0)            TrakerEfficiencyD -> beforeR -> Fill(Kbin);
	if(((int) Tup.Cutmask&11) == 11  && Tup.Beta_pre>0) 	       TrakerEfficiencyD -> afterR  -> Fill(Kbin); 	

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
        int Kbin=PRB.GetRBin (Tup.R_pre) ;
	if(((int) Tup.Cutmask&3)  == 3   && Tup.Beta_pre>0) 	      TrakerEfficiencyMCP -> beforeR -> Fill(Kbin);
        if(((int) Tup.Cutmask&11) == 11  && Tup.Beta_pre>0)	      TrakerEfficiencyMCP -> afterR  -> Fill(Kbin);
        
	return;

}


void DVSMCTrackeff_Write(){
        ECALvsR_D->Write(); 
        ECALvsR_MC->Write();
        TrakerEfficiencyMCP->Write();
	TrakerEfficiencyD  ->Write();
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
	Efficiency * TrakerEfficiencyMCP = new Efficiency(file,"TrakerEfficiencyMCP");
	Efficiency * TrakerEfficiencyD   = new Efficiency(file,"TrakerEfficiencyC"  );

	TH2F * ECALvsR_D =(TH2F*) file->Get("ECALvsR_D");
        TH2F * ECALvsR_MC =(TH2F*) file->Get("ECALvsR_MC");

	cout<<"*************** Tracker Eff: Data vs MC ***************"<<endl;	
 	TrakerEfficiencyMCP  -> Eval_Efficiency();
   	TrakerEfficiencyD    -> Eval_Efficiency();

	TH1F * TrakerEfficiencyMC    = (TH1F *)TrakerEfficiencyMCP  -> effR    ->    Clone();
	TH1F * TrakerEfficiencyData  = (TH1F *)TrakerEfficiencyD    -> effR    ->    Clone();		

	TH1F * Tracker_DvsMC_CorrectionR   = new TH1F("Tracker_DvsMC_CorrectionR  ","Tracker_DvsMC_CorrectionR  ",nbinsr  ,0,nbinsr  );
	TH1F * Tracker_DvsMC_CorrectionTOF = new TH1F("Tracker_DvsMC_CorrectionTOF","Tracker_DvsMC_CorrectionTOF",nbinsToF,0,nbinsToF);
	TH1F * Tracker_DvsMC_CorrectionNaF = new TH1F("Tracker_DvsMC_CorrectionNaF","Tracker_DvsMC_CorrectionNaF",nbinsNaF,0,nbinsNaF);	
	TH1F * Tracker_DvsMC_CorrectionAgl = new TH1F("Tracker_DvsMC_CorrectionAgl","Tracker_DvsMC_CorrectionAgl",nbinsAgl,0,nbinsAgl);

	Set_GlobalDatavsMCCorr(Tracker_DvsMC_CorrectionR  ,TrakerEfficiencyData,TrakerEfficiencyMC);
        Set_GlobalDatavsMCCorr(Tracker_DvsMC_CorrectionTOF,TrakerEfficiencyData,TrakerEfficiencyMC);
        Set_GlobalDatavsMCCorr(Tracker_DvsMC_CorrectionNaF,TrakerEfficiencyData,TrakerEfficiencyMC);
	Set_GlobalDatavsMCCorr(Tracker_DvsMC_CorrectionAgl,TrakerEfficiencyData,TrakerEfficiencyMC);
	
	cout<<"*** Updating P1 file ****"<<endl;
   	inputHistoFile->ReOpen("UPDATE");

	inputHistoFile->cd("Results");
 
         TrakerEfficiencyData -> Write("TrakerEfficiencyData");
         Tracker_DvsMC_CorrectionR     -> Write();
         Tracker_DvsMC_CorrectionTOF   -> Write();
         Tracker_DvsMC_CorrectionNaF   -> Write();
         Tracker_DvsMC_CorrectionAgl   -> Write();
	
	inputHistoFile->Write();
	
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
	TrakerEfficiencyMC    -> SetFillColor(2);
	TrakerEfficiencyData  -> SetFillColor(1);
	TrakerEfficiencyMC  ->SetBarWidth(0.5);
        TrakerEfficiencyData->SetBarWidth(0.5);
	
	TrakerEfficiencyMC    -> Draw("B");
        TrakerEfficiencyData  -> Draw("B,same");

	cout<<"*** Updating Results file ***"<<endl;
        fileFinalPlots->mkdir("DATA-driven Results/Data vs MC/Tracker Efficiency");
        fileFinalPlots->cd("DATA-driven Results/Data vs MC/Tracker Efficiency");
        c28->Write();
        c29->Write();
	fileFinalPlots->cd("Export");
	TrakerEfficiencyData -> Write("TrakerEfficiencyData");
        fileFinalPlots->Write();

}
