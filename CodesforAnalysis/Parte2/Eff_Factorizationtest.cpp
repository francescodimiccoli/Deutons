using namespace std;

Efficiency * EffFullSETselectionsMCP =  new Efficiency("EffFullSETselectionsMCP"); 

Efficiency * Eff_do_preSelMCP = new Efficiency("Eff_do_preSelMCP",3); 
Efficiency * Eff_do_preSelMCD = new Efficiency("Eff_do_preSelMCD",6,3);


Efficiency * Eff_do_DistMCP  = new Efficiency("Eff_do_DistMCP");
Efficiency * Eff_do_LikMCP   = new Efficiency("Eff_do_LikMCP");


void FluxFactorizationtest_Pre_Fill(TNtuple *ntupla, int l){
	 ntupla->GetEvent(l);
	if(Tup.Unbias!=0||Tup.Beta_pre<=0||Tup.R_pre<=0) return;
	int Rbin;
	// full set efficiency before
	if(cmask.notPassed(0)){
		Rbin=RB.GetRBin(RUsed);
		if(Massa_gen<1&&Massa_gen>0.5){
			EffFullSETselectionsMCP->beforeR->Fill(Rbin);	
		}
	}
	
	//Drop-one approach eff. calc.
	for(int iS=0;iS<3;iS++){
		Rbin=RB.GetRBin(RUsed);
		if(Massa_gen<1&&Massa_gen>0.5) {
			if(cmask.notPassed(iS)) Eff_do_preSelMCP->beforeR->Fill(Rbin,iS);
			if(cmask.passed(iS)) Eff_do_preSelMCP->afterR ->Fill(Rbin,iS);
		}				 

		if(Massa_gen>1&&Massa_gen<2) {
			if(cmask.notPassed(iS)) FillBinMGen((TH3*)Eff_do_preSelMCD->beforeR, Rbin, iS);
			if(cmask.passed(iS)) FillBinMGen((TH3*)Eff_do_preSelMCD->afterR,  Rbin, iS);

		}
	}
	////////////////////////////////
	return;
}




void FluxFactorizationtest_Qual_Fill(TNtuple *ntupla, int l){

	 ntupla->GetEvent(l);
	//cuts
	if(Tup.Beta<=0||Tup.R<=0||Tup.R<1.2*Tup.Rcutoff) return;
	//R bins
	int Kbin;
	Kbin = RB.GetRBin(RUsed);
	
	//full set efficiency after
	if(Tup.Dist5D_P<6&&Likcut)  EffFullSETselectionsMCP->afterR->Fill(Kbin);
	

	//Drop-one approach eff calc.
	//eff evaluation cuts
	if(Tup.Beta>protons->Eval(Tup.R)+0.1||Tup.Beta<protons->Eval(Tup.R)-0.1) return;
	if(!Herejcut) return;	
	
	if(Massa_gen<1&&Massa_gen>0.5) {
		Eff_do_DistMCP -> beforeR -> Fill(Kbin); 
		if(Tup.Dist5D_P<6) Eff_do_LikMCP -> beforeR -> Fill(Kbin);

		if(Tup.Dist5D_P<6){
			Eff_do_DistMCP -> afterR -> Fill(Kbin);
			if(Likcut) Eff_do_LikMCP -> afterR -> Fill(Kbin);
		}

	}
	///////////////////////////	
	return;

}


void FluxFactorizationtest_Write(){
        Eff_do_preSelMCP -> Write();
        Eff_do_preSelMCD ->Write();
	
	Eff_do_DistMCP   -> Write();
	Eff_do_LikMCP    ->Write();
	
	EffFullSETselectionsMCP ->Write();
	return; 
}


void FluxFactorizationtest(TFile * inputHistoFile){

	Efficiency * EffFullSETselectionsMCP  = new Efficiency(inputHistoFile,"EffFullSETselectionsMCP");
		
	Efficiency * Eff_do_preSelMCP = new Efficiency(inputHistoFile,"Eff_do_preSelMCP"); 
	Efficiency * Eff_do_preSelMCD = new Efficiency(inputHistoFile,"Eff_do_preSelMCD");
		
	Efficiency * Eff_do_DistMCP  = new Efficiency(inputHistoFile,"Eff_do_DistMCP"); 
        Efficiency * Eff_do_LikMCP   = new Efficiency(inputHistoFile,"Eff_do_LikMCP"); 
	
	string tagli[5]={"Matching TOF","Chi^2 R","1 Tr. Track","Distance","Likelihood"};
	string nome;

	cout<<"********* MC \"GOLDEN\" SEL. EFFICIENCIES *********"<<endl;
	
	EffFullSETselectionsMCP -> Eval_Efficiency();
	
	Eff_do_preSelMCP 	-> Eval_Efficiency();
	Eff_do_preSelMCD 	-> Eval_Efficiency();
	
	Eff_do_DistMCP  	-> Eval_Efficiency();
	Eff_do_LikMCP   	-> Eval_Efficiency();



	TH1F * Eff_FullSETMCP_R_TH1F = (TH1F *) EffFullSETselectionsMCP -> effR -> Clone();

	TH2F * Eff_do_preSelMCP_R_TH2F	= (TH2F *) Eff_do_preSelMCP -> effR -> Clone();
	
	TH1F * Eff_do_DistMCP_R_TH1F = (TH1F *) Eff_do_DistMCP -> effR -> Clone();
	TH1F * Eff_do_LikMCP_R_TH1F = (TH1F *) Eff_do_LikMCP -> effR -> Clone();


	// factorized eff. calc.
	TH1F * FactorizedEffMCP_R = (TH1F *) Eff_do_DistMCP_R_TH1F -> Clone();
	FactorizedEffMCP_R -> Multiply(Eff_do_LikMCP_R_TH1F);

	for(int iS=0;iS<Eff_do_preSelMCP_R_TH2F->GetNbinsY();iS++){
		for(int iR=0;iR<Eff_do_preSelMCP_R_TH2F->GetNbinsX();iR++)
			FactorizedEffMCP_R -> SetBinContent(iR+1,FactorizedEffMCP_R -> GetBinContent(iR+1)*Eff_do_preSelMCP_R_TH2F-> GetBinContent(iR+1,iS+1) ); 
	}

	cout<<"*** Updating P1 file ****"<<endl;
   inputHistoFile->ReOpen("UPDATE");

        inputHistoFile->cd("Results");
	inputHistoFile-> Write();
	inputHistoFile-> Close();

		
	TCanvas *c9 = new TCanvas("MC Protons Factorization Test");
	c9->cd();
	gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();

	TGraphErrors * FullsetEfficiency    = new TGraphErrors();
	TGraphErrors * FactorizedEfficiency = new TGraphErrors();

	for(int i=0; i<Eff_FullSETMCP_R_TH1F->GetNbinsX();i++) {
			FullsetEfficiency    ->SetPoint(i,RB.RigBinCent(i),Eff_FullSETMCP_R_TH1F->GetBinContent(i+1));
			FactorizedEfficiency ->SetPoint(i,RB.RigBinCent(i),FactorizedEffMCP_R->GetBinContent(i+1));		
	}

	FullsetEfficiency->SetMarkerColor(2);
        FullsetEfficiency->SetMarkerStyle(8);
        FullsetEfficiency->SetLineColor(2);
        FullsetEfficiency->SetLineWidth(2);

	FactorizedEfficiency->SetMarkerColor(2);
        FactorizedEfficiency->SetMarkerStyle(4);
        FactorizedEfficiency->SetLineColor(2);
        FactorizedEfficiency->SetLineWidth(2);



	FullsetEfficiency->Draw("APC");
	FactorizedEfficiency->Draw("PCsame");
	
	cout<<"*** Updating Results file ***"<<endl;
        string filename="./Final_plots/"+mese+".root";
        TFile *f_out=new TFile(filename.c_str(), "UPDATE");
        f_out->mkdir("MC Results/Eff. Factorization Test");
        f_out->cd("MC Results/Eff. Factorization Test");
        c9->Write();
        f_out->Write();
        f_out->Close();
 	

return;
}

