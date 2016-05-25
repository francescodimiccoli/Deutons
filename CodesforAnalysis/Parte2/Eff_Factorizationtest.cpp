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
	if(((int)Tup.Cutmask&notpassed[0])==notpassed[0]){
		Rbin=RB.GetRBin(RUsed);
		if(Massa_gen<1&&Massa_gen>0.5){
			EffFullSETselectionsMCP->beforeR->Fill(Rbin);	
		}
	}
	
	//Drop-one approach eff. calc.
	for(int iS=0;iS<3;iS++){
		Rbin=RB.GetRBin(RUsed);
		if(Massa_gen<1&&Massa_gen>0.5) {
			if(((int)Tup.Cutmask&notpassed[iS])==notpassed[iS]) Eff_do_preSelMCP->beforeR->Fill(Rbin,iS);
			if(((int)Tup.Cutmask&   passed[iS])==   passed[iS]) Eff_do_preSelMCP->afterR ->Fill(Rbin,iS);
		}				 

		if(Massa_gen>1&&Massa_gen<2) {
			if(((int)Tup.Cutmask&notpassed[iS])==notpassed[iS]) FillBinMGen((TH3*)Eff_do_preSelMCD->beforeR, Rbin, iS);
			if(((int)Tup.Cutmask&passed[iS])   ==passed[iS]   ) FillBinMGen((TH3*)Eff_do_preSelMCD->afterR,  Rbin, iS);

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


void FluxFactorizationtest(TFile * file1){

	Efficiency * EffFullSETselectionsMCP  = new Efficiency(file1,"EffFullSETselectionsMCP");
		
	Efficiency * Eff_do_preSelMCP = new Efficiency(file1,"Eff_do_preSelMCP"); 
	Efficiency * Eff_do_preSelMCD = new Efficiency(file1,"Eff_do_preSelMCD");
		
	Efficiency * Eff_do_DistMCP  = new Efficiency(file1,"Eff_do_DistMCP"); 
        Efficiency * Eff_do_LikMCP   = new Efficiency(file1,"Eff_do_LikMCP"); 
	
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
        string nomefile="../Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
        file1 =TFile::Open(nomefile.c_str(),"UPDATE");

        file1->cd("Results");
	file1-> Write();
	file1-> Close();

		
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
        nomefile="./Final_plots/"+mese+".root";
        TFile *f_out=new TFile(nomefile.c_str(), "UPDATE");
        f_out->mkdir("MC Results/Eff. Factorization Test");
        f_out->cd("MC Results/Eff. Factorization Test");
        c9->Write();
        f_out->Write();
        f_out->Close();
 	





	/*TCanvas *c9[4];
	string MCLegend[7]={"protons.B800","d.pl1.0_520_GG_Blic","d.pl1.0_520_GG_BlicDPMJet","d.pl1.0_520_GG_QMD","d.pl1.0_520_Shen_Blic","d.pl1.0_520_Shen_BlicDPMJet","d.pl1.0_520_Shen_QMD"};
	
	for(int S=0;S<3;S++){
		nome="Efficiency: "+tagli[S];
		c9[S]=new TCanvas(nome.c_str());
		c9[S]->Divide(2,1);
		
		
		c9[S]->cd(1);
		gPad->SetLogx();
		gPad->SetGridx();
		gPad->SetGridy();
		TGraph * Eff_do_preSelMCP_R = new TGraph();
		Eff_do_preSelMCP_R->SetTitle(MCLegend[0].c_str());
		for(int i=0;i<nbinsr;i++) Eff_do_preSelMCP_R->SetPoint(i,RB.RigBinCent(i),Eff_do_preSelMCP_R_TH2F->GetBinContent(i+1,S+1));
		TGraph * Eff_do_preSelMCD_R[6][3];
		Eff_do_preSelMCP_R->SetMarkerColor(2);
		Eff_do_preSelMCP_R->SetMarkerStyle(8);
		Eff_do_preSelMCP_R->SetLineColor(2);
		Eff_do_preSelMCP_R->SetLineWidth(2);
		nome="Efficiency: "+tagli[S]+ "MC (R bins)";
		Eff_do_preSelMCP_R->SetTitle(nome.c_str());
		Eff_do_preSelMCP_R->GetXaxis()->SetTitle("R [GV]");
		Eff_do_preSelMCP_R->GetYaxis()->SetTitle("Pres. Efficiency");
		Eff_do_preSelMCP_R->GetXaxis()->SetTitleSize(0.045);
		Eff_do_preSelMCP_R->GetYaxis()->SetTitleSize(0.045);
		{
			Eff_do_preSelMCP_R->Draw("ACP");
			TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
			leg->AddEntry(Eff_do_preSelMCP_R,MCLegend[0].c_str(), "ep");

			for(int h=0;h<6;h++){
				Eff_do_preSelMCD_R[h][S]= new TGraph();
				Eff_do_preSelMCD_R[h][S]->SetTitle(MCLegend[h+1].c_str());
				for(int i=1;i<nbinsr;i++) Eff_do_preSelMCD_R[h][S]->SetPoint(i,RB.RigBinCent(i),Eff_do_preSelMCD_R_TH3F->GetBinContent(i+1,h+1,S+1));
				leg->AddEntry(Eff_do_preSelMCD_R[h][S],MCLegend[h+1].c_str(), "ep");
				Eff_do_preSelMCD_R[h][S]->SetMarkerColor(4);
				Eff_do_preSelMCD_R[h][S]->SetMarkerStyle(h+3);
				Eff_do_preSelMCD_R[h][S]->SetMarkerSize(2);
				Eff_do_preSelMCD_R[h][S]->SetLineColor(4);
				Eff_do_preSelMCD_R[h][S]->SetLineWidth(2);
				Eff_do_preSelMCD_R[h][S]->Draw("Psame");
				leg->Draw();
			}
		}

		c9[S]->cd(2);
		gPad->SetLogx();
		gPad->SetGridx();
		gPad->SetGridy();
		TGraph * Eff_do_preSelMCP = new TGraph();
		for(int i=0;i<nbinsbeta;i++) Eff_do_preSelMCP->SetPoint(i,ToFPB.EkBinCent(i),Eff_do_preSelMCP_TH2F->GetBinContent(i+1,S+1));
		TGraph * Eff_do_preSelMCD[6][3];
		Eff_do_preSelMCP->SetMarkerColor(2);
		Eff_do_preSelMCP->SetMarkerStyle(8);
		Eff_do_preSelMCP->SetLineColor(2);
		Eff_do_preSelMCP->SetLineWidth(2);
		Eff_do_preSelMCP->SetTitle("Preselections Efficiency MC (Beta bins)");
		Eff_do_preSelMCP->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
		Eff_do_preSelMCP->GetYaxis()->SetTitle("Pres. Efficiency");
		Eff_do_preSelMCP->GetXaxis()->SetTitleSize(0.045);
		Eff_do_preSelMCP->GetYaxis()->SetTitleSize(0.045);
		{
			Eff_do_preSelMCP->Draw("ACP");
			TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
			leg->AddEntry(Eff_do_preSelMCP,MCLegend[0].c_str(), "ep");

			for(int h=0;h<6;h++){
				Eff_do_preSelMCD[h][S]= new TGraph();
				for(int i=0;i<nbinsbeta;i++) Eff_do_preSelMCD[h][S]->SetPoint(i,ToFPB.EkBinCent(i),Eff_do_preSelMCD_TH3F->GetBinContent(i+1,h+1,S+1));
				Eff_do_preSelMCD[h][S]->SetMarkerColor(4);
				Eff_do_preSelMCD[h][S]->SetMarkerStyle(h+3);
				leg->AddEntry(Eff_do_preSelMCD[h][S],MCLegend[h+1].c_str(), "ep");
				Eff_do_preSelMCD[h][S]->SetMarkerSize(2);
				Eff_do_preSelMCD[h][S]->SetLineColor(4);
				Eff_do_preSelMCD[h][S]->SetLineWidth(2);
				Eff_do_preSelMCD[h][S]->Draw("Psame");
				leg->Draw();
			}
		}
	}

	cout<<"*** Updating Results file ***"<<endl;
        nomefile="./Final_plots/"+mese+".root";
        TFile *f_out=new TFile(nomefile.c_str(), "UPDATE");
        f_out->mkdir("MC Results/Preselections/\"Clean-Event\" Selections");
        f_out->cd("MC Results/Preselections/\"Clean-Event\" Selections");
        for(int S=0;S<3;S++) c9[S]->Write();
	f_out->Write();
        f_out->Close();

	*/
return;
}

