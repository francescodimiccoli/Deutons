using namespace std;

Efficiency * Eff_do_preSelMCP = new Efficiency("Eff_do_preSelMCP",3); 
Efficiency * Eff_do_preSelMCD = new Efficiency("Eff_do_preSelMCD",6,3);


void MC_do_preSeleff_Fill(TNtuple *ntupla, int l){
	int k = ntupla->GetEvent(l);
	if(Unbias!=0||Beta_pre<=0||R_pre<=0) return;

	for(int S=0;S<3;S++){
		if(Massa_gen<1&&Massa_gen>0.5) {
			for(int M=0;M<nbinsr;M++) 
				if(Var<bin[M+1]&&Var>bin[M]) {
						if(((int)Cutmask&notpassed[S])==notpassed[S]) Eff_do_preSelMCP->beforeR->Fill(M,S);
						if(((int)Cutmask&passed[S])==passed[S]) Eff_do_preSelMCP->afterR->Fill(M,S);
				}
	
			for(int m=0;m<nbinsbeta;m++)  
				if(Var>BetaP[m]&&Var<=BetaP[m+1]){
						if(((int)Cutmask&notpassed[S])==notpassed[S]) Eff_do_preSelMCP->beforeTOF->Fill(m,S);
						if(((int)Cutmask&passed[S])==passed[S]) Eff_do_preSelMCP->afterTOF->Fill(m,S);	
				}
			}				 

		if(Massa_gen>1&&Massa_gen<2) {
			for(int M=0;M<nbinsr;M++) 
				if(Var<bin[M+1]&&Var>bin[M]) {
					if(((int)Cutmask&notpassed[S])==notpassed[S]) ((TH3*)(Eff_do_preSelMCD->beforeR))->Fill(1,2,3);
					if(((int)Cutmask&passed[S])==passed[S]) ((TH3*)Eff_do_preSelMCD->afterR)->Fill(M,(int)(10000*Massa_gen-18570),S);
			}
			for(int m=0;m<nbinsbeta;m++) if(Var>BetaD[m]&&Var<=BetaD[m+1]){
					if(((int)Cutmask&notpassed[S])==notpassed[S]) ((TH3*)Eff_do_preSelMCD->beforeTOF)->Fill(m,(int)(10000*Massa_gen-18570),S);
					if(((int)Cutmask&passed[S])==passed[S]) ((TH3*)Eff_do_preSelMCD->afterTOF)->Fill(m,(int)(10000*Massa_gen-18570),S);
			}
		}
	}
	return;
}


void MC_do_preSeleff_Write(){
        Eff_do_preSelMCP -> Write();
        Eff_do_preSelMCD ->Write();
	return; 
}


void MC_do_preSeleff(TFile * file1){
	
	Efficiency * Eff_do_preSelMCP = new Efficiency(file1,"Eff_do_preSelMCP"); 
	Efficiency * Eff_do_preSelMCD = new Efficiency(file1,"Eff_do_preSelMCD");

	string numero[18]={"0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17"};
	string tagli[3]={"Matching TOF","Chi^2 R","1 Tr. Track"};
	string nome;

	cout<<"********* MC \"GOLDEN\" SEL. EFFICIENCIES *********"<<endl;

	Eff_do_preSelMCP -> Eval_Efficiency();
	Eff_do_preSelMCD -> Eval_Efficiency();

	TH2F * Eff_do_preSelMCP_R_TH2F	= (TH2F *) Eff_do_preSelMCP -> effR -> Clone();
	TH2F * Eff_do_preSelMCP_TH2F	= (TH2F *) Eff_do_preSelMCP -> effTOF-> Clone();
	TH3F * Eff_do_preSelMCD_R_TH3F	= (TH3F *) Eff_do_preSelMCD -> effR -> Clone();
	TH3F * Eff_do_preSelMCD_TH3F	= (TH3F *) Eff_do_preSelMCD -> effTOF-> Clone();

	cout<<"*** Updating P1 file ****"<<endl;
        string nomefile="../Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
        file1 =TFile::Open(nomefile.c_str(),"UPDATE");

        file1->cd("Results");
	Eff_do_preSelMCP_R_TH2F	-> Write();
	Eff_do_preSelMCP_TH2F	-> Write();
	Eff_do_preSelMCD_R_TH3F	-> Write();
	Eff_do_preSelMCD_TH3F	-> Write();
	file1-> Write();
	file1-> Close();

		

	TCanvas *c9[4];
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
		for(int i=0;i<nbinsr;i++) Eff_do_preSelMCP_R->SetPoint(i,R_cent[i],Eff_do_preSelMCP_R_TH2F->GetBinContent(i+1,S+1));
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
				for(int i=1;i<nbinsr;i++) Eff_do_preSelMCD_R[h][S]->SetPoint(i,R_cent[i],Eff_do_preSelMCD_R_TH3F->GetBinContent(i+1,h+1,S+1));
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
		for(int i=0;i<nbinsbeta;i++) Eff_do_preSelMCP->SetPoint(i,Ekincent[i],Eff_do_preSelMCP_TH2F->GetBinContent(i+1,S+1));
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
				for(int i=0;i<nbinsbeta;i++) Eff_do_preSelMCD[h][S]->SetPoint(i,Ekincent[i],Eff_do_preSelMCD_TH3F->GetBinContent(i+1,h+1,S+1));
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


return;
}

