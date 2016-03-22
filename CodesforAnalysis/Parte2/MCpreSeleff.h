using namespace std;


TH2F * EffpreSelMCP1=new TH2F("EffpreSelMCP1","EffpreSelMCP1",18,0,18,3,0,3);
TH2F * EffpreSelMCP2=new TH2F("EffpreSelMCP2","EffpreSelMCP2",18,0,18,3,0,3);
TH2F * EffpreSelMCP1_R=new TH2F("EffpreSelMCP1_R","EffpreSelMCP1_R",43,0,43,3,0,3);
TH2F * EffpreSelMCP2_R=new TH2F("EffpreSelMCP2_R","EffpreSelMCP2_R",43,0,43,3,0,3);
TH3F * EffpreSelMCD1=new TH3F("EffpreSelMCD1","EffpreSelMCD1",18,0,18,6,0,6,3,0,3);
TH3F * EffpreSelMCD2=new TH3F("EffpreSelMCD2","EffpreSelMCD2",18,0,18,6,0,6,3,0,3);
TH3F * EffpreSelMCD1_R=new TH3F("EffpreSelMCD1_R","EffpreSelMCD1_R",43,0,43,6,0,6,3,0,3);
TH3F * EffpreSelMCD2_R=new TH3F("EffpreSelMCD2_R","EffpreSelMCD2_R",43,0,43,6,0,6,3,0,3);

Efficiency * EffpreSelMCP = new Efficiency("EffpreSelMCP",3); 
Efficiency * EffpreSelMCD = new Efficiency("EffpreSelMCD",6,3);



void MCpreSeleff_Fill(TNtuple *ntupla, int l){
	int k = ntupla->GetEvent(l);
	if(Unbias!=0||Beta_pre<=0||R_pre<=0) return;

	for(int S=0;S<3;S++){
		if(Massa_gen<1&&Massa_gen>0.5) {
			for(int M=0;M<43;M++) 
				if(Var<bin[M+1]&&Var>bin[M]) {
						if(((int)Cutmask&notpassed[S])==notpassed[S]) EffpreSelMCP->beforeR->Fill(M,S);
						if(((int)Cutmask&passed[S])==passed[S]) EffpreSelMCP->afterR->Fill(M,S);
				}
	
			for(int m=0;m<18;m++)  
				if(Var>BetaP[m]&&Var<=BetaP[m+1]){
						if(((int)Cutmask&notpassed[S])==notpassed[S]) EffpreSelMCP->beforeTOF->Fill(m,S);
						if(((int)Cutmask&passed[S])==passed[S]) EffpreSelMCP->afterTOF->Fill(m,S);	
				}
			}				 

		if(Massa_gen>1&&Massa_gen<2) {
			for(int M=0;M<43;M++) 
				if(Var<bin[M+1]&&Var>bin[M]) {
					if(((int)Cutmask&notpassed[S])==notpassed[S]) ((TH3*)(EffpreSelMCD->beforeR))->Fill(1,2,3);
					if(((int)Cutmask&passed[S])==passed[S]) ((TH3*)EffpreSelMCD->afterR)->Fill(M,(int)(10000*Massa_gen-18570),S);
			}
			for(int m=0;m<18;m++) if(Var>BetaD[m]&&Var<=BetaD[m+1]){
					if(((int)Cutmask&notpassed[S])==notpassed[S]) ((TH3*)EffpreSelMCD->beforeTOF)->Fill(m,(int)(10000*Massa_gen-18570),S);
					if(((int)Cutmask&passed[S])==passed[S]) ((TH3*)EffpreSelMCD->afterTOF)->Fill(m,(int)(10000*Massa_gen-18570),S);
			}
		}
	}
	return;
}


void MCpreSeleff_Write(){
        EffpreSelMCP -> Write();
        EffpreSelMCD ->Write();
	return; 
}



TCanvas *c9[4];

void MCpreSeleff(TFile * file1){
	
	Efficiency * EffpreSelMCP = new Efficiency(file1,"EffpreSelMCP"); 
	Efficiency * EffpreSelMCD = new Efficiency(file1,"EffpreSelMCD");

	string numero[18]={"0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17"};
	string tagli[3]={"Matching TOF","Chi^2 R","1 Tr. Track"};
	string nome;

	cout<<"********* MC \"GOLDEN\" SEL. EFFICIENCIES *********"<<endl;

	EffpreSelMCP -> Eval_Efficiency();
	EffpreSelMCD -> Eval_Efficiency();

	TH2F * EffPreSelMCP_R_TH2F	= (TH2F *) EffpreSelMCP -> effR -> Clone();
	TH2F * EffPreSelMCP_TH2F	= (TH2F *) EffpreSelMCP -> effTOF-> Clone();
	TH3F * EffPreSelMCD_R_TH3F	= (TH3F *) EffpreSelMCD -> effR -> Clone();
	TH3F * EffPreSelMCD_TH3F	= (TH3F *) EffpreSelMCD -> effTOF-> Clone();

	cout<<"*** Updating P1 file ****"<<endl;
        string nomefile=percorso + "/Risultati/risultati/"+mese+"_"+frac+"_P1.root";
        file1 =TFile::Open(nomefile.c_str(),"UPDATE");
        if(!file1){
                nomefile=percorso + "/Risultati/"+mese+"/"+mese+"_"+frac+"_P1.root";
                file1 =TFile::Open(nomefile.c_str(),"UPDATE");
        }

        file1->cd("Results");
	EffPreSelMCP_R_TH2F	-> Write();
	EffPreSelMCP_TH2F	-> Write();
	EffPreSelMCD_R_TH3F	-> Write();
	EffPreSelMCD_TH3F	-> Write();
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
		TGraph * EffPreSelMCP_R = new TGraph();
		EffPreSelMCP_R->SetTitle(MCLegend[0].c_str());
		for(int i=0;i<43;i++) EffPreSelMCP_R->SetPoint(i,R_cent[i],EffPreSelMCP_R_TH2F->GetBinContent(i+1,S+1));
		TGraph * EffPreSelMCD_R[6][3];
		EffPreSelMCP_R->SetMarkerColor(2);
		EffPreSelMCP_R->SetMarkerStyle(8);
		EffPreSelMCP_R->SetLineColor(2);
		EffPreSelMCP_R->SetLineWidth(2);
		nome="Efficiency: "+tagli[S]+ "MC (R bins)";
		EffPreSelMCP_R->SetTitle(nome.c_str());
		EffPreSelMCP_R->GetXaxis()->SetTitle("R [GV]");
		EffPreSelMCP_R->GetYaxis()->SetTitle("Pres. Efficiency");
		EffPreSelMCP_R->GetXaxis()->SetTitleSize(0.045);
		EffPreSelMCP_R->GetYaxis()->SetTitleSize(0.045);
		{
			EffPreSelMCP_R->Draw("ACP");
			TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
			leg->AddEntry(EffPreSelMCP_R,MCLegend[0].c_str(), "ep");

			for(int h=0;h<6;h++){
				EffPreSelMCD_R[h][S]= new TGraph();
				EffPreSelMCD_R[h][S]->SetTitle(MCLegend[h+1].c_str());
				for(int i=1;i<43;i++) EffPreSelMCD_R[h][S]->SetPoint(i,R_cent[i],EffPreSelMCD_R_TH3F->GetBinContent(i+1,h+1,S+1));
				leg->AddEntry(EffPreSelMCD_R[h][S],MCLegend[h+1].c_str(), "ep");
				EffPreSelMCD_R[h][S]->SetMarkerColor(4);
				EffPreSelMCD_R[h][S]->SetMarkerStyle(h+3);
				EffPreSelMCD_R[h][S]->SetMarkerSize(2);
				EffPreSelMCD_R[h][S]->SetLineColor(4);
				EffPreSelMCD_R[h][S]->SetLineWidth(2);
				EffPreSelMCD_R[h][S]->Draw("Psame");
				leg->Draw();
			}
		}

		c9[S]->cd(2);
		gPad->SetLogx();
		gPad->SetGridx();
		gPad->SetGridy();
		TGraph * EffPreSelMCP = new TGraph();
		for(int i=0;i<18;i++) EffPreSelMCP->SetPoint(i,Ekincent[i],EffPreSelMCP_TH2F->GetBinContent(i+1,S+1));
		TGraph * EffPreSelMCD[6][3];
		EffPreSelMCP->SetMarkerColor(2);
		EffPreSelMCP->SetMarkerStyle(8);
		EffPreSelMCP->SetLineColor(2);
		EffPreSelMCP->SetLineWidth(2);
		EffPreSelMCP->SetTitle("Preselections Efficiency MC (Beta bins)");
		EffPreSelMCP->GetXaxis()->SetTitle("Kin. En. / nucl. [GeV/nucl.]");
		EffPreSelMCP->GetYaxis()->SetTitle("Pres. Efficiency");
		EffPreSelMCP->GetXaxis()->SetTitleSize(0.045);
		EffPreSelMCP->GetYaxis()->SetTitleSize(0.045);
		{
			EffPreSelMCP->Draw("ACP");
			TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
			leg->AddEntry(EffPreSelMCP,MCLegend[0].c_str(), "ep");

			for(int h=0;h<6;h++){
				EffPreSelMCD[h][S]= new TGraph();
				for(int i=0;i<18;i++) EffPreSelMCD[h][S]->SetPoint(i,Ekincent[i],EffPreSelMCD_TH3F->GetBinContent(i+1,h+1,S+1));
				EffPreSelMCD[h][S]->SetMarkerColor(4);
				EffPreSelMCD[h][S]->SetMarkerStyle(h+3);
				leg->AddEntry(EffPreSelMCD[h][S],MCLegend[h+1].c_str(), "ep");
				EffPreSelMCD[h][S]->SetMarkerSize(2);
				EffPreSelMCD[h][S]->SetLineColor(4);
				EffPreSelMCD[h][S]->SetLineWidth(2);
				EffPreSelMCD[h][S]->Draw("Psame");
				leg->Draw();
			}
		}
	}
return;
}

