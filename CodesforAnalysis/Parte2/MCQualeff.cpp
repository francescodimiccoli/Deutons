using namespace std;

// eff. likelihood sel.
Efficiency * EffLikMCP  = new Efficiency ("EffLikMCP");
Efficiency * EffLikMCD  = new Efficiency ("EffLikMCD",6);
// eff. distance sel.
Efficiency * EffDistMCP = new Efficiency ("EffDistMCP");
Efficiency * EffDistMCD = new Efficiency ("EffDistMCD",6);

void MCQualeff_Fill(TNtuple *ntupla, int l) {

	 ntupla->GetEvent(l);
	if(Beta<=0||R<=0) return;
	int Kbin;
	if(Massa_gen<1) {
		//R bins
		Kbin=RB.GetRBin(R);

		EffLikMCP->beforeR->Fill(Kbin);
		EffDistMCP->beforeR->Fill(Kbin);
		if(Dist5D_P<6)          EffDistMCP->afterR->Fill(Kbin);
		if(Dist5D_P<6&&Likcut)	EffLikMCP->afterR->Fill(Kbin);


		//Beta bins
		Kbin=ToFPB.GetRBin(Var);

		EffLikMCP ->beforeTOF->Fill(Kbin);
		EffDistMCP->beforeTOF->Fill(Kbin);
		if(Distcut)             EffDistMCP ->afterTOF ->Fill(Kbin);
		if(Distcut&&Likcut)	EffLikMCP->afterTOF ->Fill(Kbin);


		if(((int)Cutmask)>>11 ==512) {
			Kbin=NaFPB.GetRBin(Var2);
			EffLikMCP  ->beforeNaF->Fill(Kbin);
			EffDistMCP ->beforeNaF->Fill(Kbin);
			if(Distcut)             EffDistMCP  ->afterNaF ->Fill(Kbin);
			if(Distcut&&Likcut)	EffLikMCP ->afterNaF ->Fill(Kbin);
		}

		if(((int)Cutmask)>>11==0) {
			Kbin=AglPB.GetRBin(Var2);
			EffLikMCP->beforeAgl->Fill(Kbin);
			EffDistMCP->beforeAgl->Fill(Kbin);
			if(Distcut) 		EffDistMCP->afterAgl->Fill(Kbin);
			if(Distcut&&Likcut)	EffLikMCP->afterAgl->Fill(Kbin);
		}


	}
	if(Massa_gen<2&&Massa_gen>1) {
		//R bins
		Kbin=RB.GetRBin(R);

		FillBinMGen(EffLikMCD ->beforeR, Kbin);
		FillBinMGen(EffDistMCD->beforeR, Kbin);
		if(Distcut)             FillBinMGen(EffDistMCD ->afterR,  Kbin);
		if(Dist5D<6&&Likcut)	FillBinMGen(EffLikMCD->afterR,  Kbin);


		//Beta bins
		Kbin=ToFDB.GetRBin(Var);
		FillBinMGen(EffLikMCD ->beforeTOF, Kbin);
		FillBinMGen(EffDistMCD->beforeTOF, Kbin);
		if(Distcut)             FillBinMGen(EffDistMCD ->afterTOF , Kbin);
		if(Distcut&&Likcut)	FillBinMGen(EffLikMCD->afterTOF , Kbin);

		if(((int)Cutmask)>>11==512) {
			Kbin=NaFDB.GetRBin(Var2);
			FillBinMGen(EffLikMCD ->beforeNaF, Kbin);
			FillBinMGen(EffDistMCD->beforeNaF, Kbin);
			if(Distcut)             FillBinMGen(EffDistMCD ->afterNaF , Kbin);
			if(Distcut&&Likcut)	FillBinMGen(EffLikMCD->afterNaF , Kbin);
		}

		if(((int)Cutmask)>>11==0) {
			Kbin=AglDB.GetRBin(Var2);
			FillBinMGen(EffLikMCD ->beforeAgl, Kbin);
			FillBinMGen(EffDistMCD->beforeAgl, Kbin);
			if(Distcut)             FillBinMGen(EffDistMCD ->afterAgl , Kbin);
			if(Distcut&&Likcut)	FillBinMGen(EffLikMCD->afterAgl , Kbin);
		}
	}

return;
}




void MCQualeff_Write() {
   EffLikMCP  -> Write();
   EffLikMCD  -> Write();
   EffDistMCP -> Write();
   EffDistMCD -> Write();
   return;
}


void MCQualeff(TFile * file1) {

   // eff. likelihood sel.
   Efficiency * EffLikMCP  = new Efficiency (file1,"EffLikMCP");
   Efficiency * EffLikMCD  = new Efficiency (file1,"EffLikMCD");
   // eff. distance sel.
   Efficiency * EffDistMCP = new Efficiency (file1,"EffDistMCP");
   Efficiency * EffDistMCD = new Efficiency (file1,"EffDistMCD");


   cout<<"******* MC QUALITY SEL. EFFICIENCIES ********"<<endl;


   EffLikMCP ->Eval_Efficiency();
   EffLikMCD ->Eval_Efficiency();

   EffDistMCP->Eval_Efficiency();
   EffDistMCD->Eval_Efficiency();


   TH1F * EffMCLikP_TH1F 		=(TH1F *)EffLikMCP ->effR  ->Clone();
   TH2F * EffMCLikD_TH2F 		=(TH2F *)EffLikMCD ->effR  ->Clone();
   TH1F * EffMCLikP_Beta_TH1F 	=(TH1F *)EffLikMCP ->effTOF->Clone();
   TH2F * EffMCLikD_Beta_TH2F 	=(TH2F *)EffLikMCD ->effTOF->Clone();
   TH1F * EffMCLikP_BetaNaF_TH1F 	=(TH1F *)EffLikMCP ->effNaF->Clone();
   TH2F * EffMCLikD_BetaNaF_TH2F 	=(TH2F *)EffLikMCD ->effNaF->Clone();
   TH1F * EffMCLikP_BetaAgl_TH1F 	=(TH1F *)EffLikMCP ->effAgl->Clone();
   TH2F * EffMCLikD_BetaAgl_TH2F 	=(TH2F *)EffLikMCD ->effAgl->Clone();

   TH1F * EffMCDistP_TH1F 		=(TH1F *)EffDistMCP->effR  ->Clone();
   TH2F * EffMCDistD_TH2F 		=(TH2F *)EffDistMCD->effR  ->Clone();
   TH1F * EffMCDistP_Beta_TH1F 	=(TH1F *)EffDistMCP->effTOF->Clone();
   TH2F * EffMCDistD_Beta_TH2F 	=(TH2F *)EffDistMCD->effTOF->Clone();
   TH1F * EffMCDistP_BetaNaF_TH1F  =(TH1F *)EffDistMCP->effNaF->Clone();
   TH2F * EffMCDistD_BetaNaF_TH2F  =(TH2F *)EffDistMCD->effNaF->Clone();
   TH1F * EffMCDistP_BetaAgl_TH1F  =(TH1F *)EffDistMCP->effAgl->Clone();
   TH2F * EffMCDistD_BetaAgl_TH2F  =(TH2F *)EffDistMCD->effAgl->Clone();


   cout<<"*** Updating P1 file ****"<<endl;
   string nomefile="../Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
   file1 =TFile::Open(nomefile.c_str(),"UPDATE");

   file1->cd("Results");

   EffMCLikP_TH1F 	  ->Write();
   EffMCLikD_TH2F 	  ->Write();
   EffMCLikP_Beta_TH1F      ->Write();
   EffMCLikD_Beta_TH2F      ->Write();
   EffMCLikP_BetaNaF_TH1F   ->Write();
   EffMCLikD_BetaNaF_TH2F   ->Write();
   EffMCLikP_BetaAgl_TH1F   ->Write();
   EffMCLikD_BetaAgl_TH2F   ->Write();

   EffMCDistP_TH1F 	  ->Write();
   EffMCDistD_TH2F 	  ->Write();
   EffMCDistP_Beta_TH1F     ->Write();
   EffMCDistD_Beta_TH2F     ->Write();
   EffMCDistP_BetaNaF_TH1F  ->Write();
   EffMCDistD_BetaNaF_TH2F  ->Write();
   EffMCDistP_BetaAgl_TH1F  ->Write();
   EffMCDistD_BetaAgl_TH2F  ->Write();
   file1->Write();
   file1->Close();



   TCanvas *c5	=new TCanvas("Likelihood Efficiency (R bins)");
   TCanvas *c6	=new TCanvas("Distance Efficiency (R bins)");
   TCanvas *c5_bis	=new TCanvas("Likelihood Efficiency (Beta bins)");
   TCanvas *c6_bis	=new TCanvas("Distance Efficiency (Beta bins)");
   c5->cd();
   string MCLegend[7]= {"protons.B800","d.pl1.0_520_GG_Blic","d.pl1.0_520_GG_BlicDPMJet","d.pl1.0_520_GG_QMD","d.pl1.0_520_Shen_Blic","d.pl1.0_520_Shen_BlicDPMJet","d.pl1.0_520_Shen_QMD"};
   gPad->SetLogx();
   gPad->SetGridx();
   gPad->SetGridy();
   TGraph *EffMCLikP= new TGraph();
   TGraph *EffMCLikD[6];
   for(int i=1; i<nbinsr; i++) EffMCLikP->SetPoint(i,R_cent[i],EffMCLikP_TH1F->GetBinContent(i+1));
   EffMCLikP->SetMarkerColor(2);
   EffMCLikP->SetMarkerStyle(8);
   EffMCLikP->SetLineColor(2);
   EffMCLikP->SetLineWidth(2);
   EffMCLikP->SetTitle("Likelihood Efficiency MC on top of Pres. (R bins)");
   EffMCLikP->GetXaxis()->SetTitle("R [GV]");
   EffMCLikP->GetYaxis()->SetTitle("Efficiency");
   EffMCLikP->GetXaxis()->SetTitleSize(0.045);
   EffMCLikP->GetYaxis()->SetTitleSize(0.045);
   {
      EffMCLikP->GetYaxis()->SetRangeUser(0,1);
      EffMCLikP->Draw("ACP");
      TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
      leg->AddEntry(EffMCLikP,MCLegend[0].c_str(), "ep");

      for(int h=0; h<6; h++) {
         EffMCLikD[h]= new TGraph();
         EffMCLikD[h]->SetTitle(MCLegend[h+1].c_str());
         for(int i=1; i<nbinsr; i++) EffMCLikD[h]->SetPoint(i,R_cent[i], EffMCLikD_TH2F->GetBinContent(i+1,h+1));
         //leg->AddEntry(EffMCLikD[h],MCLegend[h+1].c_str(), "ep");
         EffMCLikD[h]->SetMarkerColor(4);
         EffMCLikD[h]->SetMarkerStyle(h+3);
         EffMCLikD[h]->SetMarkerSize(2);
         EffMCLikD[h]->SetLineColor(4);
         EffMCLikD[h]->SetLineWidth(1);
         EffMCLikD[h]->Draw("Psame");
         leg->Draw();
      }
   }

   c5_bis->Divide(3,1);
   c5_bis->cd(1);
   gPad->SetLogx();
   gPad->SetGridx();
   gPad->SetGridy();
   TGraph *EffMCLikP_Beta= new TGraph();
   TGraph *EffMCLikD_Beta[6];
   for(int i=0; i<nbinsToF; i++)  EffMCLikP_Beta->SetPoint(i,Ekincent[i],EffMCLikP_Beta_TH1F->GetBinContent(i+1));
   EffMCLikP_Beta->SetMarkerColor(2);
   EffMCLikP_Beta->SetMarkerStyle(8);
   EffMCLikP_Beta->SetLineColor(2);
   EffMCLikP_Beta->SetLineWidth(2);
   EffMCLikP_Beta->SetTitle("Likelihood Efficiency MC on top of Pres. (Beta bins)");
   EffMCLikP_Beta->GetXaxis()->SetTitle("Kin. En.  [GeV/nucl.]");
   EffMCLikP_Beta->GetYaxis()->SetTitle("Efficiency");
   EffMCLikP_Beta->GetXaxis()->SetTitleSize(0.045);
   EffMCLikP_Beta->GetYaxis()->SetTitleSize(0.045);
   {
      EffMCLikP_Beta->GetYaxis()->SetRangeUser(0,1);
      EffMCLikP_Beta->Draw("ACP");
      TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
      leg->AddEntry(EffMCLikP_Beta,MCLegend[0].c_str(), "ep");

      for(int h=0; h<6; h++) {
         EffMCLikD_Beta[h]= new TGraph();
         EffMCLikD_Beta[h]->SetTitle(MCLegend[h+1].c_str());
         for(int i=1; i<nbinsToF; i++) EffMCLikD_Beta[h]->SetPoint(i,Ekincent[i],EffMCLikD_Beta_TH2F->GetBinContent(i+1,h+1));
         leg->AddEntry(EffMCLikD_Beta[h],MCLegend[h+1].c_str(), "ep");
         EffMCLikD_Beta[h]->SetMarkerColor(4);
         EffMCLikD_Beta[h]->SetMarkerStyle(h+3);
         EffMCLikD_Beta[h]->SetMarkerSize(2);
         EffMCLikD_Beta[h]->SetLineColor(4);
         EffMCLikD_Beta[h]->SetLineWidth(1);
         EffMCLikD_Beta[h]->Draw("Psame");
         leg->Draw();
      }
   }
   c5_bis->cd(2);
   gPad->SetLogx();
   gPad->SetGridx();
   gPad->SetGridy();
   TGraph *EffMCLikP_BetaNaF= new TGraph();
   TGraph *EffMCLikD_BetaNaF[6];
   for(int i=0; i<nbinsNaF; i++) EffMCLikP_BetaNaF->SetPoint(i,EkincentNaF[i], EffMCLikP_BetaNaF_TH1F->GetBinContent(i+1));
   EffMCLikP_BetaNaF->SetMarkerColor(2);
   EffMCLikP_BetaNaF->SetMarkerStyle(8);
   EffMCLikP_BetaNaF->SetLineColor(2);
   EffMCLikP_BetaNaF->SetLineWidth(2);
   EffMCLikP_BetaNaF->SetTitle("Likelihood Efficiency MC on top of Pres. (Beta bins)");
   EffMCLikP_BetaNaF->GetXaxis()->SetTitle("Kin. En.  [GeV/nucl.]");
   EffMCLikP_BetaNaF->GetYaxis()->SetTitle("Efficiency");
   EffMCLikP_BetaNaF->GetXaxis()->SetTitleSize(0.045);
   EffMCLikP_BetaNaF->GetYaxis()->SetTitleSize(0.045);
   {
      EffMCLikP_BetaNaF->GetYaxis()->SetRangeUser(0,1);
      EffMCLikP_BetaNaF->Draw("ACP");
      TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
      leg->AddEntry(EffMCLikP_BetaNaF,MCLegend[0].c_str(), "ep");

      for(int h=0; h<6; h++) {
         EffMCLikD_BetaNaF[h]= new TGraph();
         EffMCLikD_BetaNaF[h]->SetTitle(MCLegend[h+1].c_str());
         for(int i=1; i<nbinsNaF; i++) EffMCLikD_BetaNaF[h]->SetPoint(i,EkincentNaF[i],EffMCLikD_BetaNaF_TH2F->GetBinContent(i+1,h+1));
         leg->AddEntry(EffMCLikD_BetaNaF[h],MCLegend[h+1].c_str(), "ep");
         EffMCLikD_BetaNaF[h]->SetMarkerColor(4);
         EffMCLikD_BetaNaF[h]->SetMarkerStyle(h+3);
         EffMCLikD_BetaNaF[h]->SetMarkerSize(2);
         EffMCLikD_BetaNaF[h]->SetLineColor(4);
         EffMCLikD_BetaNaF[h]->SetLineWidth(1);
         EffMCLikD_BetaNaF[h]->Draw("Psame");
         leg->Draw();
      }
   }
   c5_bis->cd(3);
   gPad->SetLogx();
   gPad->SetGridx();
   gPad->SetGridy();
   TGraph *EffMCLikP_BetaAgl= new TGraph();
   TGraph *EffMCLikD_BetaAgl[6];
   for(int i=0; i<nbinsAgl; i++) EffMCLikP_BetaAgl->SetPoint(i,EkincentAgl[i],EffMCLikP_BetaAgl_TH1F->GetBinContent(i+1));
   EffMCLikP_BetaAgl->SetMarkerColor(2);
   EffMCLikP_BetaAgl->SetMarkerStyle(8);
   EffMCLikP_BetaAgl->SetLineColor(2);
   EffMCLikP_BetaAgl->SetLineWidth(2);
   EffMCLikP_BetaAgl->SetTitle("Likelihood Efficiency MC on top of Pres. (Beta bins)");
   EffMCLikP_BetaAgl->GetXaxis()->SetTitle("Kin. En.  [GeV/nucl.]");
   EffMCLikP_BetaAgl->GetYaxis()->SetTitle("Efficiency");
   EffMCLikP_BetaAgl->GetXaxis()->SetTitleSize(0.045);
   EffMCLikP_BetaAgl->GetYaxis()->SetTitleSize(0.045);
   {
      EffMCLikP_BetaAgl->GetYaxis()->SetRangeUser(0,1);
      EffMCLikP_BetaAgl->Draw("ACP");
      TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
      leg->AddEntry(EffMCLikP_BetaAgl,MCLegend[0].c_str(), "ep");

      for(int h=0; h<6; h++) {
         EffMCLikD_BetaAgl[h]= new TGraph();
         EffMCLikD_BetaAgl[h]->SetTitle(MCLegend[h+1].c_str());
         for(int i=1; i<nbinsAgl; i++) EffMCLikD_BetaAgl[h]->SetPoint(i,EkincentAgl[i], EffMCLikD_BetaAgl_TH2F->GetBinContent(i+1,h+1));
         leg->AddEntry(EffMCLikD_BetaAgl[h],MCLegend[h+1].c_str(), "ep");
         EffMCLikD_BetaAgl[h]->SetMarkerColor(4);
         EffMCLikD_BetaAgl[h]->SetMarkerStyle(h+3);
         EffMCLikD_BetaAgl[h]->SetMarkerSize(2);
         EffMCLikD_BetaAgl[h]->SetLineColor(4);
         EffMCLikD_BetaAgl[h]->SetLineWidth(1);
         EffMCLikD_BetaAgl[h]->Draw("Psame");
         leg->Draw();
      }
   }
   c6->cd();
   gPad->SetLogx();
   gPad->SetGridx();
   gPad->SetGridy();
   TGraph *EffMCDistP= new TGraph();
   TGraph *EffMCDistD[6];
   for(int i=1; i<nbinsr; i++) EffMCDistP->SetPoint(i,R_cent[i],EffMCDistP_TH1F->GetBinContent(i+1));
   EffMCDistP->SetMarkerColor(2);
   EffMCDistP->SetMarkerStyle(8);
   EffMCDistP->SetLineColor(2);
   EffMCDistP->SetLineWidth(2);
   EffMCDistP->SetTitle("Distance Efficiency MC on top of Pres. (R bins)");
   EffMCDistP->GetXaxis()->SetTitle("R [GV]");
   EffMCDistP->GetYaxis()->SetTitle("Efficiency");
   EffMCDistP->GetXaxis()->SetTitleSize(0.045);
   EffMCDistP->GetYaxis()->SetTitleSize(0.045);
   {
      EffMCDistP->GetYaxis()->SetRangeUser(0,1);
      EffMCDistP->Draw("ACP");
      TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
      leg->AddEntry(EffMCDistP,MCLegend[0].c_str(), "ep");

      for(int h=0; h<6; h++) {
         EffMCDistD[h]= new TGraph();
         EffMCDistD[h]->SetTitle(MCLegend[h+1].c_str());
         for(int i=1; i<nbinsr; i++) EffMCDistD[h]->SetPoint(i,R_cent[i],EffMCDistD_TH2F->GetBinContent(i+1,h+1));
         leg->AddEntry(EffMCDistD[h],MCLegend[h+1].c_str(), "ep");
         EffMCDistD[h]->SetMarkerColor(4);
         EffMCDistD[h]->SetMarkerStyle(h+3);
         EffMCDistD[h]->SetMarkerSize(2);
         EffMCDistD[h]->SetLineColor(4);
         EffMCDistD[h]->SetLineWidth(1);
         EffMCDistD[h]->Draw("Psame");
         leg->Draw();
      }
   }

   c6_bis->Divide(3,1);
   c6_bis->cd(1);
   gPad->SetLogx();
   gPad->SetGridx();
   gPad->SetGridy();
   TGraph *EffMCDistP_Beta= new TGraph();
   TGraph *EffMCDistD_Beta[6];
   for(int i=0; i<nbinsToF; i++) EffMCDistP_Beta->SetPoint(i,Ekincent[i],EffMCDistP_Beta_TH1F->GetBinContent(i+1));
   EffMCDistP_Beta->SetMarkerColor(2);
   EffMCDistP_Beta->SetMarkerStyle(8);
   EffMCDistP_Beta->SetLineColor(2);
   EffMCDistP_Beta->SetLineWidth(2);
   EffMCDistP_Beta->SetTitle("Distance Efficiency MC on top of Pres. (Beta bins)");
   EffMCDistP_Beta->GetXaxis()->SetTitle("Kin. En.  [GeV/nucl.]");
   EffMCDistP_Beta->GetYaxis()->SetTitle("Efficiency");
   EffMCDistP_Beta->GetXaxis()->SetTitleSize(0.045);
   EffMCDistP_Beta->GetYaxis()->SetTitleSize(0.045);
   {
      EffMCDistP_Beta->GetYaxis()->SetRangeUser(0,1);
      EffMCDistP_Beta->Draw("ACP");
      TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
      leg->AddEntry(EffMCDistP_Beta,MCLegend[0].c_str(), "ep");

      for(int h=0; h<6; h++) {
         EffMCDistD_Beta[h]= new TGraph();
         EffMCDistD_Beta[h]->SetTitle(MCLegend[h+1].c_str());
         for(int i=1; i<nbinsToF; i++) EffMCDistD_Beta[h]->SetPoint(i,Ekincent[i],EffMCDistD_Beta_TH2F->GetBinContent(i+1,h+1));
         leg->AddEntry(EffMCDistD_Beta[h],MCLegend[h+1].c_str(), "ep");
         EffMCDistD_Beta[h]->SetMarkerColor(4);
         EffMCDistD_Beta[h]->SetMarkerStyle(h+3);
         EffMCDistD_Beta[h]->SetMarkerSize(2);
         EffMCDistD_Beta[h]->SetLineColor(4);
         EffMCDistD_Beta[h]->SetLineWidth(1);
         EffMCDistD_Beta[h]->Draw("Psame");
         leg->Draw();
      }
   }

   c6_bis->cd(2);
   gPad->SetLogx();
   gPad->SetGridx();
   gPad->SetGridy();
   TGraph *EffMCDistP_BetaNaF= new TGraph();
   TGraph *EffMCDistD_BetaNaF[6];
   for(int i=0; i<nbinsNaF; i++) EffMCDistP_BetaNaF->SetPoint(i,EkincentNaF[i],EffMCDistP_BetaNaF_TH1F->GetBinContent(i+1));
   EffMCDistP_BetaNaF->SetMarkerColor(2);
   EffMCDistP_BetaNaF->SetMarkerStyle(8);
   EffMCDistP_BetaNaF->SetLineColor(2);
   EffMCDistP_BetaNaF->SetLineWidth(2);
   EffMCDistP_BetaNaF->SetTitle("Distelihood Efficiency MC on top of Pres. (Beta bins)");
   EffMCDistP_BetaNaF->GetXaxis()->SetTitle("Kin. En.  [GeV/nucl.]");
   EffMCDistP_BetaNaF->GetYaxis()->SetTitle("Efficiency");
   EffMCDistP_BetaNaF->GetXaxis()->SetTitleSize(0.045);
   EffMCDistP_BetaNaF->GetYaxis()->SetTitleSize(0.045);
   {
      EffMCDistP_BetaNaF->GetYaxis()->SetRangeUser(0,1);
      EffMCDistP_BetaNaF->Draw("ACP");
      TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
      leg->AddEntry(EffMCDistP_BetaNaF,MCLegend[0].c_str(), "ep");

      for(int h=0; h<6; h++) {
         EffMCDistD_BetaNaF[h]= new TGraph();
         EffMCDistD_BetaNaF[h]->SetTitle(MCLegend[h+1].c_str());
         for(int i=1; i<nbinsNaF; i++) EffMCDistD_BetaNaF[h]->SetPoint(i,EkincentNaF[i],EffMCDistD_BetaNaF_TH2F->GetBinContent(i+1,h+1));
         leg->AddEntry(EffMCDistD_BetaNaF[h],MCLegend[h+1].c_str(), "ep");
         EffMCDistD_BetaNaF[h]->SetMarkerColor(4);
         EffMCDistD_BetaNaF[h]->SetMarkerStyle(h+3);
         EffMCDistD_BetaNaF[h]->SetMarkerSize(2);
         EffMCDistD_BetaNaF[h]->SetLineColor(4);
         EffMCDistD_BetaNaF[h]->SetLineWidth(1);
         EffMCDistD_BetaNaF[h]->Draw("Psame");
         leg->Draw();
      }
   }
   c6_bis->cd(3);
   gPad->SetLogx();
   gPad->SetGridx();
   gPad->SetGridy();
   TGraph *EffMCDistP_BetaAgl= new TGraph();
   TGraph *EffMCDistD_BetaAgl[6];
   for(int i=0; i<nbinsAgl; i++) EffMCDistP_BetaAgl->SetPoint(i,EkincentAgl[i],EffMCDistP_BetaAgl_TH1F->GetBinContent(i+1));
   EffMCDistP_BetaAgl->SetMarkerColor(2);
   EffMCDistP_BetaAgl->SetMarkerStyle(8);
   EffMCDistP_BetaAgl->SetLineColor(2);
   EffMCDistP_BetaAgl->SetLineWidth(2);
   EffMCDistP_BetaAgl->SetTitle("Distelihood Efficiency MC on top of Pres. (Beta bins)");
   EffMCDistP_BetaAgl->GetXaxis()->SetTitle("Kin. En.  [GeV/nucl.]");
   EffMCDistP_BetaAgl->GetYaxis()->SetTitle("Efficiency");
   EffMCDistP_BetaAgl->GetXaxis()->SetTitleSize(0.045);
   EffMCDistP_BetaAgl->GetYaxis()->SetTitleSize(0.045);
   {
      EffMCDistP_BetaAgl->GetYaxis()->SetRangeUser(0,1);
      EffMCDistP_BetaAgl->Draw("ACP");
      TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
      leg->AddEntry(EffMCDistP_BetaAgl,MCLegend[0].c_str(), "ep");

      for(int h=0; h<6; h++) {
         EffMCDistD_BetaAgl[h]= new TGraph();
         EffMCDistD_BetaAgl[h]->SetTitle(MCLegend[h+1].c_str());
         for(int i=1; i<nbinsAgl; i++) EffMCDistD_BetaAgl[h]->SetPoint(i,EkincentAgl[i],EffMCDistD_BetaAgl_TH2F->GetBinContent(i+1,h+1));
         leg->AddEntry(EffMCDistD_BetaAgl[h],MCLegend[h+1].c_str(), "ep");
         EffMCDistD_BetaAgl[h]->SetMarkerColor(4);
         EffMCDistD_BetaAgl[h]->SetMarkerStyle(h+3);
         EffMCDistD_BetaAgl[h]->SetMarkerSize(2);
         EffMCDistD_BetaAgl[h]->SetLineColor(4);
         EffMCDistD_BetaAgl[h]->SetLineWidth(1);
         EffMCDistD_BetaAgl[h]->Draw("Psame");
         leg->Draw();
      }
   }

   cout<<"*** Updating Results file ***"<<endl;
   nomefile="./Final_plots/"+mese+".root";
   TFile *f_out=new TFile(nomefile.c_str(), "UPDATE");
   f_out->mkdir("MC Results/Quality");
   f_out->cd("MC Results/Quality");
   c5	->Write();
   c6      ->Write();
   c5_bis  ->Write();
   c6_bis  ->Write();
   f_out->Write();
   f_out->Close();
}



