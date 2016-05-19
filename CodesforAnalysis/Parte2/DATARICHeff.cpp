using namespace std;

LATcorr * LATrichDATA_NaF   = new LATcorr("LATrichDATA_NaF");
LATcorr * LATrichDATA_Agl   = new LATcorr("LATrichDATA_Agl");


void DATARICHeff_Fill(TNtuple *ntupla, int l,int zona) {

	 ntupla->GetEvent(l);
	//cuts
	if(Beta<=0||R<=0||R<1.2*Rcutoff||Beta>protons->Eval(R)+0.1||Beta<protons->Eval(R)-0.1) return;
	if(!((R>Rcut[zona]&&zona<10)||(zona==10)))  return;
	if(!Herejcut) return;

	int Kbin=RB.GetRBin(R);

	LATrichDATA_NaF -> beforeR -> Fill(Kbin,zona);
	LATrichDATA_Agl	-> beforeR -> Fill(Kbin,zona);
	
	//NaF
	if(((int)Cutmask>>11)==512) LATrichDATA_NaF -> afterR -> Fill(Kbin,zona); 
   	
   	//Agl
	if(((int)Cutmask>>11)==0)   LATrichDATA_Agl -> afterR -> Fill(Kbin,zona); 

   	return;
}

void DATARICHeff_Write() {
	LATrichDATA_NaF -> Write();
	LATrichDATA_Agl -> Write();  
	return;
}




void DATARICHeff(TFile * file1) {
   LATcorr * LATrichDATA_NaF = new LATcorr(file1,"LATrichDATA_NaF");
   LATcorr * LATrichDATA_Agl = new LATcorr(file1,"LATrichDATA_Agl");


   cout<<"****************************** DATA RICH SEL. EFFICIENCIES **************************************"<<endl;


   LATrichDATA_NaF -> Eval_Efficiency();
   LATrichDATA_Agl -> Eval_Efficiency();


   TH2F *LATrichDATANaF = (TH2F *)  LATrichDATA_NaF  -> effR -> Clone();
   TH2F *LATrichDATAAgl = (TH2F *)  LATrichDATA_Agl  -> effR -> Clone();


   cout<<"****************************** LAT. Eff. CORRECTION *************************************************"<<endl;

   LATrichDATA_NaF ->  Eval_LATcorr(1);
   LATrichDATA_Agl ->  Eval_LATcorr(1);

   TH2F *LATrichcorr_NaF  =	(TH2F *) LATrichDATA_NaF  -> LATcorrR -> Clone();
   TH2F *LATrichcorr_Agl  =	(TH2F *) LATrichDATA_Agl  -> LATcorrR -> Clone();

   TH1F *LATrichcorr_NaF_fit  	= (TH1F *) LATrichDATA_NaF   -> LATcorrR_fit-> Clone();
   TH1F *LATrichcorr_Agl_fit	= (TH1F *) LATrichDATA_Agl   -> LATcorrR_fit-> Clone();

   cout<<"*** Updating P1 file ****"<<endl;
   string nomefile= "../Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
   file1 =TFile::Open(nomefile.c_str(),"UPDATE");

   file1->cd("Results");

   LATrichDATANaF -> Write();
   LATrichDATAAgl -> Write();

   LATrichcorr_NaF-> Write();
   LATrichcorr_Agl-> Write();
	
   LATrichcorr_NaF_fit-> Write();
   LATrichcorr_Agl_fit-> Write();

   file1->Write();
   file1->Close();



   TCanvas *c15=new TCanvas("Latitude RICH  Efficiency");

   c15->Divide(2,2);

   c15->cd(1);
   gPad->SetLogx();
   gPad->SetGridx();
   gPad->SetGridy();
   TGraphErrors *EffDATArichNaF[11];
   for(int l=0; l<11; l++) {
      EffDATArichNaF[l]=new TGraphErrors();
      int j=0;
      for(int i=1; i<nbinsr; i++) {
	      if(LATrichDATANaF->GetBinContent(i+1,l+1)>0) {
		      EffDATArichNaF[l]->SetPoint(j,R_cent[i],LATrichDATANaF->GetBinContent(i+1,l+1));
		      EffDATArichNaF[l]->SetPointError(j,0,   LATrichDATANaF->GetBinError(i+1,l+1));
		      j++;}
      }
   }
   EffDATArichNaF[10]->SetMarkerColor(1);
   EffDATArichNaF[10]->SetMarkerStyle(8);
   EffDATArichNaF[10]->SetLineColor(1);
   EffDATArichNaF[10]->SetTitle("Latitude RICH Efficiency (NaF)");
   EffDATArichNaF[10]->GetXaxis()->SetTitle("R [GV]");
   EffDATArichNaF[10]->GetYaxis()->SetTitle("Efficiency");
   EffDATArichNaF[10]->GetXaxis()->SetTitleSize(0.045);
   EffDATArichNaF[10]->GetYaxis()->SetTitleSize(0.045);
   EffDATArichNaF[10]->GetYaxis()->SetRangeUser(0.0,0.1);
   EffDATArichNaF[10]->Draw("AP");
   for(int l=0; l<10; l++) {
      EffDATArichNaF[l]->SetMarkerColor(l);
      EffDATArichNaF[l]->SetMarkerStyle(8);
      EffDATArichNaF[l]->SetLineColor(l);
      EffDATArichNaF[l]->Draw("Psame");
   }

   c15->cd(3);
   gPad->SetLogx();
   gPad->SetGridx();
   gPad->SetGridy();
   TGraphErrors *EffDATArichAgl[11];
   for(int l=0; l<11; l++) {
      EffDATArichAgl[l]=new TGraphErrors();
      int j=0;
      for(int i=1; i<nbinsr; i++) {
	      if(LATrichDATAAgl ->GetBinContent(i+1,l+1)>0){
		      EffDATArichAgl[l]->SetPoint(j,R_cent[i],LATrichDATAAgl ->GetBinContent(i+1,l+1));
		      EffDATArichAgl[l]->SetPointError(j,0,   LATrichDATAAgl ->GetBinError(i+1,l+1));
		      j++;
	      }
      }
   }
   EffDATArichAgl[10]->SetMarkerColor(1);
   EffDATArichAgl[10]->SetMarkerStyle(8);
   EffDATArichAgl[10]->SetLineColor(1);
   EffDATArichAgl[10]->SetTitle("Latitude RICH Efficiency (Agl)");
   EffDATArichAgl[10]->GetXaxis()->SetTitle("R [GV]");
   EffDATArichAgl[10]->GetYaxis()->SetTitle("Efficiency");
   EffDATArichAgl[10]->GetXaxis()->SetTitleSize(0.045);
   EffDATArichAgl[10]->GetYaxis()->SetTitleSize(0.045);
   EffDATArichAgl[10]->GetYaxis()->SetRangeUser(0.1,1.1);
   EffDATArichAgl[10]->Draw("AP");
   for(int l=0; l<10; l++) {
      EffDATArichAgl[l]->SetMarkerColor(l);
      EffDATArichAgl[l]->SetMarkerStyle(8);
      EffDATArichAgl[l]->SetLineColor(l);
      EffDATArichAgl[l]->Draw("Psame");
   }


   TGraphErrors *CorrLATrichNaF;
   c15->cd(2);
   gPad->SetGridy();
   gPad->SetGridx();
   CorrLATrichNaF=new TGraphErrors();
   CorrLATrichNaF->SetTitle("Latitude Efficiency Corr.");
   CorrLATrichNaF->GetXaxis()->SetTitle("Latitude");
   CorrLATrichNaF->GetYaxis()->SetTitle("Eff. Corr. Factor");
   CorrLATrichNaF->GetYaxis()->SetRangeUser(0.96,1.04);
   CorrLATrichNaF->SetMarkerStyle(8);
   for(int i=1; i<11; i++) {
      CorrLATrichNaF->SetPoint(i-1,geomagC[i],LATrichcorr_NaF->GetBinContent(i+1,1));
      CorrLATrichNaF->SetPointError(i-1,0,LATrichcorr_NaF->GetBinError(i+1,1));
   }
   CorrLATrichNaF->Draw("AP");

   TGraphErrors *CorrLAT_richNaF_Spl=new TGraphErrors("CorrLAT_richNaF_Spl");
   CorrLAT_richNaF_Spl->SetName("CorrLAT_richNaF_Spl");
   for(int i=1; i<11; i++) {
      CorrLAT_richNaF_Spl->SetPoint(i-1,geomagC[i],LATrichcorr_NaF_fit->GetBinContent(i+1));
      CorrLAT_richNaF_Spl->SetPointError(i-1,0,LATrichcorr_NaF_fit->GetBinError(i+1));
   }
   CorrLAT_richNaF_Spl->SetLineColor(2);
   CorrLAT_richNaF_Spl->SetMarkerColor(2);
   CorrLAT_richNaF_Spl->SetFillColor(2);
   CorrLAT_richNaF_Spl->SetFillStyle(3001);
   CorrLAT_richNaF_Spl->SetLineWidth(2);
   CorrLAT_richNaF_Spl->Draw("Csame");

   TGraphErrors *CorrLATrichAgl;
   c15->cd(4);
   gPad->SetGridy();
   gPad->SetGridx();
   CorrLATrichAgl=new TGraphErrors();
   CorrLATrichAgl->SetTitle("Latitude Efficiency Corr.");
   CorrLATrichAgl->GetXaxis()->SetTitle("Latitude");
   CorrLATrichAgl->GetYaxis()->SetTitle("Eff. Corr. Factor");
   CorrLATrichAgl->GetYaxis()->SetRangeUser(0.96,1.04);
   CorrLATrichAgl->SetMarkerStyle(8);
   for(int i=1; i<11; i++) {
      CorrLATrichAgl->SetPoint(i-1,geomagC[i],LATrichcorr_Agl->GetBinContent(i+1,1));
      CorrLATrichAgl->SetPointError(i-1,0,LATrichcorr_Agl->GetBinError(i+1,1));
   }
   CorrLATrichAgl->Draw("AP");

   TGraphErrors *CorrLAT_richAgl_Spl=new TGraphErrors("CorrLAT_richAgl_Spl");
   CorrLAT_richAgl_Spl->SetName("CorrLAT_richAgl_Spl");
   for(int i=1; i<11; i++) {
      CorrLAT_richAgl_Spl->SetPoint(i-1,geomagC[i],LATrichcorr_Agl_fit->GetBinContent(i+1));
      CorrLAT_richAgl_Spl->SetPointError(i-1,0,LATrichcorr_Agl_fit->GetBinError(i+1));
   }
   CorrLAT_richAgl_Spl->SetLineColor(2);
   CorrLAT_richAgl_Spl->SetMarkerColor(2);
   CorrLAT_richAgl_Spl->SetFillColor(2);
   CorrLAT_richAgl_Spl->SetFillStyle(3001);
   CorrLAT_richAgl_Spl->SetLineWidth(2);
   CorrLAT_richAgl_Spl->Draw("Csame");
   
   cout<<"*** Updating Results file ***"<<endl;
   nomefile="./Final_plots/"+mese+".root";
   TFile *f_out=new TFile(nomefile.c_str(), "UPDATE");
   f_out->mkdir("DATA-driven Results/Latitude effect/RICH");
   f_out->cd("DATA-driven Results/Latitude effect/RICH");
   c15->Write();
   f_out->cd("Export");
   CorrLAT_richNaF_Spl->Write("CorrLAT_richNaF_Spl");
   CorrLAT_richAgl_Spl->Write("CorrLAT_richAgl_Spl");   	
   f_out->Write();
   f_out->Close();

}



