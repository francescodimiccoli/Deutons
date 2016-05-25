using namespace std;

DatavsMC * Dist_DvsMC_P = new DatavsMC("Dist_DvsMC_P",11);
DatavsMC * Lik_DvsMC_P  = new DatavsMC("Lik_DvsMC_P" ,11);

void DVSMCQualeff2_D_Fill(TNtuple *ntupla, int l,int zona){

	 ntupla->GetEvent(l);
	//cuts
	if(Tup.Beta<=0||Tup.R<=0||Tup.R<1.2*Tup.Rcutoff||Tup.Beta>protons->Eval(Tup.R)+0.1||Tup.Beta<protons->Eval(Tup.R)-0.1) return;
	if(!((Tup.R>Rcut[zona]&&zona<10)||(zona==10)))  return;
	if(!Herejcut) return;
	//
	int Kbin;
	
	//R bins
	Kbin = RB.GetRBin(Tup.R);
	Dist_DvsMC_P -> DataEff -> beforeR -> Fill(Kbin,zona);	
	if(Tup.Dist5D_P<6) Lik_DvsMC_P  -> DataEff -> beforeR -> Fill(Kbin,zona);
	
	if(Tup.Dist5D_P<6){
		Dist_DvsMC_P -> DataEff -> afterR -> Fill(Kbin,zona);     
        	if(Likcut) Lik_DvsMC_P  -> DataEff -> afterR -> Fill(Kbin,zona);
	}

	//Beta bins
	//ToF
	Kbin=ToFDB.GetRBin(RUsed);	
	Dist_DvsMC_P -> DataEff -> beforeTOF -> Fill(Kbin,zona);
	if(Distcut) Lik_DvsMC_P  -> DataEff -> beforeTOF -> Fill(Kbin,zona);

	if(Distcut){
		Dist_DvsMC_P -> DataEff -> afterTOF -> Fill(Kbin,zona);
		if(Likcut) Lik_DvsMC_P  -> DataEff -> afterTOF -> Fill(Kbin,zona);
	}
	//NaF
	if(((int)Tup.Cutmask)>>11==512) {	
		Kbin=NaFDB.GetRBin(RUsed);
		Dist_DvsMC_P -> DataEff -> beforeNaF -> Fill(Kbin,zona);
		if(Distcut) Lik_DvsMC_P  -> DataEff -> beforeNaF -> Fill(Kbin,zona);

		if(Distcut){
			Dist_DvsMC_P -> DataEff -> afterNaF -> Fill(Kbin,zona);
			if(Likcut) Lik_DvsMC_P  -> DataEff -> afterNaF -> Fill(Kbin,zona);
		}
	}
	//Agl
	if(((int)Tup.Cutmask)>>11==0) {
		Kbin=AglDB.GetRBin(RUsed);
		Dist_DvsMC_P -> DataEff -> beforeAgl -> Fill(Kbin,zona);
		if(Distcut) Lik_DvsMC_P  -> DataEff -> beforeAgl -> Fill(Kbin,zona);

		if(Distcut){
			Dist_DvsMC_P -> DataEff -> afterAgl -> Fill(Kbin,zona);
			if(Likcut) Lik_DvsMC_P  -> DataEff -> afterAgl -> Fill(Kbin,zona);
		}
	}
	return;

}

void DVSMCQualeff2_Fill(TNtuple *ntupla, int l){

	 ntupla->GetEvent(l);
	//cuts
	if(Tup.Beta<=0||Tup.R<=0||Tup.R<1.2*Tup.Rcutoff||Tup.Beta>protons->Eval(Tup.R)+0.1||Tup.Beta<protons->Eval(Tup.R)-0.1) return;
	if(!Herejcut) return;
	//
	int Kbin;

	//R bins
	Kbin = RB.GetRBin(Tup.R);

	if(Massa_gen<1) {
		//R bins
		Kbin = RB.GetRBin(Tup.R);	
		Dist_DvsMC_P -> MCEff -> beforeR -> Fill(Kbin);
		if(Tup.Dist5D_P<6) Lik_DvsMC_P  -> MCEff -> beforeR -> Fill(Kbin);

		if(Tup.Dist5D_P<6){
			Dist_DvsMC_P -> MCEff -> afterR -> Fill(Kbin);
			if(Likcut) Lik_DvsMC_P  -> MCEff -> afterR -> Fill(Kbin);
		}
		//Beta bins

		//ToF
		Kbin=ToFDB.GetRBin(RUsed);	
		Dist_DvsMC_P -> MCEff -> beforeTOF -> Fill(Kbin);
		if(Distcut) Lik_DvsMC_P  -> MCEff -> beforeTOF -> Fill(Kbin);

		if(Distcut){
			Dist_DvsMC_P -> MCEff -> afterTOF -> Fill(Kbin);
			if(Likcut) Lik_DvsMC_P  -> MCEff -> afterTOF -> Fill(Kbin);
		}
		//NaF
		if(((int)Tup.Cutmask)>>11==512) {	
			Kbin=NaFDB.GetRBin(RUsed);	
			Dist_DvsMC_P -> MCEff -> beforeNaF -> Fill(Kbin);
			if(Distcut) Lik_DvsMC_P  -> MCEff -> beforeNaF -> Fill(Kbin);

			if(Distcut){
				Dist_DvsMC_P -> MCEff -> afterNaF -> Fill(Kbin);
				if(Likcut) Lik_DvsMC_P  -> MCEff -> afterNaF -> Fill(Kbin);
			}

		}
		//Agl
		if(((int)Tup.Cutmask)>>11==0) {	
			Kbin=AglDB.GetRBin(RUsed);
			Dist_DvsMC_P -> MCEff -> beforeAgl -> Fill(Kbin);
			if(Distcut) Lik_DvsMC_P  -> MCEff -> beforeAgl -> Fill(Kbin);

			if(Distcut){
				Dist_DvsMC_P -> MCEff -> afterAgl -> Fill(Kbin);
				if(Likcut) Lik_DvsMC_P  -> MCEff -> afterAgl -> Fill(Kbin);
			}

		}

	}                        
}

void DVSMCQualeff2_Write(){

	Dist_DvsMC_P -> Write();
	Lik_DvsMC_P  -> Write();

	return;
}


void DVSMCQualeff2(){

	string nomefile="../Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
	TFile * file1 =TFile::Open(nomefile.c_str(),"READ");

	DatavsMC * Dist_DvsMC_P = new DatavsMC(file1,"Dist_DvsMC_P");
	DatavsMC * Lik_DvsMC_P  = new DatavsMC(file1,"Lik_DvsMC_P" );

	LATcorr * LATLikelihoodDATA_TOF = new LATcorr(file1,"LATLikDATA_TOF"   	 ,"Results");
	LATcorr * LATDistanceDATA_TOF   = new LATcorr(file1,"LATDistDATA_TOF" 	 ,"Results");

	LATcorr * LATLikelihoodDATA_NaF = new LATcorr(file1,"LATLikDATA_NaF"  	 ,"Results");
	LATcorr * LATDistanceDATA_NaF   = new LATcorr(file1,"LATDistDATA_NaF" 	 ,"Results");

	LATcorr * LATLikelihoodDATA_Agl = new LATcorr(file1,"LATLikDATA_Agl"  	 ,"Results");
	LATcorr * LATDistanceDATA_Agl   = new LATcorr(file1,"LATDistDATA_Agl" 	 ,"Results");




	cout<<"******* Data vs MC: QUALITY SEL ********"<<endl;

	Dist_DvsMC_P -> Assign_LatCorr( LATDistanceDATA_TOF   ->  LATcorrR_fit , 
			LATDistanceDATA_TOF   ->  LATcorrR_fit ,
			LATDistanceDATA_NaF   ->  LATcorrR_fit ,
			LATDistanceDATA_Agl   ->  LATcorrR_fit );

	Lik_DvsMC_P  ->	Assign_LatCorr( LATLikelihoodDATA_TOF ->  LATcorrR_fit , 	
			LATLikelihoodDATA_TOF ->  LATcorrR_fit ,
			LATLikelihoodDATA_NaF ->  LATcorrR_fit ,
			LATLikelihoodDATA_Agl ->  LATcorrR_fit );



	Dist_DvsMC_P ->Eval_DandMC_Eff();  
	Lik_DvsMC_P  ->Eval_DandMC_Eff();

	Dist_DvsMC_P ->Eval_Corrections();
	Lik_DvsMC_P  ->Eval_Corrections();


	TH1F* DistP_Correction_R   =(TH1F*) Dist_DvsMC_P -> GetCorrection_R()  ;
	TH1F* DistP_Correction_TOF =(TH1F*) Dist_DvsMC_P -> GetCorrection_TOF();
	TH1F* DistP_Correction_NaF =(TH1F*) Dist_DvsMC_P -> GetCorrection_NaF();
	TH1F* DistP_Correction_Agl =(TH1F*) Dist_DvsMC_P -> GetCorrection_Agl();

	TH1F* LikP_Correction_R    =(TH1F*) Lik_DvsMC_P -> GetCorrection_R()  ;
	TH1F* LikP_Correction_TOF  =(TH1F*) Lik_DvsMC_P -> GetCorrection_TOF();
	TH1F* LikP_Correction_NaF  =(TH1F*) Lik_DvsMC_P -> GetCorrection_NaF();
	TH1F* LikP_Correction_Agl  =(TH1F*) Lik_DvsMC_P -> GetCorrection_Agl();


	cout<<"*** Updating P1 file ****"<<endl;
	nomefile="../Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
	file1 =TFile::Open(nomefile.c_str(),"UPDATE");

	file1->cd("Results");

	DistP_Correction_R    -> Write("Dist_DvsMC_P_CorrectionR"  );
	DistP_Correction_TOF  -> Write("Dist_DvsMC_P_CorrectionTOF");
	DistP_Correction_NaF  -> Write("Dist_DvsMC_P_CorrectionNaF");
	DistP_Correction_Agl  -> Write("Dist_DvsMC_P_CorrectionAgl");

	LikP_Correction_R    -> Write("Lik_DvsMC_P_CorrectionR"  );
	LikP_Correction_TOF  -> Write("Lik_DvsMC_P_CorrectionTOF");
	LikP_Correction_NaF  -> Write("Lik_DvsMC_P_CorrectionNaF");
	LikP_Correction_Agl  -> Write("Lik_DvsMC_P_CorrectionAgl");


	file1->Write();
	file1->Close();






	TCanvas *c20=new TCanvas("Data vs MC: Likelihood (R bins)");
	TCanvas *c21=new TCanvas("Data vs MC: Distance (R bins)");

	TCanvas *c20_bis=new TCanvas("Data vs MC: Likelihood (Beta bins)");
	TCanvas *c21_bis=new TCanvas("Data vs MC: Distance (Beta bins)");

	TGraphErrors *LikDVSMC_P_Graph;
	TGraphErrors *DistDVSMC_P_Graph;

	c20->cd();
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	LikDVSMC_P_Graph=new TGraphErrors();
	LikDVSMC_P_Graph->SetName("LikDVSMC_P_Graph");
	int j=0;
	for(int i=1;i<nbinsr;i++) {
		if(LikP_Correction_R -> GetBinContent(i+1)>0){
			LikDVSMC_P_Graph->SetPoint(j,R_cent[i],LikP_Correction_R -> GetBinContent(i+1));
			LikDVSMC_P_Graph->SetPointError(j,0,LikP_Correction_R -> GetBinError(i+1));
			j++;
		}
	}
	LikDVSMC_P_Graph->SetLineColor(2);
	LikDVSMC_P_Graph->SetFillColor(2);
	LikDVSMC_P_Graph->SetFillStyle(3001);
	LikDVSMC_P_Graph->SetLineWidth(4);
	LikDVSMC_P_Graph->Draw("AP4C");

	c21->cd();
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	DistDVSMC_P_Graph=new TGraphErrors();
	DistDVSMC_P_Graph->SetName("DistDVSMC_P_Graph");
	j=0;
	for(int i=1;i<nbinsr;i++) {
		if(DistP_Correction_R -> GetBinContent(i+1)>0){
			DistDVSMC_P_Graph->SetPoint(j,R_cent[i],DistP_Correction_R -> GetBinContent(i+1));
			DistDVSMC_P_Graph->SetPointError(j,0,DistP_Correction_R -> GetBinError(i+1));
			j++;
		}
	}
	DistDVSMC_P_Graph->SetLineColor(2);
	DistDVSMC_P_Graph->SetFillColor(2);
	DistDVSMC_P_Graph->SetFillStyle(3001);
	DistDVSMC_P_Graph->SetLineWidth(4);
	DistDVSMC_P_Graph->Draw("AP4C");

	c20_bis->Divide(3,1);

	c20_bis->cd(1);
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraphErrors * LikDVSMC_P_GraphTOF=new TGraphErrors();
	LikDVSMC_P_GraphTOF->SetName("LikDVSMC_P_GraphTOF");
	j=0;
	for(int i=1;i<nbinsToF;i++) {
		if(LikP_Correction_TOF -> GetBinContent(i+1)>0){
			LikDVSMC_P_GraphTOF->SetPoint(j,Ekincent[i],LikP_Correction_TOF -> GetBinContent(i+1));
			LikDVSMC_P_GraphTOF->SetPointError(j,0,LikP_Correction_TOF -> GetBinError(i+1));
			j++;
		}
	}
	LikDVSMC_P_GraphTOF->SetLineColor(2);
	LikDVSMC_P_GraphTOF->SetFillColor(2);
	LikDVSMC_P_GraphTOF->SetFillStyle(3001);
	LikDVSMC_P_GraphTOF->SetLineWidth(4);
	LikDVSMC_P_GraphTOF->Draw("AP4C");

	c20_bis->cd(2);
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraphErrors * LikDVSMC_P_GraphNaF=new TGraphErrors();
	LikDVSMC_P_GraphNaF->SetName("LikDVSMC_P_GraphNaF");
	j=0;
	for(int i=1;i<nbinsNaF;i++) {
		if(LikP_Correction_NaF -> GetBinContent(i+1)>0){
			LikDVSMC_P_GraphNaF->SetPoint(j,EkincentNaF[i],LikP_Correction_NaF -> GetBinContent(i+1));
			LikDVSMC_P_GraphNaF->SetPointError(j,0,LikP_Correction_NaF -> GetBinError(i+1));
			j++;
		}
	}
	LikDVSMC_P_GraphNaF->SetLineColor(2);
	LikDVSMC_P_GraphNaF->SetFillColor(2);
	LikDVSMC_P_GraphNaF->SetFillStyle(3001);
	LikDVSMC_P_GraphNaF->SetLineWidth(4);
	LikDVSMC_P_GraphNaF->Draw("AP4C");

	c20_bis->cd(3);
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraphErrors * LikDVSMC_P_GraphAgl=new TGraphErrors();
	LikDVSMC_P_GraphAgl->SetName("LikDVSMC_P_GraphAgl");
	j=0;
	for(int i=1;i<nbinsAgl;i++) {
		if(LikP_Correction_Agl -> GetBinContent(i+1)>0){
			LikDVSMC_P_GraphAgl->SetPoint(j,EkincentAgl[i],LikP_Correction_Agl -> GetBinContent(i+1));
			LikDVSMC_P_GraphAgl->SetPointError(j,0,LikP_Correction_Agl -> GetBinError(i+1));
			j++;
		}
	}
	LikDVSMC_P_GraphAgl->SetLineColor(2);
	LikDVSMC_P_GraphAgl->SetFillColor(2);
	LikDVSMC_P_GraphAgl->SetFillStyle(3001);
	LikDVSMC_P_GraphAgl->SetLineWidth(4);
	LikDVSMC_P_GraphAgl->Draw("AP4C");

	c21_bis->Divide(3,1);

	c21_bis->cd(1);
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraphErrors * DistDVSMC_P_GraphTOF=new TGraphErrors();
	DistDVSMC_P_GraphTOF->SetName("DistDVSMC_P_GraphTOF");
	j=0;
	for(int i=1;i<nbinsToF;i++) {
		if(DistP_Correction_TOF -> GetBinContent(i+1)>0){
			DistDVSMC_P_GraphTOF->SetPoint(j,Ekincent[i],DistP_Correction_TOF -> GetBinContent(i+1));
			DistDVSMC_P_GraphTOF->SetPointError(j,0,DistP_Correction_TOF -> GetBinError(i+1));
			j++;
		}
	}
	DistDVSMC_P_GraphTOF->SetLineColor(2);
	DistDVSMC_P_GraphTOF->SetFillColor(2);
	DistDVSMC_P_GraphTOF->SetFillStyle(3001);
	DistDVSMC_P_GraphTOF->SetLineWidth(4);
	DistDVSMC_P_GraphTOF->Draw("AP4C");

	c21_bis->cd(2);
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraphErrors * DistDVSMC_P_GraphNaF=new TGraphErrors();
	DistDVSMC_P_GraphNaF->SetName("DistDVSMC_P_GraphNaF");
	j=0;
	for(int i=1;i<nbinsNaF;i++) {
		if(DistP_Correction_NaF -> GetBinContent(i+1)>0){
			DistDVSMC_P_GraphNaF->SetPoint(j,EkincentNaF[i],DistP_Correction_NaF -> GetBinContent(i+1));
			DistDVSMC_P_GraphNaF->SetPointError(j,0,DistP_Correction_NaF -> GetBinError(i+1));
			j++;
		}
	}
	DistDVSMC_P_GraphNaF->SetLineColor(2);
	DistDVSMC_P_GraphNaF->SetFillColor(2);
	DistDVSMC_P_GraphNaF->SetFillStyle(3001);
	DistDVSMC_P_GraphNaF->SetLineWidth(4);
	DistDVSMC_P_GraphNaF->Draw("AP4C");

	c21_bis->cd(3);
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraphErrors * DistDVSMC_P_GraphAgl=new TGraphErrors();
	DistDVSMC_P_GraphAgl->SetName("DistDVSMC_P_GraphAgl");
	j=0;
	for(int i=1;i<nbinsAgl;i++) {
		if(DistP_Correction_Agl -> GetBinContent(i+1)>0){
			DistDVSMC_P_GraphAgl->SetPoint(j,EkincentAgl[i],DistP_Correction_Agl -> GetBinContent(i+1));
			DistDVSMC_P_GraphAgl->SetPointError(j,0,DistP_Correction_Agl -> GetBinError(i+1));
			j++;
		}
	}
	DistDVSMC_P_GraphAgl->SetLineColor(2);
	DistDVSMC_P_GraphAgl->SetFillColor(2);
	DistDVSMC_P_GraphAgl->SetFillStyle(3001);
	DistDVSMC_P_GraphAgl->SetLineWidth(4);
	DistDVSMC_P_GraphAgl->Draw("AP4C");





	cout<<"*** Updating Results file ***"<<endl;
	nomefile="./Final_plots/"+mese+".root";
	TFile *f_out=new TFile(nomefile.c_str(), "UPDATE");
	f_out->mkdir("DATA-driven Results/Data vs MC/Protons");
	f_out->cd("DATA-driven Results/Data vs MC/Protons");
	c20->Write();
	c21->Write();
	c20_bis->Write();
	c21_bis->Write();

	f_out->Write();
	f_out->Close();






	return;
}

