using namespace std;

DatavsMC * PreSel_DvsMC_D = new DatavsMC("PreSel_DvsMC_D",11,3,6);

void DVSMCPreSeleffD_D_Fill(TNtuple *ntupla, int l,int zona){

	 ntupla->GetEvent(l);
	//cuts
	if(Tup.Beta_pre<=0||Tup.R_pre<=0||Tup.R_pre<1.2*Tup.Rcutoff||Tup.Beta_pre>protons->Eval(Tup.R_pre)+0.1||Tup.Beta_pre<protons->Eval(Tup.R_pre)-0.1) return;
	if(!((Tup.R_pre>Rcut[zona]&&zona<10)||(zona==10)))  return;
	if(!Herejcut) return;
	if(!(Tup.EdepTOFU<EdepTOFbeta->Eval(Tup.Beta_pre)+1&&Tup.EdepTOFU>EdepTOFbeta->Eval(Tup.Beta_pre)-1)) return;
	//
	int Kbin;
	for(int S=0;S<3;S++){
		//Beta bins
		//ToF
		Kbin=ToFDB.GetRBin(RUsed);	
		if(cmask.notPassed(S)) ((TH3*)PreSel_DvsMC_D -> DataEff -> beforeTOF) -> Fill(Kbin,zona,S);
		if(cmask.passed(S))       ((TH3*)PreSel_DvsMC_D -> DataEff -> afterTOF ) -> Fill(Kbin,zona,S);
		//NaF
		if(cmask.isFromNaF()) {	
			Kbin=NaFDB.GetRBin(RUsed);
			if(cmask.notPassed(S)) ((TH3*)PreSel_DvsMC_D -> DataEff -> beforeNaF) -> Fill(Kbin,zona,S);
			if(cmask.passed(S))       ((TH3*)PreSel_DvsMC_D -> DataEff -> afterNaF ) -> Fill(Kbin,zona,S);
		}
		//Agl
		if(cmask.isFromAgl()) {
			Kbin=AglDB.GetRBin(RUsed);
			if(cmask.notPassed(S)) ((TH3*)PreSel_DvsMC_D -> DataEff -> beforeAgl) -> Fill(Kbin,zona,S);
			if(cmask.passed(S))       ((TH3*)PreSel_DvsMC_D -> DataEff -> afterAgl ) -> Fill(Kbin,zona,S);
		}
	}
	return;

}
void DVSMCPreSeleffD_Fill(TNtuple *ntupla, int l){
	 ntupla->GetEvent(l);
	//cuts
	if(Tup.Beta_pre<=0||Tup.R_pre<=0||Tup.Beta_pre>protons->Eval(Tup.R_pre)+0.1||Tup.Beta_pre<protons->Eval(Tup.R_pre)-0.1) return;
	if(!Herejcut) return;
	if(!(Tup.EdepTOFU<EdepTOFbeta->Eval(Tup.Beta_pre)+1&&Tup.EdepTOFU>EdepTOFbeta->Eval(Tup.Beta_pre)-1)) return;
	//
	int Kbin;
	for(int S=0;S<3;S++){
		if(Massa_gen>1&&Massa_gen<2) {
			//Beta bins
			//ToF
			Kbin=ToFDB.GetRBin(RUsed);	
			if(cmask.notPassed(S)) ((TH3*)PreSel_DvsMC_D -> MCEff -> beforeTOF) -> Fill(Kbin,ReturnMCGenType(),S);
			if(cmask.passed(S))       ((TH3*)PreSel_DvsMC_D -> MCEff -> afterTOF ) -> Fill(Kbin,ReturnMCGenType(),S);

			//NaF
			if(cmask.isFromNaF()) {	
				Kbin=NaFDB.GetRBin(RUsed);	
				if(cmask.notPassed(S)) ((TH3*)PreSel_DvsMC_D -> MCEff -> beforeNaF) -> Fill(Kbin,ReturnMCGenType(),S);
				if(cmask.passed(S))       ((TH3*)PreSel_DvsMC_D -> MCEff -> afterNaF ) -> Fill(Kbin,ReturnMCGenType(),S);
			}
			//Agl
			if(cmask.isFromAgl()) {	
				Kbin=AglDB.GetRBin(RUsed);
				if(cmask.notPassed(S)) ((TH3*)PreSel_DvsMC_D -> MCEff -> beforeAgl) -> Fill(Kbin,ReturnMCGenType(),S);
				if(cmask.passed(S))       ((TH3*)PreSel_DvsMC_D -> MCEff -> afterAgl ) -> Fill(Kbin,ReturnMCGenType(),S);
			}

		}

	}
	return;
}


void DVSMCPreSeleffD_Write(){

	PreSel_DvsMC_D -> Write();
	return;
}


void DVSMCPreSeleffD(){

   inputHistoFile->ReOpen("READ");

	DatavsMC * PreSel_DvsMC_D = new DatavsMC(inputHistoFile,"PreSel_DvsMC_D",6);

	LATcorr * LATpreSelDATA = new LATcorr(inputHistoFile,"LATpreSelDATA"      ,"Results");


	cout<<"******* Data vs MC:  PRESELECTIONS (D) ********"<<endl;

	PreSel_DvsMC_D -> Assign_LatCorr( LATpreSelDATA   ->  LATcorrR_fit , 
			LATpreSelDATA   ->  LATcorrR_fit ,
			LATpreSelDATA   ->  LATcorrR_fit ,
			LATpreSelDATA   ->  LATcorrR_fit );


	PreSel_DvsMC_D ->Eval_DandMC_Eff();  

	PreSel_DvsMC_D ->Eval_Corrections();


	TH3F* PreSelD_Correction_R   =(TH3F*) PreSel_DvsMC_D -> GetCorrection_R()  ;
	TH3F* PreSelD_Correction_TOF =(TH3F*) PreSel_DvsMC_D -> GetCorrection_TOF();
	TH3F* PreSelD_Correction_NaF =(TH3F*) PreSel_DvsMC_D -> GetCorrection_NaF();
	TH3F* PreSelD_Correction_Agl =(TH3F*) PreSel_DvsMC_D -> GetCorrection_Agl();


	cout<<"*** Updating P1 file ****"<<endl;
   inputHistoFile->ReOpen("UPDATE");

	inputHistoFile->cd("Results");

	PreSelD_Correction_R    -> Write("PreSel_DvsMC_D_CorrectionR"  );
	PreSelD_Correction_TOF  -> Write("PreSel_DvsMC_D_CorrectionTOF");
	PreSelD_Correction_NaF  -> Write("PreSel_DvsMC_D_CorrectionNaF");
	PreSelD_Correction_Agl  -> Write("PreSel_DvsMC_D_CorrectionAgl");

	inputHistoFile->Write();
	inputHistoFile->Close();

	string tagli[3]={"Matching TOF","Chi^2 R","1 Tr. Track"};
	string nome;

	TCanvas *c21_D[3];


	int j=0;
	TGraphErrors * PreSelD_Correction_TOF_Graph[3][6];
	TGraphErrors * PreSelD_Correction_NaF_Graph[3][6];
	TGraphErrors * PreSelD_Correction_Agl_Graph[3][6];

	string MCLegend[6]= {"d.pl1.0_520_GG_Blic","d.pl1.0_520_GG_BlicDPMJet","d.pl1.0_520_GG_QMD","d.pl1.0_520_Shen_Blic","d.pl1.0_520_Shen_BlicDPMJet","d.pl1.0_520_Shen_QMD"};

	for(int S=0;S<3;S++){

		c21_D[S] = new TCanvas(("Deutons Data vs MC: "+tagli[S] +"(Beta Bins)").c_str());
		c21_D[S] -> Divide(3,1);

		c21_D[S]->cd(1);
		gPad->SetLogx();
		gPad->SetGridx();
		gPad->SetGridy();

		for(int mc_type=0;mc_type<6;mc_type++){
			PreSelD_Correction_TOF_Graph[S][mc_type] = new TGraphErrors();
			PreSelD_Correction_NaF_Graph[S][mc_type] = new TGraphErrors();
			PreSelD_Correction_Agl_Graph[S][mc_type] = new TGraphErrors();

			j=0;
			for(int i=1;i<nbinsToF;i++) {
				if(PreSelD_Correction_TOF -> GetBinContent(i+1,mc_type+1,S+1)>0){
					PreSelD_Correction_TOF_Graph[S][mc_type]->SetPoint(j,ToFPB.EkBinCent(i),PreSelD_Correction_TOF -> GetBinContent(i+1,mc_type+1,S+1));
					PreSelD_Correction_TOF_Graph[S][mc_type]->SetPointError(j,0,PreSelD_Correction_TOF -> GetBinError(i+1,mc_type+1,S+1));
					j++;
				}
			}
			PreSelD_Correction_TOF_Graph[S][mc_type]->SetLineColor(4);
			PreSelD_Correction_TOF_Graph[S][mc_type]->SetFillColor(4);
			PreSelD_Correction_TOF_Graph[S][mc_type]->SetMarkerColor(4);	
			PreSelD_Correction_TOF_Graph[S][mc_type]->SetFillStyle(3001);
			PreSelD_Correction_TOF_Graph[S][mc_type]->SetLineWidth(1);
			PreSelD_Correction_TOF_Graph[S][mc_type]->SetMarkerStyle(mc_type+3);

		}

		PreSelD_Correction_TOF_Graph[S][0]->Draw("AP4C");
		for(int mc_type=1;mc_type<6;mc_type++){
			PreSelD_Correction_TOF_Graph[S][mc_type]->Draw("P4Csame");	
		}


		c21_D[S]->cd(2);
		gPad->SetLogx();
		gPad->SetGridx();
		gPad->SetGridy();

		for(int mc_type=0;mc_type<6;mc_type++){
			PreSelD_Correction_NaF_Graph[S][mc_type] = new TGraphErrors();
			PreSelD_Correction_NaF_Graph[S][mc_type] = new TGraphErrors();
			PreSelD_Correction_Agl_Graph[S][mc_type] = new TGraphErrors();

			j=0;
			for(int i=1;i<nbinsToF;i++) {
				if(PreSelD_Correction_NaF -> GetBinContent(i+1,mc_type+1,S+1)>0){
					PreSelD_Correction_NaF_Graph[S][mc_type]->SetPoint(j,NaFPB.EkBinCent(i),PreSelD_Correction_NaF -> GetBinContent(i+1,mc_type+1,S+1));
					PreSelD_Correction_NaF_Graph[S][mc_type]->SetPointError(j,0,PreSelD_Correction_NaF -> GetBinError(i+1,mc_type+1,S+1));
					j++;
				}
			}
			PreSelD_Correction_NaF_Graph[S][mc_type]->SetLineColor(4);
			PreSelD_Correction_NaF_Graph[S][mc_type]->SetFillColor(4);
			PreSelD_Correction_NaF_Graph[S][mc_type]->SetMarkerColor(4);
			PreSelD_Correction_NaF_Graph[S][mc_type]->SetFillStyle(3001);
			PreSelD_Correction_NaF_Graph[S][mc_type]->SetLineWidth(1);
			PreSelD_Correction_NaF_Graph[S][mc_type]->SetMarkerStyle(mc_type+3);

		}

		{
			TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
			leg->AddEntry(PreSelD_Correction_NaF_Graph[S][0],MCLegend[0].c_str(), "ep");
			PreSelD_Correction_NaF_Graph[S][0]->Draw("AP4C");
			for(int mc_type=1;mc_type<6;mc_type++) { PreSelD_Correction_NaF_Graph[S][mc_type]->Draw("P4Csame");
				leg->AddEntry(PreSelD_Correction_NaF_Graph[S][mc_type],MCLegend[mc_type].c_str(), "ep");
			}

			leg->Draw("same");
		}


		c21_D[S]->cd(3);
		gPad->SetLogx();
		gPad->SetGridx();
		gPad->SetGridy();

		for(int mc_type=0;mc_type<6;mc_type++){
			PreSelD_Correction_Agl_Graph[S][mc_type] = new TGraphErrors();
			PreSelD_Correction_Agl_Graph[S][mc_type] = new TGraphErrors();
			PreSelD_Correction_Agl_Graph[S][mc_type] = new TGraphErrors();

			j=0;
			for(int i=1;i<nbinsToF;i++) {
				if(PreSelD_Correction_Agl -> GetBinContent(i+1,mc_type+1,S+1)>0){
					PreSelD_Correction_Agl_Graph[S][mc_type]->SetPoint(j,AglPB.EkBinCent(i),PreSelD_Correction_Agl -> GetBinContent(i+1,mc_type+1,S+1));
					PreSelD_Correction_Agl_Graph[S][mc_type]->SetPointError(j,0,PreSelD_Correction_Agl -> GetBinError(i+1,mc_type+1,S+1));
					j++;
				}
			}
			PreSelD_Correction_Agl_Graph[S][mc_type]->SetLineColor(4);
			PreSelD_Correction_Agl_Graph[S][mc_type]->SetFillColor(4);
			PreSelD_Correction_Agl_Graph[S][mc_type]->SetMarkerColor(4);
			PreSelD_Correction_Agl_Graph[S][mc_type]->SetFillStyle(3001);
			PreSelD_Correction_Agl_Graph[S][mc_type]->SetLineWidth(1);
			PreSelD_Correction_Agl_Graph[S][mc_type]->SetMarkerStyle(mc_type+3);

		}

		PreSelD_Correction_Agl_Graph[S][0]->Draw("AP4C");
		for(int mc_type=1;mc_type<6;mc_type++){
			PreSelD_Correction_Agl_Graph[S][mc_type]->Draw("P4Csame");
		}

	}





	cout<<"*** Updating Results file ***"<<endl;
	string filename="./Final_plots/"+mese+".root";
	TFile *f_out=new TFile(filename.c_str(), "UPDATE");
	f_out->mkdir("DATA-driven Results/Data vs MC/Deutons");
	f_out->cd("DATA-driven Results/Data vs MC/Deutons");

	for(int S=0;S<3;S++){
		c21_D[S]->Write();
	}
	f_out->Write();
	f_out->Close();

	return;
}

