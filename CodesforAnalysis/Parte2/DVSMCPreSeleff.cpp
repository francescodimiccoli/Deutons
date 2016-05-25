using namespace std;

DatavsMC * PreSel_DvsMC_P = new DatavsMC("PreSel_DvsMC_P",11,3);

void DVSMCPreSeleff_D_Fill(TNtuple *ntupla, int l,int zona){

	 ntupla->GetEvent(l);
	//cuts
	if(Tup.Beta_pre<=0||Tup.R_pre<=0||Tup.R_pre<1.2*Tup.Rcutoff||Tup.Beta_pre>protons->Eval(Tup.R_pre)+0.1||Tup.Beta_pre<protons->Eval(Tup.R_pre)-0.1) return;
	if(!((Tup.R_pre>Rcut[zona]&&zona<10)||(zona==10)))  return;
	if(!Herejcut) return;
	if(!(Tup.EdepTOFU<EdepTOFbeta->Eval(Tup.Beta_pre)+1&&Tup.EdepTOFU>EdepTOFbeta->Eval(Tup.Beta_pre)-1)) return;
	//
	int Kbin;
	for(int S=0;S<3;S++){
		//R bins
		Kbin = RB.GetRBin(RUsed);
		if(((int)Tup.Cutmask&notpassed[S])==notpassed[S]) ((TH3*)PreSel_DvsMC_P -> DataEff -> beforeR) -> Fill(Kbin,zona,S);	
		if(((int)Tup.Cutmask&passed[S])==passed[S])	      ((TH3*)PreSel_DvsMC_P -> DataEff -> afterR ) -> Fill(Kbin,zona,S);     

		//Beta bins
		//ToF
		Kbin=ToFDB.GetRBin(RUsed);	
		if(((int)Tup.Cutmask&notpassed[S])==notpassed[S]) ((TH3*)PreSel_DvsMC_P -> DataEff -> beforeTOF) -> Fill(Kbin,zona,S);
                if(((int)Tup.Cutmask&passed[S])==passed[S])       ((TH3*)PreSel_DvsMC_P -> DataEff -> afterTOF ) -> Fill(Kbin,zona,S);
		//NaF
		if(((int)Tup.Cutmask)>>11==512) {	
			Kbin=NaFDB.GetRBin(RUsed);
			if(((int)Tup.Cutmask&notpassed[S])==notpassed[S]) ((TH3*)PreSel_DvsMC_P -> DataEff -> beforeNaF) -> Fill(Kbin,zona,S);
                	if(((int)Tup.Cutmask&passed[S])==passed[S])       ((TH3*)PreSel_DvsMC_P -> DataEff -> afterNaF ) -> Fill(Kbin,zona,S);
		}
		//Agl
		if(((int)Tup.Cutmask)>>11==0) {
			Kbin=AglDB.GetRBin(RUsed);
			if(((int)Tup.Cutmask&notpassed[S])==notpassed[S]) ((TH3*)PreSel_DvsMC_P -> DataEff -> beforeAgl) -> Fill(Kbin,zona,S);
                	if(((int)Tup.Cutmask&passed[S])==passed[S])       ((TH3*)PreSel_DvsMC_P -> DataEff -> afterAgl ) -> Fill(Kbin,zona,S);
		}
	}
	return;
	
}
void DVSMCPreSeleff_Fill(TNtuple *ntupla, int l){
	 ntupla->GetEvent(l);
	//cuts
	if(Tup.Beta_pre<=0||Tup.R_pre<=0||Tup.Beta_pre>protons->Eval(Tup.R_pre)+0.1||Tup.Beta_pre<protons->Eval(Tup.R_pre)-0.1) return;
	if(!Herejcut) return;
	if(!(Tup.EdepTOFU<EdepTOFbeta->Eval(Tup.Beta_pre)+1&&Tup.EdepTOFU>EdepTOFbeta->Eval(Tup.Beta_pre)-1)) return;
	//
	int Kbin;
	for(int S=0;S<3;S++){
		if(Massa_gen<1) {
			//R bins
			Kbin = RB.GetRBin(RUsed);	
			if(((int)Tup.Cutmask&notpassed[S])==notpassed[S]) ((TH2*)PreSel_DvsMC_P -> MCEff -> beforeR) -> Fill(Kbin,S);
			if(((int)Tup.Cutmask&passed[S])==passed[S])       ((TH2*)PreSel_DvsMC_P -> MCEff -> afterR ) -> Fill(Kbin,S);
			//Beta bins

			//ToF
			Kbin=ToFDB.GetRBin(RUsed);	
			if(((int)Tup.Cutmask&notpassed[S])==notpassed[S]) ((TH2*)PreSel_DvsMC_P -> MCEff -> beforeTOF) -> Fill(Kbin,S);
			if(((int)Tup.Cutmask&passed[S])==passed[S])       ((TH2*)PreSel_DvsMC_P -> MCEff -> afterTOF ) -> Fill(Kbin,S);

			//NaF
			if(((int)Tup.Cutmask)>>11==512) {	
				Kbin=NaFDB.GetRBin(RUsed);	
				if(((int)Tup.Cutmask&notpassed[S])==notpassed[S]) ((TH2*)PreSel_DvsMC_P -> MCEff -> beforeNaF) -> Fill(Kbin,S);
				if(((int)Tup.Cutmask&passed[S])==passed[S])       ((TH2*)PreSel_DvsMC_P -> MCEff -> afterNaF ) -> Fill(Kbin,S);
			}
			//Agl
			if(((int)Tup.Cutmask)>>11==0) {	
				Kbin=AglDB.GetRBin(RUsed);
				if(((int)Tup.Cutmask&notpassed[S])==notpassed[S]) ((TH2*)PreSel_DvsMC_P -> MCEff -> beforeAgl) -> Fill(Kbin,S);
				if(((int)Tup.Cutmask&passed[S])==passed[S])       ((TH2*)PreSel_DvsMC_P -> MCEff -> afterAgl ) -> Fill(Kbin,S);
			}

		}

	}
	return;
}


void DVSMCPreSeleff_Write(){

	PreSel_DvsMC_P -> Write();
	return;
}


void DVSMCPreSeleff(){

	string nomefile="../Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
	TFile * file1 =TFile::Open(nomefile.c_str(),"READ");

	DatavsMC * PreSel_DvsMC_P = new DatavsMC(file1,"PreSel_DvsMC_P");

	LATcorr * LATpreSelDATA = new LATcorr(file1,"LATpreSelDATA"      ,"Results");


	cout<<"******* Data vs MC:  PRESELECTIONS ********"<<endl;

	PreSel_DvsMC_P -> Assign_LatCorr( LATpreSelDATA   ->  LATcorrR_fit , 
					  LATpreSelDATA   ->  LATcorrR_fit ,
					  LATpreSelDATA   ->  LATcorrR_fit ,
					  LATpreSelDATA   ->  LATcorrR_fit );


	PreSel_DvsMC_P ->Eval_DandMC_Eff();  

	PreSel_DvsMC_P ->Eval_Corrections();


	TH2F* PreSel_Correction_R   =(TH2F*) PreSel_DvsMC_P -> GetCorrection_R()  ;
	TH2F* PreSel_Correction_TOF =(TH2F*) PreSel_DvsMC_P -> GetCorrection_TOF();
	TH2F* PreSel_Correction_NaF =(TH2F*) PreSel_DvsMC_P -> GetCorrection_NaF();
	TH2F* PreSel_Correction_Agl =(TH2F*) PreSel_DvsMC_P -> GetCorrection_Agl();


	cout<<"*** Updating P1 file ****"<<endl;
	nomefile="../Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
	file1 =TFile::Open(nomefile.c_str(),"UPDATE");

	file1->cd("Results");

	PreSel_Correction_R    -> Write("PreSel_DvsMC_P_CorrectionR"  );
	PreSel_Correction_TOF  -> Write("PreSel_DvsMC_P_CorrectionTOF");
	PreSel_Correction_NaF  -> Write("PreSel_DvsMC_P_CorrectionNaF");
	PreSel_Correction_Agl  -> Write("PreSel_DvsMC_P_CorrectionAgl");

	file1->Write();
	file1->Close();

	string tagli[3]={"Matching TOF","Chi^2 R","1 Tr. Track"};
        string nome;

	TCanvas *c21[3];
        TCanvas *c20[3];
	
	TGraphErrors * PreSel_Correction_R_Graph[3];

	TGraphErrors * PreSel_Correction_TOF_Graph[3];
        TGraphErrors * PreSel_Correction_NaF_Graph[3];
        TGraphErrors * PreSel_Correction_Agl_Graph[3];
	
	for(int S=0;S<3;S++){
		c20[S] = new TCanvas(("Data vs MC: "+tagli[S] +"(R Bins)").c_str());
		c20[S]->cd();
		gPad->SetLogx();
		gPad->SetGridx();
		gPad->SetGridy();
		PreSel_Correction_R_Graph[S] = new TGraphErrors();        
		int j=0;
		for(int i=1;i<nbinsr;i++) {
			if(PreSel_Correction_R -> GetBinContent(i+1,S+1)>0){
				PreSel_Correction_R_Graph[S]->SetPoint(j,R_cent[i],PreSel_Correction_R -> GetBinContent(i+1,S+1));
				PreSel_Correction_R_Graph[S]->SetPointError(j,0,PreSel_Correction_R -> GetBinError(i+1,S+1));
				j++;
			}
		}
		PreSel_Correction_R_Graph[S]->SetLineColor(2);
		PreSel_Correction_R_Graph[S]->SetFillColor(2);
		PreSel_Correction_R_Graph[S]->SetFillStyle(3001);
		PreSel_Correction_R_Graph[S]->SetLineWidth(4);
		PreSel_Correction_R_Graph[S]->Draw("AP4C");
                
		c21[S] = new TCanvas(("Data vs MC: "+tagli[S] +"(Beta Bins)").c_str());
		c21[S] -> Divide(3,1);
                
		c21[S]->cd(1);
                gPad->SetLogx();
                gPad->SetGridx();
                gPad->SetGridy();
	
		PreSel_Correction_TOF_Graph[S] = new TGraphErrors();
               	PreSel_Correction_NaF_Graph[S] = new TGraphErrors();
		PreSel_Correction_Agl_Graph[S] = new TGraphErrors();

		j=0;
                for(int i=1;i<nbinsToF;i++) {
                        if(PreSel_Correction_TOF -> GetBinContent(i+1,S+1)>0){
                                PreSel_Correction_TOF_Graph[S]->SetPoint(j,ToFPB.EkBinCent(i),PreSel_Correction_TOF -> GetBinContent(i+1,S+1));
                                PreSel_Correction_TOF_Graph[S]->SetPointError(j,0,PreSel_Correction_TOF -> GetBinError(i+1,S+1));
                                j++;
                        }
                }
                PreSel_Correction_TOF_Graph[S]->SetLineColor(2);
                PreSel_Correction_TOF_Graph[S]->SetFillColor(2);
                PreSel_Correction_TOF_Graph[S]->SetFillStyle(3001);
                PreSel_Correction_TOF_Graph[S]->SetLineWidth(4);
                PreSel_Correction_TOF_Graph[S]->Draw("AP4C");

		c21[S]->cd(2);
                gPad->SetLogx();
                gPad->SetGridx();
                gPad->SetGridy();
                PreSel_Correction_NaF_Graph[S] = new TGraphErrors();
                j=0;
                for(int i=1;i<nbinsNaF;i++) {
                        if(PreSel_Correction_NaF -> GetBinContent(i+1,S+1)>0){
                                PreSel_Correction_NaF_Graph[S]->SetPoint(j,ToFPB.EkBinCent(i),PreSel_Correction_NaF -> GetBinContent(i+1,S+1));
                                PreSel_Correction_NaF_Graph[S]->SetPointError(j,0,PreSel_Correction_NaF -> GetBinError(i+1,S+1));
                                j++;
                        }
                }
                PreSel_Correction_NaF_Graph[S]->SetLineColor(2);
                PreSel_Correction_NaF_Graph[S]->SetFillColor(2);
                PreSel_Correction_NaF_Graph[S]->SetFillStyle(3001);
                PreSel_Correction_NaF_Graph[S]->SetLineWidth(4);
                PreSel_Correction_NaF_Graph[S]->Draw("AP4C");
		
		c21[S]->cd(3);
                gPad->SetLogx();
                gPad->SetGridx();
                gPad->SetGridy();
                PreSel_Correction_Agl_Graph[S] = new TGraphErrors();
                j=0;
                for(int i=1;i<nbinsToF;i++) {
                        if(PreSel_Correction_Agl -> GetBinContent(i+1,S+1)>0){
                                PreSel_Correction_Agl_Graph[S]->SetPoint(j,ToFPB.EkBinCent(i),PreSel_Correction_Agl -> GetBinContent(i+1,S+1));
                                PreSel_Correction_Agl_Graph[S]->SetPointError(j,0,PreSel_Correction_Agl -> GetBinError(i+1,S+1));
                                j++;
                        }
                }
                PreSel_Correction_Agl_Graph[S]->SetLineColor(2);
                PreSel_Correction_Agl_Graph[S]->SetFillColor(2);
                PreSel_Correction_Agl_Graph[S]->SetFillStyle(3001);
                PreSel_Correction_Agl_Graph[S]->SetLineWidth(4);
                PreSel_Correction_Agl_Graph[S]->Draw("AP4C");

	
	}
	
	



	cout<<"*** Updating Results file ***"<<endl;
	nomefile="./Final_plots/"+mese+".root";
	TFile *f_out=new TFile(nomefile.c_str(), "UPDATE");
	f_out->mkdir("DATA-driven Results/Data vs MC/Protons");
	f_out->cd("DATA-driven Results/Data vs MC/Protons");
		
	for(int S=0;S<3;S++){
		c20[S]->Write();
		c21[S]->Write();
	}
	f_out->Write();
	f_out->Close();

	return;
}

