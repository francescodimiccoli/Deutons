using namespace std;

DatavsMC * RICH_DvsMC_P = new DatavsMC("RICH_DvsMC_P",11);
DatavsMC * RICH_DvsMC_D = new DatavsMC("RICH_DvsMC_D",11,1,6);


void DVSMCRICHeff_D_Fill(int zona){

	//cuts
	if(Tup.R<1.2*Tup.Rcutoff||Tup.Beta>protons->Eval(Tup.R)+0.1||Tup.Beta<protons->Eval(Tup.R)-0.1) return;
	if(!((Tup.R>Rcut[zona]&&zona<10)||(zona==10)))  return;
	if(!Herejcut) return;
	//
	int Kbin;
	
	//Beta bins
	//NaF
	Kbin=NaFDB.GetRBin(RUsed);
	RICH_DvsMC_P -> DataEff -> beforeNaF -> Fill(Kbin,zona);	
	RICH_DvsMC_D -> DataEff -> beforeNaF -> Fill(Kbin,zona);
	if(cmask.isFromNaF()) {
        RICH_DvsMC_P -> DataEff -> afterNaF -> Fill(Kbin,zona);
        RICH_DvsMC_D -> DataEff -> afterNaF -> Fill(Kbin,zona);
    }
	//Agl
	Kbin=AglDB.GetRBin(RUsed);
    RICH_DvsMC_P -> DataEff -> beforeAgl -> Fill(Kbin,zona);	
	RICH_DvsMC_D -> DataEff -> beforeAgl -> Fill(Kbin,zona);
	if(cmask.isFromAgl()) { 
        RICH_DvsMC_P -> DataEff -> afterAgl -> Fill(Kbin,zona); 
        RICH_DvsMC_D -> DataEff -> afterAgl -> Fill(Kbin,zona); 
    }
	return;
}

void DVSMCRICHeff_Fill(){

	//
	int Kbin;

	if(Massa_gen<1) {
		//Beta bins
		//NaF
		Kbin=NaFDB.GetRBin(RUsed);
		RICH_DvsMC_P -> MCEff -> beforeNaF -> Fill(Kbin);
		for(int mc_type=0;mc_type<6;mc_type++) RICH_DvsMC_D -> MCEff -> beforeNaF -> Fill(Kbin,mc_type);

		if(cmask.isFromNaF()){
			RICH_DvsMC_P -> MCEff -> afterNaF -> Fill(Kbin);
			for(int mc_type=0;mc_type<6;mc_type++) RICH_DvsMC_D -> MCEff -> afterNaF -> Fill(Kbin,mc_type);
		}
		//Agl
    	Kbin=AglDB.GetRBin(RUsed);
		RICH_DvsMC_P -> MCEff -> beforeAgl -> Fill(Kbin);
		for(int mc_type=0;mc_type<6;mc_type++) RICH_DvsMC_D -> MCEff -> beforeAgl -> Fill(Kbin,mc_type);	

		if(cmask.isFromAgl()) {
			RICH_DvsMC_P -> MCEff -> afterAgl -> Fill(Kbin); 	
			for(int mc_type=0;mc_type<6;mc_type++) 	RICH_DvsMC_D -> MCEff -> afterAgl -> Fill(Kbin,mc_type);
		}
	}
}                        


void DVSMCRICHeff_Write(){

	RICH_DvsMC_P -> Write();
	RICH_DvsMC_D -> Write();

	return;
}


void DVSMCRICHeff(){

   inputHistoFile->ReOpen("READ");

	DatavsMC * RICH_DvsMC_P = new DatavsMC(inputHistoFile,"RICH_DvsMC_P");
	DatavsMC * RICH_DvsMC_D = new DatavsMC(inputHistoFile,"RICH_DvsMC_D",6);


	LATcorr * LATrichDATA_TOF       = new LATcorr(inputHistoFile,"LATrichDATA_Agl" 	 ,"Results");
	LATcorr * LATrichDATA_NaF       = new LATcorr(inputHistoFile,"LATrichDATA_NaF" 	 ,"Results");
	LATcorr * LATrichDATA_Agl       = new LATcorr(inputHistoFile,"LATrichDATA_Agl" 	 ,"Results");	


	cout<<"******* Data vs MC: RICH ********"<<endl;

	RICH_DvsMC_P -> Assign_LatCorr( LATrichDATA_TOF   ->  LATcorrR_fit , 
			LATrichDATA_TOF   ->  LATcorrR_fit ,
			LATrichDATA_NaF   ->  LATcorrR_fit ,
			LATrichDATA_Agl   ->  LATcorrR_fit );

	RICH_DvsMC_D -> Assign_LatCorr( LATrichDATA_TOF   ->  LATcorrR_fit , 
			LATrichDATA_TOF   ->  LATcorrR_fit ,
			LATrichDATA_NaF   ->  LATcorrR_fit ,
			LATrichDATA_Agl   ->  LATcorrR_fit );

	RICH_DvsMC_P ->Eval_DandMC_Eff();  
	RICH_DvsMC_D ->Eval_DandMC_Eff();

	RICH_DvsMC_P ->Eval_Corrections();
	RICH_DvsMC_D ->Eval_Corrections();


	TH1F* RICH_Correction_P_NaF =(TH1F*) RICH_DvsMC_P -> GetCorrection_NaF();
	TH1F* RICH_Correction_P_Agl =(TH1F*) RICH_DvsMC_P -> GetCorrection_Agl();

	TH2F* RICH_Correction_D_NaF =(TH2F*) RICH_DvsMC_D -> GetCorrection_NaF();
	TH2F* RICH_Correction_D_Agl =(TH2F*) RICH_DvsMC_D -> GetCorrection_Agl();



	cout<<"*** Updating P1 file ****"<<endl;
   inputHistoFile->ReOpen("UPDATE");

	inputHistoFile->cd("Results");

	RICH_Correction_P_NaF  -> Write("RICH_DvsMC_P_CorrectionNaF");
	RICH_Correction_P_Agl  -> Write("RICH_DvsMC_P_CorrectionAgl");

	RICH_Correction_D_NaF  -> Write("RICH_DvsMC_D_CorrectionNaF");
	RICH_Correction_D_Agl  -> Write("RICH_DvsMC_D_CorrectionAgl");

	inputHistoFile->Write();


	TCanvas *c20_bis=new TCanvas("Data vs MC: RICH");


	c20_bis->Divide(2,1);

	c20_bis->cd(1);
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraphErrors * RICHDVSMC_P_GraphNaF=new TGraphErrors();
	RICHDVSMC_P_GraphNaF->SetName("RICHDVSMC_P_GraphNaF");
	int j=0;
	for(int i=1;i<nbinsNaF;i++) {
		if(RICH_Correction_P_NaF -> GetBinContent(i+1)>0){
			RICHDVSMC_P_GraphNaF->SetPoint(j,NaFPB.EkBinCent(i),RICH_Correction_P_NaF -> GetBinContent(i+1));
			RICHDVSMC_P_GraphNaF->SetPointError(j,0,RICH_Correction_P_NaF -> GetBinError(i+1));
			j++;
		}
	}
	RICHDVSMC_P_GraphNaF->SetLineColor(2);
	RICHDVSMC_P_GraphNaF->SetFillColor(2);
	RICHDVSMC_P_GraphNaF->SetFillStyle(3001);
	RICHDVSMC_P_GraphNaF->SetLineWidth(4);
	RICHDVSMC_P_GraphNaF->Draw("AP4C");

	c20_bis->cd(2);
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraphErrors * RICHDVSMC_P_GraphAgl=new TGraphErrors();
        RICHDVSMC_P_GraphAgl->SetName("RICHDVSMC_P_GraphAgl");
        j=0;
        for(int i=1;i<nbinsAgl;i++) {
                if(RICH_Correction_P_Agl -> GetBinContent(i+1)>0){
                        RICHDVSMC_P_GraphAgl->SetPoint(j,AglPB.EkBinCent(i),RICH_Correction_P_Agl -> GetBinContent(i+1));
                        RICHDVSMC_P_GraphAgl->SetPointError(j,0,RICH_Correction_P_Agl -> GetBinError(i+1));
                        j++;
                }
        }
        RICHDVSMC_P_GraphAgl->SetLineColor(2);
        RICHDVSMC_P_GraphAgl->SetFillColor(2);
        RICHDVSMC_P_GraphAgl->SetFillStyle(3001);
        RICHDVSMC_P_GraphAgl->SetLineWidth(4);
        RICHDVSMC_P_GraphAgl->Draw("AP4C");


	finalPlots.Add(c20_bis);
	finalPlots.writeObjsInFolder("DATA-driven Results/Data vs MC/RICH");
	

	return;
}
