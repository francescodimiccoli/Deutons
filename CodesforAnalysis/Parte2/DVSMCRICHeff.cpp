using namespace std;

DatavsMC * RICH_DvsMC_P = new DatavsMC("RICH_DvsMC_P",11);

void DVSMCRICHeff_D_Fill(TNtuple *ntupla, int l,int zona){

	 ntupla->GetEvent(l);
	//cuts
	if(Tup.Beta<=0||Tup.R<=0||Tup.R<1.2*Tup.Rcutoff||Tup.Beta>protons->Eval(Tup.R)+0.1||Tup.Beta<protons->Eval(Tup.R)-0.1) return;
	if(!((Tup.R>Rcut[zona]&&zona<10)||(zona==10)))  return;
	if(!Herejcut) return;
	//
	int Kbin;
	
	//Beta bins
	//NaF
	Kbin=NaFDB.GetRBin(RUsed);
	RICH_DvsMC_P -> DataEff -> beforeNaF -> Fill(Kbin,zona);	
	if(((int)Tup.Cutmask)>>11==512) RICH_DvsMC_P -> DataEff -> afterNaF -> Fill(Kbin,zona);
	//Agl
	Kbin=AglDB.GetRBin(RUsed);
        RICH_DvsMC_P -> DataEff -> beforeAgl -> Fill(Kbin,zona);	
	if(((int)Tup.Cutmask)>>11==0) RICH_DvsMC_P -> DataEff -> afterAgl -> Fill(Kbin,zona); 
	return;

}

void DVSMCRICHeff_Fill(TNtuple *ntupla, int l){

	 ntupla->GetEvent(l);
	//cuts
	if(Tup.Beta<=0||Tup.R<=0||Tup.R<1.2*Tup.Rcutoff||Tup.Beta>protons->Eval(Tup.R)+0.1||Tup.Beta<protons->Eval(Tup.R)-0.1) return;
	if(!Herejcut) return;
	//
	int Kbin;

	if(Massa_gen<1) {
		//Beta bins
		Kbin=NaFDB.GetRBin(RUsed);
                RICH_DvsMC_P -> MCEff -> beforeNaF -> Fill(Kbin);
		//NaF
		if(((int)Tup.Cutmask)>>11==512) RICH_DvsMC_P -> MCEff -> afterNaF -> Fill(Kbin);

		//Agl
		Kbin=AglDB.GetRBin(RUsed);
                RICH_DvsMC_P -> MCEff -> beforeAgl -> Fill(Kbin);	
		if(((int)Tup.Cutmask)>>11==0) RICH_DvsMC_P -> MCEff -> afterAgl -> Fill(Kbin); 	

	}                        
}

void DVSMCRICHeff_Write(){

	RICH_DvsMC_P -> Write();
	return;
}


void DVSMCRICHeff(){

	string nomefile="../Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
	TFile * file1 =TFile::Open(nomefile.c_str(),"READ");

	DatavsMC * RICH_DvsMC_P = new DatavsMC(file1,"RICH_DvsMC_P");

	LATcorr * LATrichDATA_TOF       = new LATcorr(file1,"LATrichDATA_Agl" 	 ,"Results");
	LATcorr * LATrichDATA_NaF       = new LATcorr(file1,"LATrichDATA_NaF" 	 ,"Results");
   	LATcorr * LATrichDATA_Agl       = new LATcorr(file1,"LATrichDATA_Agl" 	 ,"Results");	


	cout<<"******* Data vs MC: QUALITY SEL ********"<<endl;

	RICH_DvsMC_P -> Assign_LatCorr( LATrichDATA_TOF   ->  LATcorrR_fit , 
					LATrichDATA_TOF   ->  LATcorrR_fit ,
					LATrichDATA_NaF   ->  LATcorrR_fit ,
					LATrichDATA_Agl   ->  LATcorrR_fit );



	RICH_DvsMC_P ->Eval_DandMC_Eff();  
	
	RICH_DvsMC_P ->Eval_Corrections();


	TH1F* RICH_Correction_NaF =(TH1F*) RICH_DvsMC_P -> GetCorrection_NaF();
	TH1F* RICH_Correction_Agl =(TH1F*) RICH_DvsMC_P -> GetCorrection_Agl();


	cout<<"*** Updating P1 file ****"<<endl;
	nomefile="../Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
	file1 =TFile::Open(nomefile.c_str(),"UPDATE");

	file1->cd("Results");

	RICH_Correction_NaF  -> Write("RICH_DvsMC_P_CorrectionNaF");
	RICH_Correction_Agl  -> Write("RICH_DvsMC_P_CorrectionAgl");


	file1->Write();
	file1->Close();


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
		if(RICH_Correction_NaF -> GetBinContent(i+1)>0){
			RICHDVSMC_P_GraphNaF->SetPoint(j,NaFPB.EkBinCent(i),RICH_Correction_NaF -> GetBinContent(i+1));
			RICHDVSMC_P_GraphNaF->SetPointError(j,0,RICH_Correction_NaF -> GetBinError(i+1));
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
                if(RICH_Correction_Agl -> GetBinContent(i+1)>0){
                        RICHDVSMC_P_GraphAgl->SetPoint(j,AglPB.EkBinCent(i),RICH_Correction_Agl -> GetBinContent(i+1));
                        RICHDVSMC_P_GraphAgl->SetPointError(j,0,RICH_Correction_Agl -> GetBinError(i+1));
                        j++;
                }
        }
        RICHDVSMC_P_GraphAgl->SetLineColor(2);
        RICHDVSMC_P_GraphAgl->SetFillColor(2);
        RICHDVSMC_P_GraphAgl->SetFillStyle(3001);
        RICHDVSMC_P_GraphAgl->SetLineWidth(4);
        RICHDVSMC_P_GraphAgl->Draw("AP4C");


	cout<<"*** Updating Results file ***"<<endl;
	nomefile="./Final_plots/"+mese+".root";
	TFile *f_out=new TFile(nomefile.c_str(), "UPDATE");
	f_out->mkdir("DATA-driven Results/Data vs MC/RICH");
	f_out->cd("DATA-driven Results/Data vs MC/RICH");
	c20_bis->Write();

	f_out->Write();
	f_out->Close();






	return;
}

