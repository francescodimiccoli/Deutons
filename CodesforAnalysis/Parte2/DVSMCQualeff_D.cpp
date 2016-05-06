using namespace std;

DatavsMC * Dist_DvsMC_D = new DatavsMC("Dist_DvsMC_D",11,1,6);
DatavsMC * Lik_DvsMC_D  = new DatavsMC("Lik_DvsMC_D" ,11,1,6);

void DVSMCQualeffD_D_Fill(TNtuple *ntupla, int l,int zona){

	int k = ntupla->GetEvent(l);
	//cuts
	if(Beta<=0||R<=0||R<1.2*Rcutoff||Beta>deutons->Eval(R)+0.1||Beta<deutons->Eval(R)-0.1) return;
	if(!((R>Rcut[zona]&&zona<10)||(zona==10)))  return;
	if(!Herejcut) return;
	if(!Betastrongcut) return;
	//
	int Kbin;
	
	//Beta bins
	//ToF
	Kbin=GetArrayBin(Var, BetaD, nbinsToF);	
	Dist_DvsMC_D -> DataEff -> beforeTOF -> Fill(Kbin,zona);
	if(Distcut) Lik_DvsMC_D  -> DataEff -> beforeTOF -> Fill(Kbin,zona);

	if(Distcut){
		Dist_DvsMC_D -> DataEff -> afterTOF -> Fill(Kbin,zona);
		if(Likcut) Lik_DvsMC_D  -> DataEff -> afterTOF -> Fill(Kbin,zona);
	}
	//NaF
	if(((int)Cutmask)>>11==512) {	
		Kbin=GetArrayBin(Var, BetaNaFD, nbinsNaF);
		Dist_DvsMC_D -> DataEff -> beforeNaF -> Fill(Kbin,zona);
		if(Distcut) Lik_DvsMC_D  -> DataEff -> beforeNaF -> Fill(Kbin,zona);

		if(Distcut){
			Dist_DvsMC_D -> DataEff -> afterNaF -> Fill(Kbin,zona);
			if(Likcut) Lik_DvsMC_D  -> DataEff -> afterNaF -> Fill(Kbin,zona);
		}
	}
	//Agl
	if(((int)Cutmask)>>11==0) {
		Kbin=GetArrayBin(Var, BetaAglD, nbinsAgl);
		Dist_DvsMC_D -> DataEff -> beforeAgl -> Fill(Kbin,zona);
		if(Distcut) Lik_DvsMC_D  -> DataEff -> beforeAgl -> Fill(Kbin,zona);

		if(Distcut){
			Dist_DvsMC_D -> DataEff -> afterAgl -> Fill(Kbin,zona);
			if(Likcut) Lik_DvsMC_D  -> DataEff -> afterAgl -> Fill(Kbin,zona);
		}
	}
	return;

}

void DVSMCQualeffD_Fill(TNtuple *ntupla, int l){

	int k = ntupla->GetEvent(l);
	//cuts
	if(Beta<=0||R<=0||R<1.2*Rcutoff||Beta>protons->Eval(R)+0.1||Beta<protons->Eval(R)-0.1) return;
	if(!Herejcut) return;
	if(!Betastrongcut) return;
	//
	int Kbin;

	if(Massa_gen>1&&Massa_gen<2) {
		//Beta bins

		//ToF
		Kbin=GetArrayBin(Var, BetaD, nbinsToF);	
		Dist_DvsMC_D -> MCEff -> beforeTOF -> Fill(Kbin,ReturnMCGenType());
		if(Distcut) Lik_DvsMC_D  -> MCEff -> beforeTOF -> Fill(Kbin,ReturnMCGenType());

		if(Distcut){
			Dist_DvsMC_D -> MCEff -> afterTOF -> Fill(Kbin,ReturnMCGenType());
			if(Likcut) Lik_DvsMC_D  -> MCEff -> afterTOF -> Fill(Kbin,ReturnMCGenType());
		}
		//NaF
		if(((int)Cutmask)>>11==512) {	
			Kbin=GetArrayBin(Var, BetaNaFD, nbinsNaF);	
			Dist_DvsMC_D -> MCEff -> beforeNaF -> Fill(Kbin,ReturnMCGenType());
			if(Distcut) Lik_DvsMC_D  -> MCEff -> beforeNaF -> Fill(Kbin,ReturnMCGenType());

			if(Distcut){
				Dist_DvsMC_D -> MCEff -> afterNaF -> Fill(Kbin,ReturnMCGenType());
				if(Likcut) Lik_DvsMC_D  -> MCEff -> afterNaF -> Fill(Kbin,ReturnMCGenType());
			}

		}
		//Agl
		if(((int)Cutmask)>>11==0) {	
			Kbin=GetArrayBin(Var, BetaAglD, nbinsAgl);
			Dist_DvsMC_D -> MCEff -> beforeAgl -> Fill(Kbin,ReturnMCGenType());
			if(Distcut) Lik_DvsMC_D  -> MCEff -> beforeAgl -> Fill(Kbin,ReturnMCGenType());

			if(Distcut){
				Dist_DvsMC_D -> MCEff -> afterAgl -> Fill(Kbin,ReturnMCGenType());
				if(Likcut) Lik_DvsMC_D  -> MCEff -> afterAgl -> Fill(Kbin,ReturnMCGenType());
			}

		}

	}                        
}

void DVSMCQualeffD_Write(){

	Dist_DvsMC_D -> Write();
	Lik_DvsMC_D  -> Write();

	return;
}


void DVSMCQualeffD(){

	string nomefile="../Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
	TFile * file1 =TFile::Open(nomefile.c_str(),"READ");

	DatavsMC * Dist_DvsMC_D = new DatavsMC(file1,"Dist_DvsMC_D");
	DatavsMC * Lik_DvsMC_D  = new DatavsMC(file1,"Lik_DvsMC_D" );

	LATcorr * LATLikelihoodDATA_TOF = new LATcorr(file1,"LATLikDATA_TOF"   	 ,"Results");
	LATcorr * LATDistanceDATA_TOF   = new LATcorr(file1,"LATDistDATA_TOF" 	 ,"Results");

	LATcorr * LATLikelihoodDATA_NaF = new LATcorr(file1,"LATLikDATA_NaF"  	 ,"Results");
	LATcorr * LATDistanceDATA_NaF   = new LATcorr(file1,"LATDistDATA_NaF" 	 ,"Results");

	LATcorr * LATLikelihoodDATA_Agl = new LATcorr(file1,"LATLikDATA_Agl"  	 ,"Results");
	LATcorr * LATDistanceDATA_Agl   = new LATcorr(file1,"LATDistDATA_Agl" 	 ,"Results");




	cout<<"******* Data vs MC: QUALITY SEL ********"<<endl;

	Dist_DvsMC_D -> Assign_LatCorr( LATDistanceDATA_TOF   ->  LATcorrR_fit , 
					LATDistanceDATA_TOF   ->  LATcorrR_fit ,
					LATDistanceDATA_NaF   ->  LATcorrR_fit ,
					LATDistanceDATA_Agl   ->  LATcorrR_fit );

	Lik_DvsMC_D  ->	Assign_LatCorr( LATLikelihoodDATA_TOF ->  LATcorrR_fit , 	
					LATLikelihoodDATA_TOF ->  LATcorrR_fit ,
					LATLikelihoodDATA_NaF ->  LATcorrR_fit ,
					LATLikelihoodDATA_Agl ->  LATcorrR_fit );



	Dist_DvsMC_D ->Eval_DandMC_Eff();  
	Lik_DvsMC_D  ->Eval_DandMC_Eff();

	Dist_DvsMC_D ->Eval_Corrections();
	Lik_DvsMC_D  ->Eval_Corrections();


	TH2F* DistD_Correction_R   =(TH2F*) Dist_DvsMC_D -> GetCorrection_R()  ;
	TH2F* DistD_Correction_TOF =(TH2F*) Dist_DvsMC_D -> GetCorrection_TOF();
	TH2F* DistD_Correction_NaF =(TH2F*) Dist_DvsMC_D -> GetCorrection_NaF();
	TH2F* DistD_Correction_Agl =(TH2F*) Dist_DvsMC_D -> GetCorrection_Agl();

	TH2F* LikD_Correction_R    =(TH2F*) Lik_DvsMC_D -> GetCorrection_R()  ;
	TH2F* LikD_Correction_TOF  =(TH2F*) Lik_DvsMC_D -> GetCorrection_TOF();
	TH2F* LikD_Correction_NaF  =(TH2F*) Lik_DvsMC_D -> GetCorrection_NaF();
	TH2F* LikD_Correction_Agl  =(TH2F*) Lik_DvsMC_D -> GetCorrection_Agl();


	cout<<"*** Updating P1 file ****"<<endl;
	nomefile="../Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
	file1 =TFile::Open(nomefile.c_str(),"UPDATE");

	file1->cd("Results");

	DistD_Correction_R    -> Write("Dist_DvsMC_D_CorrectionR"  );
	DistD_Correction_TOF  -> Write("Dist_DvsMC_D_CorrectionTOF");
	DistD_Correction_NaF  -> Write("Dist_DvsMC_D_CorrectionNaF");
	DistD_Correction_Agl  -> Write("Dist_DvsMC_D_CorrectionAgl");

	LikD_Correction_R    -> Write("Lik_DvsMC_D_CorrectionR"  );
	LikD_Correction_TOF  -> Write("Lik_DvsMC_D_CorrectionTOF");
	LikD_Correction_NaF  -> Write("Lik_DvsMC_D_CorrectionNaF");
	LikD_Correction_Agl  -> Write("Lik_DvsMC_D_CorrectionAgl");


	file1->Write();
	file1->Close();



	



	TCanvas *c20_bis=new TCanvas("Deutons Data vs MC: Likelihood (Beta bins)");
	TCanvas *c21_bis=new TCanvas("Deutons Data vs MC: Distance (Beta bins)");

	string MCLegend[6]= {"d.pl1.0_520_GG_Blic","d.pl1.0_520_GG_BlicDPMJet","d.pl1.0_520_GG_QMD","d.pl1.0_520_Shen_Blic","d.pl1.0_520_Shen_BlicDPMJet","d.pl1.0_520_Shen_QMD"};

	int j=0;

	c20_bis->Divide(3,1);

	c20_bis->cd(1);
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	TGraphErrors * LikDVSMC_D_GraphTOF[6];
	for(int mc_type=0;mc_type<6;mc_type++){
		LikDVSMC_D_GraphTOF[mc_type]=new TGraphErrors();
		j=0;
		for(int i=1;i<nbinsToF;i++) {
				LikDVSMC_D_GraphTOF[mc_type]->SetPoint(j,Ekincent[i],LikD_Correction_TOF -> GetBinContent(i+1,mc_type+1));
				LikDVSMC_D_GraphTOF[mc_type]->SetPointError(j,0,LikD_Correction_TOF -> GetBinError(i+1,mc_type+1));
				j++;
		}
		LikDVSMC_D_GraphTOF[mc_type]->SetLineColor(4);
		LikDVSMC_D_GraphTOF[mc_type]->SetFillColor(4);
		LikDVSMC_D_GraphTOF[mc_type]->SetFillStyle(3001);
		LikDVSMC_D_GraphTOF[mc_type]->SetLineWidth(1);
		LikDVSMC_D_GraphTOF[mc_type]->SetMarkerColor(4);
		LikDVSMC_D_GraphTOF[mc_type]->SetMarkerStyle(mc_type+3);
	}
	LikDVSMC_D_GraphTOF[0]->Draw("AP4C");
	for(int mc_type=1;mc_type<6;mc_type++) LikDVSMC_D_GraphTOF[mc_type]->Draw("P4Csame");

	c20_bis->cd(2);
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
        TGraphErrors * LikDVSMC_D_GraphNaF[6];
        for(int mc_type=0;mc_type<6;mc_type++){
                LikDVSMC_D_GraphNaF[mc_type]=new TGraphErrors();
                j=0;
                for(int i=1;i<nbinsToF;i++) {
                        if(LikD_Correction_NaF -> GetBinContent(i+1,mc_type+1)>0){
                                LikDVSMC_D_GraphNaF[mc_type]->SetPoint(j,EkincentNaF[i],LikD_Correction_NaF -> GetBinContent(i+1,mc_type+1));
                                LikDVSMC_D_GraphNaF[mc_type]->SetPointError(j,0,LikD_Correction_NaF -> GetBinError(i+1,mc_type+1));
                                j++;
                        }
                }
		LikDVSMC_D_GraphNaF[mc_type]->SetLineColor(4);
                LikDVSMC_D_GraphNaF[mc_type]->SetFillColor(4);
                LikDVSMC_D_GraphNaF[mc_type]->SetFillStyle(3001);
                LikDVSMC_D_GraphNaF[mc_type]->SetLineWidth(1);
                LikDVSMC_D_GraphNaF[mc_type]->SetMarkerColor(4);
		LikDVSMC_D_GraphNaF[mc_type]->SetMarkerStyle(mc_type+3);
        }
        {
	TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
        leg->AddEntry(LikDVSMC_D_GraphNaF[0],MCLegend[0].c_str(), "ep");
	LikDVSMC_D_GraphNaF[0]->Draw("AP4C");
        for(int mc_type=1;mc_type<6;mc_type++) { LikDVSMC_D_GraphNaF[mc_type]->Draw("P4Csame");
						leg->AddEntry(LikDVSMC_D_GraphNaF[mc_type],MCLegend[mc_type].c_str(), "ep");
						}
	
	leg->Draw("same");
	}
	c20_bis->cd(3);
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
        TGraphErrors * LikDVSMC_D_GraphAgl[6];
        for(int mc_type=0;mc_type<6;mc_type++){
                LikDVSMC_D_GraphAgl[mc_type]=new TGraphErrors();
                j=0;
                for(int i=1;i<nbinsToF;i++) {
                        if(LikD_Correction_Agl -> GetBinContent(i+1,mc_type+1)>0){
                                LikDVSMC_D_GraphAgl[mc_type]->SetPoint(j,EkincentAgl[i],LikD_Correction_Agl -> GetBinContent(i+1,mc_type+1));
                                LikDVSMC_D_GraphAgl[mc_type]->SetPointError(j,0,LikD_Correction_Agl -> GetBinError(i+1,mc_type+1));
                                j++;
                        }
                }
                LikDVSMC_D_GraphAgl[mc_type]->SetLineColor(4);
                LikDVSMC_D_GraphAgl[mc_type]->SetFillColor(4);
                LikDVSMC_D_GraphAgl[mc_type]->SetFillStyle(3001);
                LikDVSMC_D_GraphAgl[mc_type]->SetLineWidth(1);
		LikDVSMC_D_GraphAgl[mc_type]->SetMarkerColor(4);                
		LikDVSMC_D_GraphAgl[mc_type]->SetMarkerStyle(mc_type+3);
        }
        LikDVSMC_D_GraphAgl[0]->Draw("AP4C");
        for(int mc_type=1;mc_type<6;mc_type++) LikDVSMC_D_GraphAgl[mc_type]->Draw("P4Csame");

	
	        c21_bis->Divide(3,1);

        c21_bis->cd(1);
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
        TGraphErrors * DistDVSMC_D_GraphTOF[6];
        for(int mc_type=0;mc_type<6;mc_type++){
                DistDVSMC_D_GraphTOF[mc_type]=new TGraphErrors();
                j=0;
                for(int i=1;i<nbinsToF;i++) {
                                DistDVSMC_D_GraphTOF[mc_type]->SetPoint(j,Ekincent[i],DistD_Correction_TOF -> GetBinContent(i+1,mc_type+1));
                                DistDVSMC_D_GraphTOF[mc_type]->SetPointError(j,0,DistD_Correction_TOF -> GetBinError(i+1,mc_type+1));
                                j++;
                }
                DistDVSMC_D_GraphTOF[mc_type]->SetLineColor(4);
                DistDVSMC_D_GraphTOF[mc_type]->SetFillColor(4);
                DistDVSMC_D_GraphTOF[mc_type]->SetFillStyle(3001);
                DistDVSMC_D_GraphTOF[mc_type]->SetLineWidth(1);
                DistDVSMC_D_GraphTOF[mc_type]->SetMarkerColor(4);
                DistDVSMC_D_GraphTOF[mc_type]->SetMarkerStyle(mc_type+3);
        }
        DistDVSMC_D_GraphTOF[0]->Draw("AP4C");
        for(int mc_type=1;mc_type<6;mc_type++) DistDVSMC_D_GraphTOF[mc_type]->Draw("P4Csame");

        c21_bis->cd(2);
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
        TGraphErrors * DistDVSMC_D_GraphNaF[6];
        for(int mc_type=0;mc_type<6;mc_type++){
                DistDVSMC_D_GraphNaF[mc_type]=new TGraphErrors();
                j=0;
                for(int i=1;i<nbinsToF;i++) {
                        if(DistD_Correction_NaF -> GetBinContent(i+1,mc_type+1)>0){
                                DistDVSMC_D_GraphNaF[mc_type]->SetPoint(j,EkincentNaF[i],DistD_Correction_NaF -> GetBinContent(i+1,mc_type+1));
                                DistDVSMC_D_GraphNaF[mc_type]->SetPointError(j,0,DistD_Correction_NaF -> GetBinError(i+1,mc_type+1));
                                j++;
                        }
                }
                DistDVSMC_D_GraphNaF[mc_type]->SetLineColor(4);
                DistDVSMC_D_GraphNaF[mc_type]->SetFillColor(4);
                DistDVSMC_D_GraphNaF[mc_type]->SetFillStyle(3001);
                DistDVSMC_D_GraphNaF[mc_type]->SetLineWidth(1);
                DistDVSMC_D_GraphNaF[mc_type]->SetMarkerColor(4);
                DistDVSMC_D_GraphNaF[mc_type]->SetMarkerStyle(mc_type+3);
        }
	{
        TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
        leg->AddEntry(LikDVSMC_D_GraphAgl[0],MCLegend[0].c_str(), "ep");
        LikDVSMC_D_GraphAgl[0]->Draw("AP4C");
        for(int mc_type=1;mc_type<6;mc_type++) { LikDVSMC_D_GraphAgl[mc_type]->Draw("P4Csame");
                                                leg->AddEntry(LikDVSMC_D_GraphAgl[mc_type],MCLegend[mc_type].c_str(), "ep");
                                                }
        
        leg->Draw("same");
        }


        c21_bis->cd(3);
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
        TGraphErrors * DistDVSMC_D_GraphAgl[6];
	for(int mc_type=0;mc_type<6;mc_type++){
                DistDVSMC_D_GraphAgl[mc_type]=new TGraphErrors();
                j=0;
                for(int i=1;i<nbinsToF;i++) {
                        if(DistD_Correction_Agl -> GetBinContent(i+1,mc_type+1)>0){
				DistDVSMC_D_GraphAgl[mc_type]->SetPoint(j,EkincentAgl[i],DistD_Correction_Agl -> GetBinContent(i+1,mc_type+1));
                                DistDVSMC_D_GraphAgl[mc_type]->SetPointError(j,0,DistD_Correction_Agl -> GetBinError(i+1,mc_type+1));
                                j++;
                        }
                }
		DistDVSMC_D_GraphAgl[mc_type]->SetLineColor(4);
                DistDVSMC_D_GraphAgl[mc_type]->SetFillColor(4);
                DistDVSMC_D_GraphAgl[mc_type]->SetFillStyle(3001);
                DistDVSMC_D_GraphAgl[mc_type]->SetLineWidth(1);
                DistDVSMC_D_GraphAgl[mc_type]->SetMarkerColor(4);
                DistDVSMC_D_GraphAgl[mc_type]->SetMarkerStyle(mc_type+3);
        }
        DistDVSMC_D_GraphAgl[0]->Draw("AP4C");
        for(int mc_type=1;mc_type<6;mc_type++) DistDVSMC_D_GraphAgl[mc_type]->Draw("P4Csame");



	cout<<"*** Updating Results file ***"<<endl;
	nomefile="./Final_plots/"+mese+".root";
	TFile *f_out=new TFile(nomefile.c_str(), "UPDATE");
	f_out->mkdir("DATA-driven Results/Data vs MC/Deutons");
	f_out->cd("DATA-driven Results/Data vs MC/Deutons");
	c20_bis->Write();
	c21_bis->Write();
	f_out->Write();
	f_out->Close();






	return;
}

