using namespace std;

DatavsMC * Dist_DvsMC_P = new DatavsMC("Dist_DvsMC_P",11);
DatavsMC * Lik_DvsMC_P  = new DatavsMC("Lik_DvsMC_D" ,11);

void DVSMCQualeff2_D_Fill(TNtuple *ntupla, int l,int zona){

	int k = ntupla->GetEvent(l);
	//cuts
	if(Beta<=0||R<=0||R<1.2*Rcutoff||Beta>protons->Eval(R)+0.1||Beta<protons->Eval(R)-0.1) return;
	if(!((R>Rcut[zona]&&zona<10)||(zona==10)))  return;
	if(!Herejcut) return;
	//
	int Kbin;
	
	//R bins
	Kbin = GetRBin(R);
	Dist_DvsMC_P -> DataEff -> beforeR -> Fill(Kbin,zona);	
	if(Dist5D_P<6) Lik_DvsMC_P  -> DataEff -> beforeR -> Fill(Kbin,zona);
	
	if(Dist5D_P<6){
		Dist_DvsMC_P -> DataEff -> afterR -> Fill(Kbin,zona);     
        	if(Likcut) Lik_DvsMC_P  -> DataEff -> afterR -> Fill(Kbin,zona);
	}

	//Beta bins
	//ToF
	Kbin=GetArrayBin(Var, BetaD, nbinsToF);	
	Dist_DvsMC_P -> DataEff -> beforeTOF -> Fill(Kbin,zona);
	if(Distcut) Lik_DvsMC_P  -> DataEff -> beforeTOF -> Fill(Kbin,zona);

	if(Distcut){
		Dist_DvsMC_P -> DataEff -> afterTOF -> Fill(Kbin,zona);
		if(Likcut) Lik_DvsMC_P  -> DataEff -> afterTOF -> Fill(Kbin,zona);
	}
	//NaF
	if(((int)Cutmask)>>11==512) {	
		Kbin=GetArrayBin(Var, BetaNaFD, nbinsNaF);
		Dist_DvsMC_P -> DataEff -> beforeNaF -> Fill(Kbin,zona);
		if(Distcut) Lik_DvsMC_P  -> DataEff -> beforeNaF -> Fill(Kbin,zona);

		if(Distcut){
			Dist_DvsMC_P -> DataEff -> afterNaF -> Fill(Kbin,zona);
			if(Likcut) Lik_DvsMC_P  -> DataEff -> afterNaF -> Fill(Kbin,zona);
		}
	}
	//Agl
	if(((int)Cutmask)>>11==0) {
		Kbin=GetArrayBin(Var, BetaAglD, nbinsAgl);
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

	int k = ntupla->GetEvent(l);
	//cuts
	if(Beta<=0||R<=0||R<1.2*Rcutoff||Beta>protons->Eval(R)+0.1||Beta<protons->Eval(R)-0.1) return;
	if(!Herejcut) return;
	//
	int Kbin;

	//R bins
	Kbin = GetRBin(R);

	if(Massa_gen<1) {
		//R bins
		Kbin = GetRBin(R);	
		Dist_DvsMC_P -> MCEff -> beforeR -> Fill(Kbin);
		if(Dist5D_P<6) Lik_DvsMC_P  -> MCEff -> beforeR -> Fill(Kbin);

		if(Dist5D_P<6){
			Dist_DvsMC_P -> MCEff -> afterR -> Fill(Kbin);
			if(Likcut) Lik_DvsMC_P  -> MCEff -> afterR -> Fill(Kbin);
		}
		//Beta bins

		//ToF
		Kbin=GetArrayBin(Var, BetaD, nbinsToF);	
		Dist_DvsMC_P -> MCEff -> beforeTOF -> Fill(Kbin);
		if(Distcut) Lik_DvsMC_P  -> MCEff -> beforeTOF -> Fill(Kbin);

		if(Distcut){
			Dist_DvsMC_P -> MCEff -> afterTOF -> Fill(Kbin);
			if(Likcut) Lik_DvsMC_P  -> MCEff -> afterTOF -> Fill(Kbin);
		}
		//NaF
		if(((int)Cutmask)>>11==512) {	
			Kbin=GetArrayBin(Var, BetaNaFD, nbinsNaF);	
			Dist_DvsMC_P -> MCEff -> beforeNaF -> Fill(Kbin);
			if(Distcut) Lik_DvsMC_P  -> MCEff -> beforeNaF -> Fill(Kbin);

			if(Distcut){
				Dist_DvsMC_P -> MCEff -> afterNaF -> Fill(Kbin);
				if(Likcut) Lik_DvsMC_P  -> MCEff -> afterNaF -> Fill(Kbin);
			}

		}
		//Agl
		if(((int)Cutmask)>>11==0) {	
			Kbin=GetArrayBin(Var, BetaAglD, nbinsAgl);
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
	DatavsMC * Lik_DvsMC_P  = new DatavsMC(file1,"Lik_DvsMC_D" );

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
		LikDVSMC_P_Graph->SetPoint(j,R_cent[i],LikP_Correction_R -> GetBinContent(i+1));
		LikDVSMC_P_Graph->SetPointError(j,0,LikP_Correction_R -> GetBinError(i+1));
		j++;
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
		DistDVSMC_P_Graph->SetPoint(j,R_cent[i],DistP_Correction_R -> GetBinContent(i+1));
		DistDVSMC_P_Graph->SetPointError(j,0,DistP_Correction_R -> GetBinError(i+1));
		j++;
	}
	DistDVSMC_P_Graph->SetLineColor(2);
        DistDVSMC_P_Graph->SetFillColor(2);
        DistDVSMC_P_Graph->SetFillStyle(3001);
        DistDVSMC_P_Graph->SetLineWidth(4);
	DistDVSMC_P_Graph->Draw("AP4C");

	cout<<"*** Updating Results file ***"<<endl;
	nomefile="./Final_plots/"+mese+".root";
	TFile *f_out=new TFile(nomefile.c_str(), "UPDATE");
	f_out->mkdir("DATA-driven Results/Data vs MC/Protons");
   	f_out->cd("DATA-driven Results/Data vs MC/Protons");
	c20->Write();
	c21->Write();
	f_out->Write();
	f_out->Close();


	



	return;
}

/*

TCanvas *c20=new TCanvas("Data vs MC: Likelihood");
TCanvas *c21=new TCanvas("Data vs MC: Distance");
TGraphErrors *LikDVSMC_P_Graph;
TGraphErrors *DistDVSMC_P_Graph;
TH2F * LikDVSMC_P_graph=new TH2F("LikDVSMC_P_graph","LikDVSMC_P_graph",nbinsr,0,nbinsr,2,0,2);
TH2F * DistDVSMC_P_graph=new TH2F("DistDVSMC_P_graph","DistDVSMC_P_graph",nbinsr,0,nbinsr,2,0,2);


void DVSMCQualeff2(TFile * file1){
	

	c20->cd();
	gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
	LikDVSMC_P_Graph=new TGraphErrors();
	LikDVSMC_P_Graph->SetName("LikDVSMC_P_Graph");
	j=0;
	
	for(int i=1;i<nbinsr;i++) {
		if(EffLik2MCvsDP_MC[i]>0){
		LikDVSMC_P_Graph->SetPoint(j,R_cent[i],ratioL_smooth[i]);
		LikDVSMC_P_graph->SetBinContent(i+1,1,ratioL_smooth[i]);
		LikDVSMC_P_Graph->SetPointError(j,0,errorefitL);
		LikDVSMC_P_graph->SetBinContent(i+1,2,errorefitL);
	j++;
	}
	}
	LikDVSMC_P_Graph->SetLineColor(4);
	LikDVSMC_P_Graph->SetFillColor(4);
	LikDVSMC_P_Graph->SetFillStyle(3001);
        LikDVSMC_P_Graph->SetLineWidth(4);
	EffLik2MCvsDP_Mean->SetMarkerColor(2);
        EffLik2MCvsDP_Mean->SetMarkerStyle(4);
        EffLik2MCvsDP_Mean->SetLineColor(2);
        EffLik2MCvsDP_Mean->SetLineWidth(2);
        EffMCvsDLikP_MC->SetMarkerColor(2);
        EffMCvsDLikP_MC->SetMarkerStyle(8);
        EffMCvsDLikP_MC->SetLineColor(2);
        EffMCvsDLikP_MC->SetLineWidth(2);
	EffLik2MCvsDP_Mean->SetTitle("Likelihood Efficiency DATA/MC");
        EffLik2MCvsDP_Mean->GetXaxis()->SetTitle("R [GV]");
        EffLik2MCvsDP_Mean->GetYaxis()->SetTitle("Efficiency");
        EffLik2MCvsDP_Mean->GetXaxis()->SetTitleSize(0.045);
        EffLik2MCvsDP_Mean->GetYaxis()->SetTitleSize(0.045);
	EffLik2MCvsDP_Mean->GetYaxis()->SetRangeUser(0.7,1.4);
	{
                EffLik2MCvsDP_Mean->Draw("AP");
		EffMCvsDLikP_MC->Draw("CPsame");
                LikDVSMC_P_Graph->Draw("L4same");
		TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                //leg->AddEntry(EffMCLikP,MCLegend[0].c_str(), "ep");

        }

	c21->cd();
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
	DistDVSMC_P_Graph=new TGraphErrors();
        DistDVSMC_P_Graph->SetName("DistDVSMC_P_Graph");
	j=0;
        for(int i=1;i<nbinsr;i++) {
                if(EffLik2MCvsDP_MC[i]>0){
                DistDVSMC_P_Graph->SetPoint(j,R_cent[i],ratioD_smooth[i]);
               	DistDVSMC_P_graph->SetBinContent(i+1,1,ratioD_smooth[i]);
		DistDVSMC_P_Graph->SetPointError(j,0,errorefitD);
        	DistDVSMC_P_graph->SetBinContent(i+1,2,errorefitD);
	j++;
        }
        }
        DistDVSMC_P_Graph->SetLineColor(4);
        DistDVSMC_P_Graph->SetFillColor(4);
        DistDVSMC_P_Graph->SetFillStyle(3001);
        DistDVSMC_P_Graph->SetLineWidth(4);
	EffDistMCvsDP_Mean->SetMarkerColor(2);
        EffDistMCvsDP_Mean->SetMarkerStyle(4);
        EffDistMCvsDP_Mean->SetLineColor(2);
        EffDistMCvsDP_Mean->SetLineWidth(2);
	EffMCvsDDistP_MC->SetMarkerColor(2);
        EffMCvsDDistP_MC->SetMarkerStyle(8);
        EffMCvsDDistP_MC->SetLineColor(2);
        EffMCvsDDistP_MC->SetLineWidth(2);
        EffDistMCvsDP_Mean->SetTitle("Distance Efficiency DATA/MC");
        EffDistMCvsDP_Mean->GetXaxis()->SetTitle("R [GV]");
        EffDistMCvsDP_Mean->GetYaxis()->SetTitle("Efficiency");
        EffDistMCvsDP_Mean->GetXaxis()->SetTitleSize(0.045);
        EffDistMCvsDP_Mean->GetYaxis()->SetTitleSize(0.045);
	EffDistMCvsDP_Mean->GetYaxis()->SetRangeUser(0.7,1.4);
        {
                EffDistMCvsDP_Mean->Draw("AP");
		EffMCvsDDistP_MC->Draw("CPsame");
                DistDVSMC_P_Graph->Draw("L4same");
		TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                //leg->AddEntry(EffMCDistP,MCLegend[0].c_str(), "ep");

        }	


}
*/
