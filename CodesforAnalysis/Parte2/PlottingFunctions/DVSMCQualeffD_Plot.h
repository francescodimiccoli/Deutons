
void DVSMCQualeffD_Plot (TH1 *   DistD_Correction_R  ,
                         TH1 *   DistD_Correction_TOF,
                         TH1 *   DistD_Correction_NaF,
                         TH1 *   DistD_Correction_Agl,

                          TH1 *  LikD_Correction_R   ,
                          TH1 *  LikD_Correction_TOF ,
                          TH1 *  LikD_Correction_NaF ,
                          TH1 *  LikD_Correction_Agl

        ){



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
				LikDVSMC_D_GraphTOF[mc_type]->SetPoint(j,ToFPB.EkBinCent(i),LikD_Correction_TOF -> GetBinContent(i+1,mc_type+1));
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
                                LikDVSMC_D_GraphNaF[mc_type]->SetPoint(j,NaFPB.EkBinCent(i),LikD_Correction_NaF -> GetBinContent(i+1,mc_type+1));
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
                                LikDVSMC_D_GraphAgl[mc_type]->SetPoint(j,AglPB.EkBinCent(i),LikD_Correction_Agl -> GetBinContent(i+1,mc_type+1));
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
                                DistDVSMC_D_GraphTOF[mc_type]->SetPoint(j,ToFPB.EkBinCent(i),DistD_Correction_TOF -> GetBinContent(i+1,mc_type+1));
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
                                DistDVSMC_D_GraphNaF[mc_type]->SetPoint(j,NaFPB.EkBinCent(i),DistD_Correction_NaF -> GetBinContent(i+1,mc_type+1));
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
				DistDVSMC_D_GraphAgl[mc_type]->SetPoint(j,AglPB.EkBinCent(i),DistD_Correction_Agl -> GetBinContent(i+1,mc_type+1));
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

	
	finalPlots.Add(c20_bis);
	finalPlots.Add(c21_bis);
	finalPlots.writeObjsInFolder("DATA-driven Results/Data vs MC/Deutons");

	return;
}
