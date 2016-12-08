
void DVSMCQualeff_Plot(      

			     TH1 *   DistP_CorrectionFit_R      ,
			     TH1 *   DistP_CorrectionFit_TOF,
			     TH1 *   DistP_CorrectionFit_NaF,
			     TH1 *   DistP_CorrectionFit_Agl,

			     TH1 *  LikP_CorrectionFit_R   ,
			     TH1 *  LikP_CorrectionFit_TOF ,
			     TH1 *  LikP_CorrectionFit_NaF ,
			     TH1 *  LikP_CorrectionFit_Agl,

			     TH1 *   DistP_Correction_R      ,
			     TH1 *   DistP_Correction_TOF,
			     TH1 *   DistP_Correction_NaF,
			     TH1 *   DistP_Correction_Agl,

			     TH1 *  LikP_Correction_R   ,
			     TH1 *  LikP_Correction_TOF ,
			     TH1 *  LikP_Correction_NaF ,
			     TH1 *  LikP_Correction_Agl ,
			     
			     TH1 *DistP_MCEff_R  ,
			     TH1 *DistP_MCEff_TOF,
			     TH1 *DistP_MCEff_NaF,
			     TH1 *DistP_MCEff_Agl,

			     TH1 *LikP_MCEff_R   ,
			     TH1 *LikP_MCEff_TOF ,
			     TH1 *LikP_MCEff_NaF ,
			     TH1 *LikP_MCEff_Agl, 

			     TH1 *DistP_DataEff_R  ,
			     TH1 *DistP_DataEff_TOF,
			     TH1 *DistP_DataEff_NaF,
			     TH1 *DistP_DataEff_Agl,

			     TH1 *LikP_DataEff_R   ,
			     TH1 *LikP_DataEff_TOF ,
			     TH1 *LikP_DataEff_NaF ,
			     TH1 *LikP_DataEff_Agl 

				


        ){


	TCanvas *c20=new TCanvas("Likelihood (Data&MC Eff.)");
	TCanvas *c21=new TCanvas("Distance (Data&MC Eff.)");

	TCanvas *c20_bis=new TCanvas("Likelihood (Corrections)");
	TCanvas *c21_bis=new TCanvas("Distance (Corrections)");

	TGraphErrors *LikDVSMC_P_Graph;
	TGraphErrors *DistDVSMC_P_Graph;

	TGraphErrors *LikDVSMCFit_P_Graph;
	TGraphErrors *DistDVSMCFit_P_Graph;

	
	c20->Divide(2,2);

        c20->cd(1);
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
        TGraphErrors * LikDVSMC_P_MC=new TGraphErrors();
        TGraphErrors * LikDVSMC_P_Data=new TGraphErrors();

        int j=0;
        for(int i=0;i<nbinsr;i++) {
                if(LikP_MCEff_R -> GetBinContent(i+1)>0){
                        LikDVSMC_P_MC->SetPoint(j,PRB.RigBinCent(i),LikP_MCEff_R -> GetBinContent(i+1));
                        LikDVSMC_P_MC->SetPointError(j,0,LikP_MCEff_R -> GetBinError(i+1));
                        LikDVSMC_P_Data->SetPoint(j,PRB.RigBinCent(i),LikP_DataEff_R -> GetBinContent(i+1));
                        LikDVSMC_P_Data->SetPointError(j,0,LikP_DataEff_R -> GetBinError(i+1));

                        j++;
                }
        }
        LikDVSMC_P_MC->SetMarkerColor(2);
        LikDVSMC_P_Data->SetMarkerColor(1);
        LikDVSMC_P_MC->SetMarkerStyle(8);
        LikDVSMC_P_Data->SetMarkerStyle(21);
	LikDVSMC_P_MC->SetLineColor(2);
        LikDVSMC_P_Data->SetLineColor(1);
	LikDVSMC_P_MC->GetXaxis()->SetTitle("R[GV]");
        LikDVSMC_P_MC->GetYaxis()->SetTitle("Efficiency ");
        LikDVSMC_P_MC->SetTitle("R binning");
        LikDVSMC_P_MC->GetYaxis()->SetRangeUser(0,1);
	
	LikDVSMC_P_MC->GetXaxis()->SetTitleSize(0.045);
        LikDVSMC_P_MC->GetYaxis()->SetTitleSize(0.045);
        
	LikDVSMC_P_MC->Draw("AP");
        LikDVSMC_P_Data->Draw("Psame");
        {
                TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(LikDVSMC_P_MC,"Protons MC", "ep");
                leg->AddEntry(LikDVSMC_P_Data,"ISS Data", "ep");
                leg->SetLineWidth(2);
		leg->Draw("same");
        }



	c20->cd(2);
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
        TGraphErrors * LikDVSMCTOF_P_MC=new TGraphErrors();
        TGraphErrors * LikDVSMCTOF_P_Data=new TGraphErrors();

        j=0;
        for(int i=0;i<nbinsToF;i++) {
                if(LikP_MCEff_TOF -> GetBinContent(i+1)>0){
                        LikDVSMCTOF_P_MC->SetPoint(j,ToFPB.EkPerMassBinCent(i),LikP_MCEff_TOF -> GetBinContent(i+1));
                        LikDVSMCTOF_P_MC->SetPointError(j,0,LikP_MCEff_TOF -> GetBinError(i+1));
                        LikDVSMCTOF_P_Data->SetPoint(j,ToFPB.EkPerMassBinCent(i),LikP_DataEff_TOF -> GetBinContent(i+1));
                        LikDVSMCTOF_P_Data->SetPointError(j,0,LikP_DataEff_TOF -> GetBinError(i+1));

                        j++;
                }
        }
        LikDVSMCTOF_P_MC->SetMarkerColor(2);
        LikDVSMCTOF_P_Data->SetMarkerColor(1);
        LikDVSMCTOF_P_MC->SetMarkerStyle(8);
        LikDVSMCTOF_P_Data->SetMarkerStyle(21);
        LikDVSMCTOF_P_MC->SetLineColor(2);
        LikDVSMCTOF_P_Data->SetLineColor(1);
        LikDVSMCTOF_P_MC->GetXaxis()->SetTitle("Kin. En./nucl. [GeV/nucl.]");
        LikDVSMCTOF_P_MC->GetYaxis()->SetTitle("Efficiency ");
        LikDVSMCTOF_P_MC->SetTitle("TOF range");
        LikDVSMCTOF_P_MC->GetYaxis()->SetRangeUser(0,1);
	
	LikDVSMCTOF_P_MC->GetXaxis()->SetTitleSize(0.045);
        LikDVSMCTOF_P_MC->GetYaxis()->SetTitleSize(0.045);
        
	
        LikDVSMCTOF_P_MC->Draw("AP");
        LikDVSMCTOF_P_Data->Draw("Psame");
        {
                TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(LikDVSMCTOF_P_MC,"Protons MC", "ep");
                leg->AddEntry(LikDVSMCTOF_P_Data,"ISS Data", "ep");
                leg->SetLineWidth(2);
		leg->Draw("same");
        }


	c20->cd(3);
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
        TGraphErrors * LikDVSMCNaF_P_MC=new TGraphErrors();
        TGraphErrors * LikDVSMCNaF_P_Data=new TGraphErrors();

        j=0;
        for(int i=0;i<nbinsNaF;i++) {
                if(LikP_MCEff_NaF -> GetBinContent(i+1)>0){
                        LikDVSMCNaF_P_MC->SetPoint(j,NaFPB.EkPerMassBinCent(i),LikP_MCEff_NaF -> GetBinContent(i+1));
                        LikDVSMCNaF_P_MC->SetPointError(j,0,LikP_MCEff_NaF -> GetBinError(i+1));
                        LikDVSMCNaF_P_Data->SetPoint(j,NaFPB.EkPerMassBinCent(i),LikP_DataEff_NaF -> GetBinContent(i+1));
                        LikDVSMCNaF_P_Data->SetPointError(j,0,LikP_DataEff_NaF -> GetBinError(i+1));

                        j++;
                }
        }
        LikDVSMCNaF_P_MC->SetMarkerColor(2);
        LikDVSMCNaF_P_Data->SetMarkerColor(1);
        LikDVSMCNaF_P_MC->SetMarkerStyle(8);
        LikDVSMCNaF_P_Data->SetMarkerStyle(21);
        LikDVSMCNaF_P_MC->SetLineColor(2);
        LikDVSMCNaF_P_Data->SetLineColor(1);
        LikDVSMCNaF_P_MC->GetXaxis()->SetTitle("Kin. En./nucl. [GeV/nucl.]");
        LikDVSMCNaF_P_MC->GetYaxis()->SetTitle("Efficiency ");
        LikDVSMCNaF_P_MC->SetTitle("NaF range");
        LikDVSMCNaF_P_MC->GetYaxis()->SetRangeUser(0,1);
	
	LikDVSMCNaF_P_MC->GetXaxis()->SetTitleSize(0.045);
        LikDVSMCNaF_P_MC->GetYaxis()->SetTitleSize(0.045);
        
	
        LikDVSMCNaF_P_MC->Draw("AP");
        LikDVSMCNaF_P_Data->Draw("Psame");
        {
                TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(LikDVSMCNaF_P_MC,"Protons MC", "ep");
                leg->AddEntry(LikDVSMCNaF_P_Data,"ISS Data", "ep");
               leg->SetLineWidth(2);
		leg->Draw("same");
        }

	c20->cd(4);
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
        TGraphErrors * LikDVSMCAgl_P_MC=new TGraphErrors();
        TGraphErrors * LikDVSMCAgl_P_Data=new TGraphErrors();

        j=0;
        for(int i=0;i<nbinsAgl;i++) {
                if(LikP_MCEff_Agl -> GetBinContent(i+1)>0){
                        LikDVSMCAgl_P_MC->SetPoint(j,AglPB.EkPerMassBinCent(i),LikP_MCEff_Agl -> GetBinContent(i+1));
                        LikDVSMCAgl_P_MC->SetPointError(j,0,LikP_MCEff_Agl -> GetBinError(i+1));
                        LikDVSMCAgl_P_Data->SetPoint(j,AglPB.EkPerMassBinCent(i),LikP_DataEff_Agl -> GetBinContent(i+1));
                        LikDVSMCAgl_P_Data->SetPointError(j,0,LikP_DataEff_Agl -> GetBinError(i+1));

                        j++;
                }
        }
        LikDVSMCAgl_P_MC->SetMarkerColor(2);
        LikDVSMCAgl_P_Data->SetMarkerColor(1);
        LikDVSMCAgl_P_MC->SetMarkerStyle(8);
        LikDVSMCAgl_P_Data->SetMarkerStyle(21);
        LikDVSMCAgl_P_MC->SetLineColor(2);
        LikDVSMCAgl_P_Data->SetLineColor(1);
        LikDVSMCAgl_P_MC->GetXaxis()->SetTitle("Kin. En./nucl. [GeV/nucl.]");
        LikDVSMCAgl_P_MC->GetYaxis()->SetTitle("Efficiency ");
        LikDVSMCAgl_P_MC->SetTitle("Agl range");
        LikDVSMCAgl_P_MC->GetYaxis()->SetRangeUser(0,1);

	LikDVSMCAgl_P_MC->GetXaxis()->SetTitleSize(0.045);
        LikDVSMCAgl_P_MC->GetYaxis()->SetTitleSize(0.045);

        LikDVSMCAgl_P_MC->Draw("AP");
        LikDVSMCAgl_P_Data->Draw("Psame");
        {
                TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(LikDVSMCAgl_P_MC,"Protons MC", "ep");
                leg->AddEntry(LikDVSMCAgl_P_Data,"ISS Data", "ep");
                leg->SetLineWidth(2);
		leg->Draw("same");
        }



	 c21->Divide(2,2);

        c21->cd(1);
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
        TGraphErrors * DistDVSMC_P_MC=new TGraphErrors();
        TGraphErrors * DistDVSMC_P_Data=new TGraphErrors();

        j=0;
        for(int i=0;i<nbinsr;i++) {
                if(DistP_MCEff_R -> GetBinContent(i+1)>0){
                        DistDVSMC_P_MC->SetPoint(j,PRB.RigBinCent(i),DistP_MCEff_R -> GetBinContent(i+1));
                        DistDVSMC_P_MC->SetPointError(j,0,DistP_MCEff_R -> GetBinError(i+1));
                        DistDVSMC_P_Data->SetPoint(j,PRB.RigBinCent(i),DistP_DataEff_R -> GetBinContent(i+1));
                        DistDVSMC_P_Data->SetPointError(j,0,DistP_DataEff_R -> GetBinError(i+1));

                        j++;
                }
        }
        DistDVSMC_P_MC->SetMarkerColor(2);
        DistDVSMC_P_Data->SetMarkerColor(1);
        DistDVSMC_P_MC->SetMarkerStyle(8);
        DistDVSMC_P_Data->SetMarkerStyle(21);
        DistDVSMC_P_MC->SetLineColor(2);
        DistDVSMC_P_Data->SetLineColor(1);
        DistDVSMC_P_MC->GetXaxis()->SetTitle("R[GV]");
        DistDVSMC_P_MC->GetYaxis()->SetTitle("Efficiency ");
        DistDVSMC_P_MC->SetTitle("R binning");
        DistDVSMC_P_MC->GetYaxis()->SetRangeUser(0,1);

	DistDVSMC_P_MC->GetXaxis()->SetTitleSize(0.045);
        DistDVSMC_P_MC->GetYaxis()->SetTitleSize(0.045);

        DistDVSMC_P_MC->Draw("AP");
        DistDVSMC_P_Data->Draw("Psame");
        {
                TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(DistDVSMC_P_MC,"Protons MC", "ep");
                leg->AddEntry(DistDVSMC_P_Data,"ISS Data", "ep");
                leg->SetLineWidth(2);
		leg->Draw("same");
        }



        c21->cd(2);
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
        TGraphErrors * DistDVSMCTOF_P_MC=new TGraphErrors();
        TGraphErrors * DistDVSMCTOF_P_Data=new TGraphErrors();

        j=0;
        for(int i=0;i<nbinsToF;i++) {
                if(DistP_MCEff_TOF -> GetBinContent(i+1)>0){
                        DistDVSMCTOF_P_MC->SetPoint(j,ToFPB.EkPerMassBinCent(i),DistP_MCEff_TOF -> GetBinContent(i+1));
                        DistDVSMCTOF_P_MC->SetPointError(j,0,DistP_MCEff_TOF -> GetBinError(i+1));
                        DistDVSMCTOF_P_Data->SetPoint(j,ToFPB.EkPerMassBinCent(i),DistP_DataEff_TOF -> GetBinContent(i+1));
                        DistDVSMCTOF_P_Data->SetPointError(j,0,DistP_DataEff_TOF -> GetBinError(i+1));

                        j++;
                }
        }
        DistDVSMCTOF_P_MC->SetMarkerColor(2);
        DistDVSMCTOF_P_Data->SetMarkerColor(1);
        DistDVSMCTOF_P_MC->SetMarkerStyle(8);
        DistDVSMCTOF_P_Data->SetMarkerStyle(21);
        DistDVSMCTOF_P_MC->SetLineColor(2);
        DistDVSMCTOF_P_Data->SetLineColor(1);
        DistDVSMCTOF_P_MC->GetXaxis()->SetTitle("n. En./nucl. [GeV/nucl.]");
        DistDVSMCTOF_P_MC->GetYaxis()->SetTitle("Efficiency ");
        DistDVSMCTOF_P_MC->SetTitle("TOF range");
        DistDVSMCTOF_P_MC->GetYaxis()->SetRangeUser(0,1);

	DistDVSMCTOF_P_MC->GetXaxis()->SetTitleSize(0.045);
        DistDVSMCTOF_P_MC->GetYaxis()->SetTitleSize(0.045);

        DistDVSMCTOF_P_MC->Draw("AP");
        DistDVSMCTOF_P_Data->Draw("Psame");
        {
                TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(DistDVSMCTOF_P_MC,"Protons MC", "ep");
                leg->AddEntry(DistDVSMCTOF_P_Data,"ISS Data", "ep");
                leg->SetLineWidth(2);
		leg->Draw("same");
        }


        c21->cd(3);
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
        TGraphErrors * DistDVSMCNaF_P_MC=new TGraphErrors();
        TGraphErrors * DistDVSMCNaF_P_Data=new TGraphErrors();

        j=0;
        for(int i=0;i<nbinsNaF;i++) {
                if(DistP_MCEff_NaF -> GetBinContent(i+1)>0){
                        DistDVSMCNaF_P_MC->SetPoint(j,NaFPB.EkPerMassBinCent(i),DistP_MCEff_NaF -> GetBinContent(i+1));
                        DistDVSMCNaF_P_MC->SetPointError(j,0,DistP_MCEff_NaF -> GetBinError(i+1));
                        DistDVSMCNaF_P_Data->SetPoint(j,NaFPB.EkPerMassBinCent(i),DistP_DataEff_NaF -> GetBinContent(i+1));
                        DistDVSMCNaF_P_Data->SetPointError(j,0,DistP_DataEff_NaF -> GetBinError(i+1));

                        j++;
                }
        }
        DistDVSMCNaF_P_MC->SetMarkerColor(2);
        DistDVSMCNaF_P_Data->SetMarkerColor(1);
        DistDVSMCNaF_P_MC->SetMarkerStyle(8);
        DistDVSMCNaF_P_Data->SetMarkerStyle(21);
        DistDVSMCNaF_P_MC->SetLineColor(2);
        DistDVSMCNaF_P_Data->SetLineColor(1);
        DistDVSMCNaF_P_MC->GetXaxis()->SetTitle("Kin. En./nucl. [GeV/nucl.]");
        DistDVSMCNaF_P_MC->GetYaxis()->SetTitle("Efficiency ");
        DistDVSMCNaF_P_MC->SetTitle("NaF range");
        DistDVSMCNaF_P_MC->GetYaxis()->SetRangeUser(0,1);

	DistDVSMCNaF_P_MC->GetXaxis()->SetTitleSize(0.045);
        DistDVSMCNaF_P_MC->GetYaxis()->SetTitleSize(0.045);

        DistDVSMCNaF_P_MC->Draw("AP");
        DistDVSMCNaF_P_Data->Draw("Psame");
        {
                TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(DistDVSMCNaF_P_MC,"Protons MC", "ep");
                leg->AddEntry(DistDVSMCNaF_P_Data,"ISS Data", "ep");
                leg->SetLineWidth(2);
		leg->Draw("same");
        }

        c21->cd(4);
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
        TGraphErrors * DistDVSMCAgl_P_MC=new TGraphErrors();
        TGraphErrors * DistDVSMCAgl_P_Data=new TGraphErrors();

        j=0;
        for(int i=0;i<nbinsAgl;i++) {
                if(DistP_MCEff_Agl -> GetBinContent(i+1)>0){
                        DistDVSMCAgl_P_MC->SetPoint(j,AglPB.EkPerMassBinCent(i),DistP_MCEff_Agl -> GetBinContent(i+1));
                        DistDVSMCAgl_P_MC->SetPointError(j,0,DistP_MCEff_Agl -> GetBinError(i+1));
                        DistDVSMCAgl_P_Data->SetPoint(j,AglPB.EkPerMassBinCent(i),DistP_DataEff_Agl -> GetBinContent(i+1));
                        DistDVSMCAgl_P_Data->SetPointError(j,0,DistP_DataEff_Agl -> GetBinError(i+1));

                        j++;
                }
        }
        DistDVSMCAgl_P_MC->SetMarkerColor(2);
        DistDVSMCAgl_P_Data->SetMarkerColor(1);
        DistDVSMCAgl_P_MC->SetMarkerStyle(8);
        DistDVSMCAgl_P_Data->SetMarkerStyle(21);
        DistDVSMCAgl_P_MC->SetLineColor(2);
        DistDVSMCAgl_P_Data->SetLineColor(1);
        DistDVSMCAgl_P_MC->GetXaxis()->SetTitle("Kin. En./nucl. [GeV/nucl.]");
        DistDVSMCAgl_P_MC->GetYaxis()->SetTitle("Efficiency ");
        DistDVSMCAgl_P_MC->SetTitle("Agl range");
        DistDVSMCAgl_P_MC->GetYaxis()->SetRangeUser(0,1);

	DistDVSMCAgl_P_MC->GetXaxis()->SetTitleSize(0.045);
        DistDVSMCAgl_P_MC->GetYaxis()->SetTitleSize(0.045);

        DistDVSMCAgl_P_MC->Draw("AP");
        DistDVSMCAgl_P_Data->Draw("Psame");
        {
                TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(DistDVSMCAgl_P_MC,"Protons MC", "ep");
                leg->AddEntry(DistDVSMCAgl_P_Data,"ISS Data", "ep");
                leg->SetLineWidth(2);
		leg->Draw("same");
        }


	c20_bis->Divide(2,2);	

	c20_bis->cd(1);
	gPad->SetLogx();
	gPad->SetGridx();
	gPad->SetGridy();
	LikDVSMC_P_Graph=new TGraphErrors();
	LikDVSMC_P_Graph->SetName("LikDVSMCFit_P_Graph");
	LikDVSMCFit_P_Graph=new TGraphErrors();
	LikDVSMCFit_P_Graph->SetName("LikDVSMCFit_P_Graph");
	
	j=0;
	for(int i=0;i<nbinsr;i++) {
		if(LikP_CorrectionFit_R -> GetBinContent(i+1)>0){
			LikDVSMC_P_Graph->SetPoint(j,PRB.RigBinCent(i),LikP_Correction_R -> GetBinContent(i+1));
			LikDVSMC_P_Graph->SetPointError(j,0,LikP_Correction_R -> GetBinError(i+1));
			LikDVSMCFit_P_Graph->SetPoint(j,PRB.RigBinCent(i),LikP_CorrectionFit_R -> GetBinContent(i+1));
			LikDVSMCFit_P_Graph->SetPointError(j,0,LikP_CorrectionFit_R -> GetBinError(i+1));

			j++;
		}
	}
	LikDVSMC_P_Graph->SetLineColor(2);
	LikDVSMCFit_P_Graph->SetFillColor(2);
	LikDVSMCFit_P_Graph->SetFillStyle(3001);
	LikDVSMCFit_P_Graph->SetLineWidth(4);
	LikDVSMC_P_Graph->SetMarkerColor(2);
	LikDVSMC_P_Graph->SetMarkerStyle(8);
	LikDVSMC_P_Graph->GetXaxis()->SetTitle("R[GV]");
	LikDVSMC_P_Graph->GetYaxis()->SetTitle("Efficiency (Data/MC)");
	LikDVSMC_P_Graph->SetTitle("R binning");
	LikDVSMC_P_Graph->GetYaxis()->SetRangeUser(0.9,1.2);

	LikDVSMC_P_Graph->GetXaxis()->SetTitleSize(0.045);
        LikDVSMC_P_Graph->GetYaxis()->SetTitleSize(0.045);

	LikDVSMC_P_Graph->Draw("AP");
	LikDVSMCFit_P_Graph->Draw("P4Csame");
	{
		TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
		leg->AddEntry(LikDVSMC_P_Graph,"Efficiency correction", "ep");
		leg->AddEntry(LikDVSMCFit_P_Graph,"Param.", "l");
		leg->SetLineWidth(2);
		leg->Draw("same");
	}
	



	c20_bis->cd(2);
	gPad->SetGridx();
        gPad->SetGridy();
        TGraphErrors * LikDVSMC_P_GraphTOF=new TGraphErrors();
        LikDVSMC_P_GraphTOF->SetName("LikDVSMCFit_P_GraphTOF");
        TGraphErrors *  LikDVSMCFit_P_GraphTOF=new TGraphErrors();
        LikDVSMCFit_P_GraphTOF->SetName("LikDVSMCFit_P_GraphTOF");

        j=0;
        for(int i=0;i<nbinsToF;i++) {
                if(LikP_CorrectionFit_TOF -> GetBinContent(i+1)>0){
                        LikDVSMC_P_GraphTOF->SetPoint(j,ToFPB.EkPerMassBinCent(i),LikP_Correction_TOF -> GetBinContent(i+1));
                        LikDVSMC_P_GraphTOF->SetPointError(j,0,LikP_Correction_TOF -> GetBinError(i+1));
                        LikDVSMCFit_P_GraphTOF->SetPoint(j,ToFPB.EkPerMassBinCent(i),LikP_CorrectionFit_TOF -> GetBinContent(i+1));
                        LikDVSMCFit_P_GraphTOF->SetPointError(j,0,LikP_CorrectionFit_TOF -> GetBinError(i+1));

                        j++;
                }
        }
        LikDVSMC_P_GraphTOF->SetLineColor(2);
        LikDVSMCFit_P_GraphTOF->SetFillColor(2);
        LikDVSMCFit_P_GraphTOF->SetFillStyle(3001);
        LikDVSMCFit_P_GraphTOF->SetLineWidth(4);
        LikDVSMC_P_GraphTOF->SetMarkerColor(2);
        LikDVSMC_P_GraphTOF->SetMarkerStyle(8);
        LikDVSMC_P_GraphTOF->GetXaxis()->SetTitle("Kin. En / nucl. [GeV/nucl.]");
        LikDVSMC_P_GraphTOF->GetYaxis()->SetTitle("Efficiency (Data/MC)");
        LikDVSMC_P_GraphTOF->GetYaxis()->SetRangeUser(0.9,1.2);
	LikDVSMC_P_GraphTOF->SetTitle("TOF range");

	LikDVSMC_P_GraphTOF->GetXaxis()->SetTitleSize(0.045);
        LikDVSMC_P_GraphTOF->GetYaxis()->SetTitleSize(0.045);

	LikDVSMC_P_GraphTOF->Draw("AP");
        LikDVSMCFit_P_GraphTOF->Draw("P4Csame");

	{
                TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(LikDVSMC_P_GraphTOF,"Efficiency correction", "ep");
                leg->AddEntry(LikDVSMCFit_P_GraphTOF,"Param.", "l");
                leg->SetLineWidth(2);
		leg->Draw("same");
        }



	c20_bis->cd(3);
	gPad->SetGridx();
        gPad->SetGridy();
        TGraphErrors * LikDVSMC_P_GraphNaF=new TGraphErrors();
        LikDVSMC_P_GraphNaF->SetName("LikDVSMCFit_P_GraphNaF");
        TGraphErrors *  LikDVSMCFit_P_GraphNaF=new TGraphErrors();
        LikDVSMCFit_P_GraphNaF->SetName("LikDVSMCFit_P_GraphNaF");

        j=0;
        for(int i=0;i<nbinsNaF;i++) {
                if(LikP_CorrectionFit_NaF -> GetBinContent(i+1)>0){
                        LikDVSMC_P_GraphNaF->SetPoint(j,NaFPB.EkPerMassBinCent(i),LikP_Correction_NaF -> GetBinContent(i+1));
                        LikDVSMC_P_GraphNaF->SetPointError(j,0,LikP_Correction_NaF -> GetBinError(i+1));
                        LikDVSMCFit_P_GraphNaF->SetPoint(j,NaFPB.EkPerMassBinCent(i),LikP_CorrectionFit_NaF -> GetBinContent(i+1));
                        LikDVSMCFit_P_GraphNaF->SetPointError(j,0,LikP_CorrectionFit_NaF -> GetBinError(i+1));

                        j++;
                }
        }
        LikDVSMC_P_GraphNaF->SetLineColor(2);
        LikDVSMCFit_P_GraphNaF->SetFillColor(2);
        LikDVSMCFit_P_GraphNaF->SetFillStyle(3001);
        LikDVSMCFit_P_GraphNaF->SetLineWidth(4);
        LikDVSMC_P_GraphNaF->SetMarkerColor(2);
        LikDVSMC_P_GraphNaF->SetMarkerStyle(8);
        LikDVSMC_P_GraphNaF->GetXaxis()->SetTitle("Kin. En / nucl. [GeV/nucl.]");
        LikDVSMC_P_GraphNaF->GetYaxis()->SetTitle("Efficiency (Data/MC)");
        LikDVSMC_P_GraphNaF->GetYaxis()->SetRangeUser(0.9,1.2);
	LikDVSMC_P_GraphNaF->SetTitle("NaF range");

	LikDVSMC_P_GraphNaF->GetXaxis()->SetTitleSize(0.045);
        LikDVSMC_P_GraphNaF->GetYaxis()->SetTitleSize(0.045);
	
	LikDVSMC_P_GraphNaF->Draw("AP");
        LikDVSMCFit_P_GraphNaF->Draw("P4Csame");
	{
                TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(LikDVSMC_P_GraphNaF,"Efficiency correction", "ep");
                leg->AddEntry(LikDVSMCFit_P_GraphNaF,"Param.", "l");
                leg->SetLineWidth(2);
		leg->Draw("same");
        }





	c20_bis->cd(4);
	gPad->SetGridx();
        gPad->SetGridy();
        TGraphErrors * LikDVSMC_P_GraphAgl=new TGraphErrors();
        LikDVSMC_P_GraphAgl->SetName("LikDVSMCFit_P_GraphAgl");
        TGraphErrors *  LikDVSMCFit_P_GraphAgl=new TGraphErrors();
        LikDVSMCFit_P_GraphAgl->SetName("LikDVSMCFit_P_GraphAgl");

        j=0;
        for(int i=0;i<nbinsToF;i++) {
                if(LikP_CorrectionFit_Agl -> GetBinContent(i+1)>0){
                        LikDVSMC_P_GraphAgl->SetPoint(j,AglPB.EkPerMassBinCent(i),LikP_Correction_Agl -> GetBinContent(i+1));
                        LikDVSMC_P_GraphAgl->SetPointError(j,0,LikP_Correction_Agl -> GetBinError(i+1));
                        LikDVSMCFit_P_GraphAgl->SetPoint(j,AglPB.EkPerMassBinCent(i),LikP_CorrectionFit_Agl -> GetBinContent(i+1));
                        LikDVSMCFit_P_GraphAgl->SetPointError(j,0,LikP_CorrectionFit_Agl -> GetBinError(i+1));

                        j++;
                }
        }
        LikDVSMC_P_GraphAgl->SetLineColor(2);
        LikDVSMCFit_P_GraphAgl->SetFillColor(2);
        LikDVSMCFit_P_GraphAgl->SetFillStyle(3001);
        LikDVSMCFit_P_GraphAgl->SetLineWidth(4);
        LikDVSMC_P_GraphAgl->SetMarkerColor(2);
        LikDVSMC_P_GraphAgl->SetMarkerStyle(8);
        LikDVSMC_P_GraphAgl->GetXaxis()->SetTitle("Kin. En / nucl. [GeV/nucl.]");
        LikDVSMC_P_GraphAgl->GetYaxis()->SetTitle("Efficiency (Data/MC)");
        LikDVSMC_P_GraphAgl->GetYaxis()->SetRangeUser(0.9,1.2);
	LikDVSMC_P_GraphAgl->SetTitle("Agl range");

	LikDVSMC_P_GraphAgl->GetXaxis()->SetTitleSize(0.045);
        LikDVSMC_P_GraphAgl->GetYaxis()->SetTitleSize(0.045);

	LikDVSMC_P_GraphAgl->Draw("AP");
        LikDVSMCFit_P_GraphAgl->Draw("P4Csame");
	{
                TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(LikDVSMC_P_GraphAgl,"Efficiency correction", "ep");
                leg->AddEntry(LikDVSMCFit_P_GraphAgl,"Param.", "l");
                leg->SetLineWidth(2);
		leg->Draw("same");
        }

	 c21_bis->Divide(2,2);

        c21_bis->cd(1);
        gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
        DistDVSMC_P_Graph=new TGraphErrors();
        DistDVSMC_P_Graph->SetName("DistDVSMCFit_P_Graph");
        DistDVSMCFit_P_Graph=new TGraphErrors();
        DistDVSMCFit_P_Graph->SetName("DistDVSMCFit_P_Graph");

        j=0;
        for(int i=0;i<nbinsr;i++) {
                if(DistP_CorrectionFit_R -> GetBinContent(i+1)>0){
                        DistDVSMC_P_Graph->SetPoint(j,PRB.RigBinCent(i),DistP_Correction_R -> GetBinContent(i+1));
                        DistDVSMC_P_Graph->SetPointError(j,0,DistP_Correction_R -> GetBinError(i+1));
                        DistDVSMCFit_P_Graph->SetPoint(j,PRB.RigBinCent(i),DistP_CorrectionFit_R -> GetBinContent(i+1));
                        DistDVSMCFit_P_Graph->SetPointError(j,0,DistP_CorrectionFit_R -> GetBinError(i+1));

                        j++;
                }
        }
        DistDVSMC_P_Graph->SetLineColor(2);
        DistDVSMCFit_P_Graph->SetFillColor(2);
        DistDVSMCFit_P_Graph->SetFillStyle(3001);
        DistDVSMCFit_P_Graph->SetLineWidth(4);
        DistDVSMC_P_Graph->SetMarkerColor(2);
        DistDVSMC_P_Graph->SetMarkerStyle(8);
        DistDVSMC_P_Graph->GetXaxis()->SetTitle("R[GV]");
        DistDVSMC_P_Graph->GetYaxis()->SetTitle("Efficiency (Data/MC)");
        DistDVSMC_P_Graph->SetTitle("R binning");
        DistDVSMC_P_Graph->GetYaxis()->SetRangeUser(0.9,1.2);

	DistDVSMC_P_Graph->GetXaxis()->SetTitleSize(0.045);
        DistDVSMC_P_Graph->GetYaxis()->SetTitleSize(0.045);

        DistDVSMC_P_Graph->Draw("AP");
        DistDVSMCFit_P_Graph->Draw("P4Csame");
        {
                TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(DistDVSMC_P_Graph,"Efficiency correction", "ep");
                leg->AddEntry(DistDVSMCFit_P_Graph,"Param.", "l");
                leg->SetLineWidth(2);
		leg->Draw("same");
        }




        c21_bis->cd(2);
        gPad->SetGridx();
        gPad->SetGridy();
        TGraphErrors * DistDVSMC_P_GraphTOF=new TGraphErrors();
        DistDVSMC_P_GraphTOF->SetName("DistDVSMCFit_P_GraphTOF");
        TGraphErrors *  DistDVSMCFit_P_GraphTOF=new TGraphErrors();
        DistDVSMCFit_P_GraphTOF->SetName("DistDVSMCFit_P_GraphTOF");

        j=0;
        for(int i=0;i<nbinsToF;i++) {
                if(DistP_CorrectionFit_TOF -> GetBinContent(i+1)>0){
                        DistDVSMC_P_GraphTOF->SetPoint(j,ToFPB.EkPerMassBinCent(i),DistP_Correction_TOF -> GetBinContent(i+1));
                        DistDVSMC_P_GraphTOF->SetPointError(j,0,DistP_Correction_TOF -> GetBinError(i+1));
                        DistDVSMCFit_P_GraphTOF->SetPoint(j,ToFPB.EkPerMassBinCent(i),DistP_CorrectionFit_TOF -> GetBinContent(i+1));
                        DistDVSMCFit_P_GraphTOF->SetPointError(j,0,DistP_CorrectionFit_TOF -> GetBinError(i+1));

                        j++;
                }
        }
        DistDVSMC_P_GraphTOF->SetLineColor(2);
        DistDVSMCFit_P_GraphTOF->SetFillColor(2);
        DistDVSMCFit_P_GraphTOF->SetFillStyle(3001);
        DistDVSMCFit_P_GraphTOF->SetLineWidth(4);
        DistDVSMC_P_GraphTOF->SetMarkerColor(2);
        DistDVSMC_P_GraphTOF->SetMarkerStyle(8);
        DistDVSMC_P_GraphTOF->GetXaxis()->SetTitle("Kin. En / nucl. [GeV/nucl.]");
        DistDVSMC_P_GraphTOF->GetYaxis()->SetTitle("Efficiency (Data/MC)");
        DistDVSMC_P_GraphTOF->GetYaxis()->SetRangeUser(0.9,1.2);
        DistDVSMC_P_GraphTOF->SetTitle("TOF range");

	DistDVSMC_P_GraphTOF->GetXaxis()->SetTitleSize(0.045);
        DistDVSMC_P_GraphTOF->GetYaxis()->SetTitleSize(0.045);

        DistDVSMC_P_GraphTOF->Draw("AP");
        DistDVSMCFit_P_GraphTOF->Draw("P4Csame");

        {
                TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(DistDVSMC_P_GraphTOF,"Efficiency correction", "ep");
                leg->AddEntry(DistDVSMCFit_P_GraphTOF,"Param.", "l");
                leg->SetLineWidth(2);
		leg->Draw("same");
        }



        c21_bis->cd(3);
        gPad->SetGridx();
        gPad->SetGridy();
        TGraphErrors * DistDVSMC_P_GraphNaF=new TGraphErrors();
        DistDVSMC_P_GraphNaF->SetName("DistDVSMCFit_P_GraphNaF");
        TGraphErrors *  DistDVSMCFit_P_GraphNaF=new TGraphErrors();
        DistDVSMCFit_P_GraphNaF->SetName("DistDVSMCFit_P_GraphNaF");

        j=0;
        for(int i=0;i<nbinsNaF;i++) {
                if(DistP_CorrectionFit_NaF -> GetBinContent(i+1)>0){
                        DistDVSMC_P_GraphNaF->SetPoint(j,NaFPB.EkPerMassBinCent(i),DistP_Correction_NaF -> GetBinContent(i+1));
                        DistDVSMC_P_GraphNaF->SetPointError(j,0,DistP_Correction_NaF -> GetBinError(i+1));
                        DistDVSMCFit_P_GraphNaF->SetPoint(j,NaFPB.EkPerMassBinCent(i),DistP_CorrectionFit_NaF -> GetBinContent(i+1));
                        DistDVSMCFit_P_GraphNaF->SetPointError(j,0,DistP_CorrectionFit_NaF -> GetBinError(i+1));

                        j++;
                }
        }
        DistDVSMC_P_GraphNaF->SetLineColor(2);
        DistDVSMCFit_P_GraphNaF->SetFillColor(2);
        DistDVSMCFit_P_GraphNaF->SetFillStyle(3001);
        DistDVSMCFit_P_GraphNaF->SetLineWidth(4);
        DistDVSMC_P_GraphNaF->SetMarkerColor(2);
        DistDVSMC_P_GraphNaF->SetMarkerStyle(8);
        DistDVSMC_P_GraphNaF->GetXaxis()->SetTitle("Kin. En / nucl. [GeV/nucl.]");
        DistDVSMC_P_GraphNaF->GetYaxis()->SetTitle("Efficiency (Data/MC)");
        DistDVSMC_P_GraphNaF->GetYaxis()->SetRangeUser(0.9,1.2);
        DistDVSMC_P_GraphNaF->SetTitle("NaF range");

	DistDVSMC_P_GraphNaF->GetXaxis()->SetTitleSize(0.045);
        DistDVSMC_P_GraphNaF->GetYaxis()->SetTitleSize(0.045);

        DistDVSMC_P_GraphNaF->Draw("AP");
        DistDVSMCFit_P_GraphNaF->Draw("P4Csame");
        {
                TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(DistDVSMC_P_GraphNaF,"Efficiency correction", "ep");
                leg->AddEntry(DistDVSMCFit_P_GraphNaF,"Param.", "l");
                leg->SetLineWidth(2);
		leg->Draw("same");
        }





        c21_bis->cd(4);
        gPad->SetGridx();
        gPad->SetGridy();
        TGraphErrors * DistDVSMC_P_GraphAgl=new TGraphErrors();
        DistDVSMC_P_GraphAgl->SetName("DistDVSMCFit_P_GraphAgl");
        TGraphErrors *  DistDVSMCFit_P_GraphAgl=new TGraphErrors();
        DistDVSMCFit_P_GraphAgl->SetName("DistDVSMCFit_P_GraphAgl");

        j=0;
        for(int i=0;i<nbinsToF;i++) {
                if(DistP_CorrectionFit_Agl -> GetBinContent(i+1)>0){
                        DistDVSMC_P_GraphAgl->SetPoint(j,AglPB.EkPerMassBinCent(i),DistP_Correction_Agl -> GetBinContent(i+1));
                        DistDVSMC_P_GraphAgl->SetPointError(j,0,DistP_Correction_Agl -> GetBinError(i+1));
                        DistDVSMCFit_P_GraphAgl->SetPoint(j,AglPB.EkPerMassBinCent(i),DistP_CorrectionFit_Agl -> GetBinContent(i+1));
                        DistDVSMCFit_P_GraphAgl->SetPointError(j,0,DistP_CorrectionFit_Agl -> GetBinError(i+1));

                        j++;
                }
        }
        DistDVSMC_P_GraphAgl->SetLineColor(2);
        DistDVSMCFit_P_GraphAgl->SetFillColor(2);
        DistDVSMCFit_P_GraphAgl->SetFillStyle(3001);
        DistDVSMCFit_P_GraphAgl->SetLineWidth(4);
        DistDVSMC_P_GraphAgl->SetMarkerColor(2);
        DistDVSMC_P_GraphAgl->SetMarkerStyle(8);
        DistDVSMC_P_GraphAgl->GetXaxis()->SetTitle("Kin. En / nucl. [GeV/nucl.]");
        DistDVSMC_P_GraphAgl->GetYaxis()->SetTitle("Efficiency (Data/MC)");
        DistDVSMC_P_GraphAgl->GetYaxis()->SetRangeUser(0.9,1.2);
        DistDVSMC_P_GraphAgl->SetTitle("Agl range");

	DistDVSMC_P_GraphAgl->GetXaxis()->SetTitleSize(0.045);
        DistDVSMC_P_GraphAgl->GetYaxis()->SetTitleSize(0.045);

        DistDVSMC_P_GraphAgl->Draw("AP");
        DistDVSMCFit_P_GraphAgl->Draw("P4Csame");
        {
                TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
                leg->AddEntry(DistDVSMC_P_GraphAgl,"Efficiency correction", "ep");
                leg->AddEntry(DistDVSMCFit_P_GraphAgl,"Param.", "l");
                leg->SetLineWidth(2);
		leg->Draw("same");
        }



			



	finalPlots.Add(c20);
	finalPlots.Add(c21);
	finalPlots.writeObjsInFolder("DATA-driven Results/Data vs MC/Selections/Eff.");

	finalPlots.Add(c20_bis);
	finalPlots.Add(c21_bis);
	finalPlots.writeObjsInFolder("DATA-driven Results/Data vs MC/Selections/Corrections");

        
	finalPlots.Add(DistDVSMCFit_P_Graph   );
	finalPlots.Add(DistDVSMCFit_P_GraphTOF);
	finalPlots.Add(DistDVSMCFit_P_GraphNaF);
	finalPlots.Add(DistDVSMCFit_P_GraphAgl);
	finalPlots.Add(LikDVSMCFit_P_Graph    );
	finalPlots.Add(LikDVSMCFit_P_GraphTOF );
	finalPlots.Add(LikDVSMCFit_P_GraphNaF );
	finalPlots.Add(LikDVSMCFit_P_GraphAgl );
	finalPlots.writeObjsInFolder("Export/DvsMC");

	return;
}
