

void DVSMCPreSeleffD_Plot(   	TH1 * PreSelD_Correction_R  ,
                               	TH1 * PreSelD_Correction_TOF,
                               	TH1 * PreSelD_Correction_NaF,
                               	TH1 * PreSelD_Correction_Agl

        ){


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
					PreSelD_Correction_TOF_Graph[S][mc_type]->SetPoint(j,ToFPB.EkPerMassBinCent(i),PreSelD_Correction_TOF -> GetBinContent(i+1,mc_type+1,S+1));
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
					PreSelD_Correction_NaF_Graph[S][mc_type]->SetPoint(j,NaFPB.EkPerMassBinCent(i),PreSelD_Correction_NaF -> GetBinContent(i+1,mc_type+1,S+1));
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
					PreSelD_Correction_Agl_Graph[S][mc_type]->SetPoint(j,AglPB.EkPerMassBinCent(i),PreSelD_Correction_Agl -> GetBinContent(i+1,mc_type+1,S+1));
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

	for(int S=0;S<3;S++){
		finalPlots.Add(c21_D[S]);
	}
	finalPlots.writeObjsInFolder("DATA-driven Results/Data vs MC/Deutons");

	return;
}
