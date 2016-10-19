void MCQcheck_Plot(
		 	TGraph *Pbckgnd_ROC_QTOF,
			TGraph *Pbckgnd_ROC_QNaF,   
			TGraph *Pbckgnd_ROC_QAgl,  
			TGraph *Pbckgnd_ROC_DistTOF,
			TGraph *Pbckgnd_ROC_DistNaF,   
			TGraph *Pbckgnd_ROC_DistAgl,  
			TGraph *Hebckgnd_ROC_QTOF,    				
			TGraph *Hebckgnd_ROC_QNaF,   
	                TGraph *Hebckgnd_ROC_QAgl,
                        TGraph *Hebckgnd_ROC_DistTOF,
                        TGraph *Hebckgnd_ROC_DistNaF,
                        TGraph *Hebckgnd_ROC_DistAgl, 
			TH3F * QdistrTOF,
			TH3F * QdistrNaF,
			TH3F * QdistrAgl,
			TH3F * QdistrTOFD,
			TH3F * QdistrNaFD,
			TH3F * QdistrAglD,
			TH3F * QdistrTOFHe,
			TH3F * QdistrNaFHe,
			TH3F * QdistrAglHe


		     )	
{


	TH2F * QUtofvsLTOF_TOF  =  (TH2F *) QdistrTOF ->Project3D("zx");
	TH2F * QUtofvsTrack_TOF =  (TH2F *) QdistrTOF ->Project3D("yx");
	TH2F * QLtofvsTrack_TOF =  (TH2F *) QdistrTOF ->Project3D("yz");
	TH2F * QUtofvsLTOF_NaF  =  (TH2F *) QdistrNaF ->Project3D("zx");
	TH2F * QUtofvsTrack_NaF =  (TH2F *) QdistrNaF ->Project3D("yx");
	TH2F * QLtofvsTrack_NaF =  (TH2F *) QdistrNaF ->Project3D("yz");
	TH2F * QUtofvsLTOF_Agl  =  (TH2F *) QdistrAgl ->Project3D("zx");
	TH2F * QUtofvsTrack_Agl =  (TH2F *) QdistrAgl ->Project3D("yx");
	TH2F * QLtofvsTrack_Agl =  (TH2F *) QdistrAgl ->Project3D("yz");

	TH2F * QUtofvsLTOFD_TOF  =  (TH2F *) QdistrTOFD ->Project3D("zx");
	TH2F * QUtofvsTrackD_TOF =  (TH2F *) QdistrTOFD ->Project3D("yx");
	TH2F * QLtofvsTrackD_TOF =  (TH2F *) QdistrTOFD ->Project3D("yz");
	TH2F * QUtofvsLTOFD_NaF  =  (TH2F *) QdistrNaFD ->Project3D("zx");
	TH2F * QUtofvsTrackD_NaF =  (TH2F *) QdistrNaFD ->Project3D("yx");
	TH2F * QLtofvsTrackD_NaF =  (TH2F *) QdistrNaFD ->Project3D("yz");
	TH2F * QUtofvsLTOFD_Agl  =  (TH2F *) QdistrAglD ->Project3D("zx");
	TH2F * QUtofvsTrackD_Agl =  (TH2F *) QdistrAglD ->Project3D("yx");
	TH2F * QLtofvsTrackD_Agl =  (TH2F *) QdistrAglD ->Project3D("yz");

	TH2F * QUtofvsLTOFHe_TOF  =  (TH2F *) QdistrTOFHe ->Project3D("zx");
	TH2F * QUtofvsTrackHe_TOF =  (TH2F *) QdistrTOFHe ->Project3D("yx");
	TH2F * QLtofvsTrackHe_TOF =  (TH2F *) QdistrTOFHe ->Project3D("yz");
	TH2F * QUtofvsLTOFHe_NaF  =  (TH2F *) QdistrNaFHe ->Project3D("zx");
	TH2F * QUtofvsTrackHe_NaF =  (TH2F *) QdistrNaFHe ->Project3D("yx");
	TH2F * QLtofvsTrackHe_NaF =  (TH2F *) QdistrNaFHe ->Project3D("yz");
	TH2F * QUtofvsLTOFHe_Agl  =  (TH2F *) QdistrAglHe ->Project3D("zx");
	TH2F * QUtofvsTrackHe_Agl =  (TH2F *) QdistrAglHe ->Project3D("yx");
	TH2F * QLtofvsTrackHe_Agl =  (TH2F *) QdistrAglHe ->Project3D("yz");






	TCanvas *c1 = new TCanvas("Bad Protons rejection");
	TCanvas *d  = new TCanvas("Upper TOF vs Lower TOF Q");
	TCanvas *d1 = new TCanvas("Upper TOF vs Inner Tracker Q");
	TCanvas *d2 = new TCanvas("Lower TOF vs Inner Tracker Q");		
	TCanvas *e  = new TCanvas("Upper TOF Q");
	TCanvas *e2 = new TCanvas("Inner Tracker Q");
	TCanvas *e1 = new TCanvas("Lower TOF Q");		






	c1->Divide(3,1);
	c1->cd(1);
	gPad->SetGridx();
   		gPad->SetGridy();
		Pbckgnd_ROC_QTOF->SetMarkerColor(2);
		Pbckgnd_ROC_QTOF->SetMarkerStyle(8);
		Pbckgnd_ROC_QTOF->SetLineColor(2);
		Pbckgnd_ROC_DistTOF->SetMarkerColor(4);
		Pbckgnd_ROC_DistTOF->SetMarkerStyle(8);
		Pbckgnd_ROC_DistTOF->SetLineColor(4);
			
		Pbckgnd_ROC_QTOF->GetXaxis()->SetTitle("Efficiency on signal");
		Pbckgnd_ROC_QTOF->GetYaxis()->SetTitle("Background rejection");
		Pbckgnd_ROC_QTOF->GetYaxis()->SetRangeUser(0,1);
		Pbckgnd_ROC_QTOF->Draw("APC");
		Pbckgnd_ROC_DistTOF->Draw("PCsame");

                c1->cd(2);
                gPad->SetGridx();
                gPad->SetGridy();
                Pbckgnd_ROC_QNaF->SetMarkerColor(2);
                Pbckgnd_ROC_QNaF->SetMarkerStyle(8);
                Pbckgnd_ROC_QNaF->SetLineColor(2);  
               	Pbckgnd_ROC_DistNaF->SetMarkerColor(4);
		Pbckgnd_ROC_DistNaF->SetMarkerStyle(8);
		Pbckgnd_ROC_DistNaF->SetLineColor(4);
		
		Pbckgnd_ROC_QNaF->GetXaxis()->SetTitle("Efficiency on signal");
                Pbckgnd_ROC_QNaF->GetYaxis()->SetTitle("Background rejection");
                Pbckgnd_ROC_QNaF->GetYaxis()->SetRangeUser(0,1);
		Pbckgnd_ROC_QNaF->Draw("APC");
		Pbckgnd_ROC_DistNaF->Draw("PCsame");

                c1->cd(3);
                gPad->SetGridx();
                gPad->SetGridy();
                Pbckgnd_ROC_QAgl->SetMarkerColor(2);
                Pbckgnd_ROC_QAgl->SetMarkerStyle(8);
                Pbckgnd_ROC_QAgl->SetLineColor(2);
               	Pbckgnd_ROC_DistAgl->SetMarkerColor(4);
		Pbckgnd_ROC_DistAgl->SetMarkerStyle(8);
		Pbckgnd_ROC_DistAgl->SetLineColor(4);

		Pbckgnd_ROC_QAgl->GetXaxis()->SetTitle("Efficiency on signal");
                Pbckgnd_ROC_QAgl->GetYaxis()->SetTitle("Background rejection");
                Pbckgnd_ROC_QAgl->GetYaxis()->SetRangeUser(0,1);
		Pbckgnd_ROC_QAgl->Draw("APC");
		Pbckgnd_ROC_DistAgl->Draw("PCsame");


		TCanvas *c2 = new TCanvas("He rejection");
		
		c2->Divide(3,1);
		c2->cd(1);
	   	gPad->SetGridx();
   		gPad->SetGridy();
		Hebckgnd_ROC_QTOF->SetMarkerColor(2);
		Hebckgnd_ROC_QTOF->SetMarkerStyle(8);
		Hebckgnd_ROC_QTOF->SetLineColor(2);
		Hebckgnd_ROC_DistTOF->SetMarkerColor(4);
		Hebckgnd_ROC_DistTOF->SetMarkerStyle(8);
		Hebckgnd_ROC_DistTOF->SetLineColor(4);
			
		Hebckgnd_ROC_QTOF->GetXaxis()->SetTitle("Efficiency on signal");
		Hebckgnd_ROC_QTOF->GetYaxis()->SetTitle("Background rejection");
		Hebckgnd_ROC_QTOF->GetYaxis()->SetRangeUser(0,1);
		Hebckgnd_ROC_QTOF->Draw("APC");
		Hebckgnd_ROC_DistTOF->Draw("PCsame");

                c2->cd(2);
                gPad->SetGridx();
                gPad->SetGridy();
                Hebckgnd_ROC_QNaF->SetMarkerColor(2);
                Hebckgnd_ROC_QNaF->SetMarkerStyle(8);
                Hebckgnd_ROC_QNaF->SetLineColor(2);  
               	Hebckgnd_ROC_DistNaF->SetMarkerColor(4);
		Hebckgnd_ROC_DistNaF->SetMarkerStyle(8);
		Hebckgnd_ROC_DistNaF->SetLineColor(4);
		
		Hebckgnd_ROC_QNaF->GetXaxis()->SetTitle("Efficiency on signal");
                Hebckgnd_ROC_QNaF->GetYaxis()->SetTitle("Background rejection");
                Hebckgnd_ROC_QNaF->GetYaxis()->SetRangeUser(0,1);
		Hebckgnd_ROC_QNaF->Draw("APC");
		Hebckgnd_ROC_DistNaF->Draw("PCsame");

                c2->cd(3);
                gPad->SetGridx();
                gPad->SetGridy();
                Hebckgnd_ROC_QAgl->SetMarkerColor(2);
                Hebckgnd_ROC_QAgl->SetMarkerStyle(8);
                Hebckgnd_ROC_QAgl->SetLineColor(2);
               	Hebckgnd_ROC_DistAgl->SetMarkerColor(4);
		Hebckgnd_ROC_DistAgl->SetMarkerStyle(8);
		Hebckgnd_ROC_DistAgl->SetLineColor(4);
                
		Hebckgnd_ROC_QAgl->GetXaxis()->SetTitle("Efficiency on signal");
                Hebckgnd_ROC_QAgl->GetYaxis()->SetTitle("Background rejection");
                Hebckgnd_ROC_QAgl->GetYaxis()->SetRangeUser(0,1);
		Hebckgnd_ROC_QAgl->Draw("APC");
		Hebckgnd_ROC_DistAgl->Draw("PCsame");

		d->Divide(3,1);
                d->cd(1);
		QUtofvsLTOF_TOF->SetMarkerColor(2);
		QUtofvsLTOFD_TOF->SetMarkerColor(4);
		QUtofvsLTOFD_TOF->SetTitle("TOF Range");
		QUtofvsLTOFD_TOF->GetXaxis()->SetTitle("Q (Upper TOF)");
		QUtofvsLTOFD_TOF->GetYaxis()->SetTitle("Q (Lower TOF)");
		QUtofvsLTOFD_TOF->Draw("CONT3");
		QUtofvsLTOF_TOF->Draw("same");
		
		d->cd(2);
                QUtofvsLTOF_NaF->SetMarkerColor(2);
		QUtofvsLTOFD_NaF->SetMarkerColor(4);
		QUtofvsLTOFD_NaF->SetTitle("NaF Range");
		QUtofvsLTOFD_NaF->GetXaxis()->SetTitle("Q (Upper TOF)");
                QUtofvsLTOFD_NaF->GetYaxis()->SetTitle("Q (Lower TOF)");
                QUtofvsLTOFD_NaF->Draw("CONT3");
		QUtofvsLTOF_NaF->Draw("same");	

		d->cd(3);
                QUtofvsLTOF_Agl->SetMarkerColor(2);
                QUtofvsLTOFD_Agl->SetMarkerColor(4);
		QUtofvsLTOFD_Agl->SetTitle("Agl Range");
		QUtofvsLTOFD_Agl->GetXaxis()->SetTitle("Q (Upper TOF)");
                QUtofvsLTOFD_Agl->GetYaxis()->SetTitle("Q (Lower TOF)");
		QUtofvsLTOFD_Agl->Draw("CONT3");
		QUtofvsLTOF_Agl->Draw("same");		

		d1->Divide(3,1);
                d1->cd(1);
		QUtofvsTrack_TOF->SetMarkerColor(2);
		QUtofvsTrackD_TOF->SetMarkerColor(4);
		QUtofvsTrackD_TOF->SetTitle("TOF Range");
		QUtofvsTrackD_TOF->GetXaxis()->SetTitle("Q (Upper TOF)");
		QUtofvsTrackD_TOF->GetYaxis()->SetTitle("Q (Inner Tracker)");
		QUtofvsTrackD_TOF->Draw("CONT3");
		QUtofvsTrack_TOF->Draw("same");
		
		d1->cd(2);
                QUtofvsTrack_NaF->SetMarkerColor(2);
		QUtofvsTrackD_NaF->SetMarkerColor(4);
		QUtofvsTrackD_NaF->SetTitle("NaF Range");
		QUtofvsTrackD_NaF->GetXaxis()->SetTitle("Q (Upper TOF)");
                QUtofvsTrackD_NaF->GetYaxis()->SetTitle("Q (Inner Tracker)");
                QUtofvsTrackD_NaF->Draw("CONT3");
		QUtofvsTrack_NaF->Draw("same");	

		d1->cd(3);
                QUtofvsTrack_Agl->SetMarkerColor(2);
                QUtofvsTrackD_Agl->SetMarkerColor(4);
		QUtofvsTrackD_Agl->SetTitle("Agl Range");
		QUtofvsTrackD_Agl->GetXaxis()->SetTitle("Q (Upper TOF)");
                QUtofvsTrackD_Agl->GetYaxis()->SetTitle("Q (Inner Tracker)");
		QUtofvsTrackD_Agl->Draw("CONT3");
		QUtofvsTrack_Agl->Draw("same");		

		d2->Divide(3,1);
                d2->cd(1);
		QLtofvsTrack_TOF->SetMarkerColor(2);
		QLtofvsTrackD_TOF->SetMarkerColor(4);
		QLtofvsTrackD_TOF->SetTitle("TOF Range");
		QLtofvsTrackD_TOF->GetXaxis()->SetTitle("Q (Lower TOF)");
		QLtofvsTrackD_TOF->GetYaxis()->SetTitle("Q (Inner Tracker)");
		QLtofvsTrackD_TOF->Draw("CONT3");
		QLtofvsTrack_TOF->Draw("same");
		
		d2->cd(2);
                QLtofvsTrack_NaF->SetMarkerColor(2);
		QLtofvsTrackD_NaF->SetMarkerColor(4);
		QLtofvsTrackD_NaF->SetTitle("NaF Range");
		QLtofvsTrackD_NaF->GetXaxis()->SetTitle("Q (Lower TOF)");
                QLtofvsTrackD_NaF->GetYaxis()->SetTitle("Q (Inner Tracker)");
                QLtofvsTrackD_NaF->Draw("CONT3");
		QLtofvsTrack_NaF->Draw("same");	

		d2->cd(3);
                QLtofvsTrack_Agl->SetMarkerColor(2);
                QLtofvsTrackD_Agl->SetMarkerColor(4);
		QLtofvsTrackD_Agl->SetTitle("Agl Range");
		QLtofvsTrackD_Agl->GetXaxis()->SetTitle("Q (Lower TOF)");
                QLtofvsTrackD_Agl->GetYaxis()->SetTitle("Q (Inner Tracker)");
		QLtofvsTrackD_Agl->Draw("CONT3");
		QLtofvsTrack_Agl->Draw("same");		


		TH1F *  QUtof_TOF = ProjectionXtoTH1F(QUtofvsLTOF_TOF , "QUtof_TOF",0,QUtofvsLTOF_TOF ->GetNbinsY()+1);
		TH1F *  QUtof_NaF = ProjectionXtoTH1F(QUtofvsLTOF_NaF , "QUtof_NaF",0,QUtofvsLTOF_NaF ->GetNbinsY()+1);
		TH1F *  QUtof_Agl = ProjectionXtoTH1F(QUtofvsLTOF_Agl , "QUtof_Agl",0,QUtofvsLTOF_Agl ->GetNbinsY()+1);		
		TH1F *  QUtofD_TOF = ProjectionXtoTH1F(QUtofvsLTOFD_TOF , "QUtofD_TOF",0,QUtofvsLTOFD_TOF ->GetNbinsY()+1);
		TH1F *  QUtofD_NaF = ProjectionXtoTH1F(QUtofvsLTOFD_NaF , "QUtofD_NaF",0,QUtofvsLTOFD_NaF ->GetNbinsY()+1);
		TH1F *  QUtofD_Agl = ProjectionXtoTH1F(QUtofvsLTOFD_Agl , "QUtofD_Agl",0,QUtofvsLTOFD_Agl ->GetNbinsY()+1);		
		TH1F *  QUtofHe_TOF = ProjectionXtoTH1F(QUtofvsLTOFHe_TOF , "QUtofHe_TOF",0,QUtofvsLTOFHe_TOF ->GetNbinsY()+1);
		TH1F *  QUtofHe_NaF = ProjectionXtoTH1F(QUtofvsLTOFHe_NaF , "QUtofHe_NaF",0,QUtofvsLTOFHe_NaF ->GetNbinsY()+1);
		TH1F *  QUtofHe_Agl = ProjectionXtoTH1F(QUtofvsLTOFHe_Agl , "QUtofHe_Agl",0,QUtofvsLTOFHe_Agl ->GetNbinsY()+1);		


		TH1F *  QLtof_TOF = ProjectionYtoTH1F(QUtofvsLTOF_TOF , "QLtof_TOF",0,QUtofvsLTOF_TOF ->GetNbinsX()+1);
		TH1F *  QLtof_NaF = ProjectionYtoTH1F(QUtofvsLTOF_NaF , "QLtof_NaF",0,QUtofvsLTOF_NaF ->GetNbinsX()+1);
		TH1F *  QLtof_Agl = ProjectionYtoTH1F(QUtofvsLTOF_Agl , "QLtof_Agl",0,QUtofvsLTOF_Agl ->GetNbinsX()+1);		
		TH1F *  QLtofD_TOF = ProjectionYtoTH1F(QUtofvsLTOFD_TOF , "QLtofD_TOF",0,QUtofvsLTOFD_TOF ->GetNbinsX()+1);
		TH1F *  QLtofD_NaF = ProjectionYtoTH1F(QUtofvsLTOFD_NaF , "QLtofD_NaF",0,QUtofvsLTOFD_NaF ->GetNbinsX()+1);
		TH1F *  QLtofD_Agl = ProjectionYtoTH1F(QUtofvsLTOFD_Agl , "QLtofD_Agl",0,QUtofvsLTOFD_Agl ->GetNbinsX()+1);		
		TH1F *  QLtofHe_TOF = ProjectionYtoTH1F(QUtofvsLTOFHe_TOF , "QLtofHe_TOF",0,QUtofvsLTOFHe_TOF ->GetNbinsX()+1);
		TH1F *  QLtofHe_NaF = ProjectionYtoTH1F(QUtofvsLTOFHe_NaF , "QLtofHe_NaF",0,QUtofvsLTOFHe_NaF ->GetNbinsX()+1);
		TH1F *  QLtofHe_Agl = ProjectionYtoTH1F(QUtofvsLTOFHe_Agl , "QLtofHe_Agl",0,QUtofvsLTOFHe_Agl ->GetNbinsX()+1);		

		
		TH1F *  QTrack_TOF = ProjectionYtoTH1F(  QUtofvsTrack_TOF , "QTrack_TOF",0,	QUtofvsTrack_TOF ->GetNbinsX()+1);
		TH1F *  QTrack_NaF = ProjectionYtoTH1F(  QUtofvsTrack_NaF , "QTrack_NaF",0,	QUtofvsTrack_NaF ->GetNbinsX()+1);
		TH1F *  QTrack_Agl = ProjectionYtoTH1F(  QUtofvsTrack_Agl , "QTrack_Agl",0,	QUtofvsTrack_Agl ->GetNbinsX()+1);		
		TH1F *  QTrackD_TOF = ProjectionYtoTH1F( QUtofvsTrackD_TOF , "QTrackD_TOF",0,	QUtofvsTrackD_TOF ->GetNbinsX()+1);
		TH1F *  QTrackD_NaF = ProjectionYtoTH1F( QUtofvsTrackD_NaF , "QTrackD_NaF",0,	QUtofvsTrackD_NaF ->GetNbinsX()+1);
		TH1F *  QTrackD_Agl = ProjectionYtoTH1F( QUtofvsTrackD_Agl , "QTrackD_Agl",0,	QUtofvsTrackD_Agl ->GetNbinsX()+1);		
		TH1F *  QTrackHe_TOF = ProjectionYtoTH1F(QUtofvsTrackHe_TOF , "QTrackHe_TOF",0, QUtofvsTrackHe_TOF ->GetNbinsX()+1);
		TH1F *  QTrackHe_NaF = ProjectionYtoTH1F(QUtofvsTrackHe_NaF , "QTrackHe_NaF",0, QUtofvsTrackHe_NaF ->GetNbinsX()+1);
		TH1F *  QTrackHe_Agl = ProjectionYtoTH1F(QUtofvsTrackHe_Agl , "QTrackHe_Agl",0, QUtofvsTrackHe_Agl ->GetNbinsX()+1);		






		e->Divide(3,1);
                e->cd(1);
			
		QUtof_TOF->SetLineColor(2);
		QUtof_TOF->SetLineWidth(2);
		QUtofD_TOF->SetLineColor(4);
                QUtofD_TOF->SetLineWidth(2);
		QUtofHe_TOF->SetLineColor(3);
                QUtofHe_TOF->SetLineWidth(2);
		QUtofD_TOF->Scale(1/QUtofD_TOF->GetEntries());
		QUtof_TOF->Scale(1/QUtof_TOF->GetEntries());
		QUtofHe_TOF->Scale(1/QUtofHe_TOF->GetEntries());
		QUtofD_TOF->Draw();
		QUtof_TOF->Draw("same");
		QUtofHe_TOF->Draw("same");		

		e->cd(2);

                QUtof_NaF->SetLineColor(2);
                QUtof_NaF->SetLineWidth(2);
		QUtofD_NaF->SetLineColor(4);
                QUtofD_NaF->SetLineWidth(2);
		QUtofHe_NaF->SetLineColor(3);
                QUtofHe_NaF->SetLineWidth(2);
                QUtofD_NaF->Scale(1/QUtofD_NaF->GetEntries());
                QUtof_NaF->Scale(1/QUtof_NaF->GetEntries());
		QUtofHe_NaF->Scale(1/QUtofHe_NaF->GetEntries());
		QUtofD_NaF->Draw();
		QUtof_NaF->Draw("same");
		QUtofHe_NaF->Draw("same");
		
		e->cd(3);

                QUtof_Agl->SetLineColor(2);
                QUtof_Agl->SetLineWidth(2);
		QUtofD_Agl->SetLineColor(4);
                QUtofD_Agl->SetLineWidth(2);
		QUtofHe_Agl->SetLineColor(3);
                QUtofHe_Agl->SetLineWidth(2);
		QUtofD_Agl->Scale(1/QUtofD_Agl->GetEntries());
                QUtof_Agl->Scale(1/QUtof_Agl->GetEntries());
                QUtofHe_Agl->Scale(1/QUtofHe_Agl->GetEntries());
		QUtofD_Agl->Draw();
		QUtof_Agl->Draw("same");
		QUtofHe_Agl->Draw("same");


		e1->Divide(3,1);
                e1->cd(1);
			
		QLtof_TOF->SetLineColor(2);
		QLtof_TOF->SetLineWidth(2);
		QLtofD_TOF->SetLineColor(4);
                QLtofD_TOF->SetLineWidth(2);
		QLtofHe_TOF->SetLineColor(3);
                QLtofHe_TOF->SetLineWidth(2);
		QLtofD_TOF->Scale(1/QLtofD_TOF->GetEntries());
		QLtof_TOF->Scale(1/QLtof_TOF->GetEntries());
		QLtofHe_TOF->Scale(1/QLtofHe_TOF->GetEntries());
		QLtofD_TOF->Draw();
		QLtof_TOF->Draw("same");
		QLtofHe_TOF->Draw("same");		

		e1->cd(2);

                QLtof_NaF->SetLineColor(2);
                QLtof_NaF->SetLineWidth(2);
		QLtofD_NaF->SetLineColor(4);
                QLtofD_NaF->SetLineWidth(2);
		QLtofHe_NaF->SetLineColor(3);
                QLtofHe_NaF->SetLineWidth(2);
                QLtofD_NaF->Scale(1/QLtofD_NaF->GetEntries());
                QLtof_NaF->Scale(1/QLtof_NaF->GetEntries());
		QLtofHe_NaF->Scale(1/QUtofHe_NaF->GetEntries());
		QLtofD_NaF->Draw();
		QLtof_NaF->Draw("same");
		QLtofHe_NaF->Draw("same");
		
		e1->cd(3);

                QLtof_Agl->SetLineColor(2);
                QLtof_Agl->SetLineWidth(2);
		QLtofD_Agl->SetLineColor(4);
                QLtofD_Agl->SetLineWidth(2);
		QLtofHe_Agl->SetLineColor(3);
                QLtofHe_Agl->SetLineWidth(2);
		QLtofD_Agl->Scale(1/QLtofD_Agl->GetEntries());
                QLtof_Agl->Scale(1/QLtof_Agl->GetEntries());
                QLtofHe_Agl->Scale(1/QLtofHe_Agl->GetEntries());
		QLtofD_Agl->Draw();
		QLtof_Agl->Draw("same");
		QLtofHe_Agl->Draw("same");

		e2->Divide(3,1);
                e2->cd(1);
			
		QTrack_TOF->SetLineColor(2);
		QTrack_TOF->SetLineWidth(2);
		QTrackD_TOF->SetLineColor(4);
                QTrackD_TOF->SetLineWidth(2);
		QTrackHe_TOF->SetLineColor(3);
                QTrackHe_TOF->SetLineWidth(2);
		QTrackD_TOF->Scale(1/QTrackD_TOF->GetEntries());
		QTrack_TOF->Scale(1/QTrack_TOF->GetEntries());
		QTrackHe_TOF->Scale(1/QTrackHe_TOF->GetEntries());
		QTrackD_TOF->Draw();
		QTrack_TOF->Draw("same");
		QTrackHe_TOF->Draw("same");		

		e2->cd(2);

                QTrack_NaF->SetLineColor(2);
                QTrack_NaF->SetLineWidth(2);
		QTrackD_NaF->SetLineColor(4);
                QTrackD_NaF->SetLineWidth(2);
		QTrackHe_NaF->SetLineColor(3);
                QTrackHe_NaF->SetLineWidth(2);
                QTrackD_NaF->Scale(1/QTrackD_NaF->GetEntries());
                QTrack_NaF->Scale(1/QTrack_NaF->GetEntries());
		QTrackHe_NaF->Scale(1/QUtofHe_NaF->GetEntries());
		QTrackD_NaF->Draw();
		QTrack_NaF->Draw("same");
		QTrackHe_NaF->Draw("same");
		
		e2->cd(3);

                QTrack_Agl->SetLineColor(2);
                QTrack_Agl->SetLineWidth(2);
		QTrackD_Agl->SetLineColor(4);
                QTrackD_Agl->SetLineWidth(2);
		QTrackHe_Agl->SetLineColor(3);
                QTrackHe_Agl->SetLineWidth(2);
		QTrackD_Agl->Scale(1/QTrackD_Agl->GetEntries());
                QTrack_Agl->Scale(1/QTrack_Agl->GetEntries());
                QTrackHe_Agl->Scale(1/QTrackHe_Agl->GetEntries());
		QTrackD_Agl->Draw();
		QTrack_Agl->Draw("same");
		QTrackHe_Agl->Draw("same");

	


	









	
		finalPlots.Add(d);
		finalPlots.Add(d1);
                finalPlots.Add(d2);
		finalPlots.Add(c1);
		finalPlots.Add(c2);
		finalPlots.Add(e);
		finalPlots.Add(e1);
		finalPlots.Add(e2);
		finalPlots.writeObjsInFolder("MC Results/Q Check");
			
}

