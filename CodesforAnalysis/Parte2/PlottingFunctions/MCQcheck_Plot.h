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
                        TGraph *Hebckgnd_ROC_DistAgl 


		     )	
{

		TCanvas *c1 = new TCanvas("Bad Protons rejection");
		
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


		finalPlots.Add(c1);
		finalPlots.Add(c2);
		finalPlots.writeObjsInFolder("MC Results/Q Check");
			
}

