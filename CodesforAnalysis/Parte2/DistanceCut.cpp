
TH2F * Dist5D_PdistrP_TOF  = new TH2F("Dist5D_PdistrP_TOF","Dist5D_PdistrP_TOF"  ,1000,0,100,1000,0,100);
TH2F * Dist5D_PdistrD_TOF  = new TH2F("Dist5D_PdistrD_TOF","Dist5D_PdistrD_TOF"  ,1000,0,100,1000,0,100);
TH2F * Dist5D_PdistrHe_TOF = new TH2F("Dist5D_PdistrHe_TOF","Dist5D_PdistrHe_TOF",1000,0,100,1000,0,100);

TH2F * Dist5DdistrP_TOF  = new TH2F("Dist5DdistrP_TOF","Dist5DdistrP_TOF"  ,1000,0,100,1000,0,100);
TH2F * Dist5DdistrD_TOF  = new TH2F("Dist5DdistrD_TOF","Dist5DdistrD_TOF"  ,1000,0,100,1000,0,100);
TH2F * Dist5DdistrHe_TOF = new TH2F("Dist5DdistrHe_TOF","Dist5DdistrHe_TOF",1000,0,100,1000,0,100);


TH2F * Dist5D_PdistrP_NaF  = new TH2F("Dist5D_PdistrP_NaF","Dist5D_PdistrP_NaF"  ,1000,0,100,1000,0,100);
TH2F * Dist5D_PdistrD_NaF  = new TH2F("Dist5D_PdistrD_NaF","Dist5D_PdistrD_NaF"  ,1000,0,100,1000,0,100);
TH2F * Dist5D_PdistrHe_NaF = new TH2F("Dist5D_PdistrHe_NaF","Dist5D_PdistrHe_NaF",1000,0,100,1000,0,100);

TH2F * Dist5DdistrP_NaF  = new TH2F("Dist5DdistrP_NaF","Dist5DdistrP_NaF"  ,1000,0,100,1000,0,100);
TH2F * Dist5DdistrD_NaF  = new TH2F("Dist5DdistrD_NaF","Dist5DdistrD_NaF"  ,1000,0,100,1000,0,100);
TH2F * Dist5DdistrHe_NaF = new TH2F("Dist5DdistrHe_NaF","Dist5DdistrHe_NaF",1000,0,100,1000,0,100);


TH2F * Dist5D_PdistrP_Agl  = new TH2F("Dist5D_PdistrP_Agl","Dist5D_PdistrP_Agl"  ,1000,0,100,1000,0,100);
TH2F * Dist5D_PdistrD_Agl  = new TH2F("Dist5D_PdistrD_Agl","Dist5D_PdistrD_Agl"  ,1000,0,100,1000,0,100);
TH2F * Dist5D_PdistrHe_Agl = new TH2F("Dist5D_PdistrHe_Agl","Dist5D_PdistrHe_Agl",1000,0,100,1000,0,100);

TH2F * Dist5DdistrP_Agl  = new TH2F("Dist5DdistrP_Agl","Dist5DdistrP_Agl"  ,1000,0,100,1000,0,100);
TH2F * Dist5DdistrD_Agl  = new TH2F("Dist5DdistrD_Agl","Dist5DdistrD_Agl"  ,1000,0,100,1000,0,100);
TH2F * Dist5DdistrHe_Agl = new TH2F("Dist5DdistrHe_Agl","Dist5DdistrHe_Agl",1000,0,100,1000,0,100);


float Eval_CutEff(TH2F * Histo,float cut);
TGraph * Plot_CutEff(TH2F * Histo);

float Eval_Herej(TH1F *HistoHe);
TGraph * Plot_Herej(TH1F * Histo);

void DistanceCut_Fill(TNtuple *ntupla, int l) {

	int k = ntupla->GetEvent(l);
	if(Beta<=0||R<=0) return;

	if(Betastrongcut){
		if(Massa_gen<1)	{
			if(((int)Cutmask)>>11!=0&&((int)Cutmask)>>11!=512) Dist5D_PdistrP_TOF -> Fill(Dist5D,Dist5D_P);
			if(((int)Cutmask)>>11==512)			   Dist5D_PdistrP_NaF -> Fill(Dist5D,Dist5D_P);
			if(((int)Cutmask)>>11==0)  			   Dist5D_PdistrP_Agl -> Fill(Dist5D,Dist5D_P);			
		}
		if(Massa_gen>1&&Massa_gen<2){
			if(((int)Cutmask)>>11!=0&&((int)Cutmask)>>11!=512) Dist5D_PdistrD_TOF -> Fill(Dist5D,Dist5D_P);
			if(((int)Cutmask)>>11==512)			   Dist5D_PdistrD_NaF -> Fill(Dist5D,Dist5D_P);
			if(((int)Cutmask)>>11==0)  			   Dist5D_PdistrD_Agl -> Fill(Dist5D,Dist5D_P);						
		}
		if(Massa_gen>3&&Massa_gen<4){
			if(((int)Cutmask)>>11!=0&&((int)Cutmask)>>11!=512) Dist5D_PdistrHe_TOF -> Fill(Dist5D,Dist5D_P);
			if(((int)Cutmask)>>11==512)			   Dist5D_PdistrHe_NaF -> Fill(Dist5D,Dist5D_P);
			if(((int)Cutmask)>>11==0)  			   Dist5D_PdistrHe_Agl -> Fill(Dist5D,Dist5D_P);                      
		}

	}

	return;
}


void DistanceCut_Write(){

	Dist5D_PdistrP_TOF  ->Write();
	Dist5D_PdistrD_TOF  ->Write();
	Dist5D_PdistrHe_TOF ->Write();

	Dist5D_PdistrP_NaF  ->Write();
	Dist5D_PdistrD_NaF  ->Write();
	Dist5D_PdistrHe_NaF ->Write();

	Dist5D_PdistrP_Agl  ->Write();
	Dist5D_PdistrD_Agl  ->Write();
	Dist5D_PdistrHe_Agl ->Write();


	return;


}

void DistanceCut(TFile * file1){

	TH2F * Dist5D_PdistrP_TOF =  (TH2F *)file1->Get("Dist5D_PdistrP_TOF" 	);
	TH2F * Dist5D_PdistrD_TOF =  (TH2F *)file1->Get("Dist5D_PdistrD_TOF" 	);
	TH2F * Dist5D_PdistrHe_TOF=  (TH2F *)file1->Get("Dist5D_PdistrHe_TOF"	);

	TH2F * Dist5DdistrP_TOF   =  (TH2F *)file1->Get("Dist5DdistrP_TOF"   	);
	TH2F * Dist5DdistrD_TOF   =  (TH2F *)file1->Get("Dist5DdistrD_TOF"   	);
	TH2F * Dist5DdistrHe_TOF  =  (TH2F *)file1->Get("Dist5DdistrHe_TOF"  	);


	TH2F * Dist5D_PdistrP_NaF =  (TH2F *)file1->Get("Dist5D_PdistrP_NaF" 	);
	TH2F * Dist5D_PdistrD_NaF =  (TH2F *)file1->Get("Dist5D_PdistrD_NaF" 	);
	TH2F * Dist5D_PdistrHe_NaF=  (TH2F *)file1->Get("Dist5D_PdistrHe_NaF"	);

	TH2F * Dist5DdistrP_NaF   =  (TH2F *)file1->Get("Dist5DdistrP_NaF"   	);
	TH2F * Dist5DdistrD_NaF   =  (TH2F *)file1->Get("Dist5DdistrD_NaF"   	);
	TH2F * Dist5DdistrHe_NaF  =  (TH2F *)file1->Get("Dist5DdistrHe_NaF"  	);


	TH2F * Dist5D_PdistrP_Agl =  (TH2F *)file1->Get("Dist5D_PdistrP_Agl" 	);
	TH2F * Dist5D_PdistrD_Agl =  (TH2F *)file1->Get("Dist5D_PdistrD_Agl" 	);
	TH2F * Dist5D_PdistrHe_Agl=  (TH2F *)file1->Get("Dist5D_PdistrHe_Agl"	);

	TH2F * Dist5DdistrP_Agl   =  (TH2F *)file1->Get("Dist5DdistrP_Agl"   	);
	TH2F * Dist5DdistrD_Agl   =  (TH2F *)file1->Get("Dist5DdistrD_Agl"   	);
	TH2F * Dist5DdistrHe_Agl  =  (TH2F *)file1->Get("Dist5DdistrHe_Agl"  	);


	TH1F * DistanceHe_TOF 	= (TH1F *)Dist5D_PdistrHe_TOF -> ProjectionX("DistanceHe_TOF",0,Dist5D_PdistrHe_TOF->GetNbinsY()) -> Clone();
	TH1F * DistanceHe_NaF   = (TH1F *)Dist5D_PdistrHe_NaF -> ProjectionX("DistanceHe_NaF",0,Dist5D_PdistrHe_NaF->GetNbinsY()) -> Clone();
	TH1F * DistanceHe_Agl   = (TH1F *)Dist5D_PdistrHe_Agl -> ProjectionX("DistanceHe_Agl",0,Dist5D_PdistrHe_Agl->GetNbinsY()) -> Clone();



	TH2F * Sum_TOF = (TH2F *) Dist5D_PdistrP_TOF -> Clone();
	Sum_TOF -> Add((TH2F *)Dist5D_PdistrD_TOF,10);
	Sum_TOF -> Add((TH2F *)Dist5D_PdistrHe_TOF,10);	

	TH2F * Sum_NaF = (TH2F *) Dist5D_PdistrP_NaF -> Clone();
	Sum_NaF -> Add((TH2F *)Dist5D_PdistrD_NaF->Clone(),10);
	Sum_NaF -> Add((TH2F *)Dist5D_PdistrHe_NaF->Clone(),10);	

	TH2F * Sum_Agl = (TH2F *) Dist5D_PdistrP_Agl -> Clone();
	Sum_Agl -> Add((TH2F *)Dist5D_PdistrD_Agl->Clone(),10);
	Sum_Agl -> Add((TH2F *)Dist5D_PdistrHe_Agl->Clone(),10);	


	cout<<"******** Distance cut global eff. *************"<<endl;

	cout<<"**TOF**"<<endl;
	TGraph * P_TOF_Efficiency = Plot_CutEff(Dist5D_PdistrP_TOF);
	TGraph * D_TOF_Efficiency = Plot_CutEff(Dist5D_PdistrD_TOF);	
	cout<<"**NaF**"<<endl;
	TGraph * P_NaF_Efficiency = Plot_CutEff(Dist5D_PdistrP_NaF);
        TGraph * D_NaF_Efficiency = Plot_CutEff(Dist5D_PdistrD_NaF);
	cout<<"**Agl**"<<endl;
        TGraph * P_Agl_Efficiency = Plot_CutEff(Dist5D_PdistrP_Agl);
        TGraph * D_Agl_Efficiency = Plot_CutEff(Dist5D_PdistrD_Agl);
	
	cout<<"******** Distance cut He rej. *************"<<endl;

	cout<<"**TOF**"<<endl;
	TGraph * Herej_TOF = Plot_Herej(DistanceHe_TOF);
	cout<<"**NaF**"<<endl;
	TGraph * Herej_NaF = Plot_Herej(DistanceHe_NaF);
	cout<<"**Agl**"<<endl;
	TGraph * Herej_Agl = Plot_Herej(DistanceHe_Agl);
	
	TCanvas * c1 = new TCanvas("Distance Distributions TOF");
	c1->cd();
	gPad->SetLogy();
	gPad->SetLogx();
	gPad->SetGridy();
        gPad->SetGridx();
	Sum_TOF -> GetXaxis() -> SetTitle("Distance from D");
	Sum_TOF -> GetYaxis() -> SetTitle("Distance from P");
	Sum_TOF -> SetTitle("Distance Distribution");
	Sum_TOF -> Draw("col");

	TCanvas * c2 = new TCanvas("Distance Distributions NaF");
        c2->cd();
        gPad->SetLogy();
        gPad->SetLogx();
        gPad->SetGridy();
        gPad->SetGridx();
        Sum_NaF -> GetXaxis() -> SetTitle("Distance from D");
        Sum_NaF -> GetYaxis() -> SetTitle("Distance from P");
        Sum_NaF -> SetTitle("Distance Distribution");
        Sum_NaF -> Draw("col");

	TCanvas * c3 = new TCanvas("Distance Distributions Agl");
        c3->cd();
        gPad->SetLogy();
        gPad->SetLogx();
        gPad->SetGridy();
        gPad->SetGridx();
        Sum_Agl -> GetXaxis() -> SetTitle("Distance from D");
        Sum_Agl -> GetYaxis() -> SetTitle("Distance from P");
        Sum_Agl -> SetTitle("Distance Distribution");
        Sum_Agl -> Draw("col");

	TCanvas * c4 = new TCanvas("Distance cut Eff.");
	c4->Divide(3,1);
	c4->cd(1);
	gPad->SetLogx();
	gPad->SetGridy();
        gPad->SetGridx();
	P_TOF_Efficiency -> SetLineColor(2);
	P_TOF_Efficiency -> SetMarkerColor(2);	
	P_TOF_Efficiency -> SetLineWidth(4);
	P_TOF_Efficiency ->SetTitle("TOF");
	P_TOF_Efficiency ->GetXaxis()->SetTitle("cut value");
	P_TOF_Efficiency ->GetYaxis()->SetTitle("Efficiency");
	P_TOF_Efficiency ->Draw("APC");
	D_TOF_Efficiency -> SetLineColor(4);
        D_TOF_Efficiency -> SetMarkerColor(4);
        D_TOF_Efficiency -> SetLineWidth(4);
        D_TOF_Efficiency ->Draw("PCsame");

	c4->cd(2);
        gPad->SetLogx();
        gPad->SetGridy();
        gPad->SetGridx();
        P_NaF_Efficiency -> SetLineColor(2);
        P_NaF_Efficiency -> SetMarkerColor(2);
        P_NaF_Efficiency -> SetLineWidth(4);
        P_NaF_Efficiency ->SetTitle("NaF");
        P_NaF_Efficiency ->GetXaxis()->SetTitle("cut value");
        P_NaF_Efficiency ->GetYaxis()->SetTitle("Efficiency");
        P_NaF_Efficiency ->Draw("APC");
        D_NaF_Efficiency -> SetLineColor(4);
        D_NaF_Efficiency -> SetMarkerColor(4);
        D_NaF_Efficiency -> SetLineWidth(4);
        D_TOF_Efficiency ->Draw("PCsame");

	c4->cd(3);
        gPad->SetLogx();
        gPad->SetGridy();
        gPad->SetGridx();
        P_Agl_Efficiency -> SetLineColor(2);
        P_Agl_Efficiency -> SetMarkerColor(2);
        P_Agl_Efficiency -> SetLineWidth(4);
        P_Agl_Efficiency ->SetTitle("Agl");
        P_Agl_Efficiency ->GetXaxis()->SetTitle("cut value");
        P_Agl_Efficiency ->GetYaxis()->SetTitle("Efficiency");
        P_Agl_Efficiency ->Draw("APC");
        D_Agl_Efficiency -> SetLineColor(4);
        D_Agl_Efficiency -> SetMarkerColor(4);
        D_Agl_Efficiency -> SetLineWidth(4);
        D_Agl_Efficiency ->Draw("PCsame");
	
	TCanvas * c5 = new TCanvas("Distance cut He rej.");
        c5->Divide(3,1);
        c5->cd(1);
        gPad->SetLogx();
        gPad->SetGridy();
        gPad->SetGridx();
        Herej_TOF -> SetLineColor(3);
        Herej_TOF -> SetMarkerColor(3);
        Herej_TOF -> SetLineWidth(4);
        Herej_TOF ->SetTitle("TOF");
        Herej_TOF ->GetXaxis()->SetTitle("cut value");
        Herej_TOF ->GetYaxis()->SetTitle("He rejection");
        Herej_TOF ->Draw("APC");

	c5->cd(2);
        gPad->SetLogx();
        gPad->SetGridy();
        gPad->SetGridx();
        Herej_NaF -> SetLineColor(3);
        Herej_NaF -> SetMarkerColor(3);
        Herej_NaF -> SetLineWidth(4);
        Herej_NaF ->SetTitle("TOF");
        Herej_NaF ->GetXaxis()->SetTitle("cut value");
        Herej_NaF ->GetYaxis()->SetTitle("He rejection");
        Herej_NaF ->Draw("APC");

	c5->cd(3);
        gPad->SetLogx();
        gPad->SetGridy();
        gPad->SetGridx();
        Herej_Agl -> SetLineColor(3);
        Herej_Agl -> SetMarkerColor(3);
        Herej_Agl -> SetLineWidth(4);
        Herej_Agl ->SetTitle("TOF");
        Herej_Agl ->GetXaxis()->SetTitle("cut value");
        Herej_Agl ->GetYaxis()->SetTitle("He rejection");
        Herej_Agl ->Draw("APC");








	cout<<"*** Updating Results file ***"<<endl;
	string nomefile="./Final_plots/"+mese+".root";
	TFile *f_out=new TFile(nomefile.c_str(), "UPDATE");
	f_out->mkdir("MC Results/Distance Cut");
	f_out->cd("MC Results/Distance Cut");
	c1   ->Write();
	c2   ->Write();
	c3   ->Write();
	c4   ->Write();
	c5   ->Write();
	f_out->Write();
	f_out->Close();

	return;
}


float Eval_CutEff(TH2F * Histo,float cut){
	float counts_passed_cut = 0;
	for (int x = 0; x < Histo ->GetNbinsX(); x++)
                for (int y = 0; y < Histo ->GetNbinsY(); y++){
                        float X = Histo -> GetXaxis() -> GetBinLowEdge(x);
			float Y = Histo -> GetYaxis() -> GetBinLowEdge(y);
			if(X<cut||Y<cut) counts_passed_cut += Histo -> GetBinContent(x+1,y+1);
                }
	return counts_passed_cut/(float)Histo->GetEntries();
}

float Eval_Herej(TH1F *HistoHe,float cut){
	float counts_passed_cut_He = 0;
	for (int x = 0; x < HistoHe ->GetNbinsX(); x++){
		float X = HistoHe -> GetXaxis() -> GetBinLowEdge(x);
		if(X<cut) counts_passed_cut_He += HistoHe-> GetBinContent(x+1);
			
	}
	float Heeff= counts_passed_cut_He/(float)HistoHe->GetEntries();
	return 1 - Heeff;
}


TGraph * Plot_CutEff(TH2F * Histo){
	TGraph * Efficiency = new TGraph();
        float cut = 0;
	float point =0;
        for(int x = 0; x <Histo->GetNbinsX();x+=10){
                        cut = Histo -> GetXaxis() -> GetBinLowEdge(x);
                        Efficiency -> SetPoint(point,cut,Eval_CutEff(Histo,cut));
			point++;
	}
	return Efficiency;
}

TGraph * Plot_Herej(TH1F * Histo){
        TGraph * Herej = new TGraph();
        float cut = 0;
        float point =0;
        for(int x = 0; x <Histo->GetNbinsX();x+=10){
                        cut = Histo -> GetXaxis() -> GetBinLowEdge(x);
                        Herej -> SetPoint(point,cut,Eval_Herej(Histo,cut));
                        point++;
        }
        return Herej;
}



