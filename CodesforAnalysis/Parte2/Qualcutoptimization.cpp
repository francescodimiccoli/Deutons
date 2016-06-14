// Helium rej.
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


// Quality cuts optimization

TH2F * DistvsLikTOF_P = new TH2F("DistvsLikTOF_P","DistvsLikTOF_P",100,0,5,1000,0,100);
TH2F * DistvsLikTOF_D = new TH2F("DistvsLikTOF_D","DistvsLikTOF_D",100,0,5,1000,0,100);	

TH2F * DistvsLikNaF_P = new TH2F("DistvsLikNaF_P","DistvsLikNaF_P",100,0,5,1000,0,100);
TH2F * DistvsLikNaF_D = new TH2F("DistvsLikNaF_D","DistvsLikNaF_D",100,0,5,1000,0,100);

TH2F * DistvsLikAgl_P = new TH2F("DistvsLikAgl_P","DistvsLikAgl_P",100,0,5,1000,0,100);
TH2F * DistvsLikAgl_D = new TH2F("DistvsLikAgl_D","DistvsLikAgl_D",100,0,5,1000,0,100);



float Eval_CutEff(TH1 * Histo,float cut);
TGraph * Plot_CutEff(TH2F * Histo);

float Eval_Herej(TH1F *HistoHe);
TGraph * Plot_Herej(TH1F * Histo);

TGraph * Plot_BadPrej(TH1F * HistoP,TH1F * HistoD,bool reverse=false);


void DistanceCut_Fill(TNtuple *ntupla, int l) {

	 ntupla->GetEvent(l);
	if(Tup.Beta<=0||Tup.R<=0) return;
	float mass=0;
	if(Betastrongcut){
		// Helium rej.
		if(Massa_gen<1)	{
			if(((int)Tup.Cutmask)>>11!=0&&((int)Tup.Cutmask)>>11!=512) Dist5D_PdistrP_TOF -> Fill(Tup.Dist5D,Tup.Dist5D_P);
			if(cmask.isRichMeasureFromNaF())			   Dist5D_PdistrP_NaF -> Fill(Tup.Dist5D,Tup.Dist5D_P);
			if(cmask.isRichMeasureFromAgl())  			   Dist5D_PdistrP_Agl -> Fill(Tup.Dist5D,Tup.Dist5D_P);			
		}
		if(Massa_gen>1&&Massa_gen<2){
			if(((int)Tup.Cutmask)>>11!=0&&((int)Tup.Cutmask)>>11!=512) Dist5D_PdistrD_TOF -> Fill(Tup.Dist5D,Tup.Dist5D_P);
			if(cmask.isRichMeasureFromNaF())			   Dist5D_PdistrD_NaF -> Fill(Tup.Dist5D,Tup.Dist5D_P);
			if(cmask.isRichMeasureFromAgl())  			   Dist5D_PdistrD_Agl -> Fill(Tup.Dist5D,Tup.Dist5D_P);						
		}
		if(Massa_gen>3&&Massa_gen<4){
			if(((int)Tup.Cutmask)>>11!=0&&((int)Tup.Cutmask)>>11!=512) Dist5D_PdistrHe_TOF -> Fill(Tup.Dist5D,Tup.Dist5D_P);
			if(cmask.isRichMeasureFromNaF())			   Dist5D_PdistrHe_NaF -> Fill(Tup.Dist5D,Tup.Dist5D_P);
			if(cmask.isRichMeasureFromAgl())  			   Dist5D_PdistrHe_Agl -> Fill(Tup.Dist5D,Tup.Dist5D_P);                      
		}

		//Qual. cuts optimization
		if(((int)Tup.Cutmask)>>11!=0&&((int)Tup.Cutmask)>>11!=512) {
			mass = (Tup.R/Tup.Beta)*pow(1-pow(Tup.Beta,2),0.5);
			if(mass>1.87){
				if(Massa_gen<1) 			   DistvsLikTOF_P -> Fill(-log(1-Tup.LDiscriminant),Tup.Dist5D);
				if(Massa_gen>1&&Massa_gen<2)		   DistvsLikTOF_D -> Fill(-log(1-Tup.LDiscriminant),Tup.Dist5D);
			}
		}	
		if(cmask.isRichMeasureFromNaF()){
                        mass = (Tup.R/Tup.BetaRICH)*pow(1-pow(Tup.BetaRICH,2),0.5);
                        if(mass>1.87){
                                if(Massa_gen<1)                            DistvsLikNaF_P -> Fill(-log(1-Tup.LDiscriminant),Tup.Dist5D);
                                if(Massa_gen>1&&Massa_gen<2)               DistvsLikNaF_D -> Fill(-log(1-Tup.LDiscriminant),Tup.Dist5D);
                        }
                }       

		if(cmask.isRichMeasureFromAgl()){
                        mass = (Tup.R/Tup.BetaRICH)*pow(1-pow(Tup.BetaRICH,2),0.5);
                        if(mass>1.87){
                                if(Massa_gen<1)                            DistvsLikAgl_P -> Fill(-log(1-Tup.LDiscriminant),Tup.Dist5D);
                                if(Massa_gen>1&&Massa_gen<2)               DistvsLikAgl_D -> Fill(-log(1-Tup.LDiscriminant),Tup.Dist5D);
                        }
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

	DistvsLikTOF_P	    ->Write();
	DistvsLikTOF_D	    ->Write();
              
	DistvsLikNaF_P 	    ->Write();
	DistvsLikNaF_D	    ->Write();
              
	DistvsLikAgl_P	    ->Write();
	DistvsLikAgl_D      ->Write();

	return;
}

void DistanceCut(TFile * file1){

	TH2F * Dist5D_PdistrP_TOF =  (TH2F *)file1->Get("Dist5D_PdistrP_TOF" 	);
	TH2F * Dist5D_PdistrD_TOF =  (TH2F *)file1->Get("Dist5D_PdistrD_TOF" 	);
	TH2F * Dist5D_PdistrHe_TOF=  (TH2F *)file1->Get("Dist5D_PdistrHe_TOF"	);

	TH2F * Dist5D_PdistrP_NaF =  (TH2F *)file1->Get("Dist5D_PdistrP_NaF" 	);
	TH2F * Dist5D_PdistrD_NaF =  (TH2F *)file1->Get("Dist5D_PdistrD_NaF" 	);
	TH2F * Dist5D_PdistrHe_NaF=  (TH2F *)file1->Get("Dist5D_PdistrHe_NaF"	);

	TH2F * Dist5D_PdistrP_Agl =  (TH2F *)file1->Get("Dist5D_PdistrP_Agl" 	);
	TH2F * Dist5D_PdistrD_Agl =  (TH2F *)file1->Get("Dist5D_PdistrD_Agl" 	);
	TH2F * Dist5D_PdistrHe_Agl=  (TH2F *)file1->Get("Dist5D_PdistrHe_Agl"	);



	TH1F * DistanceHe_TOF 	= (TH1F *)Dist5D_PdistrHe_TOF -> ProjectionX("DistanceHe_TOF",0,Dist5D_PdistrHe_TOF->GetNbinsY()) -> Clone();
	TH1F * DistanceHe_NaF   = (TH1F *)Dist5D_PdistrHe_NaF -> ProjectionX("DistanceHe_NaF",0,Dist5D_PdistrHe_NaF->GetNbinsY()) -> Clone();
	TH1F * DistanceHe_Agl   = (TH1F *)Dist5D_PdistrHe_Agl -> ProjectionX("DistanceHe_Agl",0,Dist5D_PdistrHe_Agl->GetNbinsY()) -> Clone();

	TH2F * DistvsLikTOF_P = (TH2F*)file1->Get(	"DistvsLikTOF_P"		);
        TH2F * DistvsLikTOF_D = (TH2F*)file1->Get(	"DistvsLikTOF_D"		);
                                                                      
        TH2F * DistvsLikNaF_P = (TH2F*)file1->Get(	"DistvsLikNaF_P"		);
	TH2F * DistvsLikNaF_D = (TH2F*)file1->Get(	"DistvsLikNaF_D"		);
                                                                      
        TH2F * DistvsLikAgl_P = (TH2F*)file1->Get(	"DistvsLikAgl_P"		);
        TH2F * DistvsLikAgl_D = (TH2F*)file1->Get(	"DistvsLikAgl_D"		);

	TH1F * Distance_goodD_TOF = (TH1F *)DistvsLikTOF_D -> ProjectionY("Distance_goodD_TOF",0,DistvsLikTOF_D->GetNbinsY()) -> Clone();
	TH1F * Distance_badP_TOF = (TH1F *)DistvsLikTOF_P -> ProjectionY( "Distance_badP_TOF ",0,DistvsLikTOF_P->GetNbinsY()) -> Clone();
	TH1F * Distance_goodD_NaF = (TH1F *)DistvsLikNaF_D -> ProjectionY("Distance_goodD_NaF",0,DistvsLikNaF_D->GetNbinsY()) -> Clone();
        TH1F * Distance_badP_NaF = (TH1F *)DistvsLikNaF_P -> ProjectionY( "Distance_badP_NaF ",0,DistvsLikNaF_P->GetNbinsY()) -> Clone();
	TH1F * Distance_goodD_Agl = (TH1F *)DistvsLikAgl_D -> ProjectionY("Distance_goodD_Agl",0,DistvsLikAgl_D->GetNbinsY()) -> Clone();
        TH1F * Distance_badP_Agl = (TH1F *)DistvsLikAgl_P -> ProjectionY( "Distance_badP_Agl ",0,DistvsLikAgl_P->GetNbinsY()) -> Clone();

	TH1F * Lik_goodD_TOF = (TH1F *)DistvsLikTOF_D -> ProjectionX("Lik_goodD_TOF",0,DistvsLikTOF_D->GetNbinsX()) -> Clone();
        TH1F * Lik_badP_TOF = (TH1F *)DistvsLikTOF_P -> ProjectionX( "Lik_badP_TOF ",0,DistvsLikTOF_P->GetNbinsX()) -> Clone();
        TH1F * Lik_goodD_NaF = (TH1F *)DistvsLikNaF_D -> ProjectionX("Lik_goodD_NaF",0,DistvsLikNaF_D->GetNbinsX()) -> Clone();
        TH1F * Lik_badP_NaF = (TH1F *)DistvsLikNaF_P -> ProjectionX( "Lik_badP_NaF ",0,DistvsLikNaF_P->GetNbinsX()) -> Clone();
        TH1F * Lik_goodD_Agl = (TH1F *)DistvsLikAgl_D -> ProjectionX("Lik_goodD_Agl",0,DistvsLikAgl_D->GetNbinsX()) -> Clone();
        TH1F * Lik_badP_Agl = (TH1F *)DistvsLikAgl_P -> ProjectionX( "Lik_badP_Agl ",0,DistvsLikAgl_P->GetNbinsX()) -> Clone();

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
	
	cout<<"******** Distance cut Bad P rej. *************"<<endl;
        
	cout<<"**TOF**"<<endl;
	TGraph * BadPrej_TOF = Plot_BadPrej(Distance_badP_TOF,Distance_goodD_TOF);

	cout<<"**NaF**"<<endl;
	TGraph * BadPrej_NaF = Plot_BadPrej(Distance_badP_NaF,Distance_goodD_NaF);

        cout<<"**Agl**"<<endl;

	TGraph * BadPrej_Agl = Plot_BadPrej(Distance_badP_Agl,Distance_goodD_Agl);
	

	cout<<"******** Likelihood cut Bad P rej. *************"<<endl;
	bool reverse = true;
        cout<<"**TOF**"<<endl;
        
	TGraph * BadPrejLik_TOF = Plot_BadPrej(Lik_badP_TOF,Lik_goodD_TOF,reverse);

        cout<<"**NaF**"<<endl;
        TGraph * BadPrejLik_NaF = Plot_BadPrej(Lik_badP_NaF,Lik_goodD_NaF,reverse);

        cout<<"**Agl**"<<endl;
        TGraph * BadPrejLik_Agl = Plot_BadPrej(Lik_badP_Agl,Lik_goodD_Agl,reverse);


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

	
	TCanvas * c6 = new TCanvas("Distance vs Likelihood");
        c6->Divide(3,1);
        c6->cd(1);
        gPad->SetLogy();
        gPad->SetGridy();
        gPad->SetGridx();
	DistvsLikTOF_P->SetMarkerColor(2);
	DistvsLikTOF_D->SetMarkerColor(4);
	DistvsLikTOF_D->SetMarkerStyle(8);
	DistvsLikTOF_P->SetMarkerStyle(8);
	DistvsLikTOF_D->SetMarkerSize(0.2);
        DistvsLikTOF_P->SetMarkerSize(0.3);
	DistvsLikTOF_D->SetTitle("Distance vs Likelihood: TOF");
	DistvsLikTOF_D->GetXaxis()->SetRangeUser(0,2.6);
	DistvsLikTOF_D->GetXaxis()->SetTitle("-log(1-Tup.LDiscriminant)");
	DistvsLikTOF_D->GetYaxis()->SetTitle("Distance from D");
	DistvsLikTOF_D->Draw();
	DistvsLikTOF_P->Draw("same");

        c6->cd(2);
        gPad->SetLogy();
        gPad->SetGridy();
        gPad->SetGridx();
	DistvsLikNaF_P->SetMarkerColor(2);
        DistvsLikNaF_D->SetMarkerColor(4);
	DistvsLikNaF_D->SetMarkerStyle(8);
        DistvsLikNaF_P->SetMarkerStyle(8);
        DistvsLikNaF_D->SetMarkerSize(0.3);
        DistvsLikNaF_P->SetMarkerSize(0.3);
	DistvsLikNaF_D->SetTitle("Distance vs Likelihood: NaF");
	DistvsLikNaF_D->GetXaxis()->SetTitle("-log(1-Tup.LDiscriminant)");
        DistvsLikNaF_D->GetYaxis()->SetTitle("Distance from D");
	DistvsLikNaF_D->Draw();
        DistvsLikNaF_P->Draw("same");


        c6->cd(3);
        gPad->SetLogy();
        gPad->SetGridy();
        gPad->SetGridx();
	DistvsLikAgl_P->SetMarkerColor(2);
        DistvsLikAgl_D->SetMarkerColor(4);
        DistvsLikAgl_D->SetMarkerStyle(8);
        DistvsLikAgl_P->SetMarkerStyle(8);
        DistvsLikAgl_D->SetMarkerSize(0.3);
        DistvsLikAgl_P->SetMarkerSize(0.3);
	DistvsLikAgl_D->SetTitle("Distance vs Likelihood: Agl");
	DistvsLikAgl_D->GetXaxis()->SetTitle("-log(1-Tup.LDiscriminant)");
        DistvsLikAgl_D->GetYaxis()->SetTitle("Distance from D");
	DistvsLikAgl_D->Draw();
        DistvsLikAgl_P->Draw("same");

	TCanvas * c7 = new TCanvas("Distance Bad P optimization");
        c7->Divide(3,1);
        c7->cd(1);
        gPad->SetLogx();
        gPad->SetGridy();
        gPad->SetGridx();
        BadPrej_TOF -> SetLineColor(3);
        BadPrej_TOF -> SetMarkerColor(3);
        BadPrej_TOF -> SetLineWidth(4);
        BadPrej_TOF ->SetTitle("TOF");
        BadPrej_TOF ->GetXaxis()->SetTitle("cut value");
        BadPrej_TOF ->GetYaxis()->SetTitle("Bad P rejection");
        BadPrej_TOF ->Draw("APC");

        c7->cd(2);
        gPad->SetLogx();
        gPad->SetGridy();
        gPad->SetGridx();
        BadPrej_NaF  -> SetLineColor(3);
        BadPrej_NaF  -> SetMarkerColor(3);
        BadPrej_NaF  -> SetLineWidth(4);
        BadPrej_NaF  ->SetTitle("NaF");
        BadPrej_NaF  ->GetXaxis()->SetTitle("cut value");
        BadPrej_NaF  ->GetYaxis()->SetTitle("Bad P rejection");
        BadPrej_NaF  ->Draw("APC");

        c7->cd(3);
        gPad->SetLogx();
        gPad->SetGridy();
        gPad->SetGridx();
        BadPrej_Agl -> SetLineColor(3);
        BadPrej_Agl -> SetMarkerColor(3);
        BadPrej_Agl -> SetLineWidth(4);
        BadPrej_Agl ->SetTitle("Agl");
        BadPrej_Agl ->GetXaxis()->SetTitle("cut value");
        BadPrej_Agl ->GetYaxis()->SetTitle("Bad P rejection");
        BadPrej_Agl ->Draw("APC");

	TCanvas * c8 = new TCanvas("Likelihood Bad P optimization");
        c8->Divide(3,1);
        c8->cd(1);
        gPad->SetLogx();
        gPad->SetGridy();
        gPad->SetGridx();
        BadPrejLik_TOF -> SetLineColor(3);
        BadPrejLik_TOF -> SetMarkerColor(3);
        BadPrejLik_TOF -> SetLineWidth(4);
        BadPrejLik_TOF ->SetTitle("TOF");
        BadPrejLik_TOF ->GetXaxis()->SetTitle("cut value");
        BadPrejLik_TOF ->GetYaxis()->SetTitle("Bad P rejection");
        BadPrejLik_TOF ->Draw("APC");

        c8->cd(2);
        gPad->SetLogx();
        gPad->SetGridy();
        gPad->SetGridx();
        BadPrejLik_NaF  -> SetLineColor(3);
        BadPrejLik_NaF  -> SetMarkerColor(3);
        BadPrejLik_NaF  -> SetLineWidth(4);
        BadPrejLik_NaF  ->SetTitle("NaF");
        BadPrejLik_NaF  ->GetXaxis()->SetTitle("cut value");
        BadPrejLik_NaF  ->GetYaxis()->SetTitle("Bad P rejection");
        BadPrejLik_NaF  ->Draw("APC");

        c8->cd(3);
        gPad->SetLogx();
        gPad->SetGridy();
        gPad->SetGridx();
        BadPrejLik_Agl -> SetLineColor(3);
        BadPrejLik_Agl -> SetMarkerColor(3);
        BadPrejLik_Agl -> SetLineWidth(4);
        BadPrejLik_Agl ->SetTitle("Agl");
        BadPrejLik_Agl ->GetXaxis()->SetTitle("cut value");
        BadPrejLik_Agl ->GetYaxis()->SetTitle("Bad P rejection");
        BadPrejLik_Agl ->Draw("APC");


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
	c6   ->Write();
	c7   ->Write();
	c8   ->Write();
	f_out->Write();
	f_out->Close();

	return;
}


float Eval_CutEff(TH1 * Histo,float cut){
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


TGraph * Plot_BadPrej(TH1F * HistoP,TH1F * HistoD,bool reverse){
	TGraph * BadPrej = new TGraph();
	float cut = 0;
	float point =0;
	if(!reverse)
		for(int x = 0; x <HistoD->GetNbinsX();x++){
			cut = HistoP -> GetXaxis() -> GetBinLowEdge(x);
			float D_eff = HistoD ->Integral(0,x)/(float)HistoD -> GetEntries();
			float P_eff = HistoP ->Integral(0,x)/(float)HistoP -> GetEntries();
			float opt = 0;
			if(P_eff>0) opt = D_eff/pow(P_eff,0.5);
			BadPrej -> SetPoint(point,cut,opt); 
			point++;
		}
	else
		for(int x = HistoD->GetNbinsX(); x > 0;x--){
                        cut = HistoP -> GetXaxis() -> GetBinLowEdge(x);
                        float D_eff = HistoD ->Integral(x,HistoP -> GetNbinsX())/(float)HistoD -> GetEntries();
                        float P_eff = HistoP ->Integral(x,HistoP -> GetNbinsX())/(float)HistoP -> GetEntries();
                        float opt = 0;
                        if(P_eff>0) opt = D_eff/pow(P_eff,0.5);
                        BadPrej -> SetPoint(point,cut,opt);
                        point++;
                }
		
	return BadPrej;
}
