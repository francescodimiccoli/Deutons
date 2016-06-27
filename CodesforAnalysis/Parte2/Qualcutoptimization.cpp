#include "PlottingFunctions/Qualcutoptimization.cpp"


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


void DistanceCut_Fill() {
	float mass=0;
	if(Betastrongcut){
		// Helium rej.
		if(Massa_gen<1)	{
			if(cmask.isOnlyFromToF()) Dist5D_PdistrP_TOF -> Fill(Tup.Dist5D,Tup.Dist5D_P);
			if(cmask.isFromNaF())	  Dist5D_PdistrP_NaF -> Fill(Tup.Dist5D,Tup.Dist5D_P);
			if(cmask.isFromAgl())  	  Dist5D_PdistrP_Agl -> Fill(Tup.Dist5D,Tup.Dist5D_P);			
		}
		if(Massa_gen>1&&Massa_gen<2){
			if(cmask.isOnlyFromToF()) Dist5D_PdistrD_TOF -> Fill(Tup.Dist5D,Tup.Dist5D_P);
			if(cmask.isFromNaF())	  Dist5D_PdistrD_NaF -> Fill(Tup.Dist5D,Tup.Dist5D_P);
			if(cmask.isFromAgl())  	  Dist5D_PdistrD_Agl -> Fill(Tup.Dist5D,Tup.Dist5D_P);						
		}
		if(Massa_gen>3&&Massa_gen<4){
			if(cmask.isOnlyFromToF()) Dist5D_PdistrHe_TOF -> Fill(Tup.Dist5D,Tup.Dist5D_P);
			if(cmask.isFromNaF())	  Dist5D_PdistrHe_NaF -> Fill(Tup.Dist5D,Tup.Dist5D_P);
			if(cmask.isFromAgl())  	  Dist5D_PdistrHe_Agl -> Fill(Tup.Dist5D,Tup.Dist5D_P);                      
		}

		//Qual. cuts optimization
		if(cmask.isOnlyFromToF()) {
			mass = (Tup.R/Tup.Beta)*pow(1-pow(Tup.Beta,2),0.5);
			if(mass>1.87){
				if(Massa_gen<1) 		DistvsLikTOF_P -> Fill(-log(1-Tup.LDiscriminant),Tup.Dist5D);
				if(Massa_gen>1&&Massa_gen<2)	DistvsLikTOF_D -> Fill(-log(1-Tup.LDiscriminant),Tup.Dist5D);
			}
		}	
		if(cmask.isFromNaF()){
                        mass = (Tup.R/Tup.BetaRICH)*pow(1-pow(Tup.BetaRICH,2),0.5);
                        if(mass>1.87){
                                if(Massa_gen<1)                 DistvsLikNaF_P -> Fill(-log(1-Tup.LDiscriminant),Tup.Dist5D);
                                if(Massa_gen>1&&Massa_gen<2)    DistvsLikNaF_D -> Fill(-log(1-Tup.LDiscriminant),Tup.Dist5D);
                        }
                }       

		if(cmask.isFromAgl()){
                        mass = (Tup.R/Tup.BetaRICH)*pow(1-pow(Tup.BetaRICH,2),0.5);
                        if(mass>1.87){
                                if(Massa_gen<1)                 DistvsLikAgl_P -> Fill(-log(1-Tup.LDiscriminant),Tup.Dist5D);
                                if(Massa_gen>1&&Massa_gen<2)    DistvsLikAgl_D -> Fill(-log(1-Tup.LDiscriminant),Tup.Dist5D);
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

void DistanceCut(string filename){


	cout<<"*** Reading  P1 file ****"<<endl;
        TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");


	TH2F * Dist5D_PdistrP_TOF =  (TH2F *)inputHistoFile->Get("Dist5D_PdistrP_TOF" 	);
	TH2F * Dist5D_PdistrD_TOF =  (TH2F *)inputHistoFile->Get("Dist5D_PdistrD_TOF" 	);
	TH2F * Dist5D_PdistrHe_TOF=  (TH2F *)inputHistoFile->Get("Dist5D_PdistrHe_TOF"	);

	TH2F * Dist5D_PdistrP_NaF =  (TH2F *)inputHistoFile->Get("Dist5D_PdistrP_NaF" 	);
	TH2F * Dist5D_PdistrD_NaF =  (TH2F *)inputHistoFile->Get("Dist5D_PdistrD_NaF" 	);
	TH2F * Dist5D_PdistrHe_NaF=  (TH2F *)inputHistoFile->Get("Dist5D_PdistrHe_NaF"	);

	TH2F * Dist5D_PdistrP_Agl =  (TH2F *)inputHistoFile->Get("Dist5D_PdistrP_Agl" 	);
	TH2F * Dist5D_PdistrD_Agl =  (TH2F *)inputHistoFile->Get("Dist5D_PdistrD_Agl" 	);
	TH2F * Dist5D_PdistrHe_Agl=  (TH2F *)inputHistoFile->Get("Dist5D_PdistrHe_Agl"	);



	TH1F * DistanceHe_TOF 	= ProjectionXtoTH1F(Dist5D_PdistrHe_TOF , "DistanceHe_TOF",0,Dist5D_PdistrHe_TOF->GetNbinsY());
	TH1F * DistanceHe_NaF   = ProjectionXtoTH1F(Dist5D_PdistrHe_NaF , "DistanceHe_NaF",0,Dist5D_PdistrHe_NaF->GetNbinsY());
	TH1F * DistanceHe_Agl   = ProjectionXtoTH1F(Dist5D_PdistrHe_Agl , "DistanceHe_Agl",0,Dist5D_PdistrHe_Agl->GetNbinsY());

	TH2F * DistvsLikTOF_P = (TH2F*)inputHistoFile->Get( "DistvsLikTOF_P" );
        TH2F * DistvsLikTOF_D = (TH2F*)inputHistoFile->Get( "DistvsLikTOF_D" );
                                                                    
        TH2F * DistvsLikNaF_P = (TH2F*)inputHistoFile->Get( "DistvsLikNaF_P" );
	TH2F * DistvsLikNaF_D = (TH2F*)inputHistoFile->Get( "DistvsLikNaF_D" );
                                                                    
        TH2F * DistvsLikAgl_P = (TH2F*)inputHistoFile->Get( "DistvsLikAgl_P" );
        TH2F * DistvsLikAgl_D = (TH2F*)inputHistoFile->Get( "DistvsLikAgl_D" );

	TH1F * Distance_goodD_TOF = ProjectionYtoTH1F(DistvsLikTOF_D , "Distance_goodD_TOF",0,DistvsLikTOF_D->GetNbinsY());
	TH1F * Distance_badP_TOF =  ProjectionYtoTH1F(DistvsLikTOF_P , "Distance_badP_TOF ",0,DistvsLikTOF_P->GetNbinsY());
	TH1F * Distance_goodD_NaF = ProjectionYtoTH1F(DistvsLikNaF_D , "Distance_goodD_NaF",0,DistvsLikNaF_D->GetNbinsY());
        TH1F * Distance_badP_NaF =  ProjectionYtoTH1F(DistvsLikNaF_P , "Distance_badP_NaF ",0,DistvsLikNaF_P->GetNbinsY());
	TH1F * Distance_goodD_Agl = ProjectionYtoTH1F(DistvsLikAgl_D , "Distance_goodD_Agl",0,DistvsLikAgl_D->GetNbinsY());
        TH1F * Distance_badP_Agl =  ProjectionYtoTH1F(DistvsLikAgl_P , "Distance_badP_Agl ",0,DistvsLikAgl_P->GetNbinsY());

	TH1F * Lik_goodD_TOF = ProjectionXtoTH1F(DistvsLikTOF_D , "Lik_goodD_TOF",0,DistvsLikTOF_D->GetNbinsX());
        TH1F * Lik_badP_TOF = ProjectionXtoTH1F(DistvsLikTOF_P ,  "Lik_badP_TOF ",0,DistvsLikTOF_P->GetNbinsX());
        TH1F * Lik_goodD_NaF = ProjectionXtoTH1F(DistvsLikNaF_D , "Lik_goodD_NaF",0,DistvsLikNaF_D->GetNbinsX());
        TH1F * Lik_badP_NaF = ProjectionXtoTH1F(DistvsLikNaF_P ,  "Lik_badP_NaF ",0,DistvsLikNaF_P->GetNbinsX());
        TH1F * Lik_goodD_Agl = ProjectionXtoTH1F(DistvsLikAgl_D , "Lik_goodD_Agl",0,DistvsLikAgl_D->GetNbinsX());
        TH1F * Lik_badP_Agl = ProjectionXtoTH1F(DistvsLikAgl_P ,  "Lik_badP_Agl ",0,DistvsLikAgl_P->GetNbinsX());

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


	DistanceCut_Plot(
			 Sum_TOF,
			Sum_NaF,
			Sum_Agl,
			P_TOF_Efficiency,
			D_TOF_Efficiency,
			P_NaF_Efficiency,
			D_NaF_Efficiency, 	 
			P_Agl_Efficiency,
			D_Agl_Efficiency, 	 
			Herej_TOF,
			Herej_NaF, 
			Herej_Agl, 
			DistvsLikTOF_P ,   
			DistvsLikTOF_D ,   
			DistvsLikNaF_P ,   
			DistvsLikNaF_D ,   
			DistvsLikAgl_P ,   
			DistvsLikAgl_D ,   
			BadPrej_TOF,
			BadPrej_NaF,
			BadPrej_Agl,   		
			BadPrejLik_TOF,
			BadPrejLik_NaF,
			BadPrejLik_Agl
			);

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
