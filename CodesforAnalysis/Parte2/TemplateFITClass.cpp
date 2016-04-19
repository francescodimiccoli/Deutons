#include "TemplateFITClass.h"

void TemplateFIT::Write(){

	if(TemplateP  -> GetEntries() > 0)  TemplateP  -> Write();
	if(TemplateD  -> GetEntries() > 0)  TemplateD  -> Write();
	TemplateHe -> Write();

	if(DATA  -> GetEntries() > 0)  DATA  -> Write();
	return;
}

TH1F * TemplateFIT::Extract_Bin(TH1 * Histo, int bin,int lat){
	TH1F * Slice;
	Slice = new TH1F("","",Histo->GetNbinsX(),0,Histo->GetXaxis()->GetBinLowEdge(101));
                for(int i = 0; i< Histo->GetNbinsX();i++)
                        Slice->SetBinContent(i+1,Histo->GetBinContent(i+1,bin+1));
                return Slice;
	if(lat!=0){	
		Slice = (TH1F *)((TH3F*)Histo) -> ProjectionX ("",bin+1,bin+1,lat+1,lat+1) -> Clone();
                return Slice;
	}

}


void TemplateFIT::Do_TemplateFIT(TFit * Fit,int lat){
	TObjArray *Tpl;
	Tpl = new TObjArray(3);
	Tpl -> Add( Fit ->  Templ_P );
	Tpl -> Add( Fit ->  Templ_D );
	Tpl -> Add( Fit ->  Templ_He);

	Fit -> Tfit = new TFractionFitter(Fit -> Data, Tpl ,"q");
	Fit -> Tfit_outcome = 1;//fit -> Fit();
	fits[lat].push_back(Fit);

	return;
}

double TemplateFIT::GetFitWheights(int par, int bin,int lat){
	if(GetFitOutcome(bin,lat)==-1) return  1;
	if(GetFitOutcome(bin,lat)>0)   return  1;
	if(GetFitOutcome(bin,lat)==0){
		double w1,e1=0;
		fits[bin][lat]-> Tfit ->GetResult(par,w1,e1);
		TH1F * Result = (TH1F*)fits[bin][lat] -> Tfit -> GetPlot();
		float itot= Result->Integral();
		float i1;
		if(par == 0) i1 = fits[bin][lat]-> Templ_P ->Integral();
		if(par == 1) i1 = fits[bin][lat]-> Templ_D ->Integral();
		if(par == 2) i1 = fits[bin][lat]-> Templ_He ->Integral();
		return w1/(i1*itot); 
	}	
}

double TemplateFIT::GetFitErrors(int par,int bin,int lat){
	if(GetFitOutcome(bin,lat)==-1) return  0;
	if(GetFitOutcome(bin,lat)>0)   return  0;
	if(GetFitOutcome(bin,lat)==0){
		double w1,e1=0;
		double w2,e2=0;
		double w3,e3=0;
		fits[bin][lat]-> Tfit ->GetResult(0,w1,e1);
		fits[bin][lat]-> Tfit ->GetResult(1,w2,e2);
		fits[bin][lat]-> Tfit ->GetResult(2,w3,e3);

		float Cov01=fits[bin][lat]-> Tfit->GetFitter()->GetCovarianceMatrixElement(0,1);
		float Cov02=fits[bin][lat]-> Tfit->GetFitter()->GetCovarianceMatrixElement(0,2);
		float Cov12=fits[bin][lat]-> Tfit->GetFitter()->GetCovarianceMatrixElement(1,2);

		float Sigma=pow((pow(w2*e2,2)+pow(w1*e1,2)+pow(w3*e3,2)
					-2*Cov01*w1*w2-2*Cov02*w1*w3
					-2*Cov12*w2*w3)/2,0.5);

		double Err = pow((Sigma/w2,2) + pow(Sigma/w1,2),0.5); //Fit relative error
	
		TH1F * ResultPlot;  
		if(par == 0)	ResultPlot = GetResult_P (bin,lat);	
		if(par == 1)    ResultPlot = GetResult_D (bin,lat);
		if(par == 2)    ResultPlot = GetResult_He(bin,lat);

		return Err * ResultPlot->Integral(); //Fit absolute error
	}
}

void TemplateFIT::TemplateFits(){

	int loops = 1;
	if(Geomag) loops = DATA -> GetNbinsY();

	for(int lat = 1; lat < loops ; lat ++)
		for(int bin=0; bin<nbins ; bin++){
			TFit * Fit = new TFit;
			Fit->Templ_P =  (TH1F *)TemplateFIT::Extract_Bin  (TemplateP, bin);	
			Fit->Templ_D =  (TH1F *)TemplateFIT::Extract_Bin  (TemplateD, bin);	
			Fit->Templ_He=  (TH1F *)TemplateFIT::Extract_Bin  (TemplateHe,bin);
			if(Geomag) Fit->Data    =  (TH1F *)TemplateFIT::Extract_Bin(DATA      ,bin, lat);
			else 	   Fit->Data    =  (TH1F *)TemplateFIT::Extract_Bin(DATA      ,bin);

			TemplateFIT::Do_TemplateFIT(Fit,lat);

			TH1F * ResultPlot_P  = GetResult_P (bin,lat);		
			TH1F * ResultPlot_D  = GetResult_D (bin,lat);
			TH1F * ResultPlot_He = GetResult_He(bin,lat);

			if(!Geomag){
				PCounts -> SetBinContent(bin+1,ResultPlot_P->Integral());
				DCounts -> SetBinContent(bin+1,ResultPlot_D->Integral());
				PCounts -> SetBinError(bin+1,GetFitErrors(0,bin));
				DCounts -> SetBinError(bin+1,GetFitErrors(1,bin));
			}

			if(Geomag){
				PCounts -> SetBinContent(bin+1,lat+1,ResultPlot_P->Integral());
				DCounts -> SetBinContent(bin+1,lat+1,ResultPlot_D->Integral());
				PCounts -> SetBinError(bin+1,lat+1,GetFitErrors(0,bin,lat));
				DCounts -> SetBinError(bin+1,lat+1,GetFitErrors(1,bin,lat));
			}
		}

	return;
}


void TemplateFIT::TemplateFitPlot(TCanvas * c, std::string var_name,int bin,int lat){
	c -> cd();
	gPad-> SetLogy();
	gPad-> SetGridx();
	gPad-> SetGridy();
	THStack *Stack=new THStack("","");

	TH1F *PMC  = GetResult_P   (bin,lat);
	TH1F *DMC  = GetResult_D   (bin,lat);
	TH1F *HeMC = GetResult_He  (bin,lat);
	TH1F *Data = GetResult_Data(bin,lat);

	if(fits[bin][lat]->Tfit_outcome==0) TH1F * Result = (TH1F*)fits[bin][lat] -> Tfit -> GetPlot();	

	PMC -> SetFillColor(2);
	DMC -> SetFillColor(4);
	HeMC-> SetFillColor(3);
	Data->SetMarkerStyle(8);

	if(fits[bin][lat]->Tfit_outcome!=0){
		PMC -> SetFillStyle(3001);
		DMC -> SetFillStyle(3001);
		HeMC-> SetFillStyle(3001);
	}
	Stack->Add(PMC);
	Stack->Add(DMC);
	Stack->Add(HeMC);
	Stack->Draw();
	Stack-> GetXaxis()->SetTitle(var_name.c_str());
	Stack-> GetYaxis()->SetTitle("Counts");
	Data->Draw("epsame");

}
