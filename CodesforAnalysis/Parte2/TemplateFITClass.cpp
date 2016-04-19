#include "TemplateFITClass.h"

void TemplateFIT::Write(){

	if(TemplateP  -> GetEntries() > 0)  TemplateP  -> Write();
	if(TemplateD  -> GetEntries() > 0)  TemplateD  -> Write();
	TemplateHe -> Write();

	if(Data_Prim  -> GetEntries() > 0)  Data_Prim  -> Write();
	if(Data_Geomag-> GetEntries() > 0)  Data_Geomag-> Write();
	return;
}


TH1F * TemplateFIT::Extract_Bin_histos(TH1 * Histo, int bin){
	TH1F * Slice = new TH1F("","",Histo->GetNbinsX(),0,Histo->GetXaxis()->GetBinLowEdge(101));
	for(int i = 0; i< Histo->GetNbinsX();i++)
		Slice->SetBinContent(i+1,Histo->GetBinContent(i+1,bin+1));
	return Slice;
}

TH1F * TemplateFIT::Extract_Bin_histos_geo(TH1 * Histo, int bin, int lat){
	TH1F * Slice = (TH1F *)((TH3F*)Histo) -> ProjectionX ("",bin+1,bin+1,lat+1,lat+1) -> Clone();
	return Slice;
}

void TemplateFIT::Do_TemplateFIT(TFit * Fit){
	TObjArray *Tpl;
	Tpl = new TObjArray(3);
	Tpl -> Add( Fit ->  Templ_P );
	Tpl -> Add( Fit ->  Templ_D );
	Tpl -> Add( Fit ->  Templ_He);

	Fit -> Tfit = new TFractionFitter(Fit -> Data, Tpl ,"q");
	Fit -> Tfit_outcome = 1;//fit -> Fit();
	fits.push_back(Fit);

	return;
}

double TemplateFIT::GetFitWheights(int par, int bin){
	if(GetFitOutcome(bin)==-1) return  1;
	if(GetFitOutcome(bin)>0)   return  1;
	if(GetFitOutcome(bin)==0){
		double w1,e1=0;
		fits[bin]-> Tfit ->GetResult(par,w1,e1);
		TH1F * Result = (TH1F*)fits[bin] -> Tfit -> GetPlot();
		float itot= Result->Integral();
		float i1;
		if(par == 0) i1 = fits[bin]-> Templ_P ->Integral();
		if(par == 1) i1 = fits[bin]-> Templ_D ->Integral();
		if(par == 2) i1 = fits[bin]-> Templ_He ->Integral();
		return w1/(i1*itot); 
	}	
}

double TemplateFIT::GetFitErrors(int par,int bin){
	if(GetFitOutcome(bin)==-1) return  0;
	if(GetFitOutcome(bin)>0)   return  0;
	if(GetFitOutcome(bin)==0){
		double w1,e1=0;
		double w2,e2=0;
		double w3,e3=0;
		fits[bin]-> Tfit ->GetResult(0,w1,e1);
		fits[bin]-> Tfit ->GetResult(1,w2,e2);
		fits[bin]-> Tfit ->GetResult(2,w3,e3);

		float Cov01=fits[bin]-> Tfit->GetFitter()->GetCovarianceMatrixElement(0,1);
		float Cov02=fits[bin]-> Tfit->GetFitter()->GetCovarianceMatrixElement(0,2);
		float Cov12=fits[bin]-> Tfit->GetFitter()->GetCovarianceMatrixElement(1,2);

		float Sigma=pow((pow(w2*e2,2)+pow(w1*e1,2)+pow(w3*e3,2)
					-2*Cov01*w1*w2-2*Cov02*w1*w3
					-2*Cov12*w2*w3)/2,0.5);

		double Err = pow((Sigma/w2,2) + pow(Sigma/w1,2),0.5); //Fit relative error
		TH1F * ResultPlot;  
		if(par == 0)	ResultPlot = GetResult_P (bin);	
		if(par == 1)    ResultPlot = GetResult_D (bin);
		if(par == 2)    ResultPlot = GetResult_He (bin);

		return Err*ResultPlot->Integral(); //Fit absolute error
	}
}

void TemplateFIT::TemplateFits(){

	for(int bin=0; bin<nbins ; bin++){
		TFit * Fit = new TFit;
		Fit->Templ_P =  (TH1F *)TemplateFIT::Extract_Bin_histos(TemplateP, bin);	
		Fit->Templ_D =  (TH1F *)TemplateFIT::Extract_Bin_histos(TemplateD, bin);	
		Fit->Templ_He=  (TH1F *)TemplateFIT::Extract_Bin_histos(TemplateHe,bin);
		Fit->Data    =  (TH1F *)TemplateFIT::Extract_Bin_histos(Data_Prim ,bin);
		TemplateFIT::Do_TemplateFIT(Fit);

		TH1F * ResultPlot_P  = GetResult_P (bin);		
		TH1F * ResultPlot_D  = GetResult_D (bin);
		TH1F * ResultPlot_He = GetResult_He(bin);

		PCounts -> SetBinContent(bin+1,ResultPlot_P->Integral());
		DCounts -> SetBinContent(bin+1,ResultPlot_D->Integral());
		
		PCounts -> SetBinError(bin+1,GetFitErrors(0,bin));
                DCounts -> SetBinError(bin+1,GetFitErrors(1,bin));
	}
	return;
}


void TemplateFIT::TemplateFitPlot(TCanvas * c, std::string var_name,int bin){
	c -> cd();
	gPad-> SetLogy();
	gPad-> SetGridx();
	gPad-> SetGridy();
	THStack *Stack=new THStack("","");
	TH1F *PMC  = GetResult_P(bin);
	TH1F *DMC  = GetResult_D(bin);
	TH1F *HeMC = GetResult_He(bin);
	TH1F *Data = GetResult_Data(bin);
	if(fits[bin]->Tfit_outcome==0) TH1F * Result = (TH1F*)fits[bin] -> Tfit -> GetPlot();	
	PMC -> SetFillColor(2);
	DMC -> SetFillColor(4);
	HeMC-> SetFillColor(3);
	Data->SetMarkerStyle(8);
	if(fits[bin]->Tfit_outcome!=0){
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
