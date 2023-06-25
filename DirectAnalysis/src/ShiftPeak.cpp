#include "ShiftPeak.h"

TSpline3 * Model_Histo(TH1F * Histo){
	int nbins =Histo->GetNbinsX();
	double X[nbins];
	double Y[nbins];

	for(int i=0;i<Histo->GetNbinsX();i++){
		X[i]=Histo->GetBinCenter(i+1);
		Y[i]=Histo->GetBinContent(i+1);
	}

	TSpline3 * Model = new TSpline3("Model",X,Y,nbins);
	return Model;
}


TF1 * WeightGausswithSpline(TSpline3 * Spline, float mean,float mu, float sigma){
	TF1 * f1 = new TF1("f1","gaus",0,10);
	f1->SetParameters(Spline->Eval(mean),mean-mu,sigma);
	return f1;
}


TH1F * InfinitesimalTermForConvolution(TH1F * Histo, TSpline3 * Spline, float mean, float mu, float sigma,int steps){
        float min = Histo->GetXaxis()->GetBinLowEdge(1);
        float max = Histo->GetXaxis()->GetBinUpEdge(Histo->GetNbinsX());
        TH1F * Term = new TH1F("","",steps, min, max);
        TF1 * f1 = WeightGausswithSpline(Spline,mean,mu,sigma);
        for (int i=0;i<steps;i++){
                Term->SetBinContent(i+1,f1->Eval(Term->GetXaxis()->GetBinCenter(i)));
        }
        return Term;
}

void AddInfinitesimalTermForConvolution(TH1F* Term, TH1F * Histo, TSpline3 * Spline, float mean, float mu, float sigma,int steps){
	float min = Histo->GetXaxis()->GetBinLowEdge(1);
        float max = Histo->GetXaxis()->GetBinUpEdge(Histo->GetNbinsX());
        TF1 * f1 = WeightGausswithSpline(Spline,mean,mu,sigma);
        for (int i=0;i<steps;i++){
                Term->SetBinContent(i+1,Term->GetBinContent(i+1)+f1->Eval(Term->GetXaxis()->GetBinCenter(i)));
        }
        return;
} 
TH1F * ConvolveWithGaus(TH1F * Histo, float mu, float sigma,int steps){

        float min = Histo->GetXaxis()->GetBinLowEdge(1);
        float max = Histo->GetXaxis()->GetBinUpEdge(Histo->GetNbinsX());
        int nbins =Histo->GetNbinsX();
        float original_area=Histo->Integral();
	steps = 20* nbins;
        TSpline3 * Model =  Model_Histo(Histo);
        TH1F * Term = InfinitesimalTermForConvolution(Histo,Model,Histo->GetBinCenter(1),mu,sigma,steps);

        for(int i=1;i<steps;i++)
		AddInfinitesimalTermForConvolution(Term,Histo,Model,Histo->GetBinCenter(i+1),mu,sigma,steps);

        Term -> Rebin(steps/float(nbins));
        Term -> Scale(original_area/Term->Integral());
        return Term;
}

TH1F * SimpleShiftHisto(TH1F * Histo, float mu,int steps){
        float min = Histo->GetXaxis()->GetBinLowEdge(1);
        float max = Histo->GetXaxis()->GetBinUpEdge(Histo->GetNbinsX());
        int nbins =Histo->GetNbinsX();
        float original_area=Histo->Integral();
	steps = 200* nbins;
	TH1F * Term = new TH1F("","",steps, min, max);
	float binwidth = (max-min)/steps;
	int binshift = mu/(binwidth);

	for (int i=0;i<steps;i++){ 
		if(mu>0)   Term->SetBinContent(i+1+binshift,Histo->GetBinContent(Histo->FindBin( Term->GetBinCenter(i+1))));
		else {
			if(i+1+binshift>1) {
				Term->SetBinContent(i+1+binshift,Histo->GetBinContent(Histo->FindBin( Term->GetBinCenter(i+1))));
					}
				}
	}
        Term -> Rebin(steps/float(nbins));
        Term -> Scale(original_area/Term->Integral());

//	Term = ConvolveWithGaus(Histo,mu,steps);

	for (int i=0;i<nbins;i++){
			if(Histo->GetBinContent(i+1)>0)
			 Term->SetBinError(i+1,(Histo->GetBinError(i+1)/Histo->GetBinContent(i+1))*Term->GetBinContent(i+1));
	}
       	Term->SetName(Histo->GetName());
	 return Term;

}


