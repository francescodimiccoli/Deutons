#include <iostream>

#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TH3.h"
#include "TF2.h"
#include "TVector3.h"
#include "TMath.h"
#include "TKey.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSpline.h"
#include "TFractionFitter.h"
#include "THStack.h"
#include "TNtuple.h"
#include "TObject.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TObjArray.h"


using namespace std;

TH1F* TH1DtoTH1F(TH1D* hd) ;

TH1F* ProjectionXtoTH1F(TH2F* h2, string title, int binmin, int binmax) ;

TSpline3 * Model_Histo(TH1F * Histo);

TF1 * WeightGausswithSpline(TSpline3 * Spline, float mean,float mu, float sigma);

TH1F * InfinitesimalTermForConvolution(TH1F * Histo, TSpline3 * Spline, float mean, float mu, float sigma,int steps=5000);

TH1F * ConvolveWithGaus(TH1F * Histo, float mu, float sigma,int steps=5000);

TH1F * MultiplyWithGaus(TH1F * Original, float mean,float sigma);

float GetSigmatoConvolve(float distort_factor, TH1F * Original);

float GetSigmatoMultiply(float distort_factor, TH1F * Original);

TH1F * Distort_Histo(TH1F * Original,float distort_factor, float shift_factor=0);
	
std::vector<TH1F *> Multiple_Distortions(TH1F * Original,float distort_factor, float shift_factor, int steps=7);

 

int System_test(){

	string inputfile = "/storage/gpfs_ams/ams/users/fdimicco/Deutons/InnerTrackerAnalysis/Histos/2011_09/2011_09_tot_P1.root";

	TFile * input = TFile::Open(inputfile.c_str());
	TH2F * TemplateP   =  (TH2F *)input->Get("FitTOF_Dbins_P");
	
	TH1F * Original = ProjectionXtoTH1F(TemplateP,"Original",7,7);

//	TH1F * Term =  Distort_Histo(Original, 0.1,-0.1);
//	TH1F * Term2 = Distort_Histo(Original,-0.1, 0.1);

	std::vector<TH1F *> Collection = Multiple_Distortions(Original,0.07,0.07);

	Original->Scale(1/Original->Integral());

	TCanvas * c1 = new TCanvas("prova1");
	gPad->SetLogy();
	gPad->SetTickx();
	gPad->SetTicky();

	Original->SetLineColor(1);
	Original->SetLineWidth(4);
	Original->SetTitle("Proton MC mass peak (+/- 7\% shift and #sigma)");
	Original->GetXaxis()->SetTitle("Mass [GeV/nucl.]");
	Original->GetYaxis()->SetRangeUser(1e-6,0.25);
	Original->Draw("hist");

	for(int i=0; i<Collection.size(); i++)
		Collection[i]->Draw("hist,same");
	
	TLegend* leg =new TLegend(0.4, 0.7,0.95,0.95);
	leg->AddEntry(Original,"Original (Proton MC)", "l");
	leg->Draw("same");

	return 0;
}

TH1F* TH1DtoTH1F(TH1D* hd) {
   TH1F*  hf=(TH1F*)hd->Clone();
   hf->SetName(hd->GetTitle());
    return hf;
}


TH1F* ProjectionXtoTH1F(TH2F* h2, string title, int binmin, int binmax) {
   TH1D* hd=h2->ProjectionX(title.data(),binmin, binmax);
   TH1F* hf=TH1DtoTH1F(hd);
   return hf;
   }


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


TH1F * InfinitesimalTermForConvolution(TH1F * Histo, TSpline3 * Spline, float mean, float mu, float sigma,int steps=5000){
	float min = Histo->GetXaxis()->GetBinLowEdge(1);
	float max = Histo->GetXaxis()->GetBinUpEdge(Histo->GetNbinsX());
	TH1F * Term = new TH1F("","",steps, min, max);
	TF1 * f1 = WeightGausswithSpline(Spline,mean,mu,sigma);	
	for (int i=0;i<steps;i++){
		Term->SetBinContent(i+1,f1->Eval(Term->GetXaxis()->GetBinCenter(i)));
	}	
	return Term;
}

TH1F * ConvolveWithGaus(TH1F * Histo, float mu, float sigma,int steps=5000){
	
	float min = Histo->GetXaxis()->GetBinLowEdge(1);
        float max = Histo->GetXaxis()->GetBinUpEdge(Histo->GetNbinsX());
	int nbins =Histo->GetNbinsX();	

	TSpline3 * Model =  Model_Histo(Histo);	
	TH1F * Term = InfinitesimalTermForConvolution(Histo,Model,Histo->GetBinCenter(1),mu,sigma,steps);
	
	for(int i=1;i<steps;i++){
		if(i%1000==0) cout<<"Convolution: "<<i/float(steps)*100<<" %"<<endl;
		Term->Add(InfinitesimalTermForConvolution(Histo,Model,(Histo->GetBinCenter(i)),mu,sigma,steps));
	}
	Term -> Rebin(steps/float(nbins));
	Term -> Scale(1/Term->Integral());
	return Term;
}

TH1F * MultiplyWithGaus(TH1F * Original, float mean,float sigma){

	TH1F * Histo = (TH1F*) Original->Clone();
	int nbins =Histo->GetNbinsX();

	TF1 * f1 = new TF1("f1","gaus",0,10);
	f1->SetParameters(1,mean,sigma);
	
	for(int i=0;i<nbins;i++)
		Histo->SetBinContent(i+1,Histo->GetBinContent(i+1)*f1->Eval(Histo->GetBinCenter(i+1)));	

	Histo -> Scale(1/Histo->Integral());	
	return Histo;
}



float GetSigmatoConvolve(float distort_factor, TH1F * Original){

	float sigma_or = Original->GetStdDev();
	float sigma_f  = sigma_or+sigma_or*distort_factor;

	float sigma_c = pow(pow(sigma_f,2) - pow(sigma_or,2) , 0.5);
	
	return sigma_c;
}

float GetSigmatoMultiply(float distort_factor, TH1F * Original){

        float sigma_or = Original->GetStdDev();
        float sigma_f  = sigma_or-sigma_or*distort_factor;

        float sigma_c = 1/pow(1/pow(sigma_f,2) - pow(1/sigma_or,2) , 0.5);

        return sigma_c;
}



TH1F * Distort_Histo(TH1F * Original,float distort_factor, float shift_factor=0){
	
	TH1F * Distorted;
	if(distort_factor>=0)
		Distorted = ConvolveWithGaus(Original,0.0,GetSigmatoConvolve(distort_factor,Original));  	
	else    Distorted = MultiplyWithGaus(Original,Original->GetMean(),GetSigmatoMultiply(fabs(distort_factor),Original));
	
	TH1F * Distorted_Shifted = (TH1F*) Distorted ->Clone();
	if(shift_factor!=0) Distorted_Shifted  = ConvolveWithGaus(Distorted,shift_factor*Original->GetMean(),0.01);


	return Distorted_Shifted;

}


std::vector<TH1F *> Multiple_Distortions(TH1F * Original,float distort_factor, float shift_factor, int steps=7){

	std::vector<TH1F *> Collection;

	distort_step = 2*distort_factor/float(steps);
	shift_step = 2*shift_factor/float(steps);
	int c=52;

	for(int i=0;i<steps;i++)
		for(int j=0;j<steps;j++){
			TH1F * Distorted = Distort_Histo(Original,-distort_factor+i*distort_step,-shift_factor+j*shift_step);
			Distorted->SetLineColor(c);
			Collection.push_back(Distorted);
			cout<<c<<endl;
			c++;
		}
	return Collection;
}




