#include <iostream>

#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TH3.h"
#include "TF2.h"
#include <vector>
#include <string>
#include <sstream>

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

TH1F * Extract_Bin(TH1 * Histo, int bin,int third_dim=0);

TSpline3 * Model_Histo(TH1F * Histo);

TF1 * WeightGausswithSpline(TSpline3 * Spline, float mean,float mu, float sigma);

TH1F * InfinitesimalTermForConvolution(TH1F * Histo, TSpline3 * Spline, float mean, float mu, float sigma,int steps=1000);

TH1F * ConvolveWithGaus(TH1F * Histo, float mu, float sigma,int steps=1000);

TH1F * MultiplyWithGaus(TH1F * Original, float mean,float sigma);

float GetSigmatoConvolve(float distort_factor, TH1F * Original);

float GetSigmatoMultiply(float distort_factor, TH1F * Original);

TH1F * Distort_Histo(TH1F * Original,float distort_factor, float shift_factor=0);
	
std::vector<std::vector<TH1F *>> Multiple_Distortions(TH1F * Original,float distort_factor, float shift_factor, int steps=5);

std::vector<std::vector<TH1F *>> CopyCollection(std::vector< std::vector<TH1F *> > Collection);

struct TFit {
   TH1F * Templ_P ;
   TH1F * Templ_D ;
   TH1F * Data;
   float wheightP,wheightD;
   TFractionFitter *Tfit;
   int Tfit_outcome;
};
 
void Do_TemplateFIT(TFit * Fit);


int System_test(){

	int STEPS=7;
	float DEF_SIGMA=0.07;
	float SHIFT=0.04;
	float BIN=7;
	std::vector<std::vector<int>> f;

	cout<<"************************ READING DATA ***************************"<<endl;

	string inputfile = "/storage/gpfs_ams/ams/users/fdimicco/Deutons/Analysis/AnalysisFiles/1314835200/Result.root";

	TFile * input = TFile::Open(inputfile.c_str());
	TH1F * TemplateP   =  (TH1F *)input->Get("TOFfits/TemplateP/TOFfits_MCP_7");
	TH1F * TemplateD   =  (TH1F *)input->Get("TOFfits/TemplateD/TOFfits_MCD_7");
	TH1F * DATA   	   =  (TH1F *)input->Get("TOFfits/Data/TOFfits_Data_7");	

	TH1F * OriginalP = (TH1F *)TemplateP ->Clone(); 
	TH1F * OriginalD = (TH1F *)TemplateD ->Clone();
	TH1F * Data      = (TH1F *)DATA      ->Clone();  
	
	cout<<"************************ DISTORTION TEST ***************************"<<endl;

		
	std::vector<std::vector<TH1F *>> Collection = Multiple_Distortions(OriginalP,DEF_SIGMA,SHIFT,STEPS);

	TCanvas * c1 = new TCanvas("Example of distortion");
	gPad->SetLogy();
	gPad->SetTickx();
	gPad->SetTicky();

	OriginalP->SetLineColor(1);
	OriginalP->SetLineWidth(4);
	OriginalP->SetTitle("Proton MC mass peak (+/- 7\% shift and #sigma)");
	OriginalP->GetXaxis()->SetTitle("Mass [GeV/nucl.]");
	OriginalP->Draw("hist");

	for(int i=0; i<Collection.size(); i++)
		for(int j=0; j<Collection[i].size(); j++)
			Collection[i][j]->Draw("hist,same");
	OriginalP->Draw("hist,same");
	
	TLegend* leg =new TLegend(0.7, 0.3,0.95,0.93);
	leg->AddEntry(OriginalP,"Original (Proton MC)", "l");
	for(int i=0; i<Collection.size(); i++)
		for(int j=0; j<Collection[i].size(); j++)
			leg->AddEntry(Collection[i][j],Collection[i][j]->GetTitle(), "l");
	leg->SetLineWidth(2);
	leg->Draw("same");
		
	 cout<<"************************ TEMPLATE FIT TEST ***************************"<<endl;

	TCanvas * c2 = new TCanvas("Template Fit");
        gPad->SetLogy();
        gPad->SetTickx();
        gPad->SetTicky();

	std::vector< std::vector<TH1F *> > CollectionFit=CopyCollection(Collection);
	

	TFit * Fits[Collection.size()][Collection[0].size()]; 
	for(int i=0; i<Collection.size(); i++)
                for(int j=0; j<Collection[i].size(); j++)
			Fits[i][j]= new TFit;
	
	for(int i=0; i<CollectionFit.size(); i++)
		for(int j=0; j<CollectionFit[i].size(); j++){
			Fits[i][j]->Templ_P = CollectionFit[i][j];
			Fits[i][j]->Templ_D = (TH1F*) OriginalD->Clone();
			Fits[i][j]->Data    = (TH1F*) Data->Clone();	
			Do_TemplateFIT(Fits[i][j]);
		}
	
	for(int i=0; i<Collection.size(); i++)
		for(int j=0; j<Collection[i].size(); j++){
			Fits[i][j]->Templ_P->SetLineColor(2);
			Fits[i][j]->Templ_P->SetLineWidth(4);
			Fits[i][j]->Templ_P->GetXaxis()->SetTitle("Mass [GeV/nucl.]");
			Fits[i][j]->Templ_P->GetYaxis()->SetRangeUser(1e-2,1e6);
			if(i==0&&j==0) Fits[i][j]->Templ_P->Draw("hist");
			else Fits[i][j]->Templ_P->Draw("hist,same");

			Fits[i][j]->Templ_D->SetLineColor(4);
			Fits[i][j]->Templ_D->SetLineWidth(4);
			Fits[i][j]->Templ_D->SetTitle("Deuton MC Template");
			Fits[i][j]->Templ_D->GetXaxis()->SetTitle("Mass [GeV/nucl.]");
			Fits[i][j]->Templ_D->Draw("hist,SAME");
		}
	
	Data->SetLineColor(1);
        Data->SetLineWidth(4);
        Data->SetTitle("Deuton MC Template");
        Data->GetXaxis()->SetTitle("Mass [GeV/nucl.]");
        Data->Draw("hist,SAME");

	cout<<"************************ PLOTTING RESULTS ***************************"<<endl;

	TH2F * DCounts = new TH2F("Deuteron Counts","Deuteron Counts",STEPS,-DEF_SIGMA*100,DEF_SIGMA*100,STEPS,-SHIFT*100,SHIFT*100);
	
	for(int i=0; i<Collection.size(); i++)
                for(int j=0; j<Collection[i].size(); j++){
			DCounts->SetBinContent(i+1,j+1,Fits[i][j]->Templ_D->Integral()/Fits[STEPS/2+1][STEPS/2+1]->Templ_D->Integral());
	}

	TH2F * Chisquare = new TH2F("T. Fit #chi^2","T. Fit #chi^2",STEPS,-DEF_SIGMA*100,DEF_SIGMA*100,STEPS,-SHIFT*100,SHIFT*100);
	
	float dev=0;
	int goodfits=0;
	for(int i=0; i<Collection.size(); i++)
                for(int j=0; j<Collection[i].size(); j++){
			if(Fits[i][j]->Tfit->GetChisquare()>0){	
				Chisquare->SetBinContent(i+1,j+1,Fits[i][j]->Tfit->GetChisquare()/50.);
				if(Fits[i][j]->Tfit->GetChisquare()/50<200) {
					dev+=pow((Fits[i][j]->Templ_D->Integral() - Fits[STEPS/2+1][STEPS/2+1]->Templ_D->Integral()),2);
					goodfits++;
				}
			}
			else Chisquare->SetBinContent(i+1,j+1,25000/50.);
		}

	dev/=goodfits;
        dev=pow(dev,0.5);
        dev/=Fits[STEPS/2+1][STEPS/2+1]->Templ_D->Integral();
        cout<<"std dev: "<<dev<<endl;

	TH1F * Error = new TH1F("Error","Error",STEPS*STEPS,0.6*Fits[STEPS/2+1][STEPS/2+1]->Templ_D->Integral(),1.4*Fits[STEPS/2+1][STEPS/2+1]->Templ_D->Integral());

	for(int i=0; i<Collection.size(); i++)
                for(int j=0; j<Collection[i].size(); j++){
			Error->Fill(Fits[i][j]->Templ_D->Integral(),1/Chisquare->GetBinContent(i+1,j+1));
	}
	float S_Error=Error->GetStdDev()/Fits[STEPS/2+1][STEPS/2+1]->Templ_D->Integral();
	string s_err = "1";//to_string (S_Error*100);
        s_err.erase ( s_err.find_last_not_of('0') + 1, std::string::npos );
	
	TCanvas * c3 = new TCanvas("Results: Deuteron Counts");		
	gPad->SetTickx();
        gPad->SetTicky();

	DCounts->SetTitle(("Deuteron Count Deviation from "+"1"/*+to_string((int) Fits[STEPS/2+1][STEPS/2+1]->Templ_D->Integral())*/).c_str() );
	DCounts->GetXaxis()->SetTitle("#sigma deformation (%)");
	DCounts->GetYaxis()->SetTitle("Peak shift (%)");
	DCounts->GetZaxis()->SetRangeUser(0.9,1.2);
	DCounts->Draw("colz");

	TCanvas * c4 = new TCanvas("Results: Fit #chi^2");		
	gPad->SetTickx();
        gPad->SetTicky();
	gPad->SetLogz();

	Chisquare->SetTitle("T. Fits reduced #chi^2");
	Chisquare->GetXaxis()->SetTitle("#sigma deformation (%)");
	Chisquare->GetYaxis()->SetTitle("Peak shift (%)");
	Chisquare->Draw("colz");


	TCanvas * c5 = new TCanvas("Results: System. error");         
        gPad->SetTickx();
        gPad->SetTicky();


	Error->SetTitle(("Syst. error estimation: "+s_err+"%").c_str());
	Error->GetXaxis()->SetTitle("Weighted counts");
	Error->SetLineWidth(4);
	Error->SetLineColor(2);
	Error->SetFillColor(2);
	Error->SetFillStyle(3001);
	Error->Draw("hist");
	

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

TH1F * Extract_Bin(TH1 * Histo, int bin,int third_dim)
{
	TH1F* Slice;
	Slice = new TH1F("","",Histo->GetNbinsX(),Histo->GetXaxis()->GetBinLowEdge(1),Histo->GetXaxis()->GetBinLowEdge(Histo->GetNbinsX()+1));
	for(int i = 0; i< Histo->GetNbinsX(); i++)
		Slice->SetBinContent(i+1,Histo->GetBinContent(i+1,bin+1,third_dim+1));
	return Slice;
}


TSpline3 * Model_Histo(TH1F * Histo){
	int nbins =Histo->GetNbinsX();
	double X[1000];
	double Y[1000];

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

TH1F * ConvolveWithGaus(TH1F * Histo, float mu, float sigma,int steps){
	
	float min = Histo->GetXaxis()->GetBinLowEdge(1);
        float max = Histo->GetXaxis()->GetBinUpEdge(Histo->GetNbinsX());
	int nbins =Histo->GetNbinsX();	
	float original_area=Histo->Integral();

	TSpline3 * Model =  Model_Histo(Histo);	
	TH1F * Term = InfinitesimalTermForConvolution(Histo,Model,Histo->GetBinCenter(1),mu,sigma,steps);
	
	for(int i=1;i<steps;i++)
		Term->Add(InfinitesimalTermForConvolution(Histo,Model,(Histo->GetBinCenter(i)),mu,sigma,steps));
	
	Term -> Rebin(steps/float(nbins));
	Term -> Scale(original_area/Term->Integral());
	return Term;
}

TH1F * MultiplyWithGaus(TH1F * Original, float mean,float sigma){

	TH1F * Histo = (TH1F*) Original->Clone();
	int nbins =Histo->GetNbinsX();
	float original_area=Histo->Integral();

	TF1 * f1 = new TF1("f1","gaus",0,10);
	f1->SetParameters(1,mean,sigma);
	
	for(int i=0;i<nbins;i++)
		Histo->SetBinContent(i+1,Histo->GetBinContent(i+1)*f1->Eval(Histo->GetBinCenter(i+1)));	

	Histo -> Scale(original_area/Histo->Integral());	
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



TH1F * Distort_Histo(TH1F * Original,float distort_factor, float shift_factor){
	
	TH1F * Distorted;
	
	if(distort_factor>0)        Distorted = ConvolveWithGaus(Original,0.0,GetSigmatoConvolve(distort_factor,Original));  	
	else  if(distort_factor<0)  Distorted = MultiplyWithGaus(Original,Original->GetMean(),GetSigmatoMultiply(fabs(distort_factor),Original));
	else 			    Distorted = (TH1F*) Original ->Clone(); 	 

	TH1F * Distorted_Shifted = (TH1F*) Distorted ->Clone();
	if(fabs(shift_factor-0)>0.005) Distorted_Shifted  = ConvolveWithGaus(Distorted,shift_factor*Original->GetMean(),0.01); 
	cout<<"shift: "<<shift_factor<<" sigma: "<<distort_factor<<endl;

	return Distorted_Shifted;

}


std::vector< std::vector<TH1F *> > Multiple_Distortions(TH1F * Original,float distort_factor, float shift_factor, int steps){

	std::vector<std::vector<TH1F *>> Collection;
	gStyle->SetPalette(1);
	float distort_step = 2*distort_factor/float(steps-1);
	float shift_step = 2*shift_factor/float(steps-1);
	int c=0;
	int color_step = (50)/steps/steps;
	for(int i=0;i<steps;i++){
		Collection.push_back(std::vector<TH1F *>());
		for(int j=0;j<steps;j++){
			TH1F * Distorted = Distort_Histo(Original,-distort_factor+i*distort_step,-shift_factor+j*shift_step);
			
			Distorted->SetLineColor(gStyle->GetColorPalette(0+c));
			string str1 = "0";//to_string ((-shift_factor+j*shift_step)*100);
			string str2 = "1";//to_string ((-distort_factor+i*distort_step)*100);
	  		str1.erase ( str1.find_last_not_of('0') + 1, std::string::npos );			
			str2.erase ( str2.find_last_not_of('0') + 1, std::string::npos );
			Distorted->SetTitle(("Shift: " +str1 + "0% , Sigma: " + str2+"0%").c_str());
			Distorted->SetLineWidth(2);
			Collection[i].push_back(Distorted);
			c+=color_step;
		}
	}
	return Collection;
}


void Do_TemplateFIT(TFit * Fit){
	 TObjArray *Tpl;
	 Tpl = new TObjArray(2);
   	 Tpl -> Add( Fit ->  Templ_P );
         Tpl -> Add( Fit ->  Templ_D );
	 Fit -> Tfit = new TFractionFitter(Fit -> Data, Tpl ,"q");	
	 float highP=1;
	 float lowP=0.6;
	 float highD=0.25;
	 float lowD=0.00001;
	 float min=0.78;
	 float max=2.5;
	 Fit -> Tfit -> SetRangeX(Fit -> Data -> FindBin(min), Fit -> Data -> FindBin(max));
	 Fit -> Tfit -> Constrain(0, lowP ,highP );
         Fit -> Tfit -> Constrain(1, lowD ,highD );

	 Fit -> Tfit_outcome = Fit -> Tfit -> Fit();

	 for(int fit_attempt=0; fit_attempt<20; fit_attempt++) {
            cout<<fit_attempt<<endl;
	    if(Fit -> Tfit_outcome == 0) break;
            else {
               Fit -> Tfit -> Constrain(0, lowP+(float)fit_attempt/1000 ,highP );
               Fit -> Tfit -> Constrain(1, lowD-(float)fit_attempt/10000 ,highD-fit_attempt/10000 );
               Fit -> Tfit_outcome = Fit -> Tfit -> Fit();
            }
         }
	 
	 if(Fit -> Tfit_outcome==0){
	 	TH1F * Result = (TH1F *) Fit-> Tfit -> GetPlot(); 	
		float itot= Result->Integral();
		double w1,e1 = 0;
		double w2,e2 = 0;
   		Fit -> Tfit ->GetResult(0,w1,e1);
		Fit -> Tfit ->GetResult(1,w2,e2);
		float i1 = Fit-> Templ_P  ->Integral(Fit->Templ_P -> FindBin(min), Fit->Templ_P -> FindBin(max));
		float i2 = Fit-> Templ_D  ->Integral(Fit->Templ_D -> FindBin(min), Fit->Templ_D -> FindBin(max));
		Fit ->wheightP= w1*itot/i1; 
		Fit ->wheightD= w2*itot/i2;
		cout<<w1<<" "<<w2<<endl; 
		
		Fit ->  Templ_P  -> Scale(Fit ->wheightD);
		Fit ->  Templ_D  -> Scale(Fit ->wheightD);	
	 }
	else{
		Fit ->wheightP= 0;
                Fit ->wheightD= 0;
	}
	
	
	return;
}


std::vector< std::vector<TH1F *> > CopyCollection(std::vector< std::vector<TH1F *> > Collection){
	std::vector< std::vector<TH1F *> > CopyCollection;
	for(int i=0; i<Collection.size(); i++){
		CopyCollection.push_back(std::vector<TH1F *>());
		for(int j=0; j<Collection[i].size(); j++){		
			CopyCollection[i].push_back((TH1F*)Collection[i][j]->Clone());
		}
	}
	return CopyCollection;
}
