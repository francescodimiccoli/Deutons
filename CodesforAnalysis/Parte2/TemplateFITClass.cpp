#include "TemplateFITClass.h"

void TemplateFIT::Write()
{

   if(TemplateP  -> GetEntries() > 0)  TemplateP  -> Write();
   if(TemplateD  -> GetEntries() > 0)  TemplateD  -> Write();
   TemplateHe -> Write();

   if(DATA  -> GetEntries() > 0)  DATA  -> Write();
   return;
}

void TemplateFIT::DisableFit()
{
   TemplateFITenabled=false;
}

TH1F * TemplateFIT::Extract_Bin(TH1 * Histo, int bin,int third_dim)
{
   TH1F * Slice = new TH1F("","",Histo->GetNbinsX(),Histo->GetXaxis()->GetBinLowEdge(1),Histo->GetXaxis()->GetBinLowEdge(101));
   for(int i = 0; i< Histo->GetNbinsX(); i++)
      Slice->SetBinContent(i+1,Histo->GetBinContent(i+1,bin+1,third_dim+1));
   return Slice;
}

void TemplateFIT::SetFitConstraints(float LowP, float HighP, float LowD, float HighD,float LowHe, float HighHe)
{
   lowP   = LowP	;
   highP  = HighP	;
   lowD   = LowD	;
   highD  = HighD	;
   lowHe  = LowHe	;
   highHe = HighHe	;
   return;
}



void TemplateFIT::Do_TemplateFIT(TFit * Fit,int lat)
{
   TObjArray *Tpl;
   Tpl = new TObjArray(3);
   Tpl -> Add( Fit ->  Templ_P );
   Tpl -> Add( Fit ->  Templ_D );
   Tpl -> Add( Fit ->  Templ_He);
   if(!TemplateFITenabled) {
      Fit -> Tfit = 0;
      Fit -> Tfit_outcome = 1;
   } else {
      if(Fit -> Data -> Integral() > 500) {
         Fit -> Tfit = new TFractionFitter(Fit -> Data, Tpl ,"q");
	
	 Fit -> Tfit -> SetRangeX(Fit -> Data -> FindBin(1.1), Fit -> Data -> GetNbinsX()+1);
         Fit -> Tfit -> Constrain(0, lowP ,highP );
         Fit -> Tfit -> Constrain(1, lowD ,highD );
         Fit -> Tfit -> Constrain(2, lowHe,highHe);
         Fit -> Tfit_outcome = Fit -> Tfit -> Fit();
         for(int fit_attempt=0; fit_attempt<20; fit_attempt++) {
            if(Fit -> Tfit_outcome == 0) break;
            else {
               Fit -> Tfit -> Constrain(0, lowP+(float)fit_attempt/1000 ,highP );
               Fit -> Tfit -> Constrain(1, lowD-(float)fit_attempt/10000 ,highD-fit_attempt/10000 );
               Fit -> Tfit -> Constrain(2, lowHe-(float)fit_attempt/100000,highHe+(float)fit_attempt/100000);
               Fit -> Tfit_outcome = Fit -> Tfit -> Fit();
            }
         }
      } else {
         Fit -> Tfit = 0;
         Fit -> Tfit_outcome = 2;
      }
   }
   fits[lat].push_back(Fit);
   fits[lat][fits[lat].size()-1]->wheightP =GetFitWheights(0,fits[lat].size()-1,lat);
   fits[lat][fits[lat].size()-1]->wheightD =GetFitWheights(1,fits[lat].size()-1,lat);
   fits[lat][fits[lat].size()-1]->wheightHe=GetFitWheights(2,fits[lat].size()-1,lat);

   return;
}

double TemplateFIT::GetFitWheights(int par, int bin,int lat)
{
   if(GetFitOutcome(bin,lat)!=0)   return  1;
   double w1,e1 = 0;
   fits[lat][bin]-> Tfit ->GetResult(par,w1,e1);
   TH1F * Result = (TH1F*)fits[lat][bin] -> Tfit -> GetPlot();
   float itot= Result->Integral();
   float i1;
   if(par == 0) i1 = fits[lat][bin]-> Templ_P  ->Integral(fits[lat][bin]->  Data -> FindBin(1.1), fits[lat][bin]->  Data -> GetNbinsX()+1);
   if(par == 1) i1 = fits[lat][bin]-> Templ_D  ->Integral(fits[lat][bin]->  Data -> FindBin(1.1), fits[lat][bin]->  Data -> GetNbinsX()+1);
   if(par == 2) i1 = fits[lat][bin]-> Templ_He ->Integral(fits[lat][bin]->  Data -> FindBin(1.1), fits[lat][bin]->  Data -> GetNbinsX()+1);
   return w1*itot/i1;
}

double TemplateFIT::GetFitFraction(int par, int bin,int lat)
{
   if(GetFitOutcome(bin,lat)!=0)   return  0;
   double w1,e1 = 0;
   fits[lat][bin]-> Tfit ->GetResult(par,w1,e1);
   fits[lat][bin] -> Tfit -> GetPlot();
   return w1;

}

double TemplateFIT::GetFitErrors(int par,int bin,int lat)
{
   if(GetFitOutcome(bin,lat)!=0) return  0;

   double w1,e1=0;
   double w2,e2=0;
   double w3,e3=0;
   fits[lat][bin]-> Tfit ->GetResult(0,w1,e1);
   fits[lat][bin]-> Tfit ->GetResult(1,w2,e2);
   fits[lat][bin]-> Tfit ->GetResult(2,w3,e3);

   float Cov01=fits[lat][bin]-> Tfit->GetFitter()->GetCovarianceMatrixElement(0,1);
   float Cov02=fits[lat][bin]-> Tfit->GetFitter()->GetCovarianceMatrixElement(0,2);
   float Cov12=fits[lat][bin]-> Tfit->GetFitter()->GetCovarianceMatrixElement(1,2);

   float Sigma=pow((pow(w2*e2,2)+pow(w1*e1,2)+pow(w3*e3,2)
                    -2*Cov01*w1*w2-2*Cov02*w1*w3
                    -2*Cov12*w2*w3)/2,0.5);

   double Err = Sigma;//pow((Sigma/w2,2) + pow(Sigma/w1,2),0.5); //Fit relative error

   TH1F * ResultPlot;
   if(par == 0)	ResultPlot = GetResult_P (bin,lat);
   if(par == 1)   ResultPlot = GetResult_D (bin,lat);
   if(par == 2)   ResultPlot = GetResult_He(bin,lat);

   return Err * ResultPlot->Integral(); //Fit absolute error

}

void TemplateFIT::TemplateFits(int mc_type)
{
   int loops = 1;
   if(Geomag) loops = DATA -> GetNbinsZ();
   for(int lat = 0; lat < loops ; lat ++)
      for(int bin=0; bin<nbins ; bin++) {
         TFit * Fit = new TFit;
         Fit->Templ_P =  (TH1F *)TemplateFIT::Extract_Bin  (TemplateP, bin);
         Fit->Templ_D =  (TH1F *)TemplateFIT::Extract_Bin  (TemplateD, bin,mc_type);
         Fit->Templ_He=  (TH1F *)TemplateFIT::Extract_Bin  (TemplateHe,bin);
         if(Geomag) Fit->Data    =  (TH1F *)TemplateFIT::Extract_Bin(DATA      ,bin, lat);
         else 	   Fit->Data    =  (TH1F *)TemplateFIT::Extract_Bin(DATA      ,bin);

         TemplateFIT::Do_TemplateFIT(Fit,lat);
         TH1F * Data          = GetResult_Data(bin,lat);

         if(!Geomag) {
            PCounts -> SetBinContent(bin+1,Data->Integral()/*GetFitFraction(0,bin)*/);
            DCounts -> SetBinContent(bin+1,GetResult_D(bin)->Integral());
            PCounts -> SetBinError(bin+1,GetFitErrors(0,bin));
            DCounts -> SetBinError(bin+1,GetFitErrors(1,bin));
         }

         if(Geomag) {
            PCounts -> SetBinContent(bin+1,lat+1,Data->Integral()/*GetFitFraction(0,bin,lat)*/);
            DCounts -> SetBinContent(bin+1,lat+1,GetResult_D(bin,lat)->Integral());
            PCounts -> SetBinError(bin+1,lat+1,GetFitErrors(0,bin,lat));
            DCounts -> SetBinError(bin+1,lat+1,GetFitErrors(1,bin,lat));
         }
      }
   return;
}


void TemplateFIT::TemplateFitPlot(TVirtualPad * c, std::string var_name,int bin,int lat)
{
   c -> cd();
   gPad-> SetLogy();
   gPad-> SetGridx();
   gPad-> SetGridy();
   THStack *Stack=new THStack((var_name + "_" + to_string(lat) + "_" + to_string(bin)).c_str(),(var_name + "_" + to_string(lat) + "_" + to_string(bin)).c_str());

   TH1F *PMC  = GetResult_P   (bin,lat);
   TH1F *DMC  = GetResult_D   (bin,lat);
   TH1F *HeMC = GetResult_He  (bin,lat);
   TH1F *Data = GetResult_Data(bin,lat);
   TH1F * Result = NULL;
   if(GetFitOutcome(bin,lat)==0) {
      Result = (TH1F*)fits[lat][bin] -> Tfit -> GetPlot();
      Result ->SetLineColor(5);
      Result ->SetLineWidth(2);
   }
   PMC -> SetFillColor(2);
   DMC -> SetFillColor(4);
   HeMC-> SetFillColor(3);
   Data->SetMarkerStyle(8);
   if(GetFitOutcome(bin,lat)!=0) {
      PMC -> SetFillStyle(3001);
      DMC -> SetFillStyle(3001);
      HeMC-> SetFillStyle(3001);
   }
   if(TemplateFITenabled) {
      if(GetFitOutcome(bin,lat)==0) {
         Stack->Add(PMC);
         Stack->Add(DMC);
         Stack->Add(HeMC);
         Stack->Draw();
         Stack-> SetTitle(("Template Fit bin " + to_string(bin)).c_str());
         Stack-> GetXaxis()->SetTitle(var_name.c_str());
         Stack-> GetYaxis()->SetTitle("Counts");
         Data->Draw("epsame");
         Result->Draw("same");
      }
      if(GetFitOutcome(bin,lat)!=0) {
         PMC ->Draw();
         DMC ->Draw("same");
         HeMC->Draw("same");
         Data->Draw("epsame");
      }
   }
   return;
}
