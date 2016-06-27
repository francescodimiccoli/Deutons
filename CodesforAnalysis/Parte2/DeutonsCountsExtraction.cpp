#include "PlottingFunctions/DeutonsCountsExtraction_Plot.h"


TemplateFIT * FitTOF_Dbins = new TemplateFIT("FitTOF_Dbins",nbinsToF,0,3,6);
TemplateFIT * FitNaF_Dbins = new TemplateFIT("FitNaF_Dbins",nbinsNaF,0,3,6);
TemplateFIT * FitAgl_Dbins = new TemplateFIT("FitAgl_Dbins",nbinsAgl,0,3,6);

TemplateFIT * FitTOFgeo_Dbins = new TemplateFIT("FitTOFgeo_Dbins",nbinsToF,0,3,6,11);
TemplateFIT * FitNaFgeo_Dbins = new TemplateFIT("FitNaFgeo_Dbins",nbinsNaF,0,3,6,11);
TemplateFIT * FitAglgeo_Dbins = new TemplateFIT("FitAglgeo_Dbins",nbinsAgl,0,3,6,11);

TemplateFIT * FitTOF_Pbins = new TemplateFIT("FitTOF_Pbins",nbinsToF,0,3,6);
TemplateFIT * FitNaF_Pbins = new TemplateFIT("FitNaF_Pbins",nbinsNaF,0,3,6);
TemplateFIT * FitAgl_Pbins = new TemplateFIT("FitAgl_Pbins",nbinsAgl,0,3,6);


void DeutonsMC_Fill()
{

   float mass = 0;
   //cuts
   if(!(Likcut&&Distcut)) return;
   //
   for(int m=0; m<nbinsToF; m++) { //TOF
      mass = ((Tup.R/Tup.Beta)*pow((1-pow(Tup.Beta,2)),0.5));
      if(RUsed>ToFDB.MomBins()[m]&&RUsed<=ToFDB.MomBins()[m+1]) {
         if(Massa_gen<1&&Massa_gen>0.5) FitTOF_Dbins -> TemplateP -> Fill(mass,m);
         if(Massa_gen<2&&Massa_gen>1.5) ((TH3*)FitTOF_Dbins -> TemplateD) -> Fill(mass,m,ReturnMCGenType());
         if(Massa_gen<4&&Massa_gen>2.5) FitTOF_Dbins -> TemplateHe-> Fill(mass,m);
      }
      if(RUsed>ToFPB.MomBins()[m]&&RUsed<=ToFPB.MomBins()[m+1]) {
         if(Massa_gen<1&&Massa_gen>0.5) FitTOF_Pbins -> TemplateP -> Fill(mass,m);
         if(Massa_gen<2&&Massa_gen>1.5) ((TH3*)FitTOF_Pbins -> TemplateD) -> Fill(mass,m,ReturnMCGenType());
         if(Massa_gen<4&&Massa_gen>2.5) FitTOF_Pbins -> TemplateHe-> Fill(mass,m);
      }
   }
   for(int m=0; m<nbinsNaF; m++) { //NaF
      if(cmask.isFromNaF()) {
         mass = ((Tup.R/Tup.BetaRICH)*pow((1-pow(Tup.BetaRICH,2)),0.5));
         if(RUsed>NaFDB.MomBins()[m]&&RUsed<=NaFDB.MomBins()[m+1]) {
            if(Massa_gen<1&&Massa_gen>0.5) FitNaF_Dbins -> TemplateP -> Fill(mass,m);
            if(Massa_gen<2&&Massa_gen>1.5) ((TH3*)FitNaF_Dbins -> TemplateD) -> Fill(mass,m,ReturnMCGenType());
            if(Massa_gen<4&&Massa_gen>2.5) FitNaF_Dbins -> TemplateHe-> Fill(mass,m);
         }
         if(RUsed>NaFPB.MomBins()[m]&&RUsed<=NaFPB.MomBins()[m+1]) {
            if(Massa_gen<1&&Massa_gen>0.5) FitNaF_Pbins -> TemplateP -> Fill(mass,m);
            if(Massa_gen<2&&Massa_gen>1.5) ((TH3*)FitNaF_Pbins -> TemplateD) -> Fill(mass,m,ReturnMCGenType());
            if(Massa_gen<4&&Massa_gen>2.5) FitNaF_Pbins -> TemplateHe-> Fill(mass,m);
         }
      }
   }
   for(int m=0; m<nbinsAgl; m++) { //Agl
      if(cmask.isFromAgl()) {
         mass = ((Tup.R/Tup.BetaRICH)*pow((1-pow(Tup.BetaRICH,2)),0.5));
         if(RUsed>AglDB.MomBins()[m]&&RUsed<=AglDB.MomBins()[m+1]) {
            if(Massa_gen<1&&Massa_gen>0.5) FitAgl_Dbins -> TemplateP -> Fill(mass,m);
            if(Massa_gen<2&&Massa_gen>1.5) ((TH3*)FitAgl_Dbins -> TemplateD) -> Fill(mass,m,ReturnMCGenType());
            if(Massa_gen<4&&Massa_gen>2.5) FitAgl_Dbins -> TemplateHe-> Fill(mass,m);
         }
         if(RUsed>AglPB.MomBins()[m]&&RUsed<=AglPB.MomBins()[m+1]) {
            if(Massa_gen<1&&Massa_gen>0.5) FitAgl_Pbins -> TemplateP -> Fill(mass,m);
            if(Massa_gen<2&&Massa_gen>1.5) ((TH3*)FitAgl_Pbins -> TemplateD) -> Fill(mass,m,ReturnMCGenType());
            if(Massa_gen<4&&Massa_gen>2.5) FitAgl_Pbins -> TemplateHe-> Fill(mass,m);
         }
      }
   }
}


void DeutonsDATA_Fill(int zona)
{

   float mass = 0;
   //cuts
   if(!(Likcut&&Distcut)) return;
   //
   for(int m=0; m<nbinsToF; m++) { //TOF
      mass = ((Tup.R/Tup.Beta)*pow((1-pow(Tup.Beta,2)),0.5));
      if(RUsed>ToFDB.MomBins()[m]&&RUsed<=ToFDB.MomBins()[m+1]) {
         if(Tup.R>1.2*Tup.Rcutoff) FitTOF_Dbins -> DATA -> Fill(mass,m);
         ((TH3*)FitTOFgeo_Dbins -> DATA) -> Fill(mass,m,zona);
      }
      if(RUsed>ToFPB.MomBins()[m]&&RUsed<=ToFPB.MomBins()[m+1]) {
         if(Tup.R>1.2*Tup.Rcutoff) FitTOF_Pbins -> DATA -> Fill(mass,m);
      }
   }
   for(int m=0; m<nbinsNaF; m++) { //NaF
      if(cmask.isFromNaF()) {
         mass = ((Tup.R/Tup.BetaRICH)*pow((1-pow(Tup.BetaRICH,2)),0.5));
         if(RUsed>NaFDB.MomBins()[m]&&RUsed<=NaFDB.MomBins()[m+1]) {
            if(Tup.R>1.2*Tup.Rcutoff) FitNaF_Dbins -> DATA -> Fill(mass,m);
            ((TH3*)FitNaFgeo_Dbins -> DATA) -> Fill(mass,m,zona);
         }
         if(RUsed>NaFPB.MomBins()[m]&&RUsed<=NaFPB.MomBins()[m+1]) {
            if(Tup.R>1.2*Tup.Rcutoff) FitNaF_Pbins -> DATA -> Fill(mass,m);
         }
      }
   }
   for(int m=0; m<nbinsAgl; m++) { //Agl
      if(cmask.isFromAgl()) {
         mass = ((Tup.R/Tup.BetaRICH)*pow((1-pow(Tup.BetaRICH,2)),0.5));
         if(RUsed>AglDB.MomBins()[m]&&RUsed<=AglDB.MomBins()[m+1]) {
            if(Tup.R>1.2*Tup.Rcutoff) FitAgl_Dbins -> DATA -> Fill(mass,m);
            ((TH3*)FitAglgeo_Dbins -> DATA) -> Fill(mass,m,zona);
         }
         if(RUsed>AglPB.MomBins()[m]&&RUsed<=AglPB.MomBins()[m+1]) {
            if(Tup.R>1.2*Tup.Rcutoff) FitAgl_Pbins -> DATA -> Fill(mass,m);
         }
      }
   }
   return;
}


void DeutonsMC_Write()
{

   FitTOF_Dbins -> Write();
   FitNaF_Dbins -> Write();
   FitAgl_Dbins -> Write();

   FitTOFgeo_Dbins -> Write();
   FitNaFgeo_Dbins -> Write();
   FitAglgeo_Dbins -> Write();


   FitTOF_Pbins -> Write();
   FitNaF_Pbins -> Write();
   FitAgl_Pbins -> Write();
   return;
}


void DeutonsTemplFits(string filename)
{

   cout<<"******************** DEUTONS TEMPlATE FITS ************************"<<endl;
   cout<<"*** Reading  P1 file ****"<<endl;
   TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");
	

   TemplateFIT * FitTOF_Dbins	= new TemplateFIT(inputHistoFile,"FitTOF_Dbins","FitTOF_Dbins");
   TemplateFIT * FitNaF_Dbins	= new TemplateFIT(inputHistoFile,"FitNaF_Dbins","FitNaF_Dbins");
   TemplateFIT * FitAgl_Dbins	= new TemplateFIT(inputHistoFile,"FitAgl_Dbins","FitAgl_Dbins");

   TemplateFIT * FitTOFgeo_Dbins	= new TemplateFIT(inputHistoFile,"FitTOF_Dbins","FitTOFgeo_Dbins",11);
   TemplateFIT * FitNaFgeo_Dbins	= new TemplateFIT(inputHistoFile,"FitNaF_Dbins","FitNaFgeo_Dbins",11);
   TemplateFIT * FitAglgeo_Dbins	= new TemplateFIT(inputHistoFile,"FitAgl_Dbins","FitAglgeo_Dbins",11);

   TemplateFIT * FitTOF_Pbins	= new TemplateFIT(inputHistoFile,"FitTOF_Pbins","FitTOF_Pbins");
   TemplateFIT * FitNaF_Pbins	= new TemplateFIT(inputHistoFile,"FitNaF_Pbins","FitNaF_Pbins");
   TemplateFIT * FitAgl_Pbins	= new TemplateFIT(inputHistoFile,"FitAgl_Pbins","FitAgl_Pbins");

   cout<<"******************** DEUTONS TEMPlATE FITS ************************"<<endl;

   FitTOF_Dbins 	-> DisableFit();//-> SetFitConstraints(0.8,1,0.00,0.2,0.0001,0.0025);
   FitNaF_Dbins 	-> DisableFit();//-> SetFitConstraints(0.8,1,0.00,0.2,0.0001,0.0015);
   FitAgl_Dbins 	-> DisableFit();//-> SetFitConstraints(0.8,1,0.00,0.2,0.0001,0.0005);
                                          
   FitTOFgeo_Dbins 	-> DisableFit();//  -> SetFitConstraints(0.8,1,0.00,0.2,0.0001,0.0025);
   FitNaFgeo_Dbins 	-> DisableFit();//  -> SetFitConstraints(0.8,1,0.00,0.2,0.0001,0.0015);
   FitAglgeo_Dbins 	-> DisableFit();//  -> SetFitConstraints(0.8,1,0.00,0.2,0.0001,0.0005);
                                          
   FitTOF_Pbins 	-> DisableFit();//-> SetFitConstraints(0.8,1,0.00001,0.02,0.0001,0.0025);
   FitNaF_Pbins 	-> DisableFit();//-> SetFitConstraints(0.8,1,0.00001,0.02,0.0005,0.0015);
   FitAgl_Pbins 	-> DisableFit();//-> SetFitConstraints(0.8,1,0.00001,0.02,0.0001,0.0005);

   cout<<"** TOF **"<<endl;
   FitTOF_Dbins 	-> TemplateFits();
   FitTOFgeo_Dbins -> TemplateFits();
   FitTOF_Pbins    -> TemplateFits();

   cout<<"** NaF **"<<endl;
   FitNaF_Dbins 	-> TemplateFits();
   FitNaFgeo_Dbins -> TemplateFits();
   FitNaF_Pbins    -> TemplateFits();

   cout<<"** Agl **"<<endl;
   FitAgl_Dbins 	-> TemplateFits();
   FitAglgeo_Dbins -> TemplateFits();
   FitAgl_Pbins 	-> TemplateFits();

   cout<<"***** TemplateFits Outcome (Mass) ******"<<endl;
   cout<<"** TOF **"<<endl;
   for(int bin =0; bin <nbinsToF; bin++)
      cout<<FitTOF_Dbins->GetFitOutcome(bin)<<" ";
   cout<<endl;
   for(int bin =0; bin <nbinsToF; bin++)
      cout<<FitTOF_Pbins->GetFitOutcome(bin)<<" ";
   cout<<endl;
   cout<<"** NaF **"<<endl;
   for(int bin =0; bin <nbinsNaF; bin++)
      cout<<FitNaF_Dbins->GetFitOutcome(bin)<<" ";
   cout<<endl;
   for(int bin =0; bin <nbinsNaF; bin++)
      cout<<FitNaF_Pbins->GetFitOutcome(bin)<<" ";
   cout<<endl;
   cout<<"** Agl **"<<endl;
   for(int bin =0; bin <nbinsAgl; bin++)
      cout<<FitAgl_Dbins->GetFitOutcome(bin)<<" ";
   cout<<endl;
   for(int bin =0; bin <nbinsAgl; bin++)
      cout<<FitAgl_Pbins->GetFitOutcome(bin)<<" ";
   cout<<endl;

   FitTOF_Dbins -> DCounts 	-> SetName ("D_FluxCounts_TOF");
   FitNaF_Dbins -> DCounts 	-> SetName ("D_FluxCounts_NaF");
   FitAgl_Dbins -> DCounts 	-> SetName ("D_FluxCounts_Agl");

   FitTOFgeo_Dbins -> DCounts   -> SetName ("D_Flux_geoCounts_TOF");
   FitNaFgeo_Dbins -> DCounts   -> SetName ("D_Flux_geoCounts_NaF");
   FitAglgeo_Dbins -> DCounts   -> SetName ("D_Flux_geoCounts_Agl");

   FitTOF_Pbins -> PCounts 	-> SetName ("P_FluxCounts_TOF");
   FitNaF_Pbins -> PCounts 	-> SetName ("P_FluxCounts_NaF");
   FitAgl_Pbins -> PCounts 	-> SetName ("P_FluxCounts_Agl");



	finalHistos.Add(FitTOF_Dbins -> DCounts 	);	
        finalHistos.Add(FitNaF_Dbins -> DCounts 	);
        finalHistos.Add(FitAgl_Dbins -> DCounts 	);
                                                     
        finalHistos.Add(FitTOFgeo_Dbins -> DCounts   );
        finalHistos.Add(FitNaFgeo_Dbins -> DCounts   );
        finalHistos.Add(FitAglgeo_Dbins -> DCounts   );
                                                     
        finalHistos.Add(FitTOF_Pbins -> PCounts 	);
        finalHistos.Add(FitNaF_Pbins -> PCounts 	);
        finalHistos.Add(FitAgl_Pbins -> PCounts 	);
	finalHistos.writeObjsInFolder("Results");

        cout<<"*** Plotting ...  ****"<<endl;

	string varname = " Mass ";
	DeutonsTemplateFits_Plot(FitTOF_Dbins, 
	                         FitNaF_Dbins ,
				 FitAgl_Dbins,       
                                                
                                 FitTOFgeo_Dbins,
                                 FitNaFgeo_Dbins,
                                 FitAglgeo_Dbins,
                                                
                                 FitTOF_Pbins,
                                 FitNaF_Pbins,
                                 FitAgl_Pbins,
				 varname 
	);	


	return;
}


