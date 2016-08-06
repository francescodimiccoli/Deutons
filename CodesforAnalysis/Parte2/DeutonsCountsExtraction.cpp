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
   if(Tup.Beta<=0||Tup.R<=0) return; 
   if(!(Likcut&&Distcut)) return;
   //
   int Kbin;
      mass = ((Tup.R/Tup.Beta)*pow((1-pow(Tup.Beta,2)),0.5));
     //TOF
      Kbin=ToFDB.GetBin(RUsed); 
         if(Massa_gen<1&&Massa_gen>0.5) ((TH2*)FitTOF_Dbins -> TemplateP) -> Fill(mass,Kbin,Tup.mcweight);
         if(Massa_gen<2&&Massa_gen>1.5) ((TH3*)FitTOF_Dbins -> TemplateD) -> Fill(mass,Kbin,ReturnMCGenType(),Tup.mcweight);
         if(Massa_gen<4&&Massa_gen>2.5) ((TH2*)FitTOF_Dbins -> TemplateHe)-> Fill(mass,Kbin);
      Kbin=ToFPB.GetBin(RUsed);      
         if(Massa_gen<1&&Massa_gen>0.5) ((TH2*)FitTOF_Pbins -> TemplateP) -> Fill(mass,Kbin,Tup.mcweight);
         if(Massa_gen<2&&Massa_gen>1.5) ((TH3*)FitTOF_Pbins -> TemplateD) -> Fill(mass,Kbin,ReturnMCGenType(),Tup.mcweight);
         if(Massa_gen<4&&Massa_gen>2.5) ((TH2*)FitTOF_Pbins -> TemplateHe)-> Fill(mass,Kbin,Tup.mcweight);
   
    //  if(cmask.isFromNaF()) {//NaF
         mass = ((Tup.R/Tup.BetaRICH)*pow((1-pow(Tup.BetaRICH,2)),0.5));
         Kbin=NaFDB.GetBin(RUsed);
	    if(Massa_gen<1&&Massa_gen>0.5)((TH2*) FitNaF_Dbins -> TemplateP) -> Fill(mass,Kbin,Tup.mcweight);
            if(Massa_gen<2&&Massa_gen>1.5) ((TH3*)FitNaF_Dbins -> TemplateD) -> Fill(mass,Kbin,ReturnMCGenType(),Tup.mcweight);
            if(Massa_gen<4&&Massa_gen>2.5)((TH2*) FitNaF_Dbins -> TemplateHe)-> Fill(mass,Kbin,Tup.mcweight);
      	 Kbin=NaFPB.GetBin(RUsed);   
            if(Massa_gen<1&&Massa_gen>0.5)((TH2*) FitNaF_Pbins -> TemplateP) -> Fill(mass,Kbin,Tup.mcweight);
            if(Massa_gen<2&&Massa_gen>1.5) ((TH3*)FitNaF_Pbins -> TemplateD) -> Fill(mass,Kbin,ReturnMCGenType(),Tup.mcweight);
            if(Massa_gen<4&&Massa_gen>2.5)((TH2*) FitNaF_Pbins -> TemplateHe)-> Fill(mass,Kbin,Tup.mcweight);
   //	}
   //   if(cmask.isFromAgl()) {//Agl
         mass = ((Tup.R/Tup.BetaRICH)*pow((1-pow(Tup.BetaRICH,2)),0.5));
         Kbin=AglDB.GetBin(RUsed);   
	    if(Massa_gen<1&&Massa_gen>0.5)((TH2*) FitAgl_Dbins -> TemplateP) -> Fill(mass,Kbin,Tup.mcweight);
            if(Massa_gen<2&&Massa_gen>1.5) ((TH3*)FitAgl_Dbins -> TemplateD) -> Fill(mass,Kbin,ReturnMCGenType(),Tup.mcweight);
            if(Massa_gen<4&&Massa_gen>2.5)((TH2*) FitAgl_Dbins -> TemplateHe)-> Fill(mass,Kbin,Tup.mcweight);
         Kbin=AglPB.GetBin(RUsed);   
	    if(Massa_gen<1&&Massa_gen>0.5)((TH2*) FitAgl_Pbins -> TemplateP) -> Fill(mass,Kbin,Tup.mcweight);
            if(Massa_gen<2&&Massa_gen>1.5) ((TH3*)FitAgl_Pbins -> TemplateD) -> Fill(mass,Kbin,ReturnMCGenType(),Tup.mcweight);
            if(Massa_gen<4&&Massa_gen>2.5) ((TH2*)FitAgl_Pbins -> TemplateHe)-> Fill(mass,Kbin,Tup.mcweight);
    //  }
}


void DeutonsDATA_Fill(int zona)
{

   float mass = 0;
   //cuts
   if(!trgpatt.IsPhysical()) return;
   if(!(Likcut&&Distcut)) return;
   //
   int Kbin;

   mass = ((Tup.R/Tup.Beta)*pow((1-pow(Tup.Beta,2)),0.5));
   Kbin=ToFDB.GetBin(RUsed);
   if(Tup.R>1.2*Tup.Rcutoff) FitTOF_Dbins -> DATA -> Fill(mass,Kbin);
   ((TH3*)FitTOFgeo_Dbins -> DATA) -> Fill(mass,Kbin,zona);
   Kbin=ToFPB.GetBin(RUsed);
   if(Tup.R>1.2*Tup.Rcutoff) FitTOF_Pbins -> DATA -> Fill(mass,Kbin);

      if(cmask.isFromNaF()) {//NaF
         mass = ((Tup.R/Tup.BetaRICH)*pow((1-pow(Tup.BetaRICH,2)),0.5));
         Kbin=NaFDB.GetBin(RUsed);
	    if(Tup.R>1.2*Tup.Rcutoff) FitNaF_Dbins -> DATA -> Fill(mass,Kbin);
            ((TH3*)FitNaFgeo_Dbins -> DATA) -> Fill(mass,Kbin,zona);
      	 Kbin=NaFPB.GetBin(RUsed);    
            if(Tup.R>1.2*Tup.Rcutoff) FitNaF_Pbins -> DATA -> Fill(mass,Kbin);
      }

      if(cmask.isFromAgl()) {//Agl
         mass = ((Tup.R/Tup.BetaRICH)*pow((1-pow(Tup.BetaRICH,2)),0.5));
         Kbin=AglDB.GetBin(RUsed); 
	   if(Tup.R>1.2*Tup.Rcutoff) FitAgl_Dbins -> DATA -> Fill(mass,Kbin);
            ((TH3*)FitAglgeo_Dbins -> DATA) -> Fill(mass,Kbin,zona);
     	  Kbin=AglPB.GetBin(RUsed); 
           if(Tup.R>1.2*Tup.Rcutoff) FitAgl_Pbins -> DATA -> Fill(mass,Kbin);
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

   FitTOF_Dbins 	-> DisableFit();// 
   FitNaF_Dbins 	-> DisableFit();// 
   FitAgl_Dbins 	-> DisableFit();// 
                                                                                              
   FitTOFgeo_Dbins 	-> DisableFit();// 
   FitNaFgeo_Dbins 	-> DisableFit();// 
   FitAglgeo_Dbins 	-> DisableFit();// 
                                                                                              
   FitTOF_Pbins 	-> DisableFit();// 
   FitNaF_Pbins 	-> DisableFit();// 
   FitAgl_Pbins 	-> DisableFit();// 

/*
   FitTOF_Dbins 	 ->  SetFitConstraints(0.8,1,0.0001,0.2,0.0001,0.0025);
   FitNaF_Dbins 	 ->  SetFitConstraints(0.8,1,0.0001,0.2,0.0001,0.0015);
   FitAgl_Dbins 	 ->  SetFitConstraints(0.8,1,0.0001,0.2,0.0001,0.0005);
                                                                                              
   FitTOFgeo_Dbins 	     ->  SetFitConstraints(0.8,1,0.0001,0.2,0.0001,0.0025);
   FitNaFgeo_Dbins 	     ->  SetFitConstraints(0.8,1,0.0001,0.2,0.0001,0.0015);
   FitAglgeo_Dbins 	     ->  SetFitConstraints(0.8,1,0.0001,0.2,0.0001,0.0005);
                                                                                              
   FitTOF_Pbins 	 ->  SetFitConstraints(0.8,1,0.0001,0.02,0.00,0.0025);
   FitNaF_Pbins 	 ->  SetFitConstraints(0.8,1,0.0001,0.02,0.00,0.0015);
   FitAgl_Pbins 	->  SetFitConstraints(0.8,1,0.0001,0.02,0.00,0.0005);
*/

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


