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
         if(Massa_gen<2&&Massa_gen>1.5) ((TH3*)FitTOF_Dbins -> TemplateD) -> Fill(mass,Kbin,ReturnMCGenType());
         if(Massa_gen<4&&Massa_gen>2.5) ((TH2*)FitTOF_Dbins -> TemplateHe)-> Fill(mass,Kbin);
      Kbin=ToFPB.GetBin(RUsed);      
         if(Massa_gen<1&&Massa_gen>0.5) ((TH2*)FitTOF_Pbins -> TemplateP) -> Fill(mass,Kbin,Tup.mcweight);
         if(Massa_gen<2&&Massa_gen>1.5) ((TH3*)FitTOF_Pbins -> TemplateD) -> Fill(mass,Kbin,ReturnMCGenType());
         if(Massa_gen<4&&Massa_gen>2.5) ((TH2*)FitTOF_Pbins -> TemplateHe)-> Fill(mass,Kbin);
   
      if(cmask.isFromNaF()) {//NaF
         mass = ((Tup.R/Tup.BetaRICH)*pow((1-pow(Tup.BetaRICH,2)),0.5));
         Kbin=NaFDB.GetBin(RUsed);
	    if(Massa_gen<1&&Massa_gen>0.5) ((TH2*) FitNaF_Dbins -> TemplateP) -> Fill(mass,Kbin,Tup.mcweight);
            if(Massa_gen<2&&Massa_gen>1.5) for(int mctype=0;mctype<6;mctype++) ((TH3*)FitNaF_Dbins -> TemplateD) -> Fill(mass,Kbin,mctype);
            if(Massa_gen<4&&Massa_gen>2.5) ((TH2*) FitNaF_Dbins -> TemplateHe)-> Fill(mass,Kbin);
      	 Kbin=NaFPB.GetBin(RUsed);   
            if(Massa_gen<1&&Massa_gen>0.5) ((TH2*) FitNaF_Pbins -> TemplateP) -> Fill(mass,Kbin,Tup.mcweight);
            if(Massa_gen<2&&Massa_gen>1.5) for(int mctype=0;mctype<6;mctype++) ((TH3*)FitNaF_Pbins -> TemplateD) -> Fill(mass,Kbin,mctype);
            if(Massa_gen<4&&Massa_gen>2.5) ((TH2*) FitNaF_Pbins -> TemplateHe)-> Fill(mass,Kbin);
   	}
      if(cmask.isFromAgl()) {//Agl
         mass = ((Tup.R/Tup.BetaRICH)*pow((1-pow(Tup.BetaRICH,2)),0.5));
         Kbin=AglDB.GetBin(RUsed);   
	    if(Massa_gen<1&&Massa_gen>0.5) ((TH2*) FitAgl_Dbins -> TemplateP) -> Fill(mass,Kbin,Tup.mcweight);
            if(Massa_gen<2&&Massa_gen>1.5) ((TH3*)FitAgl_Dbins -> TemplateD) -> Fill(mass,Kbin,ReturnMCGenType());
            if(Massa_gen<4&&Massa_gen>2.5) ((TH2*) FitAgl_Dbins -> TemplateHe)-> Fill(mass,Kbin);
         Kbin=AglPB.GetBin(RUsed);   
	    if(Massa_gen<1&&Massa_gen>0.5) ((TH2*) FitAgl_Pbins -> TemplateP) -> Fill(mass,Kbin,Tup.mcweight);
            if(Massa_gen<2&&Massa_gen>1.5) ((TH3*)FitAgl_Pbins -> TemplateD) -> Fill(mass,Kbin,ReturnMCGenType());
            if(Massa_gen<4&&Massa_gen>2.5) ((TH2*)FitAgl_Pbins -> TemplateHe)-> Fill(mass,Kbin);
      }
}


void DeutonsDATA_Fill(int zona)
{

   float mass = 0;
   //cuts
   if(!trgpatt.IsPhysical()) return;
   if(Tup.Beta<=0||Tup.R<=0) return;
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

TH1F * ExtractDCounts(TH3F * FitResults, string name, TemplateFIT * Template){
	TH1F * DCounts = new TH1F(name.c_str(),name.c_str(),FitResults->GetNbinsX(),0,FitResults->GetNbinsX());
	int successfulfits[FitResults->GetNbinsX()]   = {0};
	float meancounts[FitResults->GetNbinsX()]   = {0};
	float stddevcounts[FitResults->GetNbinsX()] = {0};
	for (int j=0;j<5;j++) 
		for (int i=0;i<5;i++){	
			for(int x=0;x<FitResults->GetNbinsX();x++){
						if(FitResults -> GetBinContent(x+1,j+1,i+1)>5){ successfulfits[x]++;
						meancounts[x]+=FitResults -> GetBinContent(x+1,j+1,i+1);}
		}		
	}

	for (int j=0;j<5;j++)
                for (int i=0;i<5;i++){
			for(int x=0;x<FitResults->GetNbinsX();x++){
				if(Template  -> GetFitOutcome(x)==0&&FitResults -> GetBinContent(x+1,j+1,i+1)>5)
					stddevcounts[x] += pow((FitResults -> GetBinContent(x+1,j+1,i+1) - Template  -> DCounts -> GetBinContent(x+1)) ,2);
				else if(FitResults -> GetBinContent(x+1,j+1,i+1)>5) 
					stddevcounts[x] += pow((FitResults -> GetBinContent(x+1,j+1,i+1) - meancounts[x]/successfulfits[x]),2);
			}
		
	}
		
	for(int x=0;x<FitResults->GetNbinsX();x++) {
		if(Template  -> GetFitOutcome(x)==0){
			DCounts -> SetBinContent(x+1, Template  -> DCounts -> GetBinContent(x+1));
			DCounts -> SetBinError(x+1,Template -> DCounts -> GetBinError(x+1) + pow(stddevcounts[x]/successfulfits[x],0.5));
			}
		else {
			DCounts -> SetBinContent(x+1, meancounts[x]/successfulfits[x]);
                        DCounts -> SetBinError(x+1, pow(stddevcounts[x]/successfulfits[x],0.5));
		}
		}
	return DCounts; 
}

void DeutonsTemplFits(string filename, string frac)
{

   cout<<"******************** DEUTONS TEMPlATE FITS ************************"<<endl;
   cout<<"*** Reading  P1 file ****"<<endl;
   TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");


   	
  TemplateFIT * FitTOF_Dbest= new TemplateFIT(inputHistoFile,"FitTOF_Dbins","FitTOF_Dbins");
  TemplateFIT * FitNaF_Dbest= new TemplateFIT(inputHistoFile,"FitNaF_Dbins","FitNaF_Dbins");
  TemplateFIT * FitAgl_Dbest= new TemplateFIT(inputHistoFile,"FitAgl_Dbins","FitAgl_Dbins");


  TemplateFIT * FitTOF_Dbins[5][5];
  TemplateFIT * FitNaF_Dbins[5][5];
  TemplateFIT * FitAgl_Dbins[5][5];

 for(int j=0;j<5;j++)
  for(int i=0;i<5;i++){
    FitTOF_Dbins[i][j]	= new TemplateFIT(inputHistoFile,"FitTOF_Dbins","FitTOF_Dbins");
    FitNaF_Dbins[i][j]	= new TemplateFIT(inputHistoFile,"FitNaF_Dbins","FitNaF_Dbins");
    FitAgl_Dbins[i][j]	= new TemplateFIT(inputHistoFile,"FitAgl_Dbins","FitAgl_Dbins");
   }
   

   TemplateFIT * FitTOFgeo_Dbins	= new TemplateFIT(inputHistoFile,"FitTOF_Dbins","FitTOFgeo_Dbins",11);
   TemplateFIT * FitNaFgeo_Dbins	= new TemplateFIT(inputHistoFile,"FitNaF_Dbins","FitNaFgeo_Dbins",11);
   TemplateFIT * FitAglgeo_Dbins	= new TemplateFIT(inputHistoFile,"FitAgl_Dbins","FitAglgeo_Dbins",11);

   TemplateFIT * FitTOF_Pbins	= new TemplateFIT(inputHistoFile,"FitTOF_Pbins","FitTOF_Pbins");
   TemplateFIT * FitNaF_Pbins	= new TemplateFIT(inputHistoFile,"FitNaF_Pbins","FitNaF_Pbins");
   TemplateFIT * FitAgl_Pbins	= new TemplateFIT(inputHistoFile,"FitAgl_Pbins","FitAgl_Pbins");

   TH1F * ContaminationTOF = (TH1F*) inputHistoFile->Get("Results/ContaminationTOF");
   TH1F * ContaminationNaF = (TH1F*) inputHistoFile->Get("Results/ContaminationNaF");
   TH1F * ContaminationAgl = (TH1F*) inputHistoFile->Get("Results/ContaminationAgl");	

   

   cout<<"******************** DEUTONS TEMPlATE FITS ************************"<<endl;


		   FitTOF_Dbest->SetTolerance(0.1);
		   FitNaF_Dbest->SetTolerance(0.1);
		   FitAgl_Dbest->SetTolerance(0.1);

		   FitTOF_Dbest->SetFitRange(1.1,3);
		   FitNaF_Dbest->SetFitRange(1.1,3);
		   FitAgl_Dbest->SetFitRange(1.1,3);

		   if(frac!="tot"){
			   FitTOF_Dbest         -> DisableFit(); 
			   FitNaF_Dbest         -> DisableFit(); 
			   FitAgl_Dbest         -> DisableFit();
		   }
		   else{
			   FitTOF_Dbest       ->  SetFitConstraints(ContaminationTOF,0.8,1,0.0001,0.2);
			   FitNaF_Dbest       ->  SetFitConstraints(ContaminationNaF,0.8,1,0.0001,0.2);
			   FitAgl_Dbest       ->  SetFitConstraints(ContaminationAgl,0.8,1,0.0001,0.2);
		   }

		   cout<<"TOF FITS"<<endl;
		   FitTOF_Dbest    -> TemplateFits();
		   cout<<"NaF FITS"<<endl;
		   FitNaF_Dbest    -> TemplateFits();
		   cout<<"Agl FITS"<<endl;
		   FitAgl_Dbest    -> TemplateFits();


   for(int j=0;j<5;j++)
	   for(int i=0;i<5;i++){
		   FitTOF_Dbins[i][j]->SetTolerance(0.1+0.2*j);
		   FitNaF_Dbins[i][j]->SetTolerance(0.1+0.2*j);
		   FitAgl_Dbins[i][j]->SetTolerance(0.1+0.2*j);

		   FitTOF_Dbins[i][j]->SetFitRange(0.9+0.07*i,3);
		   FitNaF_Dbins[i][j]->SetFitRange(0.8+0.1*i,3);
		   FitAgl_Dbins[i][j]->SetFitRange(0.8+0.1*i,3);

		   if(frac!="ctot"){
			   FitTOF_Dbins[i][j]         -> DisableFit(); 
			   FitNaF_Dbins[i][j]         -> DisableFit(); 
			   FitAgl_Dbins[i][j]         -> DisableFit();
		   }
		   else{
			   FitTOF_Dbins[i][j]       ->  SetFitConstraints(ContaminationTOF,0.8,1,0.0001,0.2);
			   FitNaF_Dbins[i][j]       ->  SetFitConstraints(ContaminationNaF,0.8,1,0.0001,0.2);
			   FitAgl_Dbins[i][j]       ->  SetFitConstraints(ContaminationAgl,0.8,1,0.0001,0.2);
		   }
		   FitTOF_Dbins[i][j]    -> TemplateFits();
		   FitNaF_Dbins[i][j]    -> TemplateFits();
		   FitAgl_Dbins[i][j]    -> TemplateFits();

	   }


  TH3F * FitResultsTOF = new TH3F("FitResultsTOF","FitResultsTOF",nbinsToF,0,nbinsToF,5,0,5,5,0,5);
  TH3F * FitResultsNaF = new TH3F("FitResultsNaF","FitResultsNaF",nbinsNaF,0,nbinsNaF,5,0,5,5,0,5);
  TH3F * FitResultsAgl = new TH3F("FitResultsAgl","FitResultsAgl",nbinsAgl,0,nbinsAgl,5,0,5,5,0,5);


	 for(int j=0;j<5;j++)
	   for(int i=0;i<5;i++){
		for(int x=0;x<nbinsToF;x++)
			if(FitTOF_Dbins[j][i]->GetFitOutcome(x)==0) FitResultsTOF -> SetBinContent(x+1,j+1,i+1,FitTOF_Dbins[j][i] -> DCounts ->GetBinContent(x+1));
		for(int x=0;x<nbinsNaF;x++)
			if(FitNaF_Dbins[j][i]->GetFitOutcome(x)==0) FitResultsNaF -> SetBinContent(x+1,j+1,i+1,FitNaF_Dbins[j][i] -> DCounts ->GetBinContent(x+1));
		for(int x=0;x<nbinsAgl;x++)
			if(FitAgl_Dbins[j][i]->GetFitOutcome(x)==0) FitResultsAgl -> SetBinContent(x+1,j+1,i+1,FitAgl_Dbins[j][i] -> DCounts ->GetBinContent(x+1));
	}



   FitTOFgeo_Dbins 	-> DisableFit();// 
   FitNaFgeo_Dbins 	-> DisableFit();// 
   FitAglgeo_Dbins 	-> DisableFit();// 

   FitTOF_Pbins 	-> DisableFit();// 
   FitNaF_Pbins 	-> DisableFit();// 
   FitAgl_Pbins 	-> DisableFit();// 


   /*FitTOFgeo_Dbins 	     ->  SetFitConstraints(0.8,1,0.0001,0.2,0.0001,0.0025);
   FitNaFgeo_Dbins 	     ->  SetFitConstraints(0.8,1,0.0001,0.2,0.0001,0.0015);
   FitAglgeo_Dbins 	     ->  SetFitConstraints(0.8,1,0.0001,0.2,0.0001,0.0005);

   FitTOF_Pbins 	 ->  SetFitConstraints(0.8,1,0.0001,0.02,0.00,0.0025);
   FitNaF_Pbins 	 ->  SetFitConstraints(0.8,1,0.0001,0.02,0.00,0.0015);
   FitAgl_Pbins 	->  SetFitConstraints(0.8,1,0.0001,0.02,0.00,0.0005);
*/	
   cout<<"** TOF **"<<endl;
   FitTOFgeo_Dbins -> TemplateFits();
   FitTOF_Pbins    -> TemplateFits();

   cout<<"** NaF **"<<endl;
   FitNaFgeo_Dbins -> TemplateFits();
   FitNaF_Pbins    -> TemplateFits();

   cout<<"** Agl **"<<endl;
   FitAglgeo_Dbins -> TemplateFits();
   FitAgl_Pbins 	-> TemplateFits();

   cout<<"***** TemplateFits Outcome (Mass) ******"<<endl;
   cout<<"** TOF **"<<endl;
   for(int bin =0; bin <nbinsToF; bin++)
      cout<<FitTOF_Dbins[2][2]->GetFitOutcome(bin)<<" ";
   cout<<endl;
   for(int bin =0; bin <nbinsToF; bin++)
      cout<<FitTOF_Pbins->GetFitOutcome(bin)<<" ";
   cout<<endl;
   cout<<"** NaF **"<<endl;
   for(int bin =0; bin <nbinsNaF; bin++)
      cout<<FitNaF_Dbins[2][2]->GetFitOutcome(bin)<<" ";
   cout<<endl;
   for(int bin =0; bin <nbinsNaF; bin++)
      cout<<FitNaF_Pbins->GetFitOutcome(bin)<<" ";
   cout<<endl;
   cout<<"** Agl **"<<endl;
   for(int bin =0; bin <nbinsAgl; bin++)
      cout<<FitAgl_Dbins[2][2]->GetFitOutcome(bin)<<" ";
   cout<<endl;
   for(int bin =0; bin <nbinsAgl; bin++)
      cout<<FitAgl_Pbins->GetFitOutcome(bin)<<" ";
   cout<<endl;


   FitTOFgeo_Dbins -> DCounts   -> SetName ("D_Flux_geoCounts_TOF");
   FitNaFgeo_Dbins -> DCounts   -> SetName ("D_Flux_geoCounts_NaF");
   FitAglgeo_Dbins -> DCounts   -> SetName ("D_Flux_geoCounts_Agl");

   FitTOF_Pbins -> PCounts 	-> SetName ("P_FluxCounts_TOF");
   FitNaF_Pbins -> PCounts 	-> SetName ("P_FluxCounts_NaF");
   FitAgl_Pbins -> PCounts 	-> SetName ("P_FluxCounts_Agl");


	
	TH1F * DCountsTOF 	=(TH1F*)  ExtractDCounts(FitResultsTOF,"D_FluxCounts_TOF",FitTOF_Dbest); 
	TH1F * DCountsNaF 	=(TH1F*)  ExtractDCounts(FitResultsNaF,"D_FluxCounts_NaF",FitNaF_Dbest); 
	TH1F * DCountsAgl 	=(TH1F*)  ExtractDCounts(FitResultsAgl,"D_FluxCounts_Agl",FitAgl_Dbest); 

	TH1F * DCountsgeoTOF 	=(TH1F*)  FitTOFgeo_Dbins -> DCounts ; 
	TH1F * DCountsgeoNaF	=(TH1F*)  FitNaFgeo_Dbins -> DCounts ; 
	TH1F * DCountsgeoAgl 	=(TH1F*)  FitAglgeo_Dbins -> DCounts ; 

	TH1F * PCountsTOF 	=(TH1F*)  FitTOF_Pbins -> PCounts ; 
	TH1F * PCountsNaF 	=(TH1F*)  FitNaF_Pbins -> PCounts ; 
	TH1F * PCountsAgl 	=(TH1F*)  FitAgl_Pbins -> PCounts ; 



	finalHistos.Add( FitResultsTOF 		);	
        finalHistos.Add( FitResultsNaF 	       );
        finalHistos.Add( FitResultsAgl 		);

	finalHistos.Add( DCountsTOF 		);	
        finalHistos.Add( DCountsNaF 	       );
        finalHistos.Add( DCountsAgl 		);
                                                                      
        finalHistos.Add( DCountsgeoTOF 	    	);
        finalHistos.Add( DCountsgeoNaF	          );
        finalHistos.Add( DCountsgeoAgl 	 	 );
                                                                      
        finalHistos.Add( PCountsTOF 		);
        finalHistos.Add( PCountsNaF 		);
        finalHistos.Add( PCountsAgl 		);
	finalHistos.writeObjsInFolder("Results");

        cout<<"*** Plotting ...  ****"<<endl;

	string varname = " Mass ";
	DeutonsTemplateFits_Plot(FitTOF_Dbest, 
	                         FitNaF_Dbest,
				 FitAgl_Dbest,       
                                                
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


