#include "PlottingFunctions/DATAQualeff_Plot.h"

using namespace std;


LATcorr * LATLikelihoodDATA_TOF = new LATcorr ("LATLikDATA_TOF");
LATcorr * LATDistanceDATA_TOF   = new LATcorr ("LATDistDATA_TOF");

LATcorr * LATLikelihoodDATA_NaF = new LATcorr ("LATLikDATA_NaF");
LATcorr * LATDistanceDATA_NaF   = new LATcorr ("LATDistDATA_NaF");

LATcorr * LATLikelihoodDATA_Agl = new LATcorr ("LATLikDATA_Agl");
LATcorr * LATDistanceDATA_Agl   = new LATcorr ("LATDistDATA_Agl");




void DATAQualeff_Fill (int zona)
{

  // if (Tup.R<=Rcut[zona]) return;
   if (! (Tup.EdepL1 > 0 && Tup.EdepL1 < EdepL1beta->Eval (Tup.Beta)+0.1 && Tup.EdepL1 > EdepL1beta->Eval (Tup.Beta)-0.1 ) ) return;

   int Kbin=PRB.GetRBin (Tup.R);


   if (cmask.isFromNaF() ) { //NaF
      LATDistanceDATA_NaF  ->beforeR->Fill (Kbin,zona);
      if (Tup.Dist5D_P<6) {
         LATDistanceDATA_NaF  ->afterR ->Fill (Kbin,zona);
         LATLikelihoodDATA_NaF->beforeR->Fill (Kbin,zona);
      }
      if (Tup.Dist5D_P<6 && Likcut)
         LATLikelihoodDATA_NaF->afterR ->Fill (Kbin,zona);
   } else if (cmask.isFromAgl() ) { // Agl
      LATDistanceDATA_Agl     ->beforeR->Fill (Kbin,zona);
      if (Tup.Dist5D_P<6) {
         LATDistanceDATA_Agl  ->afterR ->Fill (Kbin,zona);
         LATLikelihoodDATA_Agl->beforeR->Fill (Kbin,zona);
      }
      if (Tup.Dist5D_P<6 && Likcut)
         LATLikelihoodDATA_Agl->afterR ->Fill (Kbin,zona);
   }
   else // ToF
   {
      LATDistanceDATA_TOF  ->beforeR->Fill (Kbin,zona);
      if (Tup.Dist5D_P<6) {
         LATDistanceDATA_TOF  ->afterR ->Fill (Kbin,zona);
         LATLikelihoodDATA_TOF->beforeR->Fill (Kbin,zona);
      }
      if (Tup.Dist5D_P<6 && Likcut)
         LATLikelihoodDATA_TOF->afterR ->Fill (Kbin,zona);
   }

   return;
}

void DATAQualeff_Write()
{
   LATLikelihoodDATA_TOF->Write();
   LATDistanceDATA_TOF  ->Write();
   LATLikelihoodDATA_NaF->Write();
   LATDistanceDATA_NaF  ->Write();
   LATLikelihoodDATA_Agl->Write();
   LATDistanceDATA_Agl  ->Write();

   return;
}




void DATAQualeff (string filename)
{
   cout<<"****************************** DATA QUALITY SEL. EFFICIENCIES **************************************"<<endl;

    cout<<"*** Reading  P1 file ****"<<endl;
        TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");


   LATcorr * LATLikelihoodDATA_TOF = new LATcorr (inputHistoFile,"LATLikDATA_TOF" );
   LATcorr * LATDistanceDATA_TOF   = new LATcorr (inputHistoFile,"LATDistDATA_TOF");

   LATcorr * LATLikelihoodDATA_NaF = new LATcorr (inputHistoFile,"LATLikDATA_NaF" );
   LATcorr * LATDistanceDATA_NaF   = new LATcorr (inputHistoFile,"LATDistDATA_NaF");

   LATcorr * LATLikelihoodDATA_Agl = new LATcorr (inputHistoFile,"LATLikDATA_Agl" );
   LATcorr * LATDistanceDATA_Agl   = new LATcorr (inputHistoFile,"LATDistDATA_Agl");


   cout<<"****************************** DATA QUALITY SEL. EFFICIENCIES **************************************"<<endl;

   LATLikelihoodDATA_TOF -> Eval_Efficiency();
   LATDistanceDATA_TOF   -> Eval_Efficiency();

   LATLikelihoodDATA_NaF -> Eval_Efficiency();
   LATDistanceDATA_NaF   -> Eval_Efficiency();

   LATLikelihoodDATA_Agl -> Eval_Efficiency();
   LATDistanceDATA_Agl   -> Eval_Efficiency();




   TH2F *LATDistDATATOF = (TH2F *) LATDistanceDATA_TOF     -> effR -> Clone();
   TH2F *LATLikDATATOF  = (TH2F *) LATLikelihoodDATA_TOF   -> effR -> Clone();

   TH2F *LATDistDATANaF = (TH2F *) LATDistanceDATA_NaF     -> effR -> Clone();
   TH2F *LATLikDATANaF  = (TH2F *) LATLikelihoodDATA_NaF   -> effR -> Clone();

   TH2F *LATDistDATAAgl = (TH2F *) LATDistanceDATA_Agl     -> effR -> Clone();
   TH2F *LATLikDATAAgl  = (TH2F *) LATLikelihoodDATA_Agl   -> effR -> Clone();



   cout<<"****************************** LAT. Eff. CORRECTION *************************************************"<<endl;

   LATLikelihoodDATA_TOF ->  Eval_LATcorr (1);
   LATDistanceDATA_TOF   ->  Eval_LATcorr (1);

   LATLikelihoodDATA_NaF ->  Eval_LATcorr (1);
   LATDistanceDATA_NaF   ->  Eval_LATcorr (1);

   LATLikelihoodDATA_Agl ->  Eval_LATcorr (1);
   LATDistanceDATA_Agl   ->  Eval_LATcorr (1);



   TH2F *LikLATcorr_TOF  =	(TH2F *) LATLikelihoodDATA_TOF  -> LATcorrR -> Clone();
   TH2F *DistLATcorr_TOF =	(TH2F *) LATDistanceDATA_TOF    -> LATcorrR -> Clone();

   TH2F *LikLATcorr_NaF  =	(TH2F *) LATLikelihoodDATA_NaF  -> LATcorrR -> Clone();
   TH2F *DistLATcorr_NaF =	(TH2F *) LATDistanceDATA_NaF    -> LATcorrR -> Clone();

   TH2F *LikLATcorr_Agl  =	(TH2F *) LATLikelihoodDATA_Agl  -> LATcorrR -> Clone();
   TH2F *DistLATcorr_Agl =	(TH2F *) LATDistanceDATA_Agl    -> LATcorrR -> Clone();


   TH1F *LikLATcorr_TOF_fit  	= (TH1F *) LATLikelihoodDATA_TOF   -> LATcorrR_fit-> Clone();
   TH1F *DistLATcorr_TOF_fit	= (TH1F *) LATDistanceDATA_TOF     -> LATcorrR_fit-> Clone();

   TH1F *LikLATcorr_NaF_fit    	= (TH1F *) LATLikelihoodDATA_NaF   -> LATcorrR_fit-> Clone();
   TH1F *DistLATcorr_NaF_fit  	= (TH1F *) LATDistanceDATA_NaF     -> LATcorrR_fit-> Clone();

   TH1F *LikLATcorr_Agl_fit   	= (TH1F *) LATLikelihoodDATA_Agl   -> LATcorrR_fit-> Clone();
   TH1F *DistLATcorr_Agl_fit  	= (TH1F *) LATDistanceDATA_Agl     -> LATcorrR_fit-> Clone();




	finalHistos.Add( LATDistDATATOF      );
       	finalHistos.Add( LATLikDATATOF      ); 
       	finalHistos.Add( LATDistDATANaF      );
       	finalHistos.Add( LATLikDATANaF       ); 
       	finalHistos.Add( LATDistDATAAgl      );
       	finalHistos.Add( LATLikDATAAgl       ); 
                            
       	finalHistos.Add( LikLATcorr_TOF      );
       	finalHistos.Add( DistLATcorr_TOF     ); 
       	finalHistos.Add( LikLATcorr_NaF      );
       	finalHistos.Add( DistLATcorr_NaF     ); 
       	finalHistos.Add( LikLATcorr_Agl      );
       	finalHistos.Add( DistLATcorr_Agl     ); 
                            
       	finalHistos.Add( LikLATcorr_TOF_fit  );
       	finalHistos.Add( DistLATcorr_TOF_fit ); 
       	finalHistos.Add( LikLATcorr_NaF_fit  );
       	finalHistos.Add( DistLATcorr_NaF_fit ); 
       	finalHistos.Add( LikLATcorr_Agl_fit  );
       	finalHistos.Add( DistLATcorr_Agl_fit ); 


	finalHistos.writeObjsInFolder("Results");

        cout<<"*** Plotting ...  ****"<<endl;

	DATAQualeff_Plot( LATDistDATATOF      ,
                          LATLikDATATOF      ,
                          LATDistDATANaF      ,
                          LATLikDATANaF       ,
                          LATDistDATAAgl      ,
                          LATLikDATAAgl       ,
                            
                          LikLATcorr_TOF      ,
                          DistLATcorr_TOF     ,
                          LikLATcorr_NaF      ,
                          DistLATcorr_NaF     ,
                          LikLATcorr_Agl      ,
                          DistLATcorr_Agl     ,
                            
                          LikLATcorr_TOF_fit  ,
                          DistLATcorr_TOF_fit ,
                          LikLATcorr_NaF_fit  ,
                          DistLATcorr_NaF_fit ,
                          LikLATcorr_Agl_fit  ,
                          DistLATcorr_Agl_fit ); 


	return;
}


