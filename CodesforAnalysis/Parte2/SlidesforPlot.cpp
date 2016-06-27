#include "PlottingFunctions/SlidesforPlot_Plot.h"

using namespace std;

TH2F *RvsBetaTOF_P=new TH2F ("RvsBetaTOF_P","RvsBetaTOF_P",500,0,6,500,0.4,1);
TH2F *RvsBetaNaF_P=new TH2F ("RvsBetaNaF_P","RvsBetaNaF_P",500,1,10,500,0.75,1);
TH2F *RvsBetaAgl_P=new TH2F ("RvsBetaAgl_P","RvsBetaAgl_P",500,3,15,500,0.95,1);
TH2F *RvsBetaTOF_D=new TH2F ("RvsBetaTOF_D","RvsBetaTOF_D",500,0,6,500,0.4,1);
TH2F *RvsBetaNaF_D=new TH2F ("RvsBetaNaF_D","RvsBetaNaF_D",500,1,10,500,0.75,1);
TH2F *RvsBetaAgl_D=new TH2F ("RvsBetaAgl_D","RvsBetaAgl_D",500,3,15,500,0.95,1);
TH2F *RvsBetaTOF_He=new TH2F ("RvsBetaTOF_He","RvsBetaTOF_He",500,0,6,500,0.4,1);
TH2F *RvsBetaNaF_He=new TH2F ("RvsBetaNaF_He","RvsBetaNaF_He",500,1,10,500,0.75,1);
TH2F *RvsBetaAgl_He=new TH2F ("RvsBetaAgl_He","RvsBetaAgl_He",500,3,15,500,0.95,1);

TH1F *MassTOF_P=new TH1F ("MassTOF_P","MassTOF_P",500,0,4.5);
TH1F *MassTOF_D=new TH1F ("MassTOF_D","MassTOF_D",500,0,4.5);
TH1F *MassNaF_P=new TH1F ("MassNaF_P","MassNaF_P",500,0,4.5);
TH1F *MassNaF_D=new TH1F ("MassNaF_D","MassNaF_D",500,0,4.5);
TH1F *MassAgl_P=new TH1F ("MassAgl_P","MassAgl_P",500,0,4.5);
TH1F *MassAgl_D=new TH1F ("MassAgl_D","MassAgl_D",500,0,4.5);

TH1F *MassTOF_PQ=new TH1F ("MassTOF_PQ","MassTOF_PQ",500,0,4.5);
TH1F *MassTOF_DQ=new TH1F ("MassTOF_DQ","MassTOF_DQ",500,0,4.5);
TH1F *MassNaF_PQ=new TH1F ("MassNaF_PQ","MassNaF_PQ",500,0,4.5);
TH1F *MassNaF_DQ=new TH1F ("MassNaF_DQ","MassNaF_DQ",500,0,4.5);
TH1F *MassAgl_PQ=new TH1F ("MassAgl_PQ","MassAgl_PQ",500,0,4.5);
TH1F *MassAgl_DQ=new TH1F ("MassAgl_DQ","MassAgl_DQ",500,0,4.5);


TH2F *EdepUTOFvsR_P=new TH2F ("EdepUTOFvsR_P","EdepUTOFvsR_P",500,0,10,500,0,40);
TH2F *EdepUTOFvsR_D=new TH2F ("EdepUTOFvsR_D","EdepUTOFvsR_D",500,0,10,500,0,40);
TH2F *EdepUTOFvsR_He=new TH2F ("EdepUTOFvsR_He","EdepUTOFvsR_He",500,0,10,500,0,40);
TH2F *EdepLTOFvsR_P=new TH2F ("EdepLTOFvsR_P","EdepLTOFvsR_P",500,0,10,500,0,40);
TH2F *EdepLTOFvsR_D=new TH2F ("EdepLTOFvsR_D","EdepLTOFvsR_D",500,0,10,500,0,40);
TH2F *EdepLTOFvsR_He=new TH2F ("EdepLTOFvsR_He","EdepLTOFvsR_He",500,0,10,500,0,40);
TH2F *EdepTrackvsR_P=new TH2F ("EdepTrackvsR_P","EdepTrackvsR_P",500,0,10,500,0,4);
TH2F *EdepTrackvsR_D=new TH2F ("EdepTrackvsR_D","EdepTrackvsR_D",500,0,10,500,0,4);
TH2F *EdepTrackvsR_He=new TH2F ("EdepTrackvsR_He","EdepTrackvsR_He",500,0,10,500,0,4);


TH2F *RvsBetaTOF=new TH2F ("RvsBetaTOF","RvsBetaTOF",500,0,6,500,0.4,1);
TH2F *RvsBetaNaF=new TH2F ("RvsBetaNaF","RvsBetaNaF",500,1,10,500,0.75,1);
TH2F *RvsBetaAgl=new TH2F ("RvsBetaAgl","RvsBetaAgl",500,3,15,500,0.95,1);


TH1F *MassTOF=new TH1F ("MassTOF","MassTOF",500,0,4.5);
TH1F *MassNaF=new TH1F ("MassNaF","MassNaF",500,0,4.5);
TH1F *MassAgl=new TH1F ("MassAgl","MassAgl",500,0,4.5);

TH1F *MassTOFQ=new TH1F ("MassTOFQ","MassTOFQ",500,0,4.5);
TH1F *MassNaFQ=new TH1F ("MassNaFQ","MassNaFQ",500,0,4.5);
TH1F *MassAglQ=new TH1F ("MassAglQ","MassAglQ",500,0,4.5);


TH2F *LikvsDistTOF_P=new TH2F ("LikvsDistTOF_P","LikvsDistTOF_P",500,0,6,500,-1,1);
TH2F *LikvsDistNaF_P=new TH2F ("LikvsDistNaF_P","LikvsDistNaF_P",500,0,6,500,-1,1);
TH2F *LikvsDistAgl_P=new TH2F ("LikvsDistAgl_P","LikvsDistAgl_P",500,0,6,500,-1,1);
TH2F *LikvsDistTOF_D=new TH2F ("LikvsDistTOF_D","LikvsDistTOF_D",500,0,6,500,-1,1);
TH2F *LikvsDistNaF_D=new TH2F ("LikvsDistNaF_D","LikvsDistNaF_D",500,0,6,500,-1,1);
TH2F *LikvsDistAgl_D=new TH2F ("LikvsDistAgl_D","LikvsDistAgl_D",500,0,6,500,-1,1);

TH2F *RvsDistTOF_P=new TH2F ("RvsDistTOF_P","RvsDistTOF_P",500,0,6,500,-1,1);
TH2F *RvsDistNaF_P=new TH2F ("RvsDistNaF_P","RvsDistNaF_P",500,1,10,500,-1,1);
TH2F *RvsDistAgl_P=new TH2F ("RvsDistAgl_P","RvsDistAgl_P",500,2,19,500,-1,1);
TH2F *RvsDistTOF_D=new TH2F ("RvsDistTOF_D","RvsDistTOF_D",500,0,6,500,-1,1);
TH2F *RvsDistNaF_D=new TH2F ("RvsDistNaF_D","RvsDistNaF_D",500,1,10,500,-1,1);
TH2F *RvsDistAgl_D=new TH2F ("RvsDistAgl_D","RvsDistAgl_D",500,2,19,500,-1,1);
TH2F *RvsDistTOF_He=new TH2F ("RvsDistTOF_He","RvsDistTOF_He",500,0,6,500,-1,1);
TH2F *RvsDistNaF_He=new TH2F ("RvsDistNaF_He","RvsDistNaF_He",500,1,10,500,-1,1);
TH2F *RvsDistAgl_He=new TH2F ("RvsDistAgl_He","RvsDistAgl_He",500,2,19,500,-1,1);



TH1F *DistTOF_P=new TH1F ("DistTOF_P","DistTOF_P",500,-1,1);
TH1F *DistNaF_P=new TH1F ("DistNaF_P","DistNaF_P",500,-1,1);
TH1F *DistAgl_P=new TH1F ("DistAgl_P","DistAgl_P",500,-1,1);
TH1F *DistTOF_D=new TH1F ("DistTOF_D","DistTOF_D",500,-1,1);
TH1F *DistNaF_D=new TH1F ("DistNaF_D","DistNaF_D",500,-1,1);
TH1F *DistAgl_D=new TH1F ("DistAgl_D","DistAgl_D",500,-1,1);
TH1F *DistTOF_He=new TH1F ("DistTOF_He","DistTOF_He",500,-1,1);
TH1F *DistNaF_He=new TH1F ("DistNaF_He","DistNaF_He",500,-1,1);
TH1F *DistAgl_He=new TH1F ("DistAgl_He","DistAgl_He",500,-1,1);


TH2F *sigmagen_bad =new TH2F ("sigmagen_bad","sigmagen_bad",500,0,30,500,0,30);



void SlidesforPlot_Fill ()
{

   float Betagen= pow (pow (Tup.Momento_gen/Massa_gen,2) / (1+pow (Tup.Momento_gen/Massa_gen,2) ),0.5);
   if (Herejcut) {
      if (Massa_gen<1&&Massa_gen>0.5) {
         RvsBetaTOF_P->Fill (Tup.R,Tup.Beta);
         if ( cmask.isFromNaF() ) RvsBetaNaF_P->Fill (Tup.R,Tup.BetaRICH);
         if ( cmask.isFromAgl() ) RvsBetaAgl_P->Fill (Tup.R,Tup.BetaRICH);
         if (Betastrongcut&&Tup.BetaRICH<0) MassTOF_P->Fill ( (Tup.R/Tup.Beta) *pow (1-pow (Tup.Beta,2),0.5) );
         if (Betastrongcut&& cmask.isFromNaF() ) MassNaF_P->Fill ( (Tup.R/Tup.BetaRICH) *pow (1-pow (Tup.BetaRICH,2),0.5) );
         if (Betastrongcut&& cmask.isFromAgl() ) MassAgl_P->Fill ( (Tup.R/Tup.BetaRICH) *pow (1-pow (Tup.BetaRICH,2),0.5) );
         if (Likcut&&Distcut) {
            if (Betastrongcut&&Tup.BetaRICH<0) MassTOF_PQ->Fill ( (Tup.R/Tup.Beta) *pow (1-pow (Tup.Beta,2),0.5) );
            if (Betastrongcut&& cmask.isFromNaF() ) MassNaF_PQ->Fill ( (Tup.R/Tup.BetaRICH) *pow (1-pow (Tup.BetaRICH,2),0.5) );
            if (Betastrongcut&& cmask.isFromAgl() ) MassAgl_PQ->Fill ( (Tup.R/Tup.BetaRICH) *pow (1-pow (Tup.BetaRICH,2),0.5) );
         }
         if (Betastrongcut&&Tup.BetaRICH<0&& (Tup.R/Tup.Beta) *pow (1-pow (Tup.Beta,2),0.5) >2)
            sigmagen_bad->Fill (fabs (Tup.R-Tup.Momento_gen) / (pow (Tup.Momento_gen,2) *Rig->Eval (Tup.Momento_gen) ),fabs (Tup.Beta-Betagen) / (pow (Tup.Beta,2) *beta->Eval (Tup.Beta) ) );
      }
      if (Massa_gen<2&&Massa_gen>1.5) {
         RvsBetaTOF_D->Fill (Tup.R,Tup.Beta);
         if ( cmask.isFromNaF() ) RvsBetaNaF_D->Fill (Tup.R,Tup.BetaRICH);
         if ( cmask.isFromAgl() ) RvsBetaAgl_D->Fill (Tup.R,Tup.BetaRICH);
         if (Betastrongcut&&Tup.BetaRICH<0) MassTOF_D->Fill ( (Tup.R/Tup.Beta) *pow (1-pow (Tup.Beta,2),0.5) );
         if (Betastrongcut&& cmask.isFromNaF() ) MassNaF_D->Fill ( (Tup.R/Tup.BetaRICH) *pow (1-pow (Tup.BetaRICH,2),0.5) );
         if (Betastrongcut&& cmask.isFromAgl() ) MassAgl_D->Fill ( (Tup.R/Tup.BetaRICH) *pow (1-pow (Tup.BetaRICH,2),0.5) );
         if (Likcut&&Distcut) {
            if (Betastrongcut&&Tup.BetaRICH<0) MassTOF_DQ->Fill ( (Tup.R/Tup.Beta) *pow (1-pow (Tup.Beta,2),0.5) );
            if (Betastrongcut&& cmask.isFromNaF() ) MassNaF_DQ->Fill ( (Tup.R/Tup.BetaRICH) *pow (1-pow (Tup.BetaRICH,2),0.5) );
            if (Betastrongcut&& cmask.isFromAgl() ) MassAgl_DQ->Fill ( (Tup.R/Tup.BetaRICH) *pow (1-pow (Tup.BetaRICH,2),0.5) );
         }
      }
      if (Massa_gen<4.5&&Massa_gen>2.5) {
         if (Tup.BetaRICH<0) RvsBetaTOF_He->Fill (Tup.R,Tup.Beta);
         if ( cmask.isFromNaF() ) RvsBetaNaF_He->Fill (Tup.R,Tup.BetaRICH);
         if ( cmask.isFromAgl() ) RvsBetaAgl_He->Fill (Tup.R,Tup.BetaRICH);

      }
   }
   if (Massa_gen<1&&Massa_gen>0.5) {
      EdepUTOFvsR_P->Fill (Tup.R,Tup.EdepTOFU);
      EdepLTOFvsR_P->Fill (Tup.R,Tup.EdepTOFD);
      EdepTrackvsR_P->Fill (Tup.R,Tup.EdepTrack);
      if (Betastrongcut&&Tup.BetaRICH<0) DistTOF_P->Fill ( (Tup.Dist5D_P-Tup.Dist5D) / (Tup.Dist5D_P+Tup.Dist5D) );
      if (Betastrongcut&& cmask.isFromNaF() ) DistNaF_P->Fill ( (Tup.Dist5D_P-Tup.Dist5D) / (Tup.Dist5D_P+Tup.Dist5D) );
      if (Betastrongcut&& cmask.isFromAgl() ) DistAgl_P->Fill ( (Tup.Dist5D_P-Tup.Dist5D) / (Tup.Dist5D_P+Tup.Dist5D) );

      if (Tup.BetaRICH<0) RvsDistTOF_P->Fill (Tup.R, (Tup.Dist5D_P-Tup.Dist5D) / (Tup.Dist5D_P+Tup.Dist5D) );
      if ( cmask.isFromNaF() ) RvsDistNaF_P->Fill (Tup.R, (Tup.Dist5D_P-Tup.Dist5D) / (Tup.Dist5D_P+Tup.Dist5D) );
      if ( cmask.isFromAgl() ) RvsDistAgl_P->Fill (Tup.R, (Tup.Dist5D_P-Tup.Dist5D) / (Tup.Dist5D_P+Tup.Dist5D) );

   }
   if (Massa_gen<2&&Massa_gen>1.5) {
      EdepUTOFvsR_D->Fill (Tup.R,Tup.EdepTOFU);
      EdepLTOFvsR_D->Fill (Tup.R,Tup.EdepTOFD);
      EdepTrackvsR_D->Fill (Tup.R,Tup.EdepTrack);
      if (Betastrongcut&&Tup.BetaRICH<0) DistTOF_D->Fill ( (Tup.Dist5D_P-Tup.Dist5D) / (Tup.Dist5D_P+Tup.Dist5D) );
      if (Betastrongcut&& cmask.isFromNaF() ) DistNaF_D->Fill ( (Tup.Dist5D_P-Tup.Dist5D) / (Tup.Dist5D_P+Tup.Dist5D) );
      if (Betastrongcut&& cmask.isFromAgl() ) DistAgl_D->Fill ( (Tup.Dist5D_P-Tup.Dist5D) / (Tup.Dist5D_P+Tup.Dist5D) );

      if (Tup.BetaRICH<0) RvsDistTOF_D->Fill (Tup.R, (Tup.Dist5D_P-Tup.Dist5D) / (Tup.Dist5D_P+Tup.Dist5D) );
      if ( cmask.isFromNaF() ) RvsDistNaF_D->Fill (Tup.R, (Tup.Dist5D_P-Tup.Dist5D) / (Tup.Dist5D_P+Tup.Dist5D) );
      if ( cmask.isFromAgl() ) RvsDistAgl_D->Fill (Tup.R, (Tup.Dist5D_P-Tup.Dist5D) / (Tup.Dist5D_P+Tup.Dist5D) );

   }
   if (Massa_gen<4.5&&Massa_gen>2.5) {
      EdepUTOFvsR_He->Fill (Tup.R,Tup.EdepTOFU);
      EdepLTOFvsR_He->Fill (Tup.R,Tup.EdepTOFD);
      EdepTrackvsR_He->Fill (Tup.R,Tup.EdepTrack);
      if (Betastrongcut&&Tup.BetaRICH<0&&Tup.R>1) DistTOF_He->Fill ( (Tup.Dist5D_P-Tup.Dist5D) / (Tup.Dist5D_P+Tup.Dist5D) );
      if (Betastrongcut&& cmask.isFromNaF()&&Tup.R>1) DistNaF_He->Fill ( (Tup.Dist5D_P-Tup.Dist5D) / (Tup.Dist5D_P+Tup.Dist5D) );
      if (Betastrongcut&& cmask.isFromAgl()&&Tup.R>1) DistAgl_He->Fill ( (Tup.Dist5D_P-Tup.Dist5D) / (Tup.Dist5D_P+Tup.Dist5D) );

      if (Tup.BetaRICH<0) RvsDistTOF_He->Fill (Tup.R, (Tup.Dist5D_P-Tup.Dist5D) / (Tup.Dist5D_P+Tup.Dist5D) );
      if ( cmask.isFromNaF() ) RvsDistNaF_He->Fill (Tup.R, (Tup.Dist5D_P-Tup.Dist5D) / (Tup.Dist5D_P+Tup.Dist5D) );
      if ( cmask.isFromAgl() ) RvsDistAgl_He->Fill (Tup.R, (Tup.Dist5D_P-Tup.Dist5D) / (Tup.Dist5D_P+Tup.Dist5D) );
   }


}


void SlidesforPlot_D_Fill ()
{


   if (Herejcut&&Tup.Latitude>0.8) {
      RvsBetaTOF->Fill (Tup.R,Tup.Beta);
      if ( cmask.isFromNaF() ) RvsBetaNaF->Fill (Tup.R,Tup.BetaRICH);
      if ( cmask.isFromAgl() ) RvsBetaAgl->Fill (Tup.R,Tup.BetaRICH);
      if (Betastrongcut&&Tup.BetaRICH<0) MassTOF->Fill ( (Tup.R/Tup.Beta) *pow (1-pow (Tup.Beta,2),0.5) );
      if (Betastrongcut&& cmask.isFromNaF() ) MassNaF->Fill ( (Tup.R/Tup.BetaRICH) *pow (1-pow (Tup.BetaRICH,2),0.5) );
      if (Betastrongcut&& cmask.isFromAgl() ) MassAgl->Fill ( (Tup.R/Tup.BetaRICH) *pow (1-pow (Tup.BetaRICH,2),0.5) );
      if (Likcut&&Distcut) {
         if (Betastrongcut&&Tup.BetaRICH<0) MassTOFQ->Fill ( (Tup.R/Tup.Beta) *pow (1-pow (Tup.Beta,2),0.5) );
         if (Betastrongcut&& cmask.isFromNaF() ) MassNaFQ->Fill ( (Tup.R/Tup.BetaRICH) *pow (1-pow (Tup.BetaRICH,2),0.5) );
         if (Betastrongcut&& cmask.isFromAgl() ) MassAglQ->Fill ( (Tup.R/Tup.BetaRICH) *pow (1-pow (Tup.BetaRICH,2),0.5) );
      }

   }

}




void SlidesforPlot_Write()
{
   RvsBetaTOF_P->Write();
   RvsBetaNaF_P->Write();
   RvsBetaAgl_P->Write();
   RvsBetaTOF_D->Write();
   RvsBetaNaF_D->Write();
   RvsBetaAgl_D->Write();
   RvsBetaTOF_He->Write();
   RvsBetaNaF_He->Write();
   RvsBetaAgl_He->Write();
   MassTOF_P->Write();
   MassTOF_D->Write();
   MassNaF_P->Write();
   MassNaF_D->Write();
   MassAgl_P->Write();
   MassAgl_D->Write();
   MassTOF_PQ->Write();
   MassTOF_DQ->Write();
   MassNaF_PQ->Write();
   MassNaF_DQ->Write();
   MassAgl_PQ->Write();
   MassAgl_DQ->Write();
   EdepUTOFvsR_P->Write();
   EdepUTOFvsR_D->Write();
   EdepUTOFvsR_He->Write();
   EdepLTOFvsR_P->Write();
   EdepLTOFvsR_D->Write();
   EdepLTOFvsR_He->Write();
   EdepTrackvsR_P->Write();
   EdepTrackvsR_D->Write();
   EdepTrackvsR_He->Write();
   RvsBetaTOF->Write();
   RvsBetaNaF->Write();
   RvsBetaAgl->Write();
   MassTOF->Write();
   MassNaF->Write();
   MassAgl->Write();
   MassTOFQ->Write();
   MassNaFQ->Write();
   MassAglQ->Write();
   DistTOF_P->Write();
   DistNaF_P->Write();
   DistAgl_P->Write();
   DistTOF_D->Write();
   DistNaF_D->Write();
   DistAgl_D->Write();
   DistTOF_He->Write();
   DistNaF_He->Write();
   DistAgl_He->Write();
   RvsDistTOF_P->Write();
   RvsDistNaF_P->Write();
   RvsDistAgl_P->Write();
   RvsDistTOF_D->Write();
   RvsDistNaF_D->Write();
   RvsDistAgl_D->Write();
   RvsDistTOF_He->Write();
   RvsDistNaF_He->Write();
   RvsDistAgl_He->Write();
   sigmagen_bad->Write();
}

