#include "PlottingFunctions/FluxFactorizationtest_Plot.h"

using namespace std;

Efficiency * EffFullSETselectionsMCP =  new Efficiency("EffFullSETselectionsMCP");

Efficiency * Eff_do_preSelMCP = new Efficiency("Eff_do_preSelMCP",3);
Efficiency * Eff_do_preSelMCD = new Efficiency("Eff_do_preSelMCD",6,3);


Efficiency * Eff_do_DistMCP  = new Efficiency("Eff_do_DistMCP");
Efficiency * Eff_do_LikMCP   = new Efficiency("Eff_do_LikMCP");


void FluxFactorizationtest_Pre_Fill()
{
   if(!trgpatt.IsPhysical()) return;

   int Rbin;
   // full set efficiency before
   if(cmask.isMinimumBiasTrigger()&&cmask.isMinimumBiasToF3or4Layers()&&cmask.isMinimumBiasTracker()&&Tup.Beta_pre>0&&Tup.R_pre>0) {
      Rbin=PRB.GetRBin(Tup.R_pre);
      if(Massa_gen<1&&Massa_gen>0.5) {
         (EffFullSETselectionsMCP->beforeR)->Fill(Rbin);
	}
   }

	
   if(Tup.Beta_pre<=0||Tup.R_pre<=0) return;
   
   if(!ProtonsMassWindow) return;
   if(!Herejcut) return;
	
   //Drop-one approach eff. calc.
   for(int iS=0; iS<3; iS++) {
      Rbin=PRB.GetRBin(RUsed);
      if(Massa_gen<1&&Massa_gen>0.5) {
         if(cmask.notPassed(iS)) ((TH2*)Eff_do_preSelMCP->beforeR)->Fill(Rbin,iS);
         if(cmask.passed(iS)) ((TH2*)Eff_do_preSelMCP->afterR) ->Fill(Rbin,iS);
      }

      if(Massa_gen>1&&Massa_gen<2) {
         if(cmask.notPassed(iS)) FillBinMGen((TH3*)Eff_do_preSelMCD->beforeR, Rbin, iS);
         if(cmask.passed(iS)) FillBinMGen((TH3*)Eff_do_preSelMCD->afterR,  Rbin, iS);

      }
   }
   ////////////////////////////////
   return;
}




void FluxFactorizationtest_Qual_Fill()
{

  if(!trgpatt.IsPhysical()||Tup.Beta<=0||Tup.R<=0)return;
 
   //R bins
   int Kbin;
   Kbin = PRB.GetRBin(RUsed);

   //full set efficiency after
   if(Massa_gen<1&&Massa_gen>0.5) {
	if(Tup.Dist5D_P<6&&Likcut)  (EffFullSETselectionsMCP->afterR)->Fill(Kbin);
	}

   //Drop-one approach eff calc.
   //eff evaluation cuts
  if(!ProtonsMassWindow) return;
   if(!Herejcut) return;

   if(Massa_gen<1&&Massa_gen>0.5) {
    Kbin = PRB.GetRBin(RUsed);
     ( Eff_do_DistMCP -> beforeR) -> Fill(Kbin);
      if(Tup.Dist5D_P<6) (Eff_do_LikMCP -> beforeR) -> Fill(Kbin);

      if(Tup.Dist5D_P<6) {
        ( Eff_do_DistMCP -> afterR) -> Fill(Kbin);
         if(Likcut) (Eff_do_LikMCP -> afterR) -> Fill(Kbin);
      }

   }
   ///////////////////////////
   return;

}


void FluxFactorizationtest_Write()
{
   Eff_do_preSelMCP -> Write();
   Eff_do_preSelMCD ->Write();

   Eff_do_DistMCP   -> Write();
   Eff_do_LikMCP    ->Write();

   EffFullSETselectionsMCP ->Write();
   return;
}


void FluxFactorizationtest(string filename)
{
  cout<<"********* MC \"GOLDEN\" SEL. EFFICIENCIES *********"<<endl;

  cout<<"*** Reading  P1 file ****"<<endl;
  TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");


   Efficiency * EffFullSETselectionsMCP  = new Efficiency(inputHistoFile,"EffFullSETselectionsMCP");

   Efficiency * Eff_do_preSelMCP = new Efficiency(inputHistoFile,"Eff_do_preSelMCP");
   Efficiency * Eff_do_preSelMCD = new Efficiency(inputHistoFile,"Eff_do_preSelMCD");

   Efficiency * Eff_do_DistMCP  = new Efficiency(inputHistoFile,"Eff_do_DistMCP");
   Efficiency * Eff_do_LikMCP   = new Efficiency(inputHistoFile,"Eff_do_LikMCP");

   string tagli[5]= {"Matching TOF","Chi^2 R","1 Tr. Track","Distance","Likelihood"};
   string nome;

   cout<<"********* MC \"GOLDEN\" SEL. EFFICIENCIES *********"<<endl;

   EffFullSETselectionsMCP -> Eval_Efficiency();

   Eff_do_preSelMCP 	-> Eval_Efficiency();
   Eff_do_preSelMCD 	-> Eval_Efficiency();

   Eff_do_DistMCP  	-> Eval_Efficiency();
   Eff_do_LikMCP   	-> Eval_Efficiency();



   TH1F * Eff_FullSETMCP_R_TH1F = (TH1F *) EffFullSETselectionsMCP -> effR -> Clone();

   TH2F * Eff_do_preSelMCP_R_TH2F	= (TH2F *) Eff_do_preSelMCP -> effR -> Clone();

   TH1F * Eff_do_DistMCP_R_TH1F = (TH1F *) Eff_do_DistMCP -> effR -> Clone();
   TH1F * Eff_do_LikMCP_R_TH1F = (TH1F *) Eff_do_LikMCP -> effR -> Clone();


   // factorized eff. calc.
   TH1F * FactorizedEffMCP_R = (TH1F *) Eff_do_DistMCP_R_TH1F -> Clone();
   FactorizedEffMCP_R -> Multiply(Eff_do_LikMCP_R_TH1F);

   for(int iS=0; iS<Eff_do_preSelMCP_R_TH2F->GetNbinsY(); iS++) {
      for(int iR=0; iR<Eff_do_preSelMCP_R_TH2F->GetNbinsX(); iR++)
         FactorizedEffMCP_R -> SetBinContent(iR+1,FactorizedEffMCP_R -> GetBinContent(iR+1)*Eff_do_preSelMCP_R_TH2F-> GetBinContent(iR+1,iS+1) );
   }

   cout<<"*** Plotting ...  ****"<<endl;
   FluxFactorizationtest_Plot(
	Eff_FullSETMCP_R_TH1F,
	FactorizedEffMCP_R
	);  	

   return;
}

