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


#include "Parte2/cutmask.h"
#include "Parte2/TriggPatt.h"


#include "Parte2/Definitions.cpp"
#include "Parte2/FitError.cpp"
#include "Parte2/EfficiencyClass.cpp"
#include "Parte2/LATcorrClass.cpp"
#include "Parte2/ACCEPTANCEClass.cpp"
#include "Parte2/TemplateFITClass.cpp"
#include "Parte2/DatavsMCClass.h"
#include "Parte2/OptimizationCutClass.h"

#include "Parte2/ExposureTime.cpp"

#include "Parte2/MCpreeff.cpp"
#include "Parte2/MCUnbiaseff.cpp"
#include "Parte2/Hecont.cpp"
#include "Parte2/SlidesforPlot.cpp"
#include "Parte2/Qualcutoptimization.cpp"
#include "Parte2/MCQualeff.cpp"
#include "Parte2/Cuts.cpp"
#include "Parte2/MCTrackeff.cpp"
#include "Parte2/Eff_Factorizationtest.cpp"
#include "Parte2/MigrationMatrix.cpp"
#include "Parte2/MCFullSeteff.cpp"
#include "Parte2/AntiDCutOptimization.h"
#include "Parte2/AntiDPredictions.h"
#include "Parte2/DATAUnbiaseff.cpp"
#include "Parte2/CorrelazionePreselezioni.cpp"
#include "Parte2/DATApreSeleff.cpp"
#include "Parte2/DATAQualeff.cpp"
#include "Parte2/DATARICHeff.cpp"
#include "Parte2/CorrLAT.cpp"
#include "Parte2/DeutonsCountsExtraction.cpp"
#include "Parte2/DeutonsCountsExtraction_Dist.cpp"
#include "Parte2/MCMC.cpp"
#include "Parte2/DVSMCQualeff.cpp"
#include "Parte2/DVSMCRICHeff.cpp"
#include "Parte2/DVSMCQualeff_D.cpp"
#include "Parte2/DVSMCPreSeleff.cpp"
#include "Parte2/DVSMCPreSeleff_D.cpp"
#include "Parte2/Acceptance.cpp"

//#include "Parte2/FluxImplementation.cpp"
#include "Parte2/FluxClass.cpp"
#include "Parte2/ProtonFlux.cpp"
#include "Parte2/DeutonsFlux.cpp"


#include "Parte2/PlottingFunctions/OtherExperimentsComparison.cpp"
#include "Parte2/DVSMCTrackeff.cpp"



using namespace std;

#include "FillIstogram.cpp"



int main(int argc, char * argv[])
{
   cout<<"Month _ Indx _ Frac _ output"<<endl;
   cout<<argc<<endl;
   mode INDX=READ;
   string frac="35";
   string outputpath = "../";
   string mese="2012_05";

   if(argc > 1 ) {
      mese=argv[1];
   } else {
      cout<<"No Month specified: running 2012_05"<<endl;
   }

   if(argc > 2 ) {
      int index=atoi(argv[2]);
      if (index>=0 && index<=2)
         INDX=static_cast<mode>(index);
   } else {
      cout<<"No Mode specified: running Mode 2"<<endl;
   }

   if(argc > 3 ) {
      frac= argv[3];
   } else {
      cout<<"No fraction specified: running 35"<<endl;
   }

   if(argc > 4 ) {
      outputpath=argv[4];
   } else {
      cout<<"No output path specified: writing locally"<<endl;
   }


   cout<<"****************************** INPUT PAR. ***********************************"<<endl;
   cout<<endl;
   cout<<"Month: "<<mese<<endl;
   cout<<endl;
   cout<<"Mode: "<<INDX<<endl;
   cout<<endl;
   cout<<"Fraction: "<<frac<<endl;
   cout<<endl;
   cout<<"Output dir: "<<outputpath<<"Histos/"<<mese<<endl;
   cout<<endl;



   
   cout<<"****************************** BINS ***************************************"<<endl;

   DRB.setBinsFromRigidity(nbinsr, 0.5, 100); // RB did not have Ek
   PRB.setBinsFromRigidity(nbinsr, 0.5, 100); 

   float ekmin=0.1, ekmax=1;
   ToFDB.setBinsFromEkPerMass (nbinsToF, ekmin, ekmax);
   ToFPB.setBinsFromEkPerMass(nbinsToF, ekmin, ekmax);

   ekmin=0.666, ekmax=4.025;
   NaFDB.setBinsFromEkPerMass(nbinsNaF, ekmin, ekmax);
   NaFPB.setBinsFromEkPerMass(nbinsNaF, ekmin, ekmax);

   ekmin=2.57, ekmax=9.01;
   AglDB.setBinsFromEkPerMass(nbinsAgl, ekmin, ekmax);
   AglPB.setBinsFromEkPerMass(nbinsAgl, ekmin, ekmax);

   ToFDB.UseREdges();
   ToFPB.UseREdges();
   NaFDB.UseREdges();
   NaFPB.UseREdges();
   AglDB.UseREdges();
   AglPB.UseREdges();

   DRB.UseREdges();
   PRB.UseREdges();   


   cout<<endl;
   
   cout<<"************************ ISTOGRAM FILLING **************************************************************"<<endl;

   FillIstogramAndDoAnalysis(INDX,frac,mese, outputpath);
   

   cout<<"************************** OUTPUT **************************************************************"<<endl;
   return 1;
}
