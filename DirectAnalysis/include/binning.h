#ifndef BINNING_H
#define BINNING_H

#include "particle.h"
#include "printmatrix.h"

#include <fstream>
#include <vector>
#include "TH1F.h"
#include "TF1.h"

/** @brief Class used for all the binning manipulation
 * @author $Author: L. Basara$
 * @date $Date: 2016/05/19$
*/

class Binning {
   public:
      Binning () :                  particle (1) {}
      Binning (float m) :           particle (m) {}
      Binning (float m, float z) :  particle (m, z) {}
      Binning (float m, float z, float a) :  particle (m, z, a) {}
      Binning (Particle p) : particle (p) {}

      void setBinsFromEk       (int, float, float, TF1 * SlowDownModel, float alpha_slowdown, float gamma_slowdown); ///< nbins, min, max
      void setBinsFromRigidity (int, float, float, TF1 * SlowDownModel, float alpha_slowdown, float gamma_slowdown);
      void setBinsFromEkPerMass(int, float, float, TF1 * SlowDownModel, float alpha_slowdown, float gamma_slowdown);
      void setBinsFromBeta     (int, float, float, TF1 * SlowDownModel, float alpha_slowdown, float gamma_slowdown);

      void Reset();

      int size() {return ekbincent.size(); };
      
      /** @brief Returns the rigidity bin containing the variable
       *  @param float var : the variable whose location to search
       *  @return int      : the bin number in rigidity of the variable
       */
      inline void UseBetaEdges() {Betaedges=true; Redges=false;BetaTOIedges=false;RTOIedges=false;}
      inline void UseREdges()    {Betaedges=false;Redges=true ;BetaTOIedges=false;RTOIedges=false; }	
      inline void UseBetaTOIEdges() {Betaedges=false;Redges=false;BetaTOIedges=true; RTOIedges=false;}
      inline void UseRTOIEdges()    {Betaedges=false;Redges=false;BetaTOIedges=false;RTOIedges=true; }	
      
      int GetBin (float var);
      inline int GetRBin (float var);
      inline int GetBetaBin (float var);
      inline int GetRTOIBin (float var);
      inline int GetBetaTOIBin (float var);


      Particle getParticle() {return particle; }
      inline void  RFill (TH1* h, float Var); ///< Fill the histogram with var indicating the rigidity bin
      void Print(); ///< Print the content of the bins
      float GetBinLowEdge(int bin);	
      float GetBinCenter(int bin);	

      std::vector<float> EkBins  ()   {  return   ekbin;   }
      std::vector<float> EtotBins()   {  return etotbin;   }
      std::vector<float> MomBins ()   {  return  mombin;   }
      std::vector<float> RigBins ()   {  return  rigbin;   }
      std::vector<float> BetaBins()   {  return betabin;   }
      std::vector<float> EkPerMasBins  ()   {  return   ekpermassbin;  }	

      std::vector<float> EkBinsCent()   { return    ekbincent; }  ///< bin centers in log
      std::vector<float> EtotBinsCent() { return  etotbincent; }
      std::vector<float> MomBinsCent () { return   mombincent; }
      std::vector<float> RigBinsCent () { return   rigbincent; }
      std::vector<float> BetaBinsCent() { return  betabincent; }

      std::vector<float> EkTOIBins  ()   {  return   ekbin_TOI;   }
      std::vector<float> EtotTOIBins()   {  return etotbin_TOI;   }
      std::vector<float> MomTOIBins ()   {  return  mombin_TOI;   }
      std::vector<float> RigTOIBins ()   {  return  rigbin_TOI;   }
      std::vector<float> BetaTOIBins()   {  return betabin_TOI;   }
      std::vector<float> EkPerMasTOIBins  ()   {  return   ekpermassbin;  }	

      std::vector<float> EkTOIBinsCent()   { return    ekbincent_TOI; }  ///< bin centers in log
      std::vector<float> EtotTOIBinsCent() { return  etotbincent_TOI; }
      std::vector<float> MomTOIBinsCent () { return   mombincent_TOI; }
      std::vector<float> RigTOIBinsCent () { return   rigbincent_TOI; }
      std::vector<float> BetaTOIBinsCent() { return  betabincent_TOI; }


      inline std::vector<float> EtotPerMassBins();  ///< returns Etot per mass

      float EkBin  (int bin)        {  if(IsUsingTOIEdges())  return   ekbin_TOI[bin];       else  return ekbin[bin];       }
      float EtotBin  (int bin)      {  if(IsUsingTOIEdges())  return   etotbin_TOI[bin];      else return etotbin[bin];     }
      float MomBin (int bin)        {  if(IsUsingTOIEdges())  return  mombin_TOI[bin];        else return mombin[bin];       }
      float RigBin (int bin)        {  if(IsUsingTOIEdges())  return  rigbin_TOI[bin];        else return rigbin[bin];       }
      float BetaBin (int bin)       {  if(IsUsingTOIEdges())  return betabin_TOI[bin];        else return betabin[bin];       }
      float EkPerMassBin (int bin)  {  if(IsUsingTOIEdges())  return ekpermassbin_TOI[bin];   else return ekpermassbin[bin];  }

      float EkBinCent   (int bin)       { if(IsUsingTOIEdges())  return    ekbincent_TOI[bin];     else  return    ekbincent[bin];       }  ///< bin centers in log
      float EtotBinCent (int bin)       { if(IsUsingTOIEdges())  return  etotbincent_TOI[bin];     else  return  etotbincent[bin];       }
      float MomBinCent  (int bin)       { if(IsUsingTOIEdges())  return   mombincent_TOI[bin];     else  return   mombincent[bin];       }
      float RigBinCent  (int bin)       { if(IsUsingTOIEdges())  return   rigbincent_TOI[bin];     else  return   rigbincent[bin];       }
      float BetaBinCent (int bin)       { if(IsUsingTOIEdges())  return  betabincent_TOI[bin];     else  return  betabincent[bin];       }
      float EkPerMassBinCent (int bin)  {  if(IsUsingTOIEdges())  return ekpermassbincent_TOI[bin]; else  return ekpermassbincent[bin];  }

      bool IsUsingBetaEdges() {return (Betaedges || BetaTOIedges);}
      bool IsUsingTOIEdges()  {return (BetaTOIedges || RTOIedges);}	

   private:
      Particle particle;

      std::vector<float> ekbin ;
      std::vector<float> etotbin ;
      std::vector<float> mombin ;
      std::vector<float> rigbin ;
      std::vector<float> betabin ;
      std::vector<float> ekpermassbin ;

      std::vector<float> ekbincent ;
      std::vector<float> etotbincent ;
      std::vector<float> mombincent ;
      std::vector<float> rigbincent ;
      std::vector<float> betabincent ;
      std::vector<float> ekpermassbincent ;

      std::vector<float> ekbin_TOI ;
      std::vector<float> etotbin_TOI ;
      std::vector<float> mombin_TOI ;
      std::vector<float> rigbin_TOI ;
      std::vector<float> betabin_TOI ;
      std::vector<float> ekpermassbin_TOI ;

      std::vector<float> ekbincent_TOI ;
      std::vector<float> etotbincent_TOI ;
      std::vector<float> mombincent_TOI ;
      std::vector<float> rigbincent_TOI ;
      std::vector<float> betabincent_TOI ;
      std::vector<float> ekpermassbincent_TOI ;


      inline void pushBackVelocities ();
      inline void pushBackCentralVelocities ();
      inline std::vector<float> computeLogBinEdges(int nbins, float min, float max);
      inline std::vector<float> computeLogBinCenters(int nbins, float min, float max);
      inline std::vector<float> computeConstResoBinEdges(int nbins, float min, float max);
      inline std::vector<float> computeConstResoBinCenters(int nbins, float min, float max);
	   
  
      bool Betaedges=false;
      bool Redges   =false;		
      bool BetaTOIedges=false;
      bool RTOIedges   =false;		

};

#endif
