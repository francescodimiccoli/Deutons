#ifndef BINNING_H
#define BINNING_H

#include "particle.h"
#include "printmatrix.h"

#include <fstream>
#include <vector>


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

      inline void setBinsFromEk (int, float, float); ///< nbins, min, max
      inline void setBinsFromRigidity (int, float, float);
      inline void setBinsFromEkPerMass(int, float, float);
      int size() {return ekbincent.size(); };
      
      /** @brief Returns the rigidity bin containing the variable
       *  @param float var : the variable whose location to search
       *  @return int      : the bin number in rigidity of the variable
       */
      inline int GetRBin (float var);
      Particle getParticle() {return particle; }


      template<typename T>
      void  RFill (T * h, float Var)  ///< Fill the histogram with var indicating the rigidity bin
      { h->Fill(GetRBin(Var)); }
      
      inline void Print(); ///< Print the content of the bins


      std::vector<float> EkBins  ()   {  return   ekbin;   }
      std::vector<float> EtotBins()   {  return etotbin;   }
      std::vector<float> MomBins ()   {  return  mombin;   }
      std::vector<float> RigBins ()   {  return  rigbin;   }
      std::vector<float> BetaBins()   {  return betabin;   }

      std::vector<float> EkBinsCent()   { return    ekbincent; }  ///< bin centers in log
      std::vector<float> EtotBinsCent() { return  etotbincent; }
      std::vector<float> MomBinsCent () { return   mombincent; }
      std::vector<float> RigBinsCent () { return   rigbincent; }
      std::vector<float> BetaBinsCent() { return  betabincent; }

inline std::vector<float> EtotPerMassBins();  ///< returns Etot per mass

      float EkBin  (int bin)        {  return   ekbin[bin];        }
      float EtotBin  (int bin)      {  return   etotbin[bin];      }
      float MomBin (int bin)        {  return  mombin[bin];        }
      float RigBin (int bin)        {  return  rigbin[bin];        }
      float BetaBin (int bin)       {  return betabin[bin];        }
      float EkPerMassBin (int bin)  {  return ekpermassbin[bin];   }

      float EkBinCent   (int bin)       { return    ekbincent[bin];        }  ///< bin centers in log
      float EtotBinCent (int bin)       { return  etotbincent[bin];        }
      float MomBinCent  (int bin)       { return   mombincent[bin];        }
      float RigBinCent  (int bin)       { return   rigbincent[bin];        }
      float BetaBinCent (int bin)       { return  betabincent[bin];        }
      float EkPerMassBinCent (int bin)  {  return ekpermassbincent[bin];   }



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

      inline void pushBackVelocities ();
      inline void pushBackCentralVelocities ();
      inline std::vector<float> computeLogBinEdges(int nbins, float min, float max);
      inline std::vector<float> computeLogBinCenters(int nbins, float min, float max);

};


void Binning::setBinsFromEk (int nbins, float min, float max)
{
   std::vector<float> vedg=computeLogBinEdges(nbins, min, max);
   for (float binedge:vedg) {
      particle.FillFromEk (binedge);
      pushBackVelocities();
   }
   
   std::vector<float> vcen=computeLogBinCenters(nbins, min, max);
   for (float bincenter:vcen) {
      particle.FillFromEk (bincenter);
      pushBackCentralVelocities();
   }

   return;
}

void Binning::setBinsFromRigidity (int nbins, float min, float max)
{
   std::vector<float> vedg=computeLogBinEdges(nbins, min, max);
   for (float binedge:vedg) {
      particle.FillFromRig (binedge);
      pushBackVelocities();
   }
   
   std::vector<float> vcen=computeLogBinCenters(nbins, min, max);
   for (float bincenter:vcen) {
      particle.FillFromRig (bincenter);
      pushBackCentralVelocities();
   }

   return;
}


void Binning::setBinsFromEkPerMass (int nbins, float min, float max)
{
   std::vector<float> vedg=computeLogBinEdges(nbins, min, max);
   for (float binedge:vedg) {
      particle.FillFromEkPerMass (binedge);
      pushBackVelocities();
   }
   
   std::vector<float> vcen=computeLogBinCenters(nbins, min, max);
   for (float bincenter:vcen) {
      particle.FillFromEkPerMass (bincenter);
      pushBackCentralVelocities();
   }

   return;
}





std::vector<float> Binning::computeLogBinEdges(int nbins, float min, float max) {
   std::vector<float> binEdges(nbins+1);
   float logmin=log (min), logmax=log (max);
   float binstep= (logmax-logmin)  / nbins;
   for (int ibin=0; ibin<nbins; ibin++) {
      float binbeg=min*exp (ibin * binstep);
      binEdges[ibin]=binbeg;
   }
    binEdges[nbins]=max; // So it has the "right" value for sure
   return binEdges;
}

std::vector<float> Binning::computeLogBinCenters(int nbins, float min, float max) {
   std::vector<float> binCent(nbins);
   float logmin=log (min), logmax=log (max);
   float binstep= (logmax-logmin)  / nbins;
   for (int ibin=0; ibin<nbins; ibin++) {
      float bincent=min* exp ( (ibin+0.5) * binstep);
      binCent[ibin]=bincent;
   }
   return binCent;
   
}



void Binning::pushBackVelocities ()
{
   ekbin.  push_back      (particle.getEkin()        );
   etotbin.push_back      (particle.getEtot()        );
   mombin. push_back      (particle.getMom()         );
   rigbin. push_back      (particle.getRig()         );
   betabin.push_back      (particle.getBeta()        );
   ekpermassbin.push_back (particle.getEkinPerMass() );
   return;
}

void Binning::pushBackCentralVelocities ()
{
   ekbincent.  push_back      (particle.getEkin()        );
   etotbincent.push_back      (particle.getEtot()        );
   mombincent. push_back      (particle.getMom()         );
   rigbincent. push_back      (particle.getRig()         );
   betabincent.push_back      (particle.getBeta()        );
   ekpermassbincent.push_back (particle.getEkinPerMass() );
   return;
}


void Binning::Print()
{
   printMatrix::print(
      { ekbin, mombin, rigbin, betabin },
      {"Ekin", "Momentum", "Rigidity", "Beta"}
   );
   return;
}



std::vector<float> Binning::EtotPerMassBins()
{
   std::vector<float> etotpermass (0);
   if (particle.getMass()!=0)
      for (float e : ekbin)   etotpermass.push_back (e/particle.getMass());
   return etotpermass;
}



int Binning::GetRBin (float var)
{
   if (var<rigbin[0]) return -1;
   for (uint ib=0; ib<rigbin.size(); ib++)  {
      if (var>rigbin[ib] && var<=rigbin[ib+1])
         return ib;
   }
   return -1;
}

#endif
