#ifndef BINNING_H
#define BINNING_H

#include "particle.h"
#include "printmatrix.h"

#include <fstream>
#include <vector>
#include "TH1F.h"


/** @brief Class used for all the binning manipulation
 * @author $Author: L. Basara$
 * @date $Date: 2016/05/19$
 */

class Binning {
   public:
      Binning (float m) :           particle (m)       {}
      Binning (float m, float z) :  particle (m, z) {}
      Binning (float m, float z, float a) :  particle (m, z, a) {}
      Binning (Particle p) : particle (p) {}

      void setBinsFromEk (int, float, float); ///< nbins, min, max
      void setBinsFromRigidity (int, float, float);
      int size() {return ekbincent.size(); };


      /** @brief Returns the rigidity bin containing the variable
       *  @param float var : the variable whose location to search
       *  @return int      : the bin number in rigidity of the variable
       */
      int GetRBin (float var);
      Particle getParticle() {return particle; }


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

      std::vector<float> EkPerMassBins  ();  ///< returns Ek per mass
      std::vector<float> EtotPerMassBins();  ///< returns Etot per mass

      float EkBin  (int bin)   {  return   ekbin[bin];   }
      float EtotBin  (int bin) {  return   etotbin[bin]; }
      float MomBin (int bin)   {  return  mombin[bin];   }
      float RigBin (int bin)   {  return  rigbin[bin];   }
      float BetaBin (int bin)  {  return betabin[bin];   }

      float EkBinCent  (int bin) { return    ekbincent[bin]; }  ///< bin centers in log
      float EtotBinCent (int bin) { return  etotbincent[bin]; }
      float MomBinCent (int bin) { return   mombincent[bin]; }
      float RigBinCent (int bin) { return   rigbincent[bin]; }
      float BetaBinCent (int bin) { return  betabincent[bin]; }

      void  RFill (TH1* h, float Var); ///< Fill the histogram with var indicating the rigidity bin

      void Print(); ///< Print the content of the bins


   private:
      Particle particle;
      float mass;
      int Z=1;
      int A=1;

      std::vector<float>   ekbin ;
      std::vector<float> etotbin ;
      std::vector<float>  mombin ;
      std::vector<float>  rigbin ;
      std::vector<float> betabin ;

      std::vector<float>   ekbincent ;
      std::vector<float> etotbincent ;
      std::vector<float>  mombincent ;
      std::vector<float>  rigbincent ;
      std::vector<float> betabincent ;

      void pushBackVelocities ();
      void pushBackCentralVelocities ();

};


void Binning::setBinsFromEk (int nbins, float min, float max)
{
   float logmin=log (min), logmax=log (max);
   float binbeg=min;
   float binstep= (logmax-logmin)  / nbins;
   int ibin=0;
   std::vector<float> vbin (nbins+1); // bins
   std::vector<float> vcen (nbins); // centers
   Particle Pedge (particle), Pcent (particle);

   // Filling the vectors
   while (binbeg<max) {
      particle.FillFromEk (binbeg);
      pushBackVelocities();
      float bincent= exp ( logmin + (ibin+0.5) * binstep);
      particle.FillFromEk (bincent);
      pushBackCentralVelocities ();
      ibin++;
      binbeg = exp ( logmin + ibin * binstep);
   }

   particle.FillFromEk (max); // Don't forget the last edge
   pushBackVelocities();
   return;
}

void Binning::setBinsFromRigidity (int nbins, float min, float max)
{
   float logmin=log (min), logmax=log (max);
   float binbeg=min;
   float binstep= (logmax-logmin)  / nbins;
   int ibin=0;
   std::vector<float> vbin (nbins+1); // bins
   std::vector<float> vcen (nbins); // centers

   // Filling the vectors
   while (binbeg<max) {
      particle.FillFromRig (binbeg);
      pushBackVelocities();
      float bincent= exp ( logmin + (ibin+0.5) * binstep);
      particle.FillFromRig (bincent);
      pushBackCentralVelocities ();
      ibin++;
      binbeg = exp ( logmin + ibin * binstep);
   }

   particle.FillFromRig (max); // Don't forget the last edge
   pushBackVelocities();
   return;
}




void Binning::pushBackVelocities ()
{
   ekbin.  push_back (particle.getEkin() );
   etotbin.push_back (particle.getEtot() );
   mombin. push_back (particle.getMom()  );
   rigbin. push_back (particle.getRig()  );
   betabin.push_back (particle.getBeta() );
   return;
}

void Binning::pushBackCentralVelocities ()
{
   ekbincent.  push_back (particle.getEkin() );
   etotbincent.push_back (particle.getEtot() );
   mombincent. push_back (particle.getMom()  );
   rigbincent. push_back (particle.getRig()  );
   betabincent.push_back (particle.getBeta() );
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



std::vector<float> Binning::EkPerMassBins()
{
   std::vector<float> ekpermass (0);
   if (mass!=0)
      for (float e : ekbin)   ekpermass.push_back (e/mass);
   return ekpermass;
}



std::vector<float> Binning::EtotPerMassBins()
{
   std::vector<float> etotpermass (0);
   if (mass!=0)
      for (float e : ekbin)   etotpermass.push_back (e/mass);
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

void Binning::RFill (TH1* h, float var)
{
   h->Fill (GetRBin (var) );
}




#endif
