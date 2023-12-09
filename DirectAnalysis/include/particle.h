#ifndef PARTICLE_H
#define PARTICLE_H

#include <cmath>
#include <limits>
#include "TF1.h"
//#include <gsl/gsl_math.h>
#include "Math/WrappedTF1.h"
#include "TF2.h"
#include "TError.h"
#include "RConfigure.h"
#include "Math/BrentRootFinder.h"


/** @brief Classed used to convert between kinetic energy, momentum, rigidity and beta
 * @author $Author: L. Basara$
 * @date   $Date: 2016/05/19$
 */

using namespace ROOT::Math;

class Particle {
   public:
      Particle (float m) :          mass (m) { MIPthreshold = 1.2*mass; }
      Particle (float m, float z) : mass (m), Z(z) {MIPthreshold = 1.2*mass;}
      Particle (float m, float z, float a) : mass (m), Z(z), A (a)  {MIPthreshold = 1.2*mass;}

      void FillFromEk        (float);
      void FillFromEkPerMass (float);
      void FillFromRig       (float);
      void FillFromBeta      (float);

      float getMass() {return mass;}
      int   getA   () {return A   ;}
      int   getZ   () {return Z   ;}
      
      float getEkin() {return ekin;  }
      float getMom () {return mom ;  }
      float getRig () {return rig ;  }
      float getBeta() {return beta;  }
      float getEtot() {return etot;  }
      float getEkinPerNuc() {return ekin/A;}
      float getEkinPerMass(){return ekin/mass;}
     
      float getEkin_TOI() {return ekin_TOI;  }
      float getMom_TOI () {return mom_TOI ;  }
      float getRig_TOI () {return rig_TOI ;  }
      float getBeta_TOI() {return beta_TOI;  }
      float getEtot_TOI() {return etot_TOI;  }
      float getEkinPerNuc_TOI() {return ekin_TOI/A;}
      float getEkinPerMass_TOI(){return ekin_TOI/mass;}
      void  SetSlowDownModel(TF1 * SlowDownModel, float alpha_slowdown, float gamma_slowdown); 
 
      float BetaFromEk (float ek)     { return sqrt (ek*ek + 2 * ek * mass) / (ek + mass); }
      float BetaFromEtot (float e)    { return sqrt (1- (mass*mass)/(e*e) ); }
      float GammaFromEk (float ek)    { return 1 + ek/mass;                                }
      float GammaFromBeta (float beta){ return 1/(sqrt(1-beta*beta)); }
      float RigFromEk (float ek)      { return RigFromMom (MomFromEk (ek) ) ;              }
      float MomFromEk  (float ek)     { return mass * BetaFromEk (ek) * GammaFromEk (ek);  }
      float EtotfromMom(float p)      { return sqrt ( mass*mass + p*p); }
      float EtotfromEk(float ek)      { return ek + mass; }
      float EkFromMom (float p)       { return EtotfromMom(p) - mass  ;  }
      float RigFromMom (float p)      { return p/Z ;   }
      float MomFromRig (float rig)    { return rig*Z ; }
      float EkPerMassFromEk(float ek) { return ek/mass; }
    
      float BetaTOIFromBeta (float beta);	
      float BetaMeasFromBetaTOI (float beta);	
      float EkFromBeta (float beta) {return GammaFromBeta(beta)*mass - mass;}	
      float RigFromBeta(float beta)    {return RigFromMom(MomFromEk(EkFromBeta(beta))); }


   private:
      float mass;
      int Z=1;
      int A=1;
      float ekin=0;
      float etot=0;
      float mom =0;
      float rig =0;
      float beta=0;
      float ekpermass=0;
    
      float ekin_TOI=0;
      float etot_TOI=0;
      float mom_TOI =0;
      float rig_TOI =0;
      float beta_TOI=0;
      float ekpermass_TOI=0;
      float MIPthreshold;	

      TF1 * slowdownmodel;
};

#endif 
