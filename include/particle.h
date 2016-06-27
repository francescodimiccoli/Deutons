#ifndef PARTICLE_H
#define PARTICLE_H

#include <cmath>
#include <limits>


/** @brief Classed used to convert between kinetic energy, momentum, rigidity and beta
 * @author $Author: L. Basara$
 * @date   $Date: 2016/05/19$
 */

class Particle {
   public:
      Particle (float m) :          mass (m) {}
      Particle (float m, float z) : mass (m), Z(z) {}
      Particle (float m, float z, float a) : mass (m), Z(z), A (a)  {}

      inline void FillFromEk        (float);
      inline void FillFromEkPerMass (float);
      inline void FillFromRig       (float);

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

      float BetaFromEk (float ek)     { return sqrt (ek*ek + 2 * ek * mass) / (ek + mass); }
      float BetaFromEtot (float e)    { return sqrt (1- (mass*mass)/(e*e) ); }
      float GammaFromEk (float ek)    { return 1 + ek/mass;                                }
      float RigFromEk (float ek)      { return RigFromMom (MomFromEk (ek) ) ;              }
      float MomFromEk  (float ek)     { return mass * BetaFromEk (ek) * GammaFromEk (ek);  }
      float EtotfromMom(float p)      { return sqrt ( mass*mass + p*p); }
      float EtotfromEk(float ek)      { return ek + mass; }
      float EkFromMom (float p)       { return EtotfromMom(p) - mass  ;  }
      float RigFromMom (float p)      { return p/Z ;   }
      float MomFromRig (float rig)    { return rig*Z ; }
      float EkPerMassFromEk(float ek) { return ek/mass; }

};

void Particle::FillFromEk (float ek)
{
   ekin=ek;
   ekpermass=EkPerMassFromEk(ekin);
   etot=  EtotfromEk(ekin);
   beta = BetaFromEk (ekin);
   mom=MomFromEk (ekin);
   rig=RigFromMom (mom);
}

void Particle::FillFromRig ( float r)
{
   rig=r;
   mom=MomFromRig (rig);
   ekin=EkFromMom (mom);
   ekpermass=EkPerMassFromEk(ekin);
   etot=EtotfromEk(ekin);  
   beta=BetaFromEk (ekin);
}

void Particle::FillFromEkPerMass(float ekpermass)
{
   ekin=ekpermass * mass;
   FillFromEk(ekin);
}


#endif 
