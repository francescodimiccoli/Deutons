
#include "binning.h"

class ParticleBins {
    Particle particle;
    Binning RB;
    Binning ToFB;
    Binning NaFB;
    Binning AglB;
public:
    ParticleBins(Particle p, int nbinsr, int nbinsToF, int nbinsNaF, int nbinsAgl):
        particle(p), ToFB(p), NaFB(p), AglB(p), RB(p)
    { 
        RB.setBinsFromRigidity(nbinsr, 0.5, 100); 

        float ekmin=0.1, ekmax=1;
        ToFB.setBinsFromEkPerMass(nbinsToF, ekmin, ekmax);

        ekmin=0.666, ekmax=4.025;
        NaFB.setBinsFromEkPerMass(nbinsNaF, ekmin, ekmax);

        ekmin=2.57, ekmax=9.01;
        AglB.setBinsFromEkPerMass(nbinsAgl, ekmin, ekmax);
    } 

    Binning Rigidity() { return RB; } 
    Binning TOF() { return ToFB; } 
    Binning NaF() { return NaFB; } 
    Binning Agl() { return AglB; } 
};

class AnalysisBins {
    ParticleBins proton, deuton;
public:

    AnalysisBins(int nbinsr, int nbinsToF, int nbinsNaF, int nbinsAgl):
         proton(Particle(0.9382720813, 1, 1), nbinsr, nbinsToF, nbinsNaF, nbinsAgl),
         deuton(Particle(  1.8756129 , 1, 2), nbinsr, nbinsToF, nbinsNaF, nbinsAgl)
    { }

    ParticleBins Proton()  {return proton;}
    ParticleBins Deuteron(){return deuton;}
   
};
