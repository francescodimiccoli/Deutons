#ifndef CALIBRATIONS__H
#define CALIBRATIONS__H

#include <string>

#include "TF1.h"
#include "TSpline.h"

class Calibrations {

    TSpline * Rig;           
    TSpline * beta;          
    TSpline * EdepL1beta;    

    TSpline * etofu;         
    TSpline * EdepTOFbeta;   

    TSpline * etrack;        
    TSpline * EdepTrackbeta; 

    TF1 * betaNaF;
    TF1 * betaAgl;

public:
    float OffsetTrackBeta(float beta, float edep) const;
    float OffsetTOFBeta(float beta, float edep) const;
    Calibrations(std::string filename);
};

#endif
