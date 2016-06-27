#include "Calibrations.h"
#include "TFile.h"
#include <iostream>

float Calibrations::OffsetTrackBeta(float beta, float edep) const
{
    float expectedEdep = EdepTrackbeta->Eval(beta);
    float sigmaEdep    = etrack->Eval(beta);
    return fabs(expectedEdep - edep)/(pow(expectedEdep,2) * sigmaEdep);
}

float Calibrations::OffsetTOFBeta(float beta, float edep) const
{
    float expectedEdep = EdepTOFbeta->Eval(beta);
    float sigmaEdep    = etofu->Eval(beta);
    return fabs(expectedEdep - edep)/(pow(expectedEdep, 2) * sigmaEdep);
}

float Calibrations::TOFBetaInside1(float beta, float edep) const
{
    return fabs( edep - EdepTOFbeta->Eval(beta) ) < 1;
}


Calibrations::Calibrations(std::string filename) {
   std::cout<<"*********************** READING CALIBRATIONS *********************"<<std:: endl;

   TFile * calib = TFile::Open(filename.c_str());
   if(!calib) std::cout << "ERROR: MC calibration file found" << std::endl;
 
   Rig           = (TSpline3 *) calib->Get("Fit Results/Splines/Rig"           );
   beta          = (TSpline3 *) calib->Get("Fit Results/Splines/beta"          );
   etofu         = (TSpline3 *) calib->Get("Fit Results/Splines/etofu"         );
   etrack        = (TSpline3 *) calib->Get("Fit Results/Splines/etrack"        );
   EdepL1beta    = (TSpline3 *) calib->Get("Fit Results/Splines/EdepL1beta"    );
   EdepTOFbeta   = (TSpline3 *) calib->Get("Fit Results/Splines/EdepTOFbeta"   );
   EdepTrackbeta = (TSpline3 *) calib->Get("Fit Results/Splines/EdepTrackbeta" );

   betaNaF = (TF1 *) calib->Get("Fit Results/Splines/SigmaInvBetaNaF_spl");
   betaAgl = (TF1 *) calib->Get("Fit Results/Splines/SigmaInvBetaAgl_spl");
   std::cout<<"******************************"<<std::endl;
}

