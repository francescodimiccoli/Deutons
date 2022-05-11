#ifndef _CUTOFF_H_
#define _CUTOFF_H_

#include "TString.h"

#include "GeomHashes.h"

namespace rigidityCutoff {

  const TString filePath="root://eosams.cern.ch//eos/ams/group/dbar/data";

  const TString fileName[] = { 
    "rigidityCutoffPositiveUpperdif.root",
    "rigidityCutoffPositivePenumbra.root",
    "rigidityCutoffNegativeUpperdif.root",
    "rigidityCutoffNegativePenumbra.root"
  };

  const TString hashName[] = {
    "rgCutPuP", "rgCutPud", "rgCutNuN", "rgCutNud"
  };


  const int positiveUpperdif = 0;
  const int positivePenumbra = 1;
  const int negativeUpperdif = 2;
  const int negativePenumbra = 3;

  const int cutoffdif = 0;
  const int penumbra  = 1;

  const int mean = 0;
  const int meanTruncated = 1;
  const int meanCorrected = 2;
  const int meanWeightedCorrected = 3;

}


//
// Evaluate Positive/Negative Upper Rigidity Cutoff for IGRF 2012 model
float rigidityCutoffPositiveUpper(float radS, float latS, float phiS, float thetaP, float phiP);
float rigidityCutoffNegativeUpper(float radS, float latS, float phiS, float thetaP, float phiP);

//
// Evaluate Differences between IGRF Upper Rigidity and Stormer cutoffs from interpolated maps
float rigidityCutoffPositiveUpperdif(float latS, float phiS, float thetaP, float phiP);
float rigidityCutoffNegativeUpperdif(float latS, float phiS, float thetaP, float phiP);

//
// Evaluate Difference between IGRF Upper and Lower Rigidity cutoffs from interpolated maps
float rigidityCutoffPositivePenumbra(float latS, float phiS, float thetaP, float phiP);
float rigidityCutoffNegativePenumbra(float latS, float phiS, float thetaP, float phiP);

//
// Base interpolation function on cutoff maps obtained via back tracing using IGRF 2012 model
float rigidityCutoffInterpolated(float latS, float phiS, float thetaP, float phiP,
				 bool isPositive, int hashType, GeomHashNamedEnsemble *hashE, int valType);

//
// Approximate Correction function based on Stormer Cutoff
float rigidityCutoffStormerCorrection(bool isPositive, float *xf, float *xi);

#endif
