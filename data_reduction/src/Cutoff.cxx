#include <iostream>
#include <stdio.h>
#include <math.h>
#include <float.h>
//#include <stdlib.h> 

#include "TMath.h"
#include "TFile.h"

#include "Cutoff.h"
#include "Coord.h"

using namespace std;

float rigidityCutoffPositiveUpper(float radS, float latS, float phiS, float thetaP, float phiP) {

  float rgCutPuP;
  float radSGM, latSGM, phiSGM, rgCutP, rgCutN;
  int time=0;

  coord(time, radS, latS, phiS, thetaP, phiP,
	radSGM, latSGM, phiSGM, rgCutP, rgCutN);

  rgCutPuP = rigidityCutoffPositiveUpperdif(latS, phiS, thetaP, phiP);

  if (rgCutPuP == FLT_MAX) return FLT_MAX;

  return rgCutP+rgCutPuP;
}


float rigidityCutoffNegativeUpper(float radS, float latS, float phiS, float thetaP, float phiP) {

  float rgCutNuN;
  float radSGM, latSGM, phiSGM, rgCutP, rgCutN;
  int time=0;

  coord(time, radS, latS, phiS, thetaP, phiP,
	radSGM, latSGM, phiSGM, rgCutP, rgCutN);

  rgCutNuN = rigidityCutoffNegativeUpperdif(latS, phiS, thetaP, phiP);

  if (rgCutNuN == -FLT_MAX) return -FLT_MAX;

  return rgCutN+rgCutNuN;
}


float rigidityCutoffPositiveUpperdif(float latS, float phiS, float thetaP, float phiP) {

  static bool first=true;
  static GeomHashNamedEnsemble *hashE=0;
#pragma omp threadprivate(first,hashE)

  int isPositive = true;
  int hashIndex = rigidityCutoff::positiveUpperdif;
  int hashType = rigidityCutoff::cutoffdif;
  int valType = rigidityCutoff::meanWeightedCorrected;

  if (first) {
    first=false;

    TString fileName = Form("%s/%s", (const char *)rigidityCutoff::filePath, (const char *)rigidityCutoff::fileName[hashIndex]);
    TFile *fHash = TFile::Open(fileName);
    if (!fHash) {
      cout << "rigidityCutoffPositiveUpperdif-E-FileNotFound : " << rigidityCutoff::fileName[hashIndex] << endl;
      return FLT_MAX;
    }

    hashE = (GeomHashNamedEnsemble *)fHash->Get(rigidityCutoff::hashName[hashIndex]); 
    if (!hashE) {
      cout << "rigidityCutoffPositiveUpperdif-E-HashNotFound : " << rigidityCutoff::hashName[hashIndex] << endl;
      return FLT_MAX;
    }
  }

  if (!hashE) return FLT_MAX;

  return rigidityCutoffInterpolated(latS, phiS, thetaP, phiP,  
				    isPositive, hashType, hashE, valType);
}



float rigidityCutoffNegativeUpperdif(float latS, float phiS, float thetaP, float phiP) {

  static bool first=true;
  static GeomHashNamedEnsemble *hashE=0;
#pragma omp threadprivate(first,hashE)

  bool isPositive = false;
  int hashIndex = rigidityCutoff::negativeUpperdif;
  int hashType = rigidityCutoff::cutoffdif;
  int valType = rigidityCutoff::meanWeightedCorrected;

  if (first) {
    first=false;

    TString fileName = Form("%s/%s", (const char *)rigidityCutoff::filePath, (const char *)rigidityCutoff::fileName[hashIndex]);
    TFile *fHash = TFile::Open(fileName);
    if (!fHash) {
      cout << "rigidityCutoffNegativeUpperdif-E-FileNotFound : " << rigidityCutoff::fileName[hashIndex] << endl;
      return -FLT_MAX;
    }

    hashE = (GeomHashNamedEnsemble *)fHash->Get(rigidityCutoff::hashName[hashIndex]); 
    if (!hashE) {
      cout << "rigidityCutoffNegativeUpperdif-E-HashNotFound : " << rigidityCutoff::hashName[hashIndex] << endl;
      return -FLT_MAX;
    }
  }

  if (!hashE) return -FLT_MAX;

  return rigidityCutoffInterpolated(latS, phiS, thetaP, phiP,  
				    isPositive, hashType, hashE, valType);
}



float rigidityCutoffPositivePenumbra(float latS, float phiS, float thetaP, float phiP) {

  static bool first=true;
  static GeomHashNamedEnsemble *hashE=0;
#pragma omp threadprivate(first,hashE)

  bool isPositive = true;
  int hashIndex = rigidityCutoff::positivePenumbra;
  int hashType = rigidityCutoff::penumbra;
  int valType = rigidityCutoff::meanTruncated;

  if (first) {
    first=false;

    TString fileName = Form("%s/%s", (const char *)rigidityCutoff::filePath, (const char *)rigidityCutoff::fileName[hashIndex]);
    TFile *fHash = TFile::Open(fileName);
    if (!fHash) {
      cout << "rigidityCutoffPositivePenumbra-E-FileNotFound : " << rigidityCutoff::fileName[hashIndex] << endl;
      return FLT_MAX;
    }

    hashE = (GeomHashNamedEnsemble *)fHash->Get(rigidityCutoff::hashName[hashIndex]); 
    if (!hashE) {
      cout << "rigidityCutoffPositivePenumbra-E-HashNotFound : " << rigidityCutoff::hashName[hashIndex] << endl;
      return FLT_MAX;
    }
  }

  if (!hashE) return FLT_MAX;

  return rigidityCutoffInterpolated(latS, phiS, thetaP, phiP, 
				    isPositive, hashType, hashE, valType);
}



float rigidityCutoffNegativePenumbra(float latS, float phiS, float thetaP, float phiP) {

  static bool first=true;
  static GeomHashNamedEnsemble *hashE=0;
#pragma omp threadprivate(first,hashE)

  bool isPositive = false;
  int hashIndex = rigidityCutoff::negativePenumbra;
  int hashType = rigidityCutoff::penumbra;
  int valType = rigidityCutoff::meanTruncated;

  if (first) {
    first=false;

    TString fileName = Form("%s/%s", (const char *)rigidityCutoff::filePath, (const char *)rigidityCutoff::fileName[hashIndex]);
    TFile *fHash = TFile::Open(fileName);
    if (!fHash) {
      cout << "rigidityCutoffNegativePenumbra-E-FileNotFound : " << rigidityCutoff::fileName[hashIndex] << endl;
      return -FLT_MAX;
    }

    hashE = (GeomHashNamedEnsemble *)fHash->Get(rigidityCutoff::hashName[hashIndex]); 
    if (!hashE) {
      cout << "rigidityCutoffNegativePenumbra-E-HashNotFound : " << rigidityCutoff::hashName[hashIndex] << endl;
      return -FLT_MAX;
    }
  }

  if (!hashE) return -FLT_MAX;

  return rigidityCutoffInterpolated(latS, phiS, thetaP, phiP, 
				    isPositive, hashType, hashE, valType);
}





float rigidityCutoffInterpolated(float latS, float phiS, float thetaP, float phiP,
				 bool isPositive, int hashType, GeomHashNamedEnsemble *hashE, int valType) {

   const float rISS[] = { 0.030, 0.010 };
   const float rPart[]= { 0.100, 0.010 };

   const int nCuts=4;    
   const float d2DCut[nCuts] = { 1.5, 2.0, 2.5, 3.0 };

   const float cosThetaMax = 0.7841;

   int i2DCut = -1;
   int minHashes = 2;

   float mean;
   float meant=-1, meanc=-1, meanwc=-1;

   float s0[nCuts], s1[nCuts], sc1[nCuts];
   float sw0[nCuts], swc1[nCuts];
   float distanceClosestHash=1.e10;
   float valClosestHash=-1, valcClosestHash=-1;

   float value = isPositive? FLT_MAX : -FLT_MAX;

   float thetaS = TMath::Pi()/2 - latS;

   float x[4];
   x[0] = ( 1 + cos(thetaS)/cosThetaMax ) / 2;
   x[1] = ( phiS>0 ? phiS : phiS+2*TMath::Pi() ) / (2*TMath::Pi());
   x[2] = ( 1 + cos(thetaP) ) / 2;
   x[3] = ( phiP>0 ? phiP : phiP+2*TMath::Pi() ) / (2*TMath::Pi());
 
   for (int k=0; k<nCuts; k++) {
     s0[k]=0; s1[k]=0; sc1[k]=0;
     sw0[k]=0; swc1[k]=0;
   }

   if (!hashE) return value;
   hashE->Eval(x[0],x[1],x[2],x[3]);
   mean = hashE->MeanValue;

   int numHashes = hashE->numHashes();
   for (int i=0; i<numHashes; i++) {
     GeomHashNamed &h = hashE->getHash(i);
     GeomHashNamed *hash = &h;
     int index = hash->get(x[0],x[1],x[2],x[3]);
     float val = hash->getMean(index);
     float *y = hash->getTemplate(index);
     if (!y) continue;
     float dISS = sqrt( pow(y[0]-x[0],2) + pow(y[1]-x[1],2) );
     float dPart= sqrt( pow(y[2]-x[2],2) + pow(y[3]-x[3],2) );
     float d2D  = sqrt( pow(dISS/rISS[hashType],2) +
			pow(dPart/rPart[hashType],2) );
     float wght = 1. / pow(d2D>1.e-3 ? d2D : 1e-3, 2);
     float corr = rigidityCutoffStormerCorrection(isPositive,x,y);
     float valw = wght*val;
     float valc = corr*val;
     float valwc= corr*valw;
     if (d2D<distanceClosestHash) {
       distanceClosestHash = d2D;
       valClosestHash = val;
       valcClosestHash = valc;
     }
       
     for (int k=0; k<nCuts; k++) {
       if (d2D<d2DCut[k]) {
	 s0[k]+=1;
	 s1[k]+=val;
	 sc1[k]+=valc;
	 sw0[k]+=wght;	
	 swc1[k]+=valwc;
       }
     }
   }

   for (int k=0; k<nCuts; k++) {
     if (s0[k]<minHashes) continue;
     i2DCut = k;
     meant  = s1[k]/s0[k];
     meanc  = sc1[k]/s0[k];
     meanwc = swc1[k]/sw0[k];
     break;
   }

   if (i2DCut<0) {
     meant  = valClosestHash;
     meanc  = valcClosestHash;
     meanwc = valcClosestHash;
   }

   switch (valType) {
   case rigidityCutoff::mean:
     value = mean;
     break;
   case rigidityCutoff::meanTruncated:
     value = meant;
     break;
   case rigidityCutoff::meanCorrected:
     value = meanc;
     break;
   case rigidityCutoff::meanWeightedCorrected:
   default:
     value = meanwc;
     break;
   }

   return value;
}


float rigidityCutoffStormerCorrection(bool isPositive, float *xf, float *xi) {

  int time=0;
  const float radS=6771.2e5;
  const float cosThetaMax = 0.7841;

  float latSf, phiSf, thetaPf, phiPf;
  float radSGMf, latSGMf, phiSGMf, rgCutPf, rgCutNf;

  float latSi, phiSi, thetaPi, phiPi;
  float radSGMi, latSGMi, phiSGMi, rgCutPi, rgCutNi;

  latSf   = asin(cosThetaMax*(2*xf[0]-1));
  phiSf   = 2*TMath::Pi()*xf[1];
  thetaPf = acos(2*xf[2]-1);  
  phiPf   = 2*TMath::Pi()*xf[3];  
  coord(time, radS, latSf, phiSf, thetaPf, phiPf, 
	radSGMf, latSGMf, phiSGMf, rgCutPf, rgCutNf);

  latSi   = asin(cosThetaMax*(2*xi[0]-1));
  phiSi   = 2*TMath::Pi()*xi[1];
  thetaPi = acos(2*xi[2]-1);  
  phiPi   = 2*TMath::Pi()*xi[3];  
  coord(time, radS, latSi, phiSi, thetaPi, phiPi, 
	radSGMi, latSGMi, phiSGMi, rgCutPi, rgCutNi);

  return isPositive? rgCutPf/rgCutPi : rgCutNf/rgCutNi; 
}

