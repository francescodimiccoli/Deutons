#ifndef __Tools_h__
#define __Tools_h__

#include "TrCharge.h"
#include "root_RVSP.h"
#include "amschain.h"
#include "HistoMan.h"
#include "Tofrec02_ihep.h"
#include "TkSens.h"
#include "TrdKCluster.h"
#include "AntiPG.h"
#include "TrRecon.h"
#include "TrReconQ.h"
#include "bcorr.h"

#include "TFile.h"

#include <signal.h>
#include <cstdio>

class Tools {

 public:

  static float GetValueDatacard(TFile* file, const char* par_name, int par_numb = -1);
  static int   RichQC(RichRingR *prich);
  static void  Mysort(int N,float *v, int *sort);
  static bool  GetInnerNHits(TrTrackR* track, int &nxy, int &ny);
  static int   axtof(AMSEventR* pev, BetaHR* pbetah, TrTrackR* ptrtk, int *flagoverlap, float *edepoverlap, float *resmeasoverlap, float *resborderoverlap);
  static int   TkSelect(AMSEventR* pEvent, TrTrackR* pTrTrack, BetaHR* pBetaH);
  static void  TrClusterOnLayerJ(int layerJ, AMSEventR* pev, TrTrackR &trtk, int idfit, int *nclx, float *eclx, int *ncly, float* ecly, float &maxedepx, float &d2maxedepx, float &maxedepy, float &d2maxedepy);
  static int   MyEcalShowerHits(AMSEventR &event, EcalShowerR &ecal);
  static bool  MipsTagCut(AMSEventR &event, EcalShowerR &ecal, bool looseCut);
  static bool  IsInsideTRD(float TRDCoo[2][2]);
  static bool  IsInsideTkL1(float  xin, float yin );
  static void  TofUnusedHits(AMSEventR* pev, BetaHR* pBetaH, float clsdt[4], float clsed[4], short int clsn[4]);
};

#endif
