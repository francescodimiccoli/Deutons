#ifndef __FillNtp_h__
#define __FillNtp_h__

#include "Ntp.h"
#include "Tools.h"
#include "Fiducial.h"
#include "RichOccupancy.h"
#include "RichBDT.h"

#include "TrCharge.h"
#include "root_RVSP.h"
#include "amschain.h"
#include "HistoMan.h"
#include "Tofrec02_ihep.h"
#include "TrdSCalib.h"
#include "TkSens.h"
#include "TrdKCluster.h"
#include "AntiPG.h"
#include "TrRecon.h"
#include "TrReconQ.h"
#include "TrMass.h"
#include "bcorr.h"
#include "VPreScaler.h"
#include "GeomHashes.h"
#include "Cutoff.h"
#include "RichTools.h"

#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TParticle.h>
#include <TDatabasePDG.h>

#include <string>
#include <signal.h>
#include <cstdio>
#include <ctime>

using namespace std;


class Analysis : public AMSEventR {

 public:

  //! is simulation?
  bool          isMC;
  //! first entry to be processed
  Long64_t      firstentry;
  //! single event to be processed
  Long64_t      singleevent;
  //! perform compact ealuation and storage
  int           compact;   

  /////////////////////////////
  // Output
  /////////////////////////////

  TFile*        pFileOutput;

  TTree*        pTreeRTI;
  TTree*        pTreeFile;
  TTree*        pTreeEvent;
  TTree*        pTreeProc;
  TTree*        pTreeComp;

  /////////////////////////////
  // Data structures
  /////////////////////////////

  RTIInfo       rtiInfo;
  FileInfo      fileInfo;
  FileMCInfo    fileMCInfo;
  ProcInfo      procInfo;
  NtpSHeader    ntpSHeader;
  NtpHeader     ntpHeader;
  NtpMCHeader   ntpMCHeader;
  NtpTrd        ntpTrd;
  NtpTof        ntpTof;
  NtpTracker    ntpTracker;
  NtpRich       ntpRich;
  NtpEcal       ntpEcal;
  NtpAnti       ntpAnti;
  NtpStandAlone ntpStandAlone;
  NtpCompact    ntpCompact;

  /////////////////////////////
  // Selected objects
  /////////////////////////////

  BetaR*        pBeta;
  BetaHR*       pBetaH;
  TrTrackR*     pTrTrack;
  TrTrackR*     pTrTrack2;
  TrdTrackR*    pTrdTrack;
  TrdHTrackR*   pTrdHTrack;
  RichRingR*    pRichRing;
  RichRingBR*   pRichRingB;
  EcalShowerR*  pEcalShower;

  BetaHR*       pBetaH_SA;
  TrdTrackR*    pTrdTrack_SA;

  std::unique_ptr<RichBDTMgr> RichBDT;

  /////////////////////////////
  // RICH calibration
  /////////////////////////////

  static GHBManager*       pRichCorr;
  static GeomHashEnsemble* pRichUnifAgl;
  static GeomHashEnsemble* pRichUnifNaf;

  RichOccupancy* pRichOcc;

  /////////////////////////////
  // Prescaling
  /////////////////////////////

  static int    pres_strategy;
  VPreScaler    vPrescale;

  void InitPrescaler();
  bool accepted;
  bool a_track;
  bool a_betah;
  bool is_associated;
  bool is_downgoing;
  bool is_positive;
  bool is_tofh;
  bool is_rich;
  bool is_in_lay1;
  bool is_zgt1;
  bool is_zgt2;
  bool has_tkl1;
  bool is_naf;

  bool Accepted()     const { return accepted; }
  bool HasTrack()     const { return a_track; }
  bool HasBeta()      const { return a_betah; }
  bool IsAssociated() const { return is_associated; }
  bool IsDowngoing()  const { return is_downgoing; }
  bool IsNegative()   const { return !is_positive; }
  bool Is_tofh()      const { return is_tofh; }
  bool Is_rich()      const { return is_rich; }
  bool IsInsideL1()   const { return is_in_lay1; }
  bool HasTkL1()      const { return has_tkl1; }
  bool IsZgt1()       const { return is_zgt1; }
  bool IsZgt2()       const { return is_zgt2; }
  bool IsNaF()        const { return is_naf; }

  bool IsCompactAnalysis();
  bool IsCompactStandalone();

  /////////////////////////////
  // Reference refit
  /////////////////////////////

  static int   char_ref_fit; ///< charge for reference fit
  static float mass_ref_fit; ///< mass for reference fit
  static int   algo_ref_fit; ///< fitting algorithm for reference fit
  static int   patt_ref_fit; ///< patter for reference fit
  static int   refi_ref_fit; ///< refitting code for reference fit
  static int   id_ref_fit;   ///< id of the reference fit

  /////////////////////////////
  // Constants
  /////////////////////////////

  static float m_p;
  static float m_d;
  static float tracker_layers_z[9];
  static float rich_radiator_z;
  static float rich_pmt_plane_z;
  static float top_z;

  /////////////////////////////
  // Basic methods
  /////////////////////////////

  Analysis(bool is_mc);
  ~Analysis();
  void ClearData();
  void ClearSelection();

  void UBegin();           ///< Init output
  void ProcessRTI();       ///< Fill RTI info
  void ProcessFile();      ///< Fill file info
  bool UProcessCut();      ///< Cut, select and prescale
  void UProcessFill();     ///< Fill ntuple
  void UTerminate();       ///< Save
  void DefineStandAlone(); ///< Select standalone objects

  /////////////////////////////
  // Fillers
  /////////////////////////////

  void FillRTI();

  void ClearNtp();
  void FillNtpSHeader();
  void FillNtpHeader();
  void FillNtpMCHeader();
  void FillNtpTrd();
  void FillNtpTof();
  void FillNtpTracker();
  void FillNtpRich();
  void FillNtpTrackerScattering();
  void FillNtpEcal();
  void FillNtpAnti();
  void FillNtpStandAlone();

  void ClearCompact();
  void FillNtpCompact(bool use_already_calculated);
  void FillNtpCompactStandAlone(bool use_already_calculated);

  int Loc2Gl(double Theta, double Phi, double &ThetaGl, double &PhiGl);
};

#endif
