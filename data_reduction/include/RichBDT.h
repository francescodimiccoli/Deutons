#ifndef __RichBDT_h__
#define __RichBDT_h__

#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"

#include "Ntp.h"

#ifndef _NTPLIB_
#include "Analysis.h"
class Analysis;
#endif

class RichBDTData {
public:
  RichBDTData(){};
  virtual ~RichBDTData(){};

  int UnusedLayers;
  int NTofUsed;

  float Richtotused;
  float RichPhEl;
  float RICHprob;
  float RICHcollovertotal;
  float RICHLipBetaConsistency;
  float RICHTOFBetaConsistency;
  float RICHChargeConsistency;
  float RICHPmts;
  float RICHgetExpected;
  float tot_hyp_p_uncorr;
  float Bad_ClusteringRICH;
  float NSecondariesRICHrich;
  float HitHVoutdir;
  float HitHVoutrefl;

  // dummy vars. ROOT wants all the spectators to be present...
  float Spectator1;
  float Spectator2;

  void FillData(Event *event);
#ifndef _NTPLIB_
  void FillDataDirect(Analysis *event);
#endif

  void Dump();

private:
};

class RichBDTMgr {
  friend class Event;
#ifndef _NTPLIB_
  friend class Analysis;
#endif

public:
  RichBDTMgr();
  virtual ~RichBDTMgr();

  void Init();
  void SetReaderVariables(TMVA::Reader *reader);

  double GetNafBDT();
  double GetAglBDT();
  double GetNafMinusLogRarity();
  double GetAglMinusLogRarity();

  bool IsInitialized() { return initialized; };

private:
  bool initialized;

  TMVA::Reader *nafReader;
  TMVA::Reader *aglReader;

  RichBDTData *data;
};

#endif
