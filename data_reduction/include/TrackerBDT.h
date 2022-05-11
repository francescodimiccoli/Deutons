#ifndef __TrackerBDT_h__
#define __TrackerBDT_h__

#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"

#include "Ntp.h"

class TrackerBDTData {
public:
  TrackerBDTData(){};
  virtual ~TrackerBDTData(){};

  // data here
  float TrEdepL2On;
  float TrEdepL2Off;
  float TrEdepInnOff0;
  float TrEdepInnOff1;
  float TrChiSqChY;
  float TrChiSqChX;
  float TrChiSqKaY;
  float TrChiSqKaX;
  float TrHalfRUpCh;
  float TrHalfRDwCh;
  float TrHalfRUpKa;
  float TrHalfRDwKa;
  float L2ChRes_x;
  float L2ChRes_y;
  float L2KaRes_x;
  float L2KaRes_y;
  float L34Scat_x;
  float L34Scat_y;
  float L56Scat_x;
  float L56FeetDist;
  float L2FeetDist;
  float TrQMin;
  float TrQAsymm;
  float L56ScatBR_y;

  float Category;
  float Mass;
  float Rigidity;
  float Beta;
  float RichBDT;

  void FillData(Event *event);

  void Dump();

private:
};

class TrackerBDTMgr {
  friend class Event;

public:
  TrackerBDTMgr();
  virtual ~TrackerBDTMgr();

  void Init();
  void SetReaderVariables(TMVA::Reader *reader);

  double GetBDT();
  // double GetXBDT();

  bool IsInitialized() { return initialized; };

private:
  bool initialized;

  TMVA::Reader *Reader;
  // TMVA::Reader *aglReader;

  TrackerBDTData *data;
};

#endif
