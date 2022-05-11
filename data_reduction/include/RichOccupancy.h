#ifndef __RichOccupancy_h__
#define __RichOccupancy_h__

#include "TNamed.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "Math/Math.h"
#include "Math/PdfFuncMathCore.h"
#include "Math/SpecFuncMathCore.h"

#include <string>
#include <cstdio>
#include <vector>
#include <algorithm> 

class RichOccupancy : public TNamed {

 public:

  std::vector<int> bad_chn; 

  TH1D* counters; 
  TH1D* nhit_vs_chn[3];
  TH1D* npho_vs_chn[3];

  RichOccupancy();
  ~RichOccupancy();
  void FillCounters(int index) { if (counters) counters->Fill(index); } 
  void FillChannel(int index, int channel, float np) { 
    if (nhit_vs_chn[index]) nhit_vs_chn[index]->Fill(channel);
    if (npho_vs_chn[index]) npho_vs_chn[index]->Fill(channel,np);
  }
  bool IsGood(int channel) { return (std::find(bad_chn.begin(),bad_chn.end(),channel)==bad_chn.end()); }
  void AddBadChannel(int channel) { if (IsGood(channel)) bad_chn.push_back(channel); }
  int GetNBad() { return (int)bad_chn.size(); } 
  int EvalThreshold(TH1D* nhit_vs_chn, int N_evt, double nsigma);
  void CreateOccupancyTable(double nsigma = 3.);
  ClassDef(RichOccupancy,1);
};

#endif
