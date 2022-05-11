#ifndef __RichOccupancyAnalysis_h__
#define __RichOccupancyAnalysis_h__

#include "RichOccupancy.h"

#include "root.h"
#include "point.h"

#include "TNamed.h"
#include "TH1.h"
#include "TH2.h"

#include <string>
#include <cstdio>
#include <map>

class pmt_info {

 public:

  int nhit;
  float np; 
  pmt_info() : nhit(0), np(0) {}
};

class RichOccupancyAnalysis : public AMSEventR { 

 public:

  RichOccupancy* pOccupancy;

  RichOccupancyAnalysis();
  ~RichOccupancyAnalysis();

  void UBegin();       
  bool UProcessCut() { return true; } 
  void UProcessFill(); 
  void UTerminate();   
};

#endif
