// Author P. Zuccon -- MIT
#ifndef VPSArchive_h
#define VPSArchive_h
#include "TObject.h"
#include "TCollection.h"
#include "TH1D.h"
#include <cstdlib>
#include <cstdio>
#include <bitset>
#include <iostream>
#include <vector>
#include <algorithm>

const int maxbit=12;
typedef  std::bitset<maxbit>  vpbitset;

class VPSCategory: public TObject {
private:
  long int prsfactor;
  long int override;
public:
  std::string name;
  vpbitset catid;
  vpbitset catmask;
  long int Scaler;
  long int ScalerGot;
  long int ScalerLost;
  long int GetTot() const{return ScalerLost+ScalerGot;}
  double   GetFrac()const {return (GetTot()>0)? ScalerGot*1./GetTot():0.;}
  long int GetPrf() const{ return (override>=0)?override:prsfactor;}
  
  void OverridePrf(long int aa){ override=aa; return;}
  VPSCategory& operator+=(const VPSCategory& in){
    if(catid!=in.catid || catmask!=in.catmask){
      printf("VPSCategory::operator+= -E-  Cannot sum object with different catid or catmask\n");
      return *this;
    }
    Scaler+=in.Scaler;
    ScalerGot+=in.ScalerGot;
    ScalerLost+=in.ScalerLost;
    return *this;
  }
  
  VPSCategory():
    TObject(),prsfactor(1),override(-1),name("dummy"), catid(-1), catmask(-1),  Scaler(0), ScalerGot(0), ScalerLost(0) {}

  VPSCategory(std::string nn, vpbitset id, vpbitset mask, long int psf):
    TObject(),prsfactor(psf), override(-1), name(nn), catid(id), catmask(mask) ,Scaler(0), ScalerGot(0), ScalerLost(0) {}
  
  virtual void Print(Option_t* opt="") const {
    Print(0);
  }
  
  void Report() const{
    printf("ps: %6ld  Got/Tot: %6ld/%6ld (%5.2g %% ) | %s\n",GetPrf(),ScalerGot,GetTot(),GetFrac()*100.,name.c_str());
    return;
  }

  void Report2(long int totR, long int totT) const{
    printf("ps: %6ld  Got/Tot: %8ld  %5.2f%% /%8ld  %5.2f%%| %s\n",GetPrf(),ScalerGot,ScalerGot/(totR/100.),GetTot(),GetTot()/(totT/100.),name.c_str());
   return;
  }
  void Print(int level) const {
    std::cout << "Name: " << name <<std::endl;
    if (level ==0 || level== 1) std::cout << "catid: " << catid << " catmask: " << catmask ;
    printf("Prescaling: %6ld ", GetPrf());
    if( override>=0)  printf("(Overriden)");
    else              printf("           ");
    if (level == 0)    {printf("\n\n"); return;}
    else{
      if(level<2) printf(" Scaler: %5ld ",Scaler);
      printf(" ScalerGot: %6ld ScalerLost %6ld\n\n", ScalerGot, ScalerLost);
    }
    return;
  }
  virtual ~VPSCategory() {}
  
  ClassDef(VPSCategory, 1);
};


class VPSArchive: public TObject{
public:
  std::vector<VPSCategory> cat;
  VPSArchive():TObject(){}
  virtual ~VPSArchive(){}
  void Add(const VPSCategory& aa){cat.push_back(aa);}
  
  virtual void Print(Option_t *opt="") const {
    Print(0);
  }
  void Report2() const{
    long int tt=0,ttr=0;
    for (const auto &it : cat){
      tt+=it.GetTot();
      ttr+=it.ScalerGot;
    }
     
    for (const auto &it : cat)
      it.Report2(ttr,tt);
    printf("          Total      %8ld  %5.2f%%  %8ld\n",ttr,(ttr*100.)/tt,tt);
  }  
  void Report() const{
    for (const auto &it : cat)
      it.Report();
  }  
  void Print(int level) const{
    for (const auto &it : cat)
      it.Print(level);
  }  

  virtual Long64_t Merge(TCollection *list);

  TH1D* GetHisto(int flag);
  void DrawHisto();
  ClassDef(VPSArchive,1);
};
#endif
