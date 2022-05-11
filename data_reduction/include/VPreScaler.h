//Author P. Zuccon -- MIT
#ifndef VPrescaler_h
#define VPrescaler_h
#include "Rtypes.h"
#include <vector>
#include <map>
#include <string>
#include <bitset>
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "VPSArchive.h"
#include "TFile.h"

class Analysis;

#define VPSEV Analysis

class VPreScaler {
protected:

  typedef  bool (VPSEV::* vfunc) () const;

  class VPSCondition {
  public:
    std::string  name;
    vfunc func;
    int prsFactor0;
    int prsFactor1;
  public:

    VPSCondition(): name("dummy"), func(0), prsFactor0(1), prsFactor1(1) {}
    VPSCondition(const char *nam, vfunc ff, int prs0, int prs1): name(nam), func(ff), prsFactor0(prs0), prsFactor1(prs1) {}

    virtual ~VPSCondition() {}

    bool Eval(const VPSEV& ev) const {return (ev.*func)();}
    const std::string& GetName() const {return name;}
    void SetName(const char * nn) {name = nn;}
    void Print() const {printf("Cond: %s  \t\tPrsTrue: %4d PrsFalse %4d\n", name.c_str(), prsFactor1, prsFactor0);}
  };


  class node{
  public:
    node* up;
    node* dr;
    node* dl;
    std::string name;
    std::string basename;
    long int ps;
    long int pso;
    long int psdr;
    long int psdl;
  
    node( std::string nn, std::string nn2,int  pdr=1,int pdl=1):
      up(0),dr(0),dl(0),name(nn),basename(nn2),
      ps(1),pso(-1),psdr(pdr),psdl(pdl){}

    node(const char* nn,const char* nn2,int  pdr=1,int pdl=1):
      up(0),dr(0),dl(0),name(nn),basename(nn2),ps(1),pso(-1),
      psdr(pdr),psdl(pdl){}

    virtual ~node(){if(dr) delete dr; if(dl) delete dl;up=0;}

    node* AddNodeRight(node * rr){dr=rr; rr->ps=ps*psdr; rr->up=this; return rr;}

    node* AddNodeLeft(node * rr) {dl=rr ;rr->ps=ps*psdl; rr->up=this; return rr;}
    void PrintNode(int rec=0);
    void PrintfNode(FILE* ff, int flag=0);
  };





protected:

  std::vector<VPSCondition> cond;
  std::map<unsigned long, VPSCategory*> categ;
  bool Consistency;
  bool TestConsistency;

  
  std::string GetName(vpbitset _cid, vpbitset _cmas) const;
  long int  GetPrs(vpbitset _cid, vpbitset _cmas) const ;
  void PrintTree(node* top) const;


public:


  VPreScaler(): Consistency(0), TestConsistency(0) {}
  virtual ~VPreScaler() {}
  VPSCategory * FindCat(long int evpatt) const;
  
  void AddCondition(const char* name,vfunc ff, int prs0, int prs1) {
    cond.insert(cond.begin(), VPSCondition(name, ff, prs0, prs1));
  }

  void AddCategory(const char* cid, const char* cmas,long int override=-1) {
    AddCategory(vpbitset(cid), vpbitset(cmas),override);
  }

  void AddCategory(vpbitset cid, vpbitset cmas,long int override=-1);

  bool CheckConsistency(int verbose=0);

  long int PrescaleEvent(const VPSEV& ev);

  void Report(int lev = 0) const {
    long int tt=0,ttr=0;
    printf("Report of Prescaling in %3lu Categories.\n", categ.size());
    for (const auto &it : categ){
      tt+=it.second->GetTot();
      ttr+=it.second->ScalerGot;
    }
    
    for (const auto &it : categ)
      it.second->Report2(ttr,tt);
    printf("          Total      %8ld  %5.2f%%  %8ld\n",ttr,(ttr*100.)/tt,tt);
  }  


  void PrintCateg(int lev = 0) const {
    printf("Categories are : %3lu\n", categ.size());
    for (const auto &it : categ)
      it.second->Print(lev);
  }
  void PrintCond() const {
    printf("Binary Conditions are : %2lu\n", cond.size() );
    for (int ii = cond.size() - 1; ii >= 0 ; ii--)
      cond[ii].Print();
  }

  void BuildTree();
  void WriteArchive(const char * filename) const{
    TFile * ff=TFile::Open(filename,"RECREATE");
    WriteArchive(ff);
    ff->Close();
  }
  void WriteArchive(TFile* tff) const{
    tff->cd();
    VPSArchive * ar= new VPSArchive();
    for (const auto &it : categ)
      ar->Add(*(it.second));
    ar->Write();
  }
};



#endif
