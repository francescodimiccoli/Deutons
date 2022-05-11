#ifndef __Classifier_h__
#define __Classifier_h__

#include "Rtypes.h"

#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"

#include "Ntp.h"

#include "PDF2.h"
#include "PDF2DB.h"
#include "Likelihood.h" 

#include <map>
#include <algorithm> 

class ClassifierVariableTransformation {

 public:

  //! transformation parameter
  Float_t* fDependency;
  //! transformation PDF
  PDF2*    fPDF2;
  //! transformation type (0: none, 1: reference, 2: normal) 
  Int_t    fType;

  //! c-tor 
  ClassifierVariableTransformation() { Clear(); } 
  //! c-tor
  ClassifierVariableTransformation(PDF2* pdf2, Float_t* dependency, Int_t type) : fDependency(dependency), fPDF2(pdf2), fType(type) {}
  //! d-tor
  ~ClassifierVariableTransformation() { Clear(); } 
  //! clear
  void Clear() { fDependency = 0; fPDF2 = 0; fType = 0; }    
  //! eval transformation 
  Float_t Eval(double y) {
    if ( (!fPDF2)||(fDependency==0)||(fType<1)||(fType>2) ) return y; 
    return fPDF2->Transform((double)*fDependency,(double)y,PDF2::kSimple,(fType==1)?PDF1::kUser:PDF1::kNormal,PDF1::kSpline);
  }
};

class ClassifierVariable {

 public: 

  //! name
  string fName;
  //! input value
  std::map<int,Float_t> fInput;
  //! transformation list (by detector ID) 
  std::map<int,ClassifierVariableTransformation*> fTransformationList;
  //! output value
  std::map<int,Float_t> fOutput;

  //! c-tor
  ClassifierVariable() { fInput.clear(); fTransformationList.clear(); fOutput.clear(); }  
  //! d-tor
  ~ClassifierVariable() {
    for (auto it=fTransformationList.begin(); it!=fTransformationList.end(); it++) delete it->second;
    fInput.clear();
    fTransformationList.clear();
    fOutput.clear();
  }
  //! set name
  void SetName(const char* name) { fName.append(name); }     
  //! set transformation
  void AddTransformation(int index, PDF2* pdf2, Float_t* dependency, Int_t transform_type) {
    if (fTransformationList.find(index)!=fTransformationList.end()) delete fTransformationList[index];
    fInput[index] = 0; 
    fTransformationList[index] = new ClassifierVariableTransformation(pdf2,dependency,transform_type); 
    fOutput[index] = 0; 
  }
  //! set value
  void SetValue(Float_t value, int id) { fInput[id] = value; }
  //! apply transformation
  void ApplyTransformation() {
    for (auto it=fTransformationList.begin(); it!=fTransformationList.end(); it++) { 
      int index = it->first; 
      fOutput[index] = it->second->Eval(fInput[index]); 
      if ( (std::isnan(fOutput[index]))||(std::isinf(fOutput[index])) ) fOutput[index] = 0; // NaN/Inf catcher
    }
  }
  //! print
  void Print() {
    for (auto it=fTransformationList.begin(); it!=fTransformationList.end(); it++) {
      int index = it->first;
      printf("%20s ID=%2d fInput=%7.3f fOutput=%7.3f fDependency=%7.3f fType=%1d\n",fName.c_str(),index,fInput[index],fOutput[index],*(it->second->fDependency),it->second->fType);
//    if (it->second->fPDF2) it->second->fPDF2->Print();
    }
  }  
};

class ClassifierData {

 public:

  // TOF
  ClassifierVariable _tof_logchisqtn;
  ClassifierVariable _tof_logchisqcn;
  ClassifierVariable _tof_logzprob;
  ClassifierVariable _tof_logqasym;
  ClassifierVariable _tof_logdbeta;
  ClassifierVariable _tof_logisolat;
  ClassifierVariable _tof_qu_unc;
  ClassifierVariable _tof_ql_unc;

  // TRD
  ClassifierVariable _trd_logphe;
  ClassifierVariable _trd_logep;
  ClassifierVariable _trd_nhit;
  ClassifierVariable _trd_vertex;
  ClassifierVariable _trd_amplonpath;

  // Tracker 
  ClassifierVariable _trk_logchisqkax;
  ClassifierVariable _trk_logchisqkay;
  ClassifierVariable _trk_logchisqchx;
  ClassifierVariable _trk_logchisqchy;
  ClassifierVariable _trk_logscatt34x;
  ClassifierVariable _trk_logscatt34y;
  ClassifierVariable _trk_logscatt56x;
  ClassifierVariable _trk_logscatt56y;
  ClassifierVariable _trk_logdinvru;
  ClassifierVariable _trk_logdinvrd;
  ClassifierVariable _trk_logcc;
  ClassifierVariable _trk_logdinvr;
  ClassifierVariable _trk_logdinvrms;
  ClassifierVariable _trk_logdinvrck;
  ClassifierVariable _trk_nyhits; 
  ClassifierVariable _trk_nxhits; 
  ClassifierVariable _trk_n; 
  ClassifierVariable _trk_r;
  ClassifierVariable _trk_q;
  ClassifierVariable _trk_qyl1;
  ClassifierVariable _trk_qinn_unc;
  ClassifierVariable _trk_logqyasym;
  ClassifierVariable _trk_logedep_l2_offtrack;
  ClassifierVariable _trk_logedep_inn_offtrack[2];

  // RICH 
  ClassifierVariable _rich_np_exp;
  ClassifierVariable _rich_exp_res;
  ClassifierVariable _rich_exp_rms;
  ClassifierVariable _rich_lognp_min;
  ClassifierVariable _rich_logdbeta_lip;
  ClassifierVariable _rich_logqpmtcons;
  ClassifierVariable _rich_q;
  ClassifierVariable _rich_prob;
  ClassifierVariable _rich_npmt;
  ClassifierVariable _rich_logdq_lip;
  ClassifierVariable _rich_ratio;
  ClassifierVariable _rich_logdbeta;
  ClassifierVariable _rich_nhit_uncorr;
  ClassifierVariable _rich_nhit;
  ClassifierVariable _rich_nhit_refl;
  ClassifierVariable _rich_nhit_hyp[2];
  ClassifierVariable _rich_nhit_notused[2];
  ClassifierVariable _rich_tot_hit[2][5];
  ClassifierVariable _rich_logdist;
  ClassifierVariable _rich_nclus;

  //! ancillary  
  Float_t _logk[3];
  Float_t _logr[3];

  //! c-tor
  ClassifierData() {}
  //! d-tor
  virtual ~ClassifierData() {}
  //! fill variables 
  void FillData(Event* event, int id);
  //! apply transformation
  void ApplyTransformations(); 
};


class ClassifierManager {

 public:

  //! initailization status
  bool initialized;

  //! classifier map   
  std::map<int,TMVA::Reader*> readers;
  //! classifier type
  std::map<int,string> multi;
  //! classifier name
  std::map<int,string> names;  
  //! data
  ClassifierData* data;
  //! pdf2 transform database
  PDF2DB* pdf2db_transform;
  //! pdf2 mass database
  PDF2DB* pdf2db_mass;
  //! mass likelihoods
  Likelihood* likelihood_mass[3][2];

  //! do not recalculate variables and transformations 
  unsigned int current_run; 
  //! do not recalculate variables and transformations 
  int          current_event;

  //! c-tor
  ClassifierManager();
  //! d-tor
  virtual ~ClassifierManager();
  //! init readers
  void Init();
  //! has been initialized?
  bool IsInitialized() { return initialized; } 
  //! calculate and initialize classifier vars
  void SetReaderVariables(TMVA::Reader* reader, int classifier_type, bool fgflag = false);
  //! get classifier
  /* classifier_type = DCBA   
   * D) 0:My training; 1:Francesca's training.
   * C) 0:TOF region; 1:NaF region; 2:Aerogel region.
   * B) 0:TOF; 1:TRD; 2:Traker; 3:RICH.
   * A) 0:positive mass training; 1:high rigidity/beta training; 2:negative mass training.
   */    
  double GetClassifier(Event* event, int classifier_type);
  //! get -log-likelihood 
  /* ll_type = BA   
   * B) 0:TOF region; 1:NaF region; 2:Aerogel region.
   * A) 0:Proton; 1:Deuteron.
   */
  double GetMinusLogLikelihood(Event* event, int ll_type);
};

#endif	
