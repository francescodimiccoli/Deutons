#ifndef __Likelihood_h__
#define __Likelihood_h__

#include "LikelihoodVar.h"
#include "PDF2DB.h"

#include "Rtypes.h"
#include "TF1.h"
#include "TNamed.h"

#include <cmath>
#include <cstdio>
#include <vector>

//! Likelihood class.
class Likelihood : public TNamed {

public:
  //! pdf DB pointer
  PDF2DB *Pdf2DB;
  //! variables map
  std::vector<LikelihoodVar *> VarList;

  //! Debug flag
  static int Debug;
  //! Minimum allowed probability
  static double MinProb;

  //! c-tor
  Likelihood(PDF2DB *pdf2db);
  //! d-tor
  virtual ~Likelihood() { Clear(); }
  //! clear
  void Clear(Option_t *opt = "");
  //! print
  void Print(int verbose ) const;
  void Print(Option_t *opt = "") const { Print((int)0); };
  //! register variable
  void AddVariable(std::string name, Float_t *x, Float_t *y, std::string pdf_name, std::string rulex, std::string ruley,
                   int interp_type, int bell_type, int fit_type);
  //! register variable
  void AddVariable(std::string name, Float_t *x, Int_t *y, std::string pdf_name, std::string rulex, std::string ruley,
                   int interp_type, int bell_type, int fit_type);
  //! add another pdf and address of variable to be used
  void Add(LikelihoodVar *var) { VarList.push_back(var); }
  //! eval
  double EvalLog(int pid = -1);
  //! eval for a x
  double Eval(Float_t x, int pid = -1);
  //! operator used for maximization
  double operator()(double *x, double *par) { return Eval(x[0], int(par[0])); }
  //! eval with maximization
  double Eval(Float_t &x, double xmin, double xmax, int pid = -1);

  ClassDef(Likelihood, 1);
};

#endif
