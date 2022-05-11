#ifndef __LikelihoodVar_h__
#define __LikelihoodVar_h__

#include "PDF2.h"

#include "Rtypes.h"
#include "TNamed.h"

#include <cmath>
#include <cstdio>
#include <string>

//! Class to calculate the p.d.f. value in the form p(x(X,Y,PID),y(X,Y,PID))
class LikelihoodVar : public TNamed {

public:
  //! X
  Float_t *FloatX;
  //! Y
  Float_t *FloatY;
  //! Y
  Int_t *IntY;
  //! is y an int?
  bool IsYInt;
  //! rule to re-calculate x = x(X,Y,PID)
  std::string RuleX;
  //! rule to re-calculate y = y(X,Y,PID)
  std::string RuleY;
  //! pdf
  PDF2 *Pdf2;

  //! 2D Interpolation Type (as in PDF2)
  int InterpType;
  //! Bell Type (as in PDF1 and PDF2)
  int BellType;
  //! Fit Type (as in PDF1 and PDF2)
  int FitType;

  //! c-tor
  LikelihoodVar() { Clear(); }
  //! c-tor
  LikelihoodVar(std::string name, Float_t *x, Int_t *y, PDF2 *pdf2, std::string rulex, std::string ruley,
                int interp_type, int bell_type, int fit_type) {
    Clear();
    Set(name, x, y, pdf2, rulex, ruley, interp_type, bell_type, fit_type);
  }
  //! c-tor
  LikelihoodVar(std::string name, Float_t *x, Float_t *y, PDF2 *pdf2, std::string rulex, std::string ruley,
                int interp_type, int bell_type, int fit_type) {
    Clear();
    Set(name, x, y, pdf2, rulex, ruley, interp_type, bell_type, fit_type);
  }
  //! d-tor
  virtual ~LikelihoodVar() { Clear(); }
  //! set
  void Set(std::string name, Float_t *x, Int_t *y, PDF2 *pdf2, std::string rulex, std::string ruley, int interp_type,
           int bell_type, int fit_type);
  //! set
  void Set(std::string name, Float_t *x, Float_t *y, PDF2 *pdf2, std::string rulex, std::string ruley, int interp_type,
           int bell_type, int fit_type);
  //! set X
  void SetX(Float_t *x) { FloatX = x; }
  //! set Y
  void SetY(Float_t *y) {
    FloatY = y;
    IntY = 0;
  }
  //! set Y
  void SetY(Int_t *y) {
    IntY = y;
    FloatY = 0;
  }
  //! clear
  void Clear(Option_t *opt = "");
  //! print
  void Print(int verbose) const;
  void Print(Option_t *opt = "") const { Print((int)0); };
  //! get X
  double GetX() const { return double(*FloatX); }
  //! get Y
  double GetY() const {
    if (FloatY)
      return double(*FloatY);
    if (IntY)
      return double(*IntY);
    return 0;
  }
  //! eval x(X,Y,PID)
  double EvalX(int pid = -1);
  //! eval y(X,Y,PID)
  double EvalY(int pid = -1);
  //! eval p(x(X,Y,PID),y(X,Y,PID))
  double Eval(int pid = -1);

  //! mass number of Geant3 PID
  static double __A(int pid);
  //! charge number of Geant3 PID
  static double __Z(int pid);
  //! mass of Geant3 PID
  static double __M(int pid);

  ClassDef(LikelihoodVar, 1);
};

#endif
