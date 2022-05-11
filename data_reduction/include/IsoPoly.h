#ifndef __IsoPoly_h__
#define __IsoPoly_h__

#include "TNamed.h"

#include <cmath>
#include <cstdio>

//! Isotonic polynomial class
class IsoPoly : public TNamed {

public:
  //! number of isotonic 2nd degree polynomials (2 times)
  int fNK;
  //! parameters of 2nd degree polynomials
  double *fK; //[fNK]
  //! fType of 2nd degree polynomials (with b=K[2*i] and c=K[2*i+1]):
  /* Form=1: (b^2+c^2) + 2bx + x^2
   * Form=2: 1         + 2bx + (b^2+c^2)x^2
   * Form=3: b^2       + 2bx + (1+c^2)x^2
   */
  int fType;
  //! number of parameters of polynomial derivative description
  int fNP;
  //! P[0] + P[1]*x + ... + P[fNP-1]*x^{fNP-1}
  double *fP; //[fNP]
  //! offset
  double fY0;
  //! multiplicative factor
  double fAlpha;
  //! X min
  double fXMin;
  //! X max
  double fXMax;

  //! c-tor
  IsoPoly();
  //! c-tor
  IsoPoly(const IsoPoly &orig);
  //! d-tor
  virtual ~IsoPoly() { Clear(); }
  //! clear
  void Clear(Option_t *opt = "");
  //! print
  void Print(Option_t *opt = "") const;
  //! calculate derivative coefficients (in the form prod_i^n_k{k_{2*i}^2+k_{2*i+1}^2+2*k_{2*i}*x+x^2})
  void Calculate(int nk, double *k);
  //! eval polynomial
  double EvalPoly(double x) const;
  //! eval polynomial derivative
  double EvalPoly1(double x) const;
  //! eval polynomial with linear extrapolation
  double EvalPolyLinExt(double x) const;
  //! eval polynomial derivative with linear extrapolation
  double EvalPoly1LinExt(double x) const;
  //! eval
  double Eval(double x) const;
  //! eval (fXMin=par[0], fXMax=par[1], fNK=par[2], fY0=par[3], fAlpha=par[4], K[i]=par[i-5])
  double operator()(double *x, double *par);

  //! product of polynomials [N=na-1, M=nb-1, (a0+a1*x+...+aN*x^N)*(b0+b1*x+...+bM*x^M)=c0+c1*x+...+c(N+M)*x^(N+M)]
  static void _poly_product(int na, double *a, int nb, double *b, int &nc, double *&c);

  ClassDef(IsoPoly, 1);
};

#endif
