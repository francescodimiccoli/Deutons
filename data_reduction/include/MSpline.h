#ifndef __MSpline_h__
#define __MSpline_h__

#include "TGraphErrors.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

//! Monotonic cubic sp-line.
/* Implementation of M. Steffen, Astron. Astrophys. 239, 443-450 (1990), A simple method for monotonic interpolation in
 * one dimension. The linear extrapolation to angular coefficients decided by user is smooth.
 */
class MSpline : public TNamed {

public:
  //! number of knots
  int fN;
  //! number of knots+1 (only for storage)
  int fM;
  //! X knots
  double *fX; //[fN]
  //! Y knots
  double *fY; //[fN]
  //! vector fA[fN+1] of the spline
  double *fA; //[fM]
  //! vector fB[fN+1] of the spline
  double *fB; //[fM]
  //! vector fC[fN+1] of the spline
  double *fC; //[fM]
  //! vector fD[fN+1] of the spline
  double *fD; //[fM]
  //! vector X0[fN+1] of the spline
  double *fX0; //[fM]
  //! angular coefficient for linear extrapolation before first knot
  double fBlow;
  //! angular coefficient for linear extrapolation after last knot
  double fBhig;

  //! c-tor
  MSpline() { Init(); }
  //! c-tor
  MSpline(const MSpline &orig);
  //! d-tor
  virtual ~MSpline() { Clear(); }
  //! initializer
  void Init();
  //! clear
  void Clear(Option_t *option = "");
  //! print
  void Print(Option_t *option = "") const;
  //! set memory
  void Set(int nknots);
  //! set x knots
  void SetX(const double *xknots) {
    for (int i = 0; i < fN; i++)
      fX[i] = xknots[i];
  }
  //! set y knots
  void SetY(const double *yknots) {
    for (int i = 0; i < fN; i++)
      fY[i] = yknots[i];
  }
  //! set x and y knots from a graph
  void Set(TGraphErrors *graph) {
    Set(graph->GetN());
    SetX(graph->GetX());
    SetY(graph->GetY());
  }
  //! interpolate y knots from a graph given over the x knots
  void SetY(TGraphErrors *graph);
  //! set from vector of x knots
  void Set(std::vector<double> &xknots);
  //! set angular coefficients for extrapolation
  void SetBlow(double blow) { fBlow = blow; }
  //! set angular coefficients for extrapolation
  void SetBhig(double bhig) { fBhig = bhig; }
  //! calculate coefficients
  void CalculateCoefficients();
  //! find bin
  int FindBin(double x);
  //! evaluate
  double Eval(double x);
  //! evaluate derivative
  double EvalDerivative(double x);
  //! interface for fitting (free y0, dy1, ..., blow, bhig)
  double EvalIsoSpline1(double *x, double *par);
  //! interface for fitting (free y0, dy1, ..., x0, x1, ..., blow, bhig)
  double EvalIsoSpline2(double *x, double *par);
  //! bah
  double EvalIsoSpline3(double *x, double *par);

  ClassDef(MSpline, 1);
};

#endif
