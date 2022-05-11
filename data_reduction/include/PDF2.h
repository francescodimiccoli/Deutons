#ifndef __PDF2_h__
#define __PDF2_h__

#include "PDF1.h"

#include "TFile.h"
#include "TH2.h"
#include "TNamed.h"

#include <cmath>
#include <cstdio>

//! Bidimensional p.d.f. class p(x;k), k is morphing parameter between the N pdfs f_i(x) at k_i
class PDF2 : public TNamed {

public:
  //! Kind of interpolations
  enum InterpType {
    kQQ = 0,  // use QQ-plot parametrisation
    kIDF = 1, // use IDF interpolation
    kFast = 2 // just simple interpolation (unitarity not guaranteed)
  };
  //! Transform type
  enum TransformType {
    kSimple = 0, // interpolate inside bin
    kSpline = 1  // use full spline
  };

  //! number of bins
  int fN;
  //! values
  double *fX; //[fN]
  //! PDFs
  PDF1 **fPdf1; //[fN]

  //! c-tor
  PDF2() { Init(); }
  //! c-tor
  PDF2(const PDF2 &orig);
  //! c-tor
  PDF2(TH2 *h, float xmin, float xmax, double thr, TGraph *user_pdf, TGraph *user_cdf, TGraph *user_idf) {
    Init();
    Set(h, xmin, xmax, thr, user_pdf, user_cdf, user_idf);
  }
  //! c-tor
  PDF2(TH2 *h, float xmin, float xmax, double thr, TH1 *reference) {
    Init();
    Set(h, xmin, xmax, thr, reference);
  }
  //! d-tor
  virtual ~PDF2() { Clear(); }
  //! init
  void Init();
  //! clear
  void Clear(Option_t *opt = "");
  //! print
  void Print(Option_t *opt = "") const;
  //! draw
  void Draw(Option_t *opt = "");
  //! set
  void Set(TH2 *h, float xmin, float xmax, double thr, TGraph *user_pdf, TGraph *user_cdf, TGraph *user_idf);
  //! set
  void Set(TH2 *h, float xmin, float xmax, double thr, TH1 *reference);
  //! find bin
  int FindBin(double x);
  //! find good boundaries
  void GetBoundaries(double x, int &low, int &up, int bell_type, int fit_type);
  //! get PDF
  PDF1 *GetPDF1(int i) const { return ((i >= 0) && (i < fN)) ? fPdf1[i] : 0; }
  //! transform y variable using QQ-plot
  double Transform(double x, double y, int transform_type, int bell_type, int fit_type);
  //! eval
  double Eval(double x, double y, int interp_type, int bell_type, int fit_type);

  ClassDef(PDF2, 1);
};

#endif
