#ifndef __PDF1_h__
#define __PDF1_h__

#include "IsoPoly.h"
#include "MSpline.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TMath.h"
#include "TNamed.h"

#include <cmath>
#include <cstdio>
#include <string>

//! One dimensional p.d.f. class, p(x)
class PDF1 : public TNamed {

public:
  //! Number of availabe bell-shaped distributions
  static const int NBellType = 5;
  //! Number of available QQ fits
  static const int NFitType = 4;

  //! Kind of PDF to be used
  enum BellType {
    kRaw = -3,      // histogram
    kKDE = -2,      // histogram: kernel density estimator
    kAKDE = -1,     // histogram: adaptive kernel estimator
    kLogistic1 = 0, // function: derivative of logistic function
    kNormal = 1,    // function: normal distribution
    kCauchy = 2,    // function: cauchy distribution
    kGumbel = 3,    // function: gumbel distribution
    kUser = 4       // function: user-defined function (via pdf, cdf and idf)
  };
  //! Each PDF can be evaluated with a different QQ fit
  enum FitType {
    kSpline = -1,  // no fit, just spline interpolation
    kPoly = 0,     // simple polynomial fit (does not guarantee unitarity)
    kOddPoly = 1,  // odd polynomial (guarantees unitarity, but not very flexible)
    kIsoPoly = 2,  // generic isotonic polinomial (in 3 versions, guarantees unitarity, but not very flexible)
    kIsoSpline = 3 // isotonic spline
  };

  //! Perform all fits
  static bool PerformFits;
  //! Draw (for checks)
  static bool DrawChecks;
  //! Print (for checks)
  static bool PrintChecks;

  //! Original histogram
  TH1 *Histo;
  //! Lower boundary for raw PDF histogram
  float XMin;
  //! Upper boundary for raw PDF histogram
  float XMax;
  //! Probability Density Function
  TH1 *Pdf[3];
  //! Cumulative Distribution Function
  TH1 *Cdf[3];
  //! Inverse Distribution Function
  TGraphErrors *Idf[3];
  //! QQ plot for every type of bell-shaped distribution
  TGraphErrors *QQ[NBellType];
  //! QQ plot spline
  MSpline *QQSpline[NBellType];
  //! Fit to QQ plot
  TF1 *QQFit[NFitType][NBellType];
  //! Isotonic polinomial fit object
  IsoPoly *QQIsoPoly[NBellType];
  //! Monotonic spline fit (used as isotonic)
  MSpline *QQIsoSpline[NBellType];

  //! User-defined Reference PDF
  TGraph *UserPdf;
  //! User-defined Reference CDF
  TGraph *UserCdf;
  //! User-defined Reference IDF
  TGraph *UserIdf;

  //! max allowed chi2 for fitting
  static double Chi2Max;

  //! c-tor
  PDF1() { Init(); }
  //! c-tor
  PDF1(const PDF1 &orig);
  //! c-tor
  PDF1(TH1 *histo, float xmin = 0, float xmax = 0, TGraph *pdf = 0, TGraph *cdf = 0, TGraph *idf = 0) {
    Init();
    Set(histo, xmin, xmax, pdf, cdf, idf);
  }
  //! d-tor
  virtual ~PDF1() { Clear(); }
  //! Init
  void Init();
  //! clear
  void Clear(Option_t *opt = "");
  //! print
  void Print(Option_t *opt = "") const;
  //! draw
  void Draw(Option_t *opt = "");
  //! set pdf and perform all sort of fits
  void Set(TH1 *histo, float xmin, float xmax, TGraph *pdf, TGraph *cdf, TGraph *idf);
  //! is the calculation done
  bool IsDone(int bell_type, int fit_type) const;
  //! is good (not well defined for the moment)
  bool IsGood(int bell_type, int fit_type) { return IsDone(bell_type, fit_type); }
  //! get chi2/ndf of the fits
  double GetChi2(int bell_type, int fit_type) const;
  //! eval QQ-plot
  double EvalQQ(double x, int bell_type, int fit_type) const;
  //! eval derivative of QQ-plot
  double EvalQQ1(double x, int bell_type, int fit_type) const;
  //! get fitted PDF graph
  TGraph *GetGraphPDF(int bell_type, int fit_type, double xmin, double xmax, int npx = 1000) const;
  //! get fitted CDF graph
  TGraph *GetGraphCDF(int bell_type, int fit_type, double xmin, double xmax, int npx = 1000) const;
  //! get fitted IDF graph
  TGraph *GetGraphIDF(int bell_type, int fit_type, double xmin, double xmax, int npx = 1000) const;
  //! eval fitted PDF
  double EvalPDF(double x, int bell_type, int fit_type) const;
  //! eval fitted CDF
  double EvalCDF(double x, int bell_type, int fit_type) const;
  //! eval fitted IDF
  double EvalIDF(double x, int bell_type, int fit_type) const;
  //! get the QQ plot
  TGraphErrors *GetQQ(int bell_type) const;

  //! create a PDF from an histogram
  static TH1 *_PDF(TH1 *histo);
  //! create a CDF
  static TH1 *_CDF(TH1 *pdf);
  //! create a IDF
  static TGraphErrors *_IDF(TH1 *cdf);
  //! create a Kernel Density Estimator
  static TH1 *_KDE(TH1 *pdf);
  //! create an Adaptive Kernel Estimator
  static TH1 *_AKDE(TH1 *pdf, TH1 *kde);
  //! fit QQ plot with polynomial with linear extrapolation
  static TF1 *_FitPoly(TH1 *pdf, TGraphErrors *qq);
  //! fit QQ plot with add polynomial with linear extrapolation
  static TF1 *_FitOddPoly(TH1 *pdf, TGraphErrors *qq);
  //! fit QQ plot with isotonic polynomial
  static TF1 *_FitIsoPoly(IsoPoly *&isopoly, TH1 *pdf, TGraphErrors *qq, int isopoly_type);
  //! fit QQ plot with isotonic spline (adaptative approach)
  static TF1 *_FitIsoSpline(MSpline *&isospline, TH1 *pdf, TH1 *cdf, TGraphErrors *idf, TH1 *histo, TGraphErrors *qq);

  //! integral
  static void _integral(TH1 *histo, double &integral, double &eintegral, double left, double right);
  //! binomial error (calculated as propagation of errors of the uncorrelated a,b variables for which c=a+b)
  static double _binomial_error(double a, double ea, double c, double ec);
  //! logit (inverse of logistic)
  static double _logit(double y) { return -std::log(1 / y - 1); }
  //! probit (inverse of cumulative normal distribution)
  static double _probit(double y) { return TMath::ErfInverse(2 * y - 1) * std::sqrt(2); }
  //! cauchit (inverse of cumulative Cauchy distribution)
  static double _cauchit(double y) { return std::tan(M_PI * (y - 0.5)); }
  //! gumbit (inverse of Gumbel)
  static double _gumbit(double y) { return std::log(1. / std::log(1. / y)); }
  //! landit (inverse of approx. Landau)
  static double _landit(double y) { return -2 * std::log(std::sqrt(2) * TMath::ErfInverse(-y)); }
  //! polynomial (N=par[0], a_i=par[i-1], p(x) = a_0 + a_1 x + ... + a{N-1} x^{N-1})
  static double _poly_fun(double *x, double *par);
  //! derivative of polynomial
  static double _poly_der_fun(double *x, double *par);
  //! polynomial with linear extrapolationi  (xmin=par[0], xmax=par[1], order=par[2], a_i=par[i-3], p(x) = p(x) = a_0 +
  //! a_1 x + ... + a{N-1} x^{N-1})
  static double _poly_extlin_fun(double *x, double *par);
  //! derivative of polynomial with linear extrapolation
  static double _poly_extlin_der_fun(double *x, double *par);
  //! odd polynomial (N=par[0], a_i=par[i-1], p(x) = a_0 + a_1 x + a_2 x^3 ... + a_{N-1} x^{2N-1})
  static double _odd_poly_fun(double *x, double *par);
  //! derivative of odd polynomial
  static double _odd_poly_der_fun(double *x, double *par);
  //! odd polynomial with linear extrapolation
  static double _odd_poly_extlin_fun(double *x, double *par);
  //! derivative of odd polynomial with linear extrapolation
  static double _odd_poly_extlin_der_fun(double *x, double *par);
  //! interpolate between QQ-plots
  static double _InterpolateQQ(double x, double f, PDF1 *a, PDF1 *b, int bell_type, int fit_type_a, int fit_type_b);
  //! interpolate between two parametrizations
  static double _InterpolateThroughQQ(double x, double f, PDF1 *a, PDF1 *b, int bell_type, int fit_type_a,
                                      int fit_type_b);
  //! interpolate between two distributions using IDF iterations (iterative method, kind of slow)
  static double _InterpolateThroughIDF(double x, double f, PDF1 *a, PDF1 *b, int bell_type_a, int bell_type_b,
                                       int fit_type_a, int fit_type_b);
  //! interpolate between two distributions (PDF unitarity not guaranteed)
  static double _InterpolateFast(double x, double f, PDF1 *a, PDF1 *b, int bell_type_a, int bell_type_b, int fit_type_a,
                                 int fit_type_b);

  ClassDef(PDF1, 1);
};

#endif
