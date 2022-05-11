#include "MSpline.h"

ClassImp(MSpline);

MSpline::MSpline(const MSpline &orig) {
  Init();
  fN = orig.fN;
  fM = orig.fM;
  fA = new double[fN + 1];
  fB = new double[fN + 1];
  fC = new double[fN + 1];
  fD = new double[fN + 1];
  fX0 = new double[fN + 1];
  for (int i = 0; i < fN + 1; i++) {
    fA[i] = orig.fA[i];
    fB[i] = orig.fB[i];
    fC[i] = orig.fC[i];
    fD[i] = orig.fD[i];
    fX0[i] = orig.fX0[i];
  }
  fBlow = orig.fBlow;
  fBhig = orig.fBhig;
}

void MSpline::Init() {
  fN = 0;
  fM = 0;
  fX = 0;
  fY = 0;
  fA = 0;
  fB = 0;
  fC = 0;
  fD = 0;
  fX0 = 0;
  fBlow = 0;
  fBhig = 0;
}

void MSpline::Clear(Option_t *opt) {
  fN = 0;
  fM = 0;
  if (fX)  { delete [] fX;  fX = 0; }
  if (fY)  { delete [] fY;  fY = 0; }
  if (fA)  { delete [] fA;  fA = 0; }
  if (fB)  { delete [] fB;  fB = 0; }
  if (fC)  { delete [] fC;  fC = 0; }
  if (fD)  { delete [] fD;  fD = 0; }
  if (fX0) { delete [] fX0; fX0 = 0; }
  fBlow = 0;
  fBhig = 0;
}

void MSpline::Set(int nknots) {
  Clear();
  fN = nknots;
  fM = nknots + 1;
  fX = new double[fN];
  fY = new double[fN];
  for (int i = 0; i < fN; i++) {
    fX[i] = 0;
    fY[i] = 0;
  }
  fA = new double[fN + 1];
  fB = new double[fN + 1];
  fC = new double[fN + 1];
  fD = new double[fN + 1];
  fX0 = new double[fN + 1];
  for (int i = 0; i < fN + 1; i++) {
    fA[i] = 0;
    fB[i] = 0;
    fC[i] = 0;
    fD[i] = 0;
    fX0[i] = 0;
  }
  fBlow = 0;
  fBhig = 0;
}

void MSpline::SetY(TGraphErrors *graph) {
  if (!graph)
    return;
  double *y = new double[fN];
  for (int i = 0; i < fN; i++)
    y[i] = graph->Eval(fX[i]);
  SetY(y);
  delete[] y;
}

void MSpline::Set(std::vector<double> &xknots) {
  int N = (int)xknots.size();
  Set(N);
  for (int i = 0; i < fN; i++)
    fX[i] = xknots.at(i);
}

void MSpline::Print(Option_t *opt) const {
  std::printf("MSpline::Print: Nknots=%d\n", fN);
  for (int i = 0; i < fN; i++)
    std::printf("(%f,%f)\n", fX[i], fY[i]);
  for (int i = 0; i < fN + 1; i++)
    std::printf("(%f,%f,%f,%f,%f)\n", fA[i], fB[i], fC[i], fD[i], fX0[i]);
  std::printf("-------------------------------------------------------\n");
}

void MSpline::CalculateCoefficients() {
  // lower extrapolation
  fA[0] = fY[0];
  fB[0] = fBlow;
  fC[0] = 0;
  fD[0] = 0;
  fX0[0] = fX[0];
  // upper extrapolation
  fA[fN] = fY[fN - 1];
  fB[fN] = fBhig;
  fD[fN] = 0;
  fC[fN] = 0;
  fX0[fN] = fX[fN - 1];
  // loop on all other bins
  for (int ibin = 0; ibin < fN - 1; ibin++) {
    int imin = (ibin == 0) ? 1 : 0;
    int imax = (ibin == fN - 2) ? 2 : 3;
    double h[3] = {0, 0, 0};
    double s[3] = {0, 0, 0};
    double S[3] = {0, 0, 0};
    double p[3] = {0, 0, 0};
    double y1[3] = {0, 0, 0};
    for (int i = imin; i < imax; i++) {
      h[i] = fX[ibin - 1 + i + 1] - fX[ibin - 1 + i];
      s[i] = (fY[ibin - 1 + i + 1] - fY[ibin - 1 + i]) / h[i];
      S[i] = (s[i] >= 0) ? 1 : -1;
    }
    // set for first bin
    y1[1] = 0;
    y1[1] = fBlow;
    // set for last bin
    y1[2] = 0;
    y1[2] = fBhig;
    // monotonic correction
    for (int i = imin + 1; i < imax; i++) {
      p[i] = (s[i - 1] * h[i] + s[i] * h[i - 1]) / (h[i] + h[i - 1]);
      y1[i] = (S[i - 1] + S[i]) * (std::min(std::min(std::fabs(s[i - 1]), std::fabs(s[i])), 0.5 * std::fabs(p[i])));
    }
    fA[ibin + 1] = fY[ibin];
    fB[ibin + 1] = y1[1];
    fC[ibin + 1] = (3 * s[1] - 2 * y1[1] - y1[2]) / h[1];
    fD[ibin + 1] = (y1[1] + y1[2] - 2 * s[1]) / (h[1] * h[1]);
    fX0[ibin + 1] = fX[ibin];
  }
}

int MSpline::FindBin(double x) {
  int n = fN;
  if (n<0) return 0;
  if (x<fX[0]) return 0;
  if (x>=fX[n-1]) return n;
  int l = (n-1);
  int f = 0;
  // best guess from uniform binning
  int i = floor((x-fX[0])/((fX[n-1]-fX[0])/n));
  if ((i<f)||(i>l)) i = int(f+(l-f)/2);
  // fast checks
  if      (          (x>=fX[i  ])&&(x<fX[i+1]) ) return i+1;
  else if ( (i>0  )&&(x>=fX[i-1])&&(x<fX[i  ]) ) return i;
  else if ( (i<l-1)&&(x>=fX[i+1])&&(x<fX[i+2]) ) return i+2;
  // go for it
  int how_many = 0;
  while ((!((x>=fX[i])&&(x<fX[i+1])))&&(how_many<100)) {
    if (x< fX[i])   { l = i;   i = int(f+(l-f)/2); }
    if (x>=fX[i+1]) { f = i+1; i = int(f+(l-f)/2); }
    how_many++;
  }
  return i+1;
}

double MSpline::Eval(double x) {
  int i = FindBin(x);
  return fA[i] + fB[i] * (x - fX0[i]) + fC[i] * std::pow(x - fX0[i], 2.) + fD[i] * std::pow(x - fX0[i], 3.);
}

double MSpline::EvalDerivative(double x) {
  int i = FindBin(x);
  return fB[i] + 2 * fC[i] * (x - fX0[i]) + 3 * fD[i] * std::pow(x - fX0[i], 2.);
}

double MSpline::EvalIsoSpline1(double *x, double *par) {
  double Y = 0;
  for (int i = 0; i < fN; i++) {
    Y += par[i];
    fY[i] = Y;
  }
  fBlow = par[fN - 1 + 1];
  fBhig = par[fN - 1 + 2];
  CalculateCoefficients();
  return Eval(x[0]);
}

double MSpline::EvalIsoSpline2(double *x, double *par) {
  double Y = 0;
  for (int i = 0; i < fN; i++) {
    Y += par[i];
    fY[i] = Y;
    fX[i] = par[i + fN];
  }
  fBlow = par[2 * fN - 1 + 1];
  fBhig = par[2 * fN - 1 + 2];
  CalculateCoefficients();
  return Eval(x[0]);
}

double MSpline::EvalIsoSpline3(double *x, double *par) {
  double Y = 0;
  for (int i = 0; i < fN; i++) {
    Y += par[i];
    fY[i] = Y;
  }
  fBlow = par[fN - 1 + 1];
  fBhig = par[fN - 1 + 2];
  CalculateCoefficients();
  double valu = Eval(x[0]);
  double deri = EvalDerivative(x[0]);
  if ((valu > 30) || (valu < -30))
    return 0;
  return deri * std::exp(-0.5 * valu * valu) / std::sqrt(2 * M_PI);
}
