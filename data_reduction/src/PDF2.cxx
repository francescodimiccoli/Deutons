#include "PDF2.h"

ClassImp(PDF2);

PDF2::PDF2(const PDF2 &orig) {
  fN = orig.fN;
  fX = new double[fN];
  fPdf1 = new PDF1 *[fN];
  for (int ix = 0; ix < fN; ix++) {
    fX[ix] = orig.fX[ix];
    fPdf1[ix] = new PDF1(*orig.fPdf1[ix]);
  }
}

void PDF2::Init() {
  fN = 0;
  fX = 0;
  fPdf1 = 0;
}

void PDF2::Clear(Option_t *opt) {
  for (int ix = 0; ix < fN; ix++) {
    delete fPdf1[ix];
    fPdf1[ix] = 0;
  }
  delete[] fX;
  delete[] fPdf1;
  Init();
}

void PDF2::Print(Option_t *opt) const {
  std::printf("PDF2::Print %s\n", GetName());
  for (int ix = 0; ix < fN; ix++) {
    std::printf("bin %3d/%3d x=%7.3f ", ix + 1, fN, fX[ix]);
    fPdf1[ix]->Print();
  }
}

void PDF2::Draw(Option_t *opt) {
  for (int ix = 0; ix < fN; ix++)
    fPdf1[ix]->Draw();
}

void PDF2::Set(TH2 *h, float xmin, float xmax, double thr, TGraph *user_pdf, TGraph *user_cdf, TGraph *user_idf) {
  Clear();
  if (!h)
    return;
  std::printf("PDF2::Set processing %s ...\n", h->GetName());
  int nx = h->GetXaxis()->GetNbins();
  double *x = new double[nx];
  PDF1 **pdf1 = new PDF1 *[nx];
  fN = 0;
  for (int ibin = 0; ibin < nx; ibin++) {
    TH1 *proj = (TH1 *)h->ProjectionY(Form("%s_proj%d", h->GetName(), ibin + 1), ibin + 1, ibin + 1);
    if (proj->Integral() < thr) {
      delete proj;
      continue;
    }
    x[fN] = h->GetXaxis()->GetBinCenter(ibin);
    pdf1[fN] = new PDF1(proj, xmin, xmax, user_pdf, user_cdf, user_idf);
    fN++;
    delete proj;
  }
  fX = new double[fN];
  fPdf1 = new PDF1 *[fN];
  for (int ix = 0; ix < fN; ix++) {
    fX[ix] = x[ix];
    fPdf1[ix] = pdf1[ix];
  }
  delete[] x;
  delete[] pdf1;
}

void PDF2::Set(TH2 *h, float xmin, float xmax, double thr, TH1 *reference) {
  PDF1 *pdf = new PDF1(reference);
  int which = 0;
  TGraph *user_pdf = new TGraph();
  TGraph *user_cdf = new TGraph();
  TGraph *user_idf = new TGraph();
  for (int i = 1; i <= pdf->Pdf[which]->GetXaxis()->GetNbins(); i++)
    user_pdf->SetPoint(i - 1, pdf->Pdf[which]->GetXaxis()->GetBinCenter(i), pdf->Pdf[which]->GetBinContent(i));
  for (int i = 1; i <= pdf->Cdf[which]->GetXaxis()->GetNbins(); i++) {
    user_cdf->SetPoint(i - 1, pdf->Cdf[which]->GetXaxis()->GetBinUpEdge(i), pdf->Cdf[which]->GetBinContent(i));
    user_idf->SetPoint(i - 1, pdf->Cdf[which]->GetBinContent(i), pdf->Cdf[which]->GetXaxis()->GetBinUpEdge(i));
  }
  Set(h, xmin, xmax, thr, user_pdf, user_cdf, user_idf);
  delete user_pdf;
  delete user_cdf;
  delete user_idf;
  delete pdf;
}

int PDF2::FindBin(double x) {
  int n = fN;
  if (n < 1)
    return 0;
  if (x < fX[0])
    return 0;
  if (x >= fX[n - 1])
    return n;
  int l = (n - 1);
  int f = 0;
  // best guess from uniform binning
  int i = floor((x - fX[0]) / ((fX[n - 1] - fX[0]) / n));
  if ((i < f) || (i > l))
    i = int(f + (l - f) / 2);
  // fast checks
  if ((x >= fX[i]) && (x < fX[i + 1]))
    return i + 1;
  else if ((i > 0) && (x >= fX[i - 1]) && (x < fX[i]))
    return i;
  else if ((i < l - 1) && (x >= fX[i + 1]) && (x < fX[i + 2]))
    return i + 2;
  // go for it
  int how_many = 0;
  while ((!((x >= fX[i]) && (x < fX[i + 1]))) && (how_many < 100)) {
    if (x < fX[i]) {
      l = i;
      i = int(f + (l - f) / 2);
    }
    if (x >= fX[i + 1]) {
      f = i + 1;
      i = int(f + (l - f) / 2);
    }
    how_many++;
  }
  return i + 1;
}

void PDF2::GetBoundaries(double x, int &low, int &up, int bell_type, int fit_type) {
  int n = fN;
  low = 0;
  up = 0;
  if (n == 0)
    return;
  int i = FindBin(x);
  if (i == 0) {
    low = 0;
    up = 0;
  } else if (i == n) {
    low = n - 1;
    up = n - 1;
  } else {
    low = i - 1;
    up = i;
  }
  for (; low >= 1; low--)
    if (GetPDF1(low) || (GetPDF1(low)->IsGood(bell_type, fit_type)))
      break; // take the first one anyway ...
  for (; up < n - 1; up++)
    if (GetPDF1(up) || (GetPDF1(up)->IsGood(bell_type, fit_type)))
      break; // take the last one anyway ...
}

double PDF2::Transform(double x, double y, int transform_type, int bell_type, int fit_type) {
  if (transform_type == kSimple) {
    int ia, ib;
    GetBoundaries(x, ia, ib, bell_type, fit_type);
    PDF1 *a = GetPDF1(ia);
    PDF1 *b = GetPDF1(ib);
    if ((a) && (!b))
      return a->EvalQQ(y, bell_type, fit_type);
    else if ((!a) && (b))
      return b->EvalQQ(y, bell_type, fit_type);
    else if ((!a) && (!b))
      return 0;
    double xa = fX[ia];
    double xb = fX[ib];
    double frac = (ia != ib) ? (x - xa) / (xb - xa) : 0;
    return PDF1::_InterpolateQQ(y, frac, a, b, bell_type, fit_type, fit_type);
  } else if (transform_type == kSpline) {
    double *Y = new double[fN];
    for (int ix = 0; ix < fN; ix++)
      Y[ix] = fPdf1[ix]->EvalQQ(y, bell_type, fit_type);
    MSpline spline;
    spline.Set(fN);
    spline.SetX(fX);
    spline.SetY(Y);
    spline.CalculateCoefficients();
    double value = spline.Eval(x);
    delete[] Y;
    return value;
  }
  return 0;
}

double PDF2::Eval(double x, double y, int interp_type, int bell_type, int fit_type) {
  int ia, ib;
  GetBoundaries(x, ia, ib, bell_type, fit_type);
  PDF1 *a = GetPDF1(ia);
  PDF1 *b = GetPDF1(ib);
  if ((a) && (!b))
    return a->EvalPDF(y, bell_type, fit_type);
  else if ((!a) && (b))
    return b->EvalPDF(y, bell_type, fit_type);
  else if ((!a) && (!b))
    return 0;
  double xa = fX[ia];
  double xb = fX[ib];
  double frac = (ia != ib) ? (x - xa) / (xb - xa) : 0;
  if (interp_type == kQQ)
    return PDF1::_InterpolateThroughQQ(y, frac, a, b, bell_type, fit_type, fit_type);
  else if (interp_type == kIDF)
    return PDF1::_InterpolateThroughIDF(y, frac, a, b, bell_type, bell_type, fit_type, fit_type);
  else if (interp_type == kFast)
    return PDF1::_InterpolateFast(y, frac, a, b, bell_type, bell_type, fit_type, fit_type);
  return 0;
}
