#include "PDF1.h"

ClassImp(PDF1);

double PDF1::Chi2Max = 3;

const char *const bell_name[PDF1::NBellType] = {"Logis1", "Normal", "Cauchy", "Gumbel", "User"};
const char *const fit_name[PDF1::NFitType] = {"Polynom", "OddPoly", "IsoPoly", "IsoSpli"};

bool PDF1::PerformFits = true;
bool PDF1::DrawChecks = true;
bool PDF1::PrintChecks = true;

PDF1::PDF1(const PDF1 &orig) {
  Init();
  static int copy_index = 0;
  if (orig.Histo) {
    Histo = (TH1 *)orig.Histo->Clone(Form("%s_%d", Histo->GetName(), copy_index));
    Histo->SetDirectory(0);
  }
  XMin = orig.XMin;
  XMax = orig.XMax;
  for (int i = 0; i < 3; i++) {
    if (orig.Pdf[i]) {
      Pdf[i] = (TH1 *)orig.Pdf[i]->Clone(Form("%s_%d", Pdf[i]->GetName(), copy_index));
      Pdf[i]->SetDirectory(0);
    }
    if (orig.Cdf[i]) {
      Cdf[i] = (TH1 *)orig.Cdf[i]->Clone(Form("%s_%d", Cdf[i]->GetName(), copy_index));
      Cdf[i]->SetDirectory(0);
    }
    if (orig.Idf[i])
      Idf[i] = new TGraphErrors(*orig.Idf[i]);
  }
  for (int bell_type = 0; bell_type < NBellType; bell_type++) {
    if (orig.QQ[bell_type])
      QQ[bell_type] = new TGraphErrors(*orig.QQ[bell_type]);
    if (orig.QQSpline[bell_type])
      QQSpline[bell_type] = new MSpline(*orig.QQSpline[bell_type]);
    for (int fit_type = 0; fit_type < NFitType; fit_type++) {
      if (orig.QQFit[fit_type][bell_type])
        QQFit[fit_type][bell_type] = new TF1(*orig.QQFit[fit_type][bell_type]);
    }
    if (orig.QQIsoPoly[bell_type])
      QQIsoPoly[bell_type] = new IsoPoly(*orig.QQIsoPoly[bell_type]);
    if (orig.QQIsoSpline[bell_type])
      QQIsoSpline[bell_type] = new MSpline(*orig.QQIsoSpline[bell_type]);
  }
  if (orig.UserPdf)
    UserPdf = (TGraph *)orig.UserPdf->Clone(Form("%s_%d", orig.UserPdf->GetName(), copy_index));
  if (orig.UserCdf)
    UserCdf = (TGraph *)orig.UserCdf->Clone(Form("%s_%d", orig.UserCdf->GetName(), copy_index));
  if (orig.UserIdf)
    UserIdf = (TGraph *)orig.UserIdf->Clone(Form("%s_%d", orig.UserIdf->GetName(), copy_index));
  copy_index++;
}

void PDF1::Init() {
  Histo = 0;
  XMin = 0;
  XMax = 0;
  for (int i = 0; i < 3; i++) {
    Pdf[i] = 0;
    Cdf[i] = 0;
    Idf[i] = 0;
  }
  for (int bell_type = 0; bell_type < NBellType; bell_type++) {
    QQ[bell_type] = 0;
    QQSpline[bell_type] = 0;
    for (int fit_type = 0; fit_type < NFitType; fit_type++)
      QQFit[fit_type][bell_type] = 0;
    QQIsoPoly[bell_type] = 0;
    QQIsoSpline[bell_type] = 0;
  }
  UserPdf = 0;
  UserCdf = 0;
  UserIdf = 0;
}

void PDF1::Clear(Option_t *opt) {
  delete Histo;
  Histo = 0;
  XMin = 0;
  XMax = 0;
  for (int i = 0; i < 3; i++) {
    delete Pdf[i];
    Pdf[i] = 0;
    delete Cdf[i];
    Cdf[i] = 0;
    delete Idf[i];
    Idf[i] = 0;
  }
  for (int bell_type = 0; bell_type < NBellType; bell_type++) {
    delete QQ[bell_type];
    QQ[bell_type] = 0;
    delete QQSpline[bell_type];
    QQSpline[bell_type] = 0;
    for (int fit_type = 0; fit_type < NFitType; fit_type++) {
      delete QQFit[fit_type][bell_type];
      QQFit[fit_type][bell_type] = 0;
    }
    delete QQIsoPoly[bell_type];
    QQIsoPoly[bell_type] = 0;
    delete QQIsoSpline[bell_type];
    QQIsoSpline[bell_type] = 0;
  }
  delete UserPdf;
  UserPdf = 0;
  delete UserCdf;
  UserCdf = 0;
  delete UserIdf;
  UserIdf = 0;
}

void PDF1::Print(Option_t *opt) const {
  std::printf("PDF1::Print %s ", GetName());
  std::printf("Entries=%10d ", (int)Histo->GetEntries());
  std::printf("KDE=%10.3f ", GetChi2(-2, 0));
  std::printf("AKDE=%10.3f ", GetChi2(-1, 0));
  for (int bell_type = 0; bell_type < NBellType; bell_type++) {
    for (int fit_type = 0; fit_type < NFitType; fit_type++) {
      if (!QQFit[fit_type][bell_type])
        continue;
      std::printf("%s-%s=%10.3f ", bell_name[bell_type], fit_name[fit_type], GetChi2(bell_type, fit_type));
    }
  }
  std::printf("\n");
}

void PDF1::Set(TH1 *histo, float xmin, float xmax, TGraph *pdf, TGraph *cdf, TGraph *idf) {
  Clear();
  if (!histo)
    return;
  if ((pdf) && (cdf) && (idf)) {
    delete UserPdf;
    delete UserCdf;
    delete UserIdf;
    UserPdf = (TGraph *)pdf->Clone(Form("user_%s", pdf->GetName()));
    UserCdf = (TGraph *)cdf->Clone(Form("user_%s", cdf->GetName()));
    UserIdf = (TGraph *)idf->Clone(Form("user_%s", idf->GetName()));
  }
  Histo = (TH1 *)histo->Clone(Form("h_%s", histo->GetName()));
  Histo->SetDirectory(0);
  XMin = xmin;
  XMax = xmax;
  if (XMin < XMax) {
    for (int i = 0; i <= Histo->GetNbinsX() + 1; i++) {
      float x = Histo->GetXaxis()->GetBinCenter(i);
      if ((x > XMin) && (x < XMax))
        continue;
      Histo->SetBinContent(i, 0);
      Histo->SetBinError(i, 0);
    }
  }
  Pdf[0] = _PDF(Histo);
  Pdf[1] = _KDE(Pdf[0]);
  Pdf[2] = _AKDE(Pdf[0], Pdf[1]);
  for (int i = 0; i < 3; i++) {
    Pdf[i]->SetDirectory(0);
    Cdf[i] = _CDF(Pdf[i]);
    Cdf[i]->SetDirectory(0);
    Idf[i] = _IDF(Cdf[i]);
  }
  for (int bell_type = 0; bell_type < NBellType; bell_type++) {
    if ((bell_type >= 2) && (bell_type < 4))
      continue; // Cauchy has very long tails
    QQ[bell_type] = GetQQ(bell_type);
    if (QQ[bell_type]) {
      QQSpline[bell_type] = new MSpline();
      QQSpline[bell_type]->Set(QQ[bell_type]);
      QQSpline[bell_type]->CalculateCoefficients();
    }
    if ((PerformFits) && (QQ[bell_type]) && (Histo->Integral() > 0)) {
      QQFit[0][bell_type] = _FitPoly(Pdf[0], QQ[bell_type]);
      // QQFit[1][bell_type] = _FitOddPoly(Pdf[0],QQ[bell_type]); // too rigid
      // QQFit[2][bell_type] = _FitIsoPoly(QQIsoPoly[bell_type],Pdf[0],QQ[bell_type],1); // very slow
      // QQFit[2][bell_type] = _FitIsoPoly(QQIsoPoly[bell_type],Pdf[0],QQ[bell_type],2); // very slow
      // QQFit[2][bell_type] = _FitIsoPoly(QQIsoPoly[bell_type],Pdf[0],QQ[bell_type],3); // very slow
      // QQFit[3][bell_type] = _FitIsoSpline(QQIsoSpline[bell_type],Pdf[0],Cdf[0],Idf[0],Histo,QQ[bell_type]);
    }
  }
  if (DrawChecks)
    Draw();
  if (PrintChecks)
    Print();
}

void PDF1::Draw(Option_t *opt) {
  int pdf_color[3] = {kBlack, kOrange + 1, kMagenta + 1};
  int bell_type_color[5] = {kBlue + 1, kRed + 1, kGray, kGreen + 1, kAzure - 4};
  // QQ plot
  TCanvas canvas("PDF1_canvas", "PDF1::Draw", 800, 600);
  canvas.Divide(2, 2);
  canvas.cd(1);
  std::string draw("APZ");
  for (int bell_type = 0; bell_type < NBellType; bell_type++) {
    if (!QQ[bell_type])
      continue;
    QQ[bell_type]->SetMarkerStyle(20);
    QQ[bell_type]->SetMarkerSize(0.4);
    QQ[bell_type]->SetMarkerColor(bell_type_color[bell_type]);
    QQ[bell_type]->SetLineColor(bell_type_color[bell_type]);
    QQ[bell_type]->Draw(draw.c_str());
    draw.assign("PZ");
  }
  for (int bell_type = 0; bell_type < NBellType; bell_type++) {
    for (int fit_type = 0; fit_type < NFitType; fit_type++) {
      if (!QQFit[fit_type][bell_type])
        continue;
      QQFit[fit_type][bell_type]->SetLineStyle(fit_type + 1);
      QQFit[fit_type][bell_type]->SetLineColor(bell_type_color[bell_type]);
      QQFit[fit_type][bell_type]->SetMarkerColor(bell_type_color[bell_type]);
      QQFit[fit_type][bell_type]->SetLineWidth(2);
      QQFit[fit_type][bell_type]->Draw("SAME");
    }
  }
  // PDF
  TVirtualPad *pad = canvas.cd(2);
  for (int i = 0; i < 3; i++) {
    Pdf[i]->SetStats(0);
    Pdf[i]->SetLineColor(pdf_color[i]);
    Pdf[i]->SetMarkerColor(pdf_color[i]);
    Pdf[i]->SetLineWidth(2);
    Pdf[i]->Draw((i == 0) ? "H" : "H SAME");
  }
  double xmin = Pdf[0]->GetXaxis()->GetXmin();
  double xmax = Pdf[0]->GetXaxis()->GetXmax();
  TGraph *g[NFitType][NBellType] = {{0}};
  for (int bell_type = 0; bell_type < NBellType; bell_type++) {
    for (int fit_type = 0; fit_type < NFitType; fit_type++) {
      g[fit_type][bell_type] = GetGraphPDF(bell_type, fit_type, xmin, xmax);
      if (!g[fit_type][bell_type])
        continue;
      g[fit_type][bell_type]->SetLineStyle(fit_type + 1);
      g[fit_type][bell_type]->SetLineColor(bell_type_color[bell_type]);
      g[fit_type][bell_type]->SetMarkerColor(bell_type_color[bell_type]);
      g[fit_type][bell_type]->SetLineWidth(2);
      g[fit_type][bell_type]->Draw("LX");
    }
  }
  pad->SetLogy();
  // CDF
  pad = canvas.cd(3);
  for (int i = 0; i < 3; i++) {
    Cdf[i]->SetStats(0);
    Cdf[i]->SetLineColor(pdf_color[i]);
    Cdf[i]->SetMarkerColor(pdf_color[i]);
    Cdf[i]->SetLineWidth(2);
    Cdf[i]->Draw((i == 0) ? "H" : "H SAME");
  }
  TGraph *c[NFitType][NBellType] = {{0}};
  for (int bell_type = 0; bell_type < NBellType; bell_type++) {
    for (int fit_type = 0; fit_type < NFitType; fit_type++) {
      c[fit_type][bell_type] = GetGraphCDF(bell_type, fit_type, xmin, xmax);
      if (!c[fit_type][bell_type])
        continue;
      c[fit_type][bell_type]->SetLineStyle(fit_type + 1);
      c[fit_type][bell_type]->SetLineColor(bell_type_color[bell_type]);
      c[fit_type][bell_type]->SetMarkerColor(bell_type_color[bell_type]);
      c[fit_type][bell_type]->SetLineWidth(2);
      c[fit_type][bell_type]->Draw("LX");
    }
  }
  pad->SetLogy(0);
  // IDF
  pad = canvas.cd(4);
  for (int i = 0; i < 3; i++) {
    Idf[i]->SetLineColor(pdf_color[i]);
    Idf[i]->SetMarkerColor(pdf_color[i]);
    Idf[i]->SetLineWidth(2);
    Idf[i]->Draw((i == 0) ? "APZ" : "PZ");
  }
  TGraph *i[NFitType][NBellType] = {{0}};
  for (int bell_type = 0; bell_type < NBellType; bell_type++) {
    for (int fit_type = 0; fit_type < NFitType; fit_type++) {
      i[fit_type][bell_type] = GetGraphIDF(bell_type, fit_type, 0, 1);
      if (!i[fit_type][bell_type])
        continue;
      i[fit_type][bell_type]->SetLineStyle(fit_type + 1);
      i[fit_type][bell_type]->SetLineColor(bell_type_color[bell_type]);
      i[fit_type][bell_type]->SetMarkerColor(bell_type_color[bell_type]);
      i[fit_type][bell_type]->SetLineWidth(1);
      i[fit_type][bell_type]->Draw("LX");
    }
  }
  pad->SetLogy(0);
  canvas.Update();
  canvas.WaitPrimitive();
  canvas.cd(0);
  for (int bell_type = 0; bell_type < NBellType; bell_type++) {
    for (int fit_type = 0; fit_type < NFitType; fit_type++) {
      if (g[fit_type][bell_type])
        delete g[fit_type][bell_type];
      if (c[fit_type][bell_type])
        delete c[fit_type][bell_type];
      if (i[fit_type][bell_type])
        delete i[fit_type][bell_type];
    }
  }
  canvas.Clear();
}

bool PDF1::IsDone(int bell_type, int fit_type) const {
  if ((bell_type < -3) || (bell_type >= NBellType))
    return false;
  if ((fit_type < -1) || (fit_type >= NFitType))
    return false;
  if (bell_type < 0)
    return true;
  if (fit_type < 0)
    return true;
  return (QQFit[fit_type][bell_type] != 0);
}

double PDF1::GetChi2(int bell_type, int fit_type) const {
  if (!Pdf[0])
    return -1;
  if (!IsDone(bell_type, fit_type))
    return -1;
  double chi2 = 0;
  int ndof = 0;
  for (int i = 0; i <= Pdf[0]->GetNbinsX() + 1; i++) {
    double y = Pdf[0]->GetBinContent(i);
    double ey = Pdf[0]->GetBinError(i);
    double y_th = EvalPDF(Pdf[0]->GetXaxis()->GetBinCenter(i), bell_type, fit_type);
    if (ey <= 0)
      continue;
    chi2 += std::pow((y - y_th) / ey, 2);
    ndof++;
  }
  return (ndof > 0) ? chi2 / ndof : 0;
}

double PDF1::EvalQQ(double x, int bell_type, int fit_type) const {
  if (!IsDone(bell_type, fit_type))
    return 0;
  if (fit_type == kSpline)
    return (!QQSpline[bell_type]) ? 0 : QQSpline[bell_type]->Eval(x);
  else if (fit_type == kPoly)
    return (!QQFit[fit_type][bell_type]) ? 0 : _poly_extlin_fun(&x, QQFit[fit_type][bell_type]->GetParameters());
  else if (fit_type == kOddPoly)
    return (!QQFit[fit_type][bell_type]) ? 0 : _odd_poly_extlin_fun(&x, QQFit[fit_type][bell_type]->GetParameters());
  else if (fit_type == kIsoPoly)
    return (!QQIsoPoly[bell_type]) ? 0 : QQIsoPoly[bell_type]->EvalPolyLinExt(x);
  else if (fit_type == kIsoSpline)
    return (!QQIsoSpline[bell_type]) ? 0 : QQIsoSpline[bell_type]->Eval(x);
  return 0;
}

double PDF1::EvalQQ1(double x, int bell_type, int fit_type) const {
  if (!IsDone(bell_type, fit_type))
    return 0;
  double deri = 0;
  if (fit_type == kSpline)
    deri = (!QQSpline[bell_type]) ? 0 : QQSpline[bell_type]->EvalDerivative(x);
  else if (fit_type == kPoly)
    deri = (!QQFit[fit_type][bell_type]) ? 0 : _poly_extlin_der_fun(&x, QQFit[fit_type][bell_type]->GetParameters());
  else if (fit_type == kOddPoly)
    deri =
        (!QQFit[fit_type][bell_type]) ? 0 : _odd_poly_extlin_der_fun(&x, QQFit[fit_type][bell_type]->GetParameters());
  else if (fit_type == kIsoPoly)
    deri = (!QQIsoPoly[bell_type]) ? 0 : QQIsoPoly[bell_type]->EvalPoly1LinExt(x);
  else if (fit_type == kIsoSpline)
    deri = (!QQIsoSpline[bell_type]) ? 0 : QQIsoSpline[bell_type]->EvalDerivative(x);
  return (deri < 0) ? 0 : deri;
}

TGraph *PDF1::GetGraphPDF(int bell_type, int fit_type, double xmin, double xmax, int npx) const {
  if (!IsDone(bell_type, fit_type))
    return 0;
  TGraph *g = new TGraph();
  for (int i = 0; i < npx; i++) {
    double x = xmin + (0.5 + i) * (xmax - xmin) / npx;
    double y = EvalPDF(x, bell_type, fit_type);
    g->SetPoint(g->GetN(), x, y);
    if (y < 0)
      std::printf("PDF1::GetGraphPDF-W y=%f is negative for x=%f.\n", x, y);
  }
  return g;
}

TGraph *PDF1::GetGraphCDF(int bell_type, int fit_type, double xmin, double xmax, int npx) const {
  if (!IsDone(bell_type, fit_type))
    return 0;
  TGraph *g = new TGraph();
  for (int i = 0; i < npx; i++) {
    double x = xmin + (0.5 + i) * (xmax - xmin) / npx;
    double y = EvalCDF(x, bell_type, fit_type);
    g->SetPoint(g->GetN(), x, y);
  }
  return g;
}

TGraph *PDF1::GetGraphIDF(int bell_type, int fit_type, double xmin, double xmax, int npx) const {
  if (!IsDone(bell_type, fit_type))
    return 0;
  TGraph *g = new TGraph();
  for (int i = 0; i < npx; i++) {
    double x = xmin + (0.5 + i) * (xmax - xmin) / npx;
    double y = EvalIDF(x, bell_type, fit_type);
    g->SetPoint(g->GetN(), x, y);
  }
  return g;
}

double PDF1::EvalPDF(double x, int bell_type, int fit_type) const {
  if (!IsDone(bell_type, fit_type))
    return 0;
  if (bell_type < 0)
    return Pdf[bell_type + 3]->Interpolate(x);
  if (fit_type < -1)
    return 0;
  double valu = EvalQQ(x, bell_type, fit_type);
  double deri = EvalQQ1(x, bell_type, fit_type);
  if ((valu > 30) || (valu < -30))
    return 0;
  if (bell_type == kLogistic1) {
    double logi = 1 / (1 + std::exp(-valu));
    return deri * logi * (1 - logi);
  } else if (bell_type == kNormal)
    return deri * std::exp(-0.5 * valu * valu) / std::sqrt(2 * M_PI);
  else if (bell_type == kCauchy)
    return deri / (1 + valu * valu) / M_PI;
  else if (bell_type == kGumbel)
    return deri * std::exp(-(valu + std::exp(-valu)));
  else if ((bell_type == kUser) && (UserPdf)) {
    double v = deri * UserPdf->Eval(valu);
    return (v <= 0) ? 0 : v;
  }
  return 0;
}

double PDF1::EvalCDF(double x, int bell_type, int fit_type) const {
  if (!IsDone(bell_type, fit_type))
    return -1;
  if (bell_type < 0)
    return Cdf[bell_type + 3]->Interpolate(x);
  if (fit_type < -1)
    return -2;
  double valu = EvalQQ(x, bell_type, fit_type);
  if (valu > 30)
    return 1;
  if (valu < -30)
    return 0;
  if (bell_type == kLogistic1)
    return 1 / (1 + std::exp(-valu));
  else if (bell_type == kNormal)
    return 0.5 + 0.5 * TMath::Erf(valu / std::sqrt(2));
  else if (bell_type == kCauchy)
    return 0.5 + std::atan(valu) / M_PI;
  else if (bell_type == kGumbel)
    return std::exp(-std::exp(-valu));
  else if ((bell_type == kUser) && (UserCdf)) {
    double v = UserCdf->Eval(valu);
    return (v <= 0) ? 0 : v;
  }
  return 0;
}

double PDF1::EvalIDF(double x, int bell_type, int fit_type) const {
  if (!IsDone(bell_type, fit_type))
    return -1;
  if (bell_type < 0)
    return Idf[bell_type + 3]->Eval(x);
  if (fit_type < -1)
    return -2;
  double y = 0;
  double eps = 1e-6;
  if ((x >= 1 - eps) || (x <= eps))
    return -3;
  if (bell_type == kLogistic1)
    y = _logit(x);
  else if (bell_type == kNormal)
    y = _probit(x);
  else if (bell_type == kCauchy)
    y = _cauchit(x);
  else if (bell_type == kGumbel)
    y = _gumbit(x);
  else if ((bell_type == kUser) && (UserIdf))
    y = UserIdf->Eval(x);
  return QQFit[fit_type][bell_type]->GetX(y, Idf[0]->Eval(x - 0.1 * fabs(x)), Idf[0]->Eval(x + 0.1 * fabs(x)));
}

TGraphErrors *PDF1::GetQQ(int bell_type) const {
  if ((!Pdf[0]) || (!Cdf[0]))
    return 0;
  if ((bell_type < 0) || (bell_type >= NBellType))
    return 0;
  if ((bell_type == kUser) && (UserIdf == 0))
    return 0;
  TGraphErrors *qq = new TGraphErrors();
  double eps = 1e-6;
  for (int i = 0; i <= Cdf[0]->GetXaxis()->GetNbins() + 1; i++) {
    double c = Cdf[0]->GetBinContent(i);
    double ec = Cdf[0]->GetBinError(i);
    if ((c >= 1 - eps) || (c <= eps))
      continue;
    double x = Cdf[0]->GetXaxis()->GetBinUpEdge(i); // the bin value corresponds to the full integral in the bin
    double h = 0;
    double eh = 0;
    if (bell_type == kLogistic1) {
      h = _logit(c);
      eh = fabs(_logit(std::min(1 - eps, c + ec)) - _logit(std::max(eps, c - ec)));
    } else if (bell_type == kNormal) {
      h = _probit(c);
      eh = fabs(_probit(std::min(1 - eps, c + ec)) - _probit(std::max(eps, c - ec)));
    } else if (bell_type == kCauchy) {
      h = _cauchit(c);
      eh = fabs(_cauchit(std::min(1 - eps, c + ec)) - _cauchit(std::max(eps, c - ec)));
    } else if (bell_type == kGumbel) {
      h = _gumbit(c);
      eh = fabs(_gumbit(std::min(1 - eps, c + ec)) - _gumbit(std::max(eps, c - ec)));
    } else if (bell_type == kUser) {
      h = UserIdf->Eval(c);
      eh = fabs(UserIdf->Eval(std::min(1 - eps, c + ec)) - UserIdf->Eval(std::max(eps, c - ec)));
    }
    if (eh <= 0)
      continue;
    int index = qq->GetN();
    qq->SetPoint(index, x, h);
    qq->SetPointError(index, 0, eh);
  }
  if (qq->GetN() == 0) {
    delete qq;
    qq = 0;
  }
  return qq;
}

TF1 *PDF1::_FitPoly(TH1 *pdf, TGraphErrors *qq) {
  if (!qq)
    return 0;
  double xmin, xmax, dummy;
  qq->GetPoint(0, xmin, dummy);
  qq->GetPoint(qq->GetN() - 1, xmax, dummy);
  TF1 *fit_qq = 0;
  double mean = pdf->GetMean();
  double rms = pdf->GetRMS();
  for (int n = 2; n < 10; n++) {
    delete fit_qq;
    fit_qq = new TF1(Form("fit_qq_poly_%s", pdf->GetName()), &_poly_extlin_fun, xmin, xmax, n + 3);
    TF1 der_qq("der_qq", &_poly_extlin_der_fun, xmin, xmax, n + 3);
    // fit core
    fit_qq->FixParameter(0, TMath::Max(xmin, mean - 3 * rms));
    fit_qq->FixParameter(1, TMath::Min(xmax, mean + 3 * rms));
    fit_qq->FixParameter(2, n);
    qq->Fit(fit_qq, "QN", "", TMath::Max(xmin, mean - 3 * rms), TMath::Min(xmax, mean + 3 * rms));
    // extend
    qq->Fit(fit_qq, "QN");
    // leave free
    fit_qq->SetParLimits(0, xmin, TMath::Max(xmin, mean - 3 * rms));
    fit_qq->SetParLimits(1, TMath::Min(xmax, mean + 3 * rms), xmax);
    qq->Fit(fit_qq, "QN");
    for (int ipar = 0; ipar < n + 3; ipar++)
      der_qq.SetParameter(ipar, fit_qq->GetParameter(ipar));
    if ((fit_qq->GetChisquare() < Chi2Max * fit_qq->GetNDF()) && (der_qq.GetMinimum() > 0))
      break;
  }
  return fit_qq;
}

TF1 *PDF1::_FitOddPoly(TH1 *pdf, TGraphErrors *qq) {
  if (!qq)
    return 0;
  double xmin, xmax, dummy;
  qq->GetPoint(0, xmin, dummy);
  qq->GetPoint(qq->GetN() - 1, xmax, dummy);
  TF1 *fit_qq = 0;
  double mean = pdf->GetMean();
  double rms = pdf->GetRMS();
  for (int n = 2; n < 10; n++) {
    delete fit_qq;
    fit_qq = new TF1(Form("fit_qq_oddpoly_%s", pdf->GetName()), &_odd_poly_extlin_fun, xmin, xmax, n + 3);
    for (int i = 1; i < n; i++)
      fit_qq->SetParLimits(i + 3, 0, 1e+6); // positive derivative
    // fit core
    fit_qq->FixParameter(0, mean - 3 * rms);
    fit_qq->FixParameter(1, mean + 3 * rms);
    fit_qq->FixParameter(2, n);
    // extend
    qq->Fit(fit_qq, "QN");
    // leave free
    fit_qq->SetParLimits(0, xmin, mean - 3 * rms);
    fit_qq->SetParLimits(1, mean + 3 * rms, xmax);
    qq->Fit(fit_qq, "QN");
    if (fit_qq->GetChisquare() < Chi2Max * fit_qq->GetNDF())
      break;
  }
  return fit_qq;
}

TF1 *PDF1::_FitIsoPoly(IsoPoly *&isopoly, TH1 *pdf, TGraphErrors *qq, int isopoly_type) {
  if (!qq)
    return 0;
  double xmin, xmax, dummy;
  qq->GetPoint(0, xmin, dummy);
  qq->GetPoint(qq->GetN() - 1, xmax, dummy);
  isopoly = new IsoPoly();
  isopoly->fType = isopoly_type;
  TF1 *fit_qq = 0;
  double mean = pdf->GetMean();
  double rms = pdf->GetRMS();
  for (int nk = 0; nk < 7; nk++) {
    delete fit_qq;
    fit_qq = new TF1(Form("fit_qq_isopoly_%s", pdf->GetName()), isopoly, xmin, xmax, 2 * nk + 5);
    for (int i = 5; i < 2 * nk + 5; i++)
      fit_qq->SetParLimits(i, -10, 10);
    fit_qq->SetParameter(3, 0);
    fit_qq->SetParLimits(4, 0, 1e6);
    // fit core
    fit_qq->FixParameter(0, mean - 3 * rms);
    fit_qq->FixParameter(1, mean + 3 * rms);
    fit_qq->FixParameter(2, nk);
    // extend
    qq->Fit(fit_qq, "QN");
    // leave free
    fit_qq->SetParLimits(0, xmin, mean - 3 * rms);
    fit_qq->SetParLimits(1, mean + 3 * rms, xmax);
    qq->Fit(fit_qq, "QN");
    if (fit_qq->GetChisquare() < Chi2Max * fit_qq->GetNDF())
      break;
  }
  return fit_qq;
}

TF1 *PDF1::_FitIsoSpline(MSpline *&isospline, TH1 *pdf, TH1 *cdf, TGraphErrors *idf, TH1 *histo, TGraphErrors *qq) {
  if (!qq)
    return 0;
  TF1 *fit_qq = 0;
  isospline = new MSpline();
  for (int itry = 0; itry < 6; itry++) {
    double sigma = 1.7 + 0.3 * itry;
    // init
    double xmin, xmax, ymin, ymax;
    qq->GetPoint(0, xmin, ymin);
    qq->GetPoint(qq->GetN() - 1, xmax, ymax);
    ymin -= 0.1 * fabs(ymin);
    ymax += 0.1 * fabs(ymax);
    std::vector<double> xknots;
    // put points equidistant from median with number of sigma, up to a fixed amount of statistics
    xknots.push_back(idf->Eval(0.5)); // median
    for (int is = 1; is <= 10; is++) {
      double integral, eintegral;
      double nsigma = is / sigma;
      double bottom = (1 - TMath::Erf(nsigma / sqrt(2))) / 2;
      _integral(histo, integral, eintegral, histo->GetXaxis()->GetBinLowEdge(1), idf->Eval(bottom));
      if (eintegral < 0.1 * integral)
        xknots.push_back(idf->Eval(bottom));
      double top = 1 - bottom;
      _integral(histo, integral, eintegral, idf->Eval(top), histo->GetXaxis()->GetBinUpEdge(histo->GetNbinsX()));
      if (eintegral < 0.1 * integral)
        xknots.push_back(idf->Eval(top));
    }
    sort(xknots.begin(), xknots.end());
    // fit isospline
    int nknots = int(xknots.size());
    int npars = nknots + 2;
    isospline->Set(xknots);
    delete fit_qq;
    fit_qq = new TF1(Form("fit_qq_isospline_%s", pdf->GetName()), isospline, &MSpline::EvalIsoSpline1, xmin, xmax,
                     npars, "MSpline", "EvalIsoSpline1");
    for (int ipar = 0; ipar < nknots; ipar++) {
      // double x = xknots.at(ipar);
      double y = qq->Eval(xknots.at(ipar));
      if (ipar == 0) {
        fit_qq->SetParameter(ipar, TMath::Min(TMath::Max(y, ymin), ymax));
        fit_qq->SetParLimits(ipar, ymin, ymax);
      } else {
        double y_prev = qq->Eval(xknots.at(ipar - 1));
        double dy = y - y_prev;
        if (dy < 0)
          dy = 0;
        fit_qq->SetParameter(ipar, dy);
        fit_qq->SetParLimits(ipar, 0, TMath::Max(dy, ymax - ymin));
      }
    }
    for (int ipar = nknots; ipar < nknots + 2; ipar++) {
      fit_qq->SetParameter(ipar, 0);
      fit_qq->SetParLimits(ipar, 0, 1e3);
    }
    qq->Fit(fit_qq, "QNM");
    if (fit_qq->GetChisquare() < Chi2Max * fit_qq->GetNDF())
      break;
  }
  return fit_qq;
}

double PDF1::_poly_fun(double *x, double *par) {
  int n = floor(par[0] + 0.5);
  double poly = 0;
  for (int i = 0; i < n; i++)
    poly += par[i + 1] * std::pow(x[0], i);
  return poly;
}

double PDF1::_poly_der_fun(double *x, double *par) {
  int n = floor(par[0] + 0.5);
  double deri = 0;
  for (int i = 1; i < n; i++)
    deri += i * par[i + 1] * std::pow(x[0], i - 1);
  return deri;
}

double PDF1::_poly_extlin_fun(double *x, double *par) {
  double xmin = par[0];
  double xmax = par[1];
  if (x[0] <= xmin)
    return _poly_der_fun(&xmin, &par[2]) * (x[0] - xmin) + _poly_fun(&xmin, &par[2]);
  else if (x[0] >= xmax)
    return _poly_der_fun(&xmax, &par[2]) * (x[0] - xmax) + _poly_fun(&xmax, &par[2]);
  return _poly_fun(x, &par[2]);
}

double PDF1::_poly_extlin_der_fun(double *x, double *par) {
  double xmin = par[0];
  double xmax = par[1];
  if (x[0] <= xmin)
    return _poly_der_fun(&xmin, &par[2]);
  else if (x[0] >= xmax)
    return _poly_der_fun(&xmax, &par[2]);
  return _poly_der_fun(x, &par[2]);
}

double PDF1::_odd_poly_fun(double *x, double *par) {
  int n = floor(par[0] + 0.5) + 1;
  double poly = 0;
  for (int i = 0; i < n; i++) {
    int k = (i < 1) ? 0 : 2 * i - 1;
    poly += par[i + 1] * std::pow(x[0], k);
  }
  return poly;
}

double PDF1::_odd_poly_der_fun(double *x, double *par) {
  int n = floor(par[0] + 0.5) + 1;
  double deri = 0;
  for (int i = 1; i < n; i++) {
    int k = 2 * i - 1;
    deri += k * par[i + 1] * std::pow(x[0], k - 1);
  }
  return deri;
}

double PDF1::_odd_poly_extlin_fun(double *x, double *par) {
  double xmin = par[0];
  double xmax = par[1];
  if (x[0] <= xmin)
    return _odd_poly_der_fun(&xmin, &par[2]) * (x[0] - xmin) + _odd_poly_fun(&xmin, &par[2]);
  else if (x[0] >= xmax)
    return _odd_poly_der_fun(&xmax, &par[2]) * (x[0] - xmax) + _odd_poly_fun(&xmax, &par[2]);
  return _odd_poly_fun(x, &par[2]);
}

double PDF1::_odd_poly_extlin_der_fun(double *x, double *par) {
  double xmin = par[0];
  double xmax = par[1];
  if (x[0] <= xmin)
    return _odd_poly_der_fun(&xmin, &par[2]);
  else if (x[0] >= xmax)
    return _odd_poly_der_fun(&xmax, &par[2]);
  return _odd_poly_der_fun(x, &par[2]);
}

void PDF1::_integral(TH1 *histo, double &integral, double &eintegral, double left, double right) {
  integral = 0;
  eintegral = 0;
  if (!histo)
    return;
  int ibinmin = histo->GetXaxis()->FindBin(left);
  int ibinmax = histo->GetXaxis()->FindBin(right);
  for (int ibin = ibinmin; ibin <= ibinmax; ibin++) {
    double x1 = histo->GetXaxis()->GetBinLowEdge(ibin);
    double x2 = histo->GetXaxis()->GetBinUpEdge(ibin);
    if ((x1 < left) && (x2 < left))
      continue;
    if ((x1 > right) && (x2 > right))
      continue;
    double X1 = TMath::Max(x1, left);
    double X2 = TMath::Min(x2, right);
    integral += histo->GetBinContent(ibin) * (X2 - X1) / (x2 - x1);
    eintegral += pow(histo->GetBinError(ibin) * (X2 - X1) / (x2 - x1), 2);
  }
  eintegral = sqrt(eintegral);
}

TH1 *PDF1::_PDF(TH1 *histo) {
  if (!histo)
    return 0;
  TH1 *pdf = (TH1 *)histo->Clone(Form("pdf_%s", histo->GetName()));
  pdf->SetDirectory(0);
  double integral = 0;
  double eintegral = 0;
  for (int i = 0; i <= pdf->GetNbinsX() + 1; i++) {
    integral += pdf->GetBinContent(i);
    eintegral += std::pow(pdf->GetBinError(i), 2);
  }
  eintegral = std::sqrt(eintegral);
  for (int i = 0; i <= pdf->GetNbinsX() + 1; i++) {
    double a = pdf->GetBinContent(i);
    double ea = pdf->GetBinError(i);
    double width = pdf->GetXaxis()->GetBinWidth(i);
    pdf->SetBinContent(i, ((fabs(integral) > 0) && (fabs(width) > 0)) ? a / integral / width : 0);
    pdf->SetBinError(i, (fabs(width) > 0) ? _binomial_error(a, ea, integral, eintegral) / width : 0);
  }
  return pdf;
}

TH1 *PDF1::_CDF(TH1 *pdf) {
  if (!pdf)
    return 0;
  double eintegral = 0;
  for (int i = 0; i <= pdf->GetNbinsX() + 1; i++) {
    double width = pdf->GetXaxis()->GetBinWidth(i);
    eintegral += std::pow(pdf->GetBinError(i) * width, 2);
  }
  eintegral = std::sqrt(eintegral);
  TH1 *cdf = (TH1 *)pdf->Clone(Form("cdf_%s", pdf->GetName()));
  cdf->SetDirectory(0);
  cdf->Reset();
  for (int i = 0; i <= pdf->GetNbinsX() + 1; i++) {
    double width = pdf->GetXaxis()->GetBinWidth(i);
    if (i == 0) {
      cdf->SetBinContent(i, pdf->GetBinContent(i) * width);
      cdf->SetBinError(i, pdf->GetBinError(i) * width);
    } else {
      cdf->SetBinContent(i, pdf->GetBinContent(i) * width + cdf->GetBinContent(i - 1));
      cdf->SetBinError(i, std::sqrt(std::pow(pdf->GetBinError(i) * width, 2) + std::pow(cdf->GetBinError(i - 1), 2)));
    }
  }
  for (int i = 0; i <= cdf->GetNbinsX() + 1; i++)
    cdf->SetBinError(i, _binomial_error(cdf->GetBinContent(i), cdf->GetBinError(i), 1., eintegral));
  return cdf;
}

TGraphErrors *PDF1::_IDF(TH1 *cdf) {
  if (!cdf)
    return 0;
  TGraphErrors *idf = new TGraphErrors();
  for (int i = 0; i <= cdf->GetNbinsX() + 1; i++) {
    double x = cdf->GetXaxis()->GetBinCenter(i);
    double y = cdf->GetBinContent(i);
    double ey = cdf->GetBinError(i);
    idf->SetPoint(idf->GetN(), y, x);
    idf->SetPointError(idf->GetN() - 1, ey, 0);
  }
  return idf;
}

TH1 *PDF1::_KDE(TH1 *pdf) {
  if (!pdf)
    return 0;
  TH1 *kde = (TH1 *)pdf->Clone(Form("kde_%s", pdf->GetName()));
  kde->SetDirectory(0);
  kde->Reset();
  double sigma = pdf->GetRMS();
  double hh = 1.06 * sigma * std::pow(pdf->GetEntries(), -0.2); // rule-of-thumb
  if (hh <= 0)
    return kde;
  for (int i = 0; i <= pdf->GetNbinsX() + 1; i++) {
    double A = pdf->GetBinContent(i);
    double eA = pdf->GetBinError(i);
    if ((A <= 0) || (eA <= 0))
      continue;
    double mu = pdf->GetXaxis()->GetBinCenter(i);
    for (int j = 0; j <= kde->GetNbinsX() + 1; j++) {
      double t1 = (kde->GetXaxis()->GetBinLowEdge(j) - mu) / hh;
      if (t1 <= -1)
        t1 = -1;
      else if (t1 >= 1)
        t1 = 1;
      double t2 = (kde->GetXaxis()->GetBinUpEdge(j) - mu) / hh;
      if (t2 <= -1)
        t2 = -1;
      else if (t2 >= 1)
        t2 = 1;
      double kernel_integral = 0.25 * (t2 * (3 - t2 * t2) - t1 * (3 - t1 * t1)); // Epanechnikov kernel integral
      kde->SetBinContent(j, kde->GetBinContent(j) + A * kernel_integral);
      kde->SetBinError(j, std::sqrt(std::pow(kde->GetBinError(j), 2) + std::pow(eA * kernel_integral, 2)));
    }
  }
  return kde;
}

TH1 *PDF1::_AKDE(TH1 *pdf, TH1 *kde) {
  if ((!pdf) || (!kde))
    return 0;
  TH1 *akde = (TH1 *)pdf->Clone(Form("akde_%s", pdf->GetName()));
  akde->SetDirectory(0);
  akde->Reset();
  double sigma = pdf->GetRMS();
  for (int i = 0; i <= pdf->GetNbinsX() + 1; i++) {
    double A = pdf->GetBinContent(i);
    double eA = pdf->GetBinError(i);
    if ((A <= 0) || (eA <= 0))
      continue;
    double mu = pdf->GetXaxis()->GetBinCenter(i);
    if (kde->Interpolate(mu) <= 0)
      continue;
    double hstar = 1.06 * std::sqrt(sigma) * std::pow(pdf->GetEntries(), -0.2) /
                   std::sqrt(kde->Interpolate(mu)); // modified-by-me-by-hand >>> not sure
    if (hstar <= 0)
      continue;
    for (int j = 0; j <= akde->GetNbinsX() + 1; j++) {
      double t1 = (akde->GetXaxis()->GetBinLowEdge(j) - mu) / hstar;
      if (t1 <= -1)
        t1 = -1;
      else if (t1 >= 1)
        t1 = 1;
      double t2 = (akde->GetXaxis()->GetBinUpEdge(j) - mu) / hstar;
      if (t2 <= -1)
        t2 = -1;
      else if (t2 >= 1)
        t2 = 1;
      double kernel_integral = 0.25 * (t2 * (3 - t2 * t2) - t1 * (3 - t1 * t1)); // Epanechnikov kernel integral
      akde->SetBinContent(j, akde->GetBinContent(j) + A * kernel_integral);
      akde->SetBinError(j, std::sqrt(std::pow(akde->GetBinError(j), 2) + std::pow(eA * kernel_integral, 2)));
    }
  }
  return akde;
}

double PDF1::_binomial_error(double a, double ea, double c, double ec) {
  if (ec < ea)
    return 0;
  double b = c - a;
  double eb = std::sqrt(ec * ec - ea * ea);
  return (c > 0) ? std::sqrt(b * b * ea * ea + a * a * eb * eb) / (c * c) : 0;
}

double PDF1::_InterpolateQQ(double x, double f, PDF1 *a, PDF1 *b, int bell_type, int fit_type_a, int fit_type_b) {
  if ((!a) || (!b))
    return -1;
  double va = a->EvalQQ(x, bell_type, fit_type_a);
  double vb = b->EvalQQ(x, bell_type, fit_type_b);
  return f * va + (1 - f) * vb;
}

double PDF1::_InterpolateThroughQQ(double x, double f, PDF1 *a, PDF1 *b, int bell_type, int fit_type_a,
                                   int fit_type_b) {
  if ((!a) || (!b))
    return -1;
  if ((fit_type_a < 0) || (fit_type_b < 0))
    return -2; // no any valid fit
  double va = a->EvalQQ(x, bell_type, fit_type_a);
  double vb = b->EvalQQ(x, bell_type, fit_type_b);
  double da = a->EvalQQ1(x, bell_type, fit_type_a);
  double db = b->EvalQQ1(x, bell_type, fit_type_b);
  double valu = f * va + (1 - f) * vb;
  double deri = f * da + (1 - f) * db;
  if ((valu > 30) || (valu < -30))
    return 0;
  if (bell_type == kLogistic1) {
    double logi = 1 / (1 + std::exp(-valu));
    return deri * logi * (1 - logi);
  } else if (bell_type == kNormal)
    return deri * std::exp(-0.5 * valu * valu) / std::sqrt(2 * M_PI);
  else if (bell_type == kCauchy)
    return deri / (1 + valu * valu) / M_PI;
  else if (bell_type == kGumbel)
    return deri * std::exp(-(valu + std::exp(-valu)));
  // UserPdf not implemented (can be done only if UserPdf1 = UserPdf2)
  return 0;
}

double PDF1::_InterpolateThroughIDF(double x, double f, PDF1 *a, PDF1 *b, int bell_type_a, int bell_type_b,
                                    int fit_type_a, int fit_type_b) {
  if ((!a) || (!b))
    return -1;
  double eps = 0.0001;
  double ymin = 0;
  double ymax = 1;
  double dist = eps;
  double ymid = ymin + 0.5 * (ymax - ymin);
  double xa = a->EvalIDF(ymid, bell_type_a, fit_type_a);
  double xb = b->EvalIDF(ymid, bell_type_b, fit_type_b);
  int maxt = 0;
  while ((dist >= eps) && (maxt < 100)) {
    double xmid = f * xa + (1 - f) * xb;
    if (x > xmid)
      ymin = ymid;
    else
      ymax = ymid;
    dist = fabs(x - xmid);
    ymid = ymin + 0.5 * (ymax - ymin);
    if ((ymid < eps) || (ymid > 1 - eps))
      break;
    xa = a->EvalIDF(ymid, bell_type_a, fit_type_a);
    xb = b->EvalIDF(ymid, bell_type_b, fit_type_b);
    xmid = f * xa + (1 - f) * xb;
    if (maxt == 99)
      std::printf("PDF1::_InterpolateThroughIDF max number of iterations reached x=%7.2e f=%7.2e bell_type_a=%2d "
                  "bell_type_b=%2d fit_type_a=%2d fit_type_b=%2d\n",
                  x, f, bell_type_a, bell_type_b, fit_type_a, fit_type_b);
    maxt++;
  }
  double fa = a->EvalPDF(xa, bell_type_a, fit_type_a);
  double fb = b->EvalPDF(xb, bell_type_b, fit_type_b);
  double fmid = fa * fb / (f * fb + (1 - f) * fa);
  return ((f * fb + (1 - f) * fa) > 0) ? fmid : 0;
}

double PDF1::_InterpolateFast(double x, double f, PDF1 *a, PDF1 *b, int bell_type_a, int bell_type_b, int fit_type_a,
                              int fit_type_b) {
  if ((!a) || (!b))
    return -1;
  double ya = a->EvalPDF(x, bell_type_a, fit_type_a);
  double yb = b->EvalPDF(x, bell_type_b, fit_type_b);
  return f * ya + (1 - f) * yb;
}
