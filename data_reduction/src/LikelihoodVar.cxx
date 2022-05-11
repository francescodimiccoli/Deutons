#include "LikelihoodVar.h"

ClassImp(LikelihoodVar);

void LikelihoodVar::Clear(Option_t *opt) {
  FloatX = 0;
  FloatY = 0;
  IntY = 0;
  RuleX.clear();
  RuleY.clear();
  Pdf2 = 0; // no ownership
  InterpType = 0;
  BellType = 0;
  FitType = 0;
}

void LikelihoodVar::Set(std::string name, Float_t *x, Int_t *y, PDF2 *pdf2, std::string rulex, std::string ruley,
                        int interp_type, int bell_type, int fit_type) {
  Clear();
  SetName(name.c_str());
  FloatX = x;
  FloatY = 0;
  IntY = y;
  RuleX = rulex;
  RuleY = ruley;
  Pdf2 = pdf2;
  InterpType = interp_type;
  BellType = bell_type;
  FitType = fit_type;
}

void LikelihoodVar::Set(std::string name, Float_t *x, Float_t *y, PDF2 *pdf2, std::string rulex, std::string ruley,
                        int interp_type, int bell_type, int fit_type) {
  Clear();
  SetName(name.c_str());
  FloatX = x;
  FloatY = y;
  IntY = 0;
  RuleX = rulex;
  RuleY = ruley;
  Pdf2 = pdf2;
  InterpType = interp_type;
  BellType = bell_type;
  FitType = fit_type;
}

void LikelihoodVar::Print(int verbose) const {
  std::printf("LikelihoodVar::Print %s = %f, %s = %f\n", RuleX.c_str(), GetX(), RuleY.c_str(), GetY());
  if (verbose > 0)
    Pdf2->Print();
}

double LikelihoodVar::EvalX(int pid) {
  double x = GetX();
  if ((RuleX.compare("") == 0) || (pid < 0))
    return x;
  return x;
}

double LikelihoodVar::EvalY(int pid) {
  double y = GetY();
  if ((RuleY.compare("") == 0) || (pid < 0))
    return y;
  if (RuleY.compare("dinvr_vs_invr") == 0) {
    double kn = std::pow(10, GetX());
    double r_rec = y;
    double r_gen = std::sqrt(std::pow(kn * __A(pid) + __M(pid), 2) - std::pow(__M(pid), 2)) / __Z(pid);
    y = r_gen / r_rec - 1;
  } else if (RuleY.compare("dinvb_vs_invb") == 0) {
    double kn = std::pow(10, GetX());
    double gamma = 1 + kn * __A(pid) / __M(pid);
    double b_rec = y;
    double b_gen = std::sqrt(1 - 1 / gamma / gamma);
    y = b_gen / b_rec - 1;
  } else if (RuleY.compare("db_vs_b") == 0) {
    double kn = std::pow(10, GetX());
    double gamma = 1 + kn * __A(pid) / __M(pid);
    double b_rec = y;
    double b_gen = std::sqrt(1 - 1 / gamma / gamma);
    y = 1 - b_rec / b_gen;
  }
  return y;
}

double LikelihoodVar::Eval(int pid) {
  return (Pdf2) ? Pdf2->Eval(EvalX(pid), EvalY(pid), InterpType, BellType, FitType) : 0;
}

double LikelihoodVar::__Z(int pid) {
  if ((pid == 2) || (pid == 14))
    return 1;
  else if ((pid == 3) || (pid == 15))
    return -1;
  else if (pid == 45)
    return 1;
  else if (pid == 46)
    return 1;
  else if (pid == 47)
    return 2;
  return 0;
}

double LikelihoodVar::__A(int pid) {
  if ((pid == 2) || (pid == 3))
    return 1;
  else if ((pid == 14) || (pid == 15))
    return 1;
  else if (pid == 45)
    return 2;
  else if (pid == 46)
    return 3;
  else if (pid == 47)
    return 4;
  return 0;
}

double LikelihoodVar::__M(int pid) {
  if ((pid == 2) || (pid == 3))
    return 0.0005109989461;
  else if ((pid == 14) || (pid == 15))
    return 0.938272081;
  else if (pid == 45)
    return 1.875612928;
  else if (pid == 46)
    return 2.808921112;
  else if (pid == 47)
    return 3.727379378;
  return 0;
}
