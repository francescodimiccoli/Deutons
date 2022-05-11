#include "Likelihood.h"

ClassImp(Likelihood);

int Likelihood::Debug = 0;
double Likelihood::MinProb = DBL_MIN;

Likelihood::Likelihood(PDF2DB *pdf2db) {
  Pdf2DB = pdf2db;
  VarList.clear();
}

void Likelihood::Clear(Option_t *opt) {
  Pdf2DB = 0; // no ownership
  for (auto it : VarList)
    delete it;
  VarList.clear();
}

void Likelihood::Print(int verbose) const {
  for (auto it : VarList)
    it->Print(verbose);
}

void Likelihood::AddVariable(std::string name, Float_t *x, Int_t *y, std::string pdf_name, std::string rulex,
                             std::string ruley, int interp_type, int bell_type, int fit_type) {
  if (!Pdf2DB) {
    std::printf("Likelihood::AddVariable, no PDF2DB available. Variable not registered.\n");
    return;
  }
  if (Pdf2DB->Pdf2Map.find(pdf_name) == Pdf2DB->Pdf2Map.end()) {
    std::printf("Likelihood::AddVariable, no PDF2 named %s available in PDF2DB. Variable not registered.\n",
                pdf_name.c_str());
    return;
  }
  PDF2 *pdf2 = (PDF2 *)Pdf2DB->Pdf2Map[pdf_name];
  VarList.push_back(new LikelihoodVar(name, x, y, pdf2, rulex, ruley, interp_type, bell_type, fit_type));
}

void Likelihood::AddVariable(std::string name, Float_t *x, Float_t *y, std::string pdf_name, std::string rulex,
                             std::string ruley, int interp_type, int bell_type, int fit_type) {
  if (!Pdf2DB) {
    std::printf("Likelihood::AddVariable, no PDF2DB available. Variable not registered.\n");
    return;
  }
  if (Pdf2DB->Pdf2Map.find(pdf_name) == Pdf2DB->Pdf2Map.end()) {
    std::printf("Likelihood::AddVariable, no PDF2 named %s available in PDF2DB. Variable not registered.\n",
                pdf_name.c_str());
    return;
  }
  PDF2 *pdf2 = (PDF2 *)Pdf2DB->Pdf2Map[pdf_name];
  VarList.push_back(new LikelihoodVar(name, x, y, pdf2, rulex, ruley, interp_type, bell_type, fit_type));
}

double Likelihood::EvalLog(int pid) {
  double logval = 0;
  for (auto it : VarList) {
    double this_val = it->Eval(pid);
    double this_logval = (this_val < Likelihood::MinProb)
                             ? -3000
                             : log10(this_val); // enhance what is too low, to a very low value ... is this ok?
    if (Debug > 1) {
      if (this_val < Likelihood::MinProb)
        printf("Likelihood::EvalLog-Debug %s has pdf(%f,%f)=%f, use 100*log10(p_min)=%f.\n", it->GetName(),
               it->GetX(), it->GetY(), this_val, this_logval);
      else
        printf("Likelihood::EvalLog-Debug %s has pdf(%f,%f)=%f.\n", it->GetName(), it->GetX(), it->GetY(),
               this_val);
    }
    logval += this_logval;
  }
  if (Debug > 1)
    printf("Likelihood::EvalLog-Debug log10(L)=%f.\n", logval);
  return logval;
}

double Likelihood::Eval(Float_t x, int pid) {
  double val = 1;
  for (auto it : VarList) {
    it->SetX(&x);
    val *= it->Eval(pid);
  }
  return val;
}

double Likelihood::Eval(Float_t &x, double xmin, double xmax, int pid) {
  TF1 Maximizer("maximizer", this, xmin, xmax, 1, "Likelihood");
  Maximizer.SetParameter(0, pid);
  Maximizer.SetNpx(10);
  x = Maximizer.GetMaximumX();
  /*
    sostituire con minimizzazione agricola
  */
  return Eval(x, pid);
}
