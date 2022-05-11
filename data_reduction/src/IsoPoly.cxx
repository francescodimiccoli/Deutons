#include "IsoPoly.h"

ClassImp(IsoPoly);

IsoPoly::IsoPoly() {
  fNK = -1;
  fK = 0;
  fType = 0;
  fNP = -1;
  fP = 0;
  fXMin = 0;
  fXMax = 0;
  fY0 = 0;
  fAlpha = 0;
}

IsoPoly::IsoPoly(const IsoPoly &orig) {
  fNK = orig.fNK;
  fK = new double[fNK];
  for (int ik = 0; ik < fNK; ik++)
    fK[ik] = orig.fK[ik];
  fType = orig.fType;
  fNP = orig.fNP;
  fP = new double[fNP];
  for (int ip = 0; ip < fNP; ip++)
    fP[ip] = orig.fP[ip];
  fY0 = orig.fY0;
  fAlpha = orig.fAlpha;
  fXMin = orig.fXMin;
  fXMax = orig.fXMax;
}

void IsoPoly::Clear(Option_t *opt) {
  fNK = -1;
  delete[] fK;
  fK = 0;
  fNP = -1;
  delete[] fP;
  fP = 0;
  fXMin = 0;
  fXMax = 0;
  fY0 = 0;
  fAlpha = 0;
}

void IsoPoly::Print(Option_t *opt) const {
  std::printf("fType=%1d fXMin=%10.3e fXMax=%10.3e fY0=%10.3e fAlpha=%10.3e\n",fType,fXMin,fXMax,fY0,fAlpha);
  for (int ik=0; ik<fNK/2; ik++) std::printf("%2d/%2d pars={%10.3e,%10.3e}\n",ik,fNK/2,fK[2*ik],fK[2*ik+1]);
  for (int ip=0; ip<fNP; ip++) std::printf("%2d/%2d par=%10.3e\n",ip,fNP,fP[ip]);
}

void IsoPoly::Calculate(int nk, double *k) {
  if ((fType < 1) || (fType > 3)) {
    std::printf("IsoPoly::Calculate-W undefined form (fType=%1d), no polynomial constructed.\n", fType);
    return;
  }
  if ((nk == fNK / 2) && (fNK != -1)) {
    bool redo = false;
    for (int ik = 0; ik < 2 * nk; ik++) {
      if (k[ik] != fK[ik]) {
        redo = true;
        break;
      }
    }
    if (!redo)
      return;
  }
  Clear();
  if (nk == 0) {
    fNK = 0;
    fNP = 1;
    fP = new double[1];
    fP[0] = 1;
    return;
  }
  fNK = 2 * nk;
  fK = new double[fNK];
  for (int ik = 0; ik < fNK; ik++)
    fK[ik] = k[ik];
  fNP = fNK + 1;
  fP = new double[fNP];
  double *a = new double[fNP];
  int na = 1;
  a[0] = 1;
  for (int ik = 0; ik < fNK / 2; ik++) {
    int nb = 3;
    double b[3] = {0};
    if (fType == 1) {
      b[0] = std::pow(fK[ik * 2], 2) + pow(fK[1 + ik * 2], 2);
      b[1] = 2 * fK[ik * 2];
      b[2] = 1;
    } else if (fType == 2) {
      b[0] = 1;
      b[1] = 2 * fK[ik * 2];
      b[2] = std::pow(fK[ik * 2], 2) + std::pow(fK[1 + ik * 2], 2);
    } else if (fType == 3) {
      b[0] = std::pow(fK[ik * 2], 2);
      b[1] = 2 * fK[ik * 2];
      b[2] = 1 + std::pow(fK[1 + ik * 2], 2);
    }
    _poly_product(na, a, nb, b, fNP, fP);
    na = fNP;
    for (int ip = 0; ip < fNP; ip++)
      a[ip] = fP[ip];
  }
  delete[] a;
}

double IsoPoly::EvalPoly(double x) const {
  double poly = 0;
  for (int i = 0; i < fNP; i++)
    poly += fP[i] * std::pow(x, i + 1) / (i + 1);
  return fY0 + fAlpha * poly;
}

double IsoPoly::EvalPoly1(double x) const {
  double deri = 0;
  for (int i = 0; i < fNP; i++)
    deri += fP[i] * std::pow(x, i);
  return fAlpha * deri;
}

double IsoPoly::EvalPolyLinExt(double x) const {
  if (x <= fXMin)
    return EvalPoly1(fXMin) * (x - fXMin) + EvalPoly(fXMin);
  else if (x >= fXMax)
    return EvalPoly1(fXMax) * (x - fXMax) + EvalPoly(fXMax);
  return EvalPoly(x);
}

double IsoPoly::EvalPoly1LinExt(double x) const {
  if (x <= fXMin)
    return EvalPoly1(fXMin);
  else if (x >= fXMax)
    return EvalPoly1(fXMax);
  return EvalPoly1(x);
}

double IsoPoly::Eval(double x) const { return (fXMin >= fXMax) ? EvalPoly(x) : EvalPolyLinExt(x); }

double IsoPoly::operator()(double *x, double *par) {
  fXMin = par[0];
  fXMax = par[1];
  int nk = std::floor(par[2] + 0.5);
  fY0 = par[3];
  fAlpha = par[4];
  Calculate(nk, &par[5]);
  return (fXMin >= fXMax) ? EvalPoly(x[0]) : EvalPolyLinExt(x[0]);
}

void IsoPoly::_poly_product(int na, double *a, int nb, double *b, int &nc, double *&c) {
  nc = na + nb - 1;
  for (int ic = 0; ic < nc; ic++)
    c[ic] = 0;
  for (int ia = 0; ia < na; ia++) {
    for (int ib = 0; ib < nb; ib++) {
      int ic = ia + ib;
      c[ic] += a[ia] * b[ib];
    }
  }
}
