//Author P. Zuccon -- MIT
#include "VPSArchive.h"

TH1D* VPSArchive::GetHisto(int flag) {

  int nbin = cat.size();
  char hname[200];
  if (flag) sprintf(hname, "PreScaleRep");
  else  sprintf(hname, "PreScaleRepS");
  TH1D* hout = new TH1D(hname, "PreScaler Report; categories; events", nbin * 2, 0, nbin * 2);
  for (int ii = 1; ii <= nbin * 2; ii += 2) {
    std::string aa = cat[ii / 2].catid.to_string();
    if (ii % 2 == 0) aa += ".Sel";
    hout->GetXaxis()->SetBinLabel(ii, aa.c_str());
    if (flag)hout->SetBinContent(ii, cat[ii / 2].ScalerGot + cat[ii / 2].ScalerLost);
    else    hout->SetBinContent(ii + 1, cat[ii / 2].ScalerGot);
  }
  return hout;

}


Long64_t VPSArchive::Merge(TCollection* li) {
  if (!li) return 0;
  if (li->IsEmpty()) return 0;
  TIter next(li);
  while (VPSArchive* other = (VPSArchive*)next()) {
    for (int ii = 0; ii < (int)cat.size(); ii++) {
      if (
        (other->cat[ii].catid    == cat[ii].catid) &&
        (other->cat[ii].catmask  == cat[ii].catmask) &&
        (other->cat[ii].GetPrf() == cat[ii].GetPrf())
      )
        cat[ii] += other->cat[ii];
    }
  }

  return 1;
}


void VPSArchive::DrawHisto() {

  TH1D* h1 = GetHisto(0);
  TH1D* h2 = GetHisto(1);
  h1->SetLineColor(2);
  h1->Draw();
  h2->Draw("same");
  return;
}
