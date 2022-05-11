#include "RichOccupancy.h"

RichOccupancy::RichOccupancy() {
  counters = new TH1D("counters","; Index",20,-0.5,19.5);
  counters->SetDirectory(0);
  for (int i=0; i<3; i++) {
    nhit_vs_chn[i] = new TH1D(Form("nhit_vs_chan_%1d",i),"; Channel; N_{hit}",680*16,-0.5,680*16-0.5);
    nhit_vs_chn[i]->Sumw2();
    nhit_vs_chn[i]->SetDirectory(0);
    npho_vs_chn[i] = new TH1D(Form("npho_vs_chan_%1d",i),"; Channel; N_{p}",680*16,-0.5,680*16-0.5);
    npho_vs_chn[i]->Sumw2();
    npho_vs_chn[i]->SetDirectory(0);
  }
}

RichOccupancy::~RichOccupancy() {
  delete counters; counters = 0;
  for (int i=0; i<3; i++) {
    delete nhit_vs_chn[i]; nhit_vs_chn[i] = 0;
    delete npho_vs_chn[i]; npho_vs_chn[i] = 0;
  }
}

int RichOccupancy::EvalThreshold(TH1D* nhit_vs_chn, int N_evt, double nsigma) {
  int N_chn = 680*16;
  TH1D* hit = new TH1D("hit","",200,-0.5,199.5);
  for (int i=0; i<N_chn; i++) hit->Fill(nhit_vs_chn->GetBinContent(i));
  int k = std::ceil(hit->GetMean()+nsigma*hit->GetRMS());
  // bad list
  int nbad = 0;
  for (int i=0; i<N_chn; i++) {
    if (nhit_vs_chn->GetBinContent(i)<=k) continue;
    nbad++;
  }
  std::cout << "RichOccupancy::EvalThreshold-events:    " << N_evt << std::endl;
  std::cout << "RichOccupancy::EvalThreshold-entries:   " << hit->GetEntries() << " (" << N_chn << ")" << std::endl;
  std::cout << "RichOccupancy::EvalThreshold-threshold: " << k << " (mean:" << hit->GetMean() << ",rms:" << hit->GetRMS() << ")" << std::endl;
  std::cout << "RichOccupancy::EvalThreshold-bad:       " << nbad << std::endl;
  delete hit;
  return k;
}

void RichOccupancy::CreateOccupancyTable(double nsigma) {
  int N_chn = 680*16;
  int N_evt = counters->GetBinContent(15);
  int k = EvalThreshold(nhit_vs_chn[0],N_evt,nsigma);
  // bad list
  for (int i=0; i<N_chn; i++) {
    if (nhit_vs_chn[0]->GetBinContent(i)<=k) continue;
    AddBadChannel(i);
  }
  std::cout << "RichOccupancy::CreateOccupancyTable-bad channels: " << (int)bad_chn.size() << std::endl;
}
