#include "RichOccupancyAnalysis.h"

RichOccupancyAnalysis::RichOccupancyAnalysis() {
  pOccupancy = 0;
}

RichOccupancyAnalysis::~RichOccupancyAnalysis() {
  delete pOccupancy;
  pOccupancy = 0; 
}

void RichOccupancyAnalysis::UBegin() {
  std::cout << "___________________________ RichOccupancyAnalysis::UBegin() ___________________________ " << std::endl;
  if (!pOccupancy) pOccupancy = new RichOccupancy();
  pOccupancy->SetName("RichOccupancy");
}

void RichOccupancyAnalysis::UProcessFill() {
  int index = 0;
  // particle
  pOccupancy->FillCounters(index++);
  if (!pParticle(0)) return; 
  pOccupancy->FillCounters(index++);
  // beta
  if (!pParticle(0)->pBetaH()) return; 
  pOccupancy->FillCounters(index++);
  BetaHR* pBetaH = pParticle(0)->pBetaH();
  if (pBetaH->GetUseHit()!=4) return; 
  pOccupancy->FillCounters(index++);
  if (pBetaH->GetNormChi2C()>10) return; 
  pOccupancy->FillCounters(index++);
  if (pBetaH->GetNormChi2T()>10) return; 
  pOccupancy->FillCounters(index++);
  int   beta_q_nhit;
  float beta_q_err; 
  float beta_q = pBetaH->GetQ(beta_q_nhit,beta_q_err,2,TofClusterHR::DefaultQOpt,1111);
  if (beta_q_nhit<3) return; 
  pOccupancy->FillCounters(index++);
  float beta_q1 = pBetaH->GetQ(beta_q_nhit,beta_q_err,2,TofClusterHR::DefaultQOpt,11);
  float beta_q2 = pBetaH->GetQ(beta_q_nhit,beta_q_err,2,TofClusterHR::DefaultQOpt,1100);
  if ( (beta_q1<0.7)||(beta_q1>1.3) ) return; // drop high-Z and multiparticle
  if ( (beta_q2<0.7)||(beta_q2>1.3) ) return; // drop high-Z and multiparticle
  pOccupancy->FillCounters(index++);
  if (beta_q_err>0.1*beta_q) return; 
  pOccupancy->FillCounters(index++);
  float beta_beta = pBetaH->GetBeta();
  // track
  if (!pParticle(0)->pTrTrack()) return; 
  pOccupancy->FillCounters(index++);
  TrTrackR* pTrTrack = pParticle(0)->pTrTrack();
  if (pTrTrack->GetNormChisqY()>10) return; 
  pOccupancy->FillCounters(index++);
  if (pTrTrack->GetRigidity()<0) return; // drop electrons  
  pOccupancy->FillCounters(index++);
  float trk_q = pTrTrack->GetInnerQ(beta_beta);
  if ( (trk_q<0.7)||(trk_q>1.3) ) return; // drop high-Z (too many e- sec.)
  pOccupancy->FillCounters(index++);   
  // under threshold logic
  AMSPoint point;
  AMSDir dir;
  pTrTrack->Interpolate(-71.87,point,dir,0);
  float rad_x = point[0];
  float rad_y = point[1];
  float cut_aerogel_external_border = 3500.;   // Aerogel external border (r**2)
  float cut_aerogel_naf_border[2] = {17.,19.}; // NaF/Aerogel border
  int   rad_index = -1; 
  if      (max(fabs(rad_x),fabs(rad_y))<cut_aerogel_naf_border[0]) rad_index = 0;
  else if (max(fabs(rad_x),fabs(rad_y))>cut_aerogel_naf_border[1]) rad_index = 1;
  if ((rad_x*rad_x+rad_y*rad_y)>cut_aerogel_external_border) rad_index = 2;
  double beta_thr[2] = {0.75-3*0.04,0.95-3*0.04};  
  if (rad_index==-1) return; 
  pOccupancy->FillCounters(index++);
  bool outside_radiator = (rad_index==0); 
  bool under_threshold = ( (rad_index>0)&&(beta_beta>beta_thr[rad_index-1]) ); 
  // integrate on PMT-by-PMT for alternative charge center finding
  std::map<int,pmt_info> pmt_info_map;
  for (int ihit=0; ihit<nRichHit(); ihit++) {
    RichHitR* hit = pRichHit(ihit);
    int pmt = int(hit->Channel/16);
    float np = hit->Npe;
    pmt_info_map[pmt].nhit++;
    pmt_info_map[pmt].np += np;
  }
  bool no_ring = ((under_threshold)||(outside_radiator))&&(nRichRing()==0); 
  if (no_ring) pOccupancy->FillCounters(index++); 
  // loop on hits 
  for (int ihit=0; ihit<nRichHit(); ihit++) {
    RichHitR* hit = pRichHit(ihit);
    if (!hit) continue;
    if (hit->IsCrossed()) continue; // from Carlos
    int pmt = int(hit->Channel/16);
    if (pmt_info_map.find(pmt)!=pmt_info_map.end()) if (pmt_info_map[pmt].np>5.) continue; // alternative charge center finding
    bool used = ((hit->Status&0x3ff)!=0);
    if ( (no_ring)&&(!used) ) pOccupancy->FillChannel(0,hit->Channel,hit->Npe);
    if (used) pOccupancy->FillChannel(1,hit->Channel,hit->Npe);
    else pOccupancy->FillChannel(2,hit->Channel,hit->Npe);   
  }
}

void RichOccupancyAnalysis::UTerminate() {
  // perform occupancy calculation 
  pOccupancy->CreateOccupancyTable();
  std::cout << "___________________________ RichOccupancyAnalysis::Terminate() ___________________________ " << std::endl;
}

