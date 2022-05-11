#include "Analysis.h"
#include "HitMCTruth.h"
#include "TrdSimpleVertex.h"

#include "richrecOff.h"
#include "richidOff.h"
#include "richradidOff.h"
#include "RichCharge.h"
#include "TEnv.h"
#include <iostream>
#include "timeid.h"
#include "TkCoo.h"

#include <iostream>
#include <map>
#include <set>
#include <algorithm>
#include <functional>

float Analysis::m_p = TrFit::Mproton;
float Analysis::m_d = 1.875613;
float Analysis::tracker_layers_z[9] = {158.920,53.060,29.228,25.212,1.698,-2.318,-25.212,-29.228,-135.882};
float Analysis::rich_radiator_z = -71.87;
float Analysis::rich_pmt_plane_z = -121.89;
float Analysis::top_z = 160;

static bool disable = false; // al things that give problems

void Analysis::FillRTI() {
  AMSSetupR::RTI rti;
  GetRTI(rti);
  rtiInfo.isinsaa  = rti.IsInSAA();
  rtiInfo.run      = rti.run;
  rtiInfo.evno     = rti.evno;
  rtiInfo.evnol    = rti.evnol;
  rtiInfo.lf       = rti.lf;
  static int index_cf[7] = {0,1,12,22,32,3,4};
  for (int icf=0; icf<7; icf++) rtiInfo.ret_cf[icf] = rti.getcutoff(rtiInfo.cf[icf],index_cf[icf]);
  // a bug fix for IGRF-12, IGRF-12+, and IGRF-12++
  for (int icf=2; icf<5; icf++) {
    for (int i=0; i<4; i++) {
      for (int j=0; j<2; j++) {
        if (fabs(rtiInfo.cf[icf][i][j])<=0) rtiInfo.cf[icf][i][j] = rti.cfi[i][j];
      }
    }
  }
  rtiInfo.mphe     = rti.mphe;
  rtiInfo.theta    = rti.theta;
  rtiInfo.phi      = rti.phi;
  rtiInfo.r        = rti.r;
  rtiInfo.zenith   = rti.zenith;
  rtiInfo.glat     = rti.glat;
  rtiInfo.glong    = rti.glong;
  rtiInfo.nev      = rti.nev;
  rtiInfo.nerr     = rti.nerr;
  rtiInfo.ntrig    = rti.ntrig;
  rtiInfo.nhwerr   = rti.nhwerr;
  rtiInfo.npart    = rti.npart;
  for (int i=0; i<2; i++) {
    for (int j=0; j<2; j++) rtiInfo.nl1l9[i][j] = rti.nl1l9[i][j];
    for (int j=0; j<3; j++) rtiInfo.dl1l9[i][j] = rti.dl1l9[i][j];
  }
  rtiInfo.mtrdh    = rti.mtrdh;
  rtiInfo.good     = rti.good;
  rtiInfo.utime    = rti.utime;
  for (int i=0; i<2; i++) {
    rtiInfo.usec[i]    = rti.usec[i];
    rtiInfo.utctime[i] = rti.utctime[i];
  }
  rtiInfo.betasun = rti.getbetasun();
}

void Analysis::ClearNtp() {
  // save
  float weight     = ntpHeader.pres_weight;
  float trck_rig   = ntpHeader.pres_trck_rig;
  float tofh_beta  = ntpHeader.pres_tofh_beta;
  float rich_beta  = ntpHeader.pres_rich_beta;
  float trck_inn_q = ntpHeader.pres_trck_inn_q;
  int   patt       = ntpHeader.pres_patt;
  // clear
  memset(&ntpHeader,0,sizeof(NtpHeader));
  memset(&ntpSHeader,0,sizeof(NtpSHeader));
  memset(&ntpMCHeader,0,sizeof(NtpMCHeader));
  memset(&ntpTrd,0,sizeof(NtpTrd));
  memset(&ntpTof,0,sizeof(NtpTof));
  memset(&ntpTracker,0,sizeof(NtpTracker));
  memset(&ntpRich,0,sizeof(NtpRich));
  memset(&ntpEcal,0,sizeof(NtpEcal));
  memset(&ntpAnti,0,sizeof(NtpAnti));
  memset(&ntpStandAlone,0,sizeof(NtpStandAlone));
  // set
  ntpHeader.pres_weight = weight;
  ntpHeader.pres_trck_rig = trck_rig;
  ntpHeader.pres_tofh_beta = tofh_beta;
  ntpHeader.pres_rich_beta = rich_beta;
  ntpHeader.pres_trck_inn_q = trck_inn_q;
  ntpHeader.pres_patt = patt;
}


void Analysis::ClearCompact() {
  memset(&ntpSHeader,0,sizeof(NtpSHeader));
  memset(&ntpCompact,0,sizeof(NtpCompact));
}

void Analysis::FillNtpSHeader(){
  ntpSHeader.run = fHeader.Run;
  ntpSHeader.event = fHeader.Event;
  ntpSHeader.herror = fHeader.Error;
  AMSSetupR::RTI rti;
  GetRTI(rti);
  ntpSHeader.utime = rti.utime;
}

void Analysis::FillNtpHeader(){
  // fill
  UTCTime(ntpHeader.utc_time,ntpHeader.utc_time_error);
  ntpHeader.nmceventg     = nMCEventg();
  ntpHeader.nparticle     = nParticle();
  ntpHeader.nanti         = nAntiCluster();
  ntpHeader.ntof          = nTofCluster();
  ntpHeader.ntofh         = nTofClusterH();
  ntpHeader.ntrtrack      = nTrTrack();
  ntpHeader.ntrrechit     = nTrRecHit();
  ntpHeader.ntrdtrack     = nTrdTrack();
  ntpHeader.ntrdhtrack    = nTrdHTrack();
  for (int iseg=0; iseg<nTrdSegment(); iseg++)
    ntpHeader.ntrdsegment[pTrdSegment(iseg)->Orientation]++;
  ntpHeader.ntrdcluster   = nTrdCluster();
  ntpHeader.nrich         = nRichRing();
  ntpHeader.nrichhit      = nRichHit();
  ntpHeader.nrichhit_crossed = 0;
  for (int ihit=0; ihit<nRichHit(); ihit++) {
    RichHitR* hit = pRichHit(ihit);
    if (!hit) continue;
    if (hit->IsCrossed()) ntpHeader.nrichhit_crossed++;
  }
  ntpHeader.necal         = nEcalShower();
  // DAQ
  if (nDaqEvent()) {
    DaqEventR* daq = pDaqEvent(0);
    bool jinj_error     = daq->HasHWError();   // error bits (0x7F00)c on the four JINJStatus word
    bool assembly_error = daq->L3EventError(); // whatever bad REPLY at JINJ, JINF or LVL1 (ROOM, Sync, 0-replies, ...)
    bool sync_error     = daq->L3RunError();   // desync error at whatever level
    bool proc_error     = daq->L3ProcError();  // decoding procedure error
    bool glob_error     = false;               // whatever board with bad REPLY
    for (int inode=0; inode<7; inode++) if (daq->L3NodeError(inode)>0) glob_error = true;
    ntpHeader.error = 0;
    if (jinj_error)     ntpHeader.error |= 0x1;
    if (assembly_error) ntpHeader.error |= 0x2;
    if (sync_error)     ntpHeader.error |= 0x4;
    if (proc_error)     ntpHeader.error |= 0x8;
    if (glob_error)     ntpHeader.error |= 0x10;
    if (Status(30))     ntpHeader.error |= 0x20;
    ntpHeader.size = daq->Length;
  }
  ntpHeader.livetime = LiveTime();
  // trigger
  if (nLevel1()) {
    Level1R* lvl1 = pLevel1(0); // anyone
    ntpHeader.antipatt = lvl1->AntiPatt;
    ntpHeader.sublvl1  = lvl1->PhysBPatt;
    ntpHeader.trigpatt = lvl1->JMembPatt;
    if (nMCEventg()>0) lvl1->RebuildTrigPatt(ntpHeader.trigpatt,ntpHeader.sublvl1,ntpHeader.antipatt);
  }
  // orbit
  ntpHeader.rads      = fHeader.RadS;
  ntpHeader.thetas    = fHeader.ThetaS;
  ntpHeader.phis      = fHeader.PhiS;
  ntpHeader.yaw       = fHeader.Yaw;
  ntpHeader.pitch     = fHeader.Pitch;
  ntpHeader.roll      = fHeader.Roll;
  ntpHeader.velocitys = fHeader.VelocityS;
  ntpHeader.velthetas = fHeader.VelTheta;
  ntpHeader.velphis   = fHeader.VelPhi;
  ntpHeader.thetam    = fHeader.ThetaM;
  ntpHeader.phim      = fHeader.PhiM;

  // some global estimators
  ntpHeader.dedx_beta[0] = TrMass::GetBeta(this,1,pTrTrack);
  ntpHeader.dedx_beta[1] = TrMass::GetBeta(this,10,pTrTrack);
  ntpHeader.dedx_beta[2] = TrMass::GetBeta(this,100,pTrTrack);
  ntpHeader.dedx_beta[3] = TrMass::GetBeta(this,111,pTrTrack);
  ntpHeader.ntrdseg_vertex = TrMass::GetNtrdSegTrk(this);
  if (pTrTrack) {
    ntpHeader.mass_quality = TrMass::GetMQL(this,pTrTrack);
    ntpHeader.bl2 = TrMass::GetBL2(pTrTrack,1);
  }
}

void Analysis::FillNtpMCHeader() {
  if (nMCEventg()==0) return;

  // primary
  MCEventgR* mcpart = GetPrimaryMC();
  ntpMCHeader.pid = mcpart->Particle;
  ntpMCHeader.mass = mcpart->Mass;
  ntpMCHeader.charge = mcpart->Charge;
  ntpMCHeader.momentum[0] = mcpart->Momentum;
  for (int i=0; i<3; i++) {
    ntpMCHeader.dir[i] = mcpart->Dir[i];
    ntpMCHeader.coo[0][i] = mcpart->Coo[i];
  }
  for (int imc=0; imc<nMCEventg(); imc++) {
    MCEventgR* mcpart = pMCEventg(imc);
    if ( (mcpart->Nskip>-1000)||(mcpart->Nskip<-1020) ) continue;
    int ipri = 1+abs(mcpart->Nskip)-1000;
    ntpMCHeader.momentum[ipri] = mcpart->Momentum;
    for (int i=0; i<3; i++) ntpMCHeader.coo[ipri][i] = mcpart->Coo[i];
  }

  // Tracker MC hit
  if (nTrTrack()!=0) {
  for (int i=0; i<9; i++) {
      mchit hh = getJLayerGeantID(this,0,i+1);
      ntpMCHeader.hit_pid[i] = hh.pid;
      ntpMCHeader.hit_edep[i] = hh.edep;
      for (int jj=0; jj<3; jj++) {
        ntpMCHeader.hit_pos[i][jj] = hh.pos[jj];
        ntpMCHeader.hit_mom[i][jj] = hh.mom[jj];
      }
    }
  }

  // count RICH photons production
  map<int,int> photons_radiator;
  map<int,int> photons_pmtplane;
  for (int imc=0; imc<(int)NRichMCCluster(); imc++) {
    RichMCClusterR* rich_mc = pRichMCCluster(imc);
    int trkid = rich_mc->GparentID;
    if (rich_mc->Origin[2]>-100) photons_radiator[trkid]++;
    else                         photons_pmtplane[trkid]++;
  }
  vector<int> usedr;
  vector<int> usedp;
  int indr = 0;
  int indp = 0;
  for (int imc = 0; imc < nMCEventg(); imc++) {
    MCEventgR* mcpart = pMCEventg(imc);
    int trkid = mcpart->trkID;
    if (photons_radiator.find(trkid)!=photons_radiator.end()) {
      if (find(usedr.begin(),usedr.end(),trkid)!=usedr.end()) continue;
      if (indr>9) continue;
      ntpMCHeader.rich_np[0][indr] = photons_radiator[trkid];
      ntpMCHeader.rich_part[0][indr] = mcpart->Particle;
      ntpMCHeader.rich_mom[0][indr] = mcpart->Momentum;
      usedr.push_back(trkid);
      indr++;
    }
    if (photons_pmtplane.find(trkid)!=photons_pmtplane.end()) {
      if (find(usedp.begin(),usedp.end(),trkid)!=usedp.end()) continue;
      if (indp>9) continue;
      ntpMCHeader.rich_np[1][indp] = photons_pmtplane[trkid];
      ntpMCHeader.rich_part[1][indp] = mcpart->Particle;
      ntpMCHeader.rich_mom[1][indp] = mcpart->Momentum;
      usedp.push_back(trkid);
      indp++;
    }
  }
}

void Analysis::FillNtpTrd() {
  if (nTrdRawHit()==0) return;
  // a) Default
  if (pTrdTrack) {
    ntpTrd.theta = pTrdTrack->Theta;
    ntpTrd.phi   = pTrdTrack->Phi;
    ntpTrd.chisq = pTrdTrack->Chi2;
    ntpTrd.nseg  = pTrdTrack->NTrdSegment();
    for (int i=0; i<3; i++) ntpTrd.coo[i] = pTrdTrack->Coo[i];
    ntpTrd.q = pTrdTrack->Q;
  }
  // b) MIT
  int isvalid[5] = {0};
  double lh[5][3] = {{0}};
  double ll[5][3] = {{0}};
  int nh[5] = {0};
  int nhoff[5] = {0};
  float ampoff[5] = {0};
  for (int i=0; i<5; i++) {
    isvalid[i] = -1;
    for (int j=0; j<3; j++) {
      lh[i][j] = -1;
      ll[i][j] = -1;
    }
  }
  if (pTrTrack) {
    // charge
    TrdKCluster* trdK = new TrdKCluster(this,pTrTrack,id_ref_fit);
    ntpTrd.trdk_is_align_ok = (trdK->IsReadAlignmentOK!=0); // 0: Alignment not performed,  1: Static Alignment of Layer level,  2: Dynamic Alignment for entire TRD
    ntpTrd.trdk_is_calib_ok = (trdK->IsReadCalibOK!=0); // 0: Gain Calibration not performed,  1: Gain Calibration Succeeded
    for (int iopt=0; iopt<3; iopt++) { // 2: Only deltas - 1:Only dE/dX - 0:both(default)
      trdK->CalculateTRDCharge(iopt);
      ntpTrd.trdk_q[iopt]     = trdK->GetTRDCharge();
      ntpTrd.trdk_q_err[iopt] = trdK->GetTRDChargeError();
    }
    ntpTrd.trdk_q_upper       = trdK->GetTRDChargeUpper();
    ntpTrd.trdk_q_lower       = trdK->GetTRDChargeLower();
    ntpTrd.trdk_q_nhit        = trdK->GetQTRDHitCollection().size();
    ntpTrd.trdk_q_nhit_nuclei = trdK->GetQTRDHitCollectionNuclei().size();
    ntpTrd.trdk_q_nhit_refit  = trdK->GetQTRDHitCollectionRefit().size();
    ntpTrd.trdk_q_nhit_dr     = trdK->GetQTRDHitCollectionDeltaRay().size();
    // edeps (Melanie)
    AMSPoint pnt;
    AMSDir dir;
    ntpTrd.trdk_nhit = trdK->NHits();
    int x = trdK->GetTrTrackExtrapolation(pnt,dir);
    if (x>=0) {
      for (int i=0; i<TMath::Min(25,trdK->NHits()); i++){
        TrdKHit *hit = trdK->GetHit(i);
        if (!hit) continue;
        float tube_length = hit->Tube_Track_3DLength(&pnt,&dir);
        if (tube_length<=0.02) continue;
        float amplitude = hit->TRDHit_Amp;
        if (amplitude<=0) continue;
        if (amplitude>32767) amplitude = 32767;
        int ilay = hit->TRDHit_Layer;
        if (amplitude<=ntpTrd.trdk_ampl[ilay]) continue; // MaxADC
        ntpTrd.trdk_ampl[ilay] = amplitude;
        ntpTrd.trdk_path[ilay] = tube_length;
      }
    }
    // likelihoods
    float threshold = 15;
    float TRDCenter = 115;
    float total_pathlength, total_amp;
    // [0] inner track
    int fit_id_inner = pTrTrack->iTrTrackPar(1,3,refi_ref_fit,mass_ref_fit,char_ref_fit);
    if (fit_id_inner>=0) {
      trdK->SetTrTrack(pTrTrack,fit_id_inner);
      isvalid[0] = trdK->GetLikelihoodRatio_TrTrack(threshold,lh[0],ll[0],nh[0],total_pathlength,total_amp,-1,0);
      trdK->GetOffTrackHit_TrTrack(nhoff[0],ampoff[0]);
      // trick for getting the deuteron likelihood
      AMSDir Dir;
      AMSPoint P0;
      pTrTrack->Interpolate(TRDCenter,P0,Dir,fit_id_inner);
      float rigidity = pTrTrack->GetRigidity(fit_id_inner);
      trdK->SetTrTrack(&P0,&Dir,2*rigidity);
      trdK->GetLikelihoodRatio_TrTrack(threshold,lh[3],ll[3],nh[3],total_pathlength,total_amp,-1,0);
    }
    // [1] max span track
    int fit_id_maxspan = pTrTrack->iTrTrackPar(1,0,refi_ref_fit,mass_ref_fit,char_ref_fit);
    if (fit_id_maxspan>=0) {
      trdK->SetTrTrack(pTrTrack,fit_id_maxspan);
      isvalid[1] = trdK->GetLikelihoodRatio_TrTrack(threshold,lh[1],ll[1],nh[1],total_pathlength,total_amp,-1,0);
      trdK->GetOffTrackHit_TrTrack(nhoff[1],ampoff[1]);
      // trick for getting the deuteron likelihood
      AMSDir Dir;
      AMSPoint P0;
      pTrTrack->Interpolate(TRDCenter,P0,Dir,fit_id_maxspan);
      float rigidity = pTrTrack->GetRigidity(fit_id_maxspan);
      trdK->SetTrTrack(&P0,&Dir,2*rigidity);
      trdK->GetLikelihoodRatio_TrTrack(threshold,lh[4],ll[4],nh[4],total_pathlength,total_amp,-1,0);
    }
    // [2] refit TrdTrack
    isvalid[2] = trdK->GetLikelihoodRatio_TRDRefit(threshold,lh[2],ll[2],nh[2],total_pathlength,total_amp,1,1,-1,0);
    trdK->GetOffTrackHit_TRDRefit(nhoff[2],ampoff[2]);
    if (trdK) delete trdK;
  }
  for (int i=0; i<3; i++) {
    ntpTrd.trdk_like_valid[i] = isvalid[i];
    ntpTrd.trdk_like_e[i]  = (ll[i][0]>0) ? -log(ll[i][0]) : -1;
    ntpTrd.trdk_like_p[i]  = (ll[i][1]>0) ? -log(ll[i][1]) : -1;
    ntpTrd.trdk_like_he[i] = (ll[i][2]>0) ? -log(ll[i][2]) : -1;
    ntpTrd.trdk_like_nhit[i] = nh[i];
    ntpTrd.trdk_like_nhit_off[i] = nhoff[i];
    ntpTrd.trdk_like_ampl_off[i] = ampoff[i];
    if (i<2) ntpTrd.trdk_like_d[i] = (ll[i+3][1]>0) ? -log(ll[i+3][1]) : -1;
  }
  // simple vertex reconstruction
  TrdSimpleVertexRecon rec;
  rec.Build(this,2,1000,0);
  ntpTrd.vertex_nx = (int)rec.vertices2D[0].size();
  ntpTrd.vertex_ny = (int)rec.vertices2D[1].size();
  if (rec.vertices3D.size()>0) {
    ntpTrd.vertex_nsegx = rec.vertices3D.at(0)->NSegX;
    ntpTrd.vertex_nsegy = rec.vertices3D.at(0)->NSegY;
    for (int i=0; i<3; i++) ntpTrd.vertex_coo[i] = rec.vertices3D.at(0)->Coo[i];
    ntpTrd.vertex_d2 = rec.vertices3D.at(0)->Dist2;
  }
}

void Analysis::FillNtpTof() {
  if (pBeta) {
    ntpTof.evgeni_beta      = pBeta->Beta;
    ntpTof.evgeni_beta_patt = pBeta->Pattern;
  }
  if (pBetaH) {
    // flags
    // tofhflag bits: 0:!IsGoodBeta - 1:!IsTkTofMatch
    int tofhflag = 0;
    if (!pBetaH->IsGoodBeta())   tofhflag = (tofhflag | (1<<0));
    if (!pBetaH->IsTkTofMatch()) tofhflag = (tofhflag | (1<<1));
    // tofhbetapatt
    int tofhbetapatt = 0;
    for (int i=0; i<4; i++) if (!(pBetaH->TestExistHL(i))) tofhbetapatt = (tofhbetapatt | (1<<i));
    // tofhflagp[i]=-1 => No betaH Cluster in plane i
    //    tofhflagp[i]=tofhflagp[i]+100*tofhflags[0]+10000*tofhflags[1];
    //    tofhflagp[i]%100 bits:               0:!IsIsolation - 1:IsOnOverlap
    //    tofhflagp[i]/100)%100 (side 0 bits)  0:!TestExistHS - 1:!IsGoodSide - 2:!IsOneLT - 3:!IsExistHT
    //    tofhflagp[i]/10000)%100 (side 1 bits)
    int tofhflagp[4];
    for (int i=0; i<4; i++) {
      if (pBetaH->TestExistHL(i)) {
        TofClusterHR* ptofhcl = pBetaH->GetClusterHL(i);
        tofhflagp[i] = 0;
        if(!pBetaH->IsIsolationHL(i)) tofhflagp[i] = (tofhflagp[i] | (1<<0));    // 0:!IsIsolation
        AMSPoint fitcoo;
        AMSDir   fitdir;
        // double   time;
        // double pathlength = (pTrTrack) ?
        //   pTrTrack->Interpolate(ptofhcl->Coo[2],fitcoo,fitdir,id_ref_fit) :
        //   pBetaH->TInterpolate(ptofhcl->Coo[2],fitcoo,fitdir,time,false);;
        if (pBetaH->IsInOverlap(i,fitcoo[0],fitcoo[1],2)) tofhflagp[i] = (tofhflagp[i] | (1<<1)); // 1:IsOnOverlap
        int tofhflags[2];
        for (int j=0; j<2; j++) {
          tofhflags[j] = 0;
          if (!ptofhcl->TestExistHS(j)) tofhflags[j] = (tofhflags[j] | (1<<0)); // 0:!TestExistHS
          if (!ptofhcl->IsGoodSide(j))  tofhflags[j] = (tofhflags[j] | (1<<1)); // 1:!IsGoodSide
          if (!ptofhcl->IsOneLT(j))     tofhflags[j] = (tofhflags[j] | (1<<2)); // 2:!IsOneLT
          if (!ptofhcl->IsExistHT(j))   tofhflags[j] = (tofhflags[j] | (1<<3)); // 3:!IsExistHT
        }
        tofhflagp[i] = tofhflagp[i]+100*tofhflags[0]+10000*tofhflags[1];
      }
      else tofhflagp[i] = -1;                                                    // -1: No BetaH in plane
    }
    // Edep [plane][5] 0:Used Cluster - 1:nearby bar - 2:other bars in the plane - 3:Maximum Cluster not in nearby bar
    //                                - 4:overlap cluster not in betah
    int nclp[4][3];
    float eclp[4][5];
    for (int i=0;i<4;i++) {
      for (int j=0;j<3;j++) nclp[i][j] = 0;
      for (int j=0;j<5;j++) eclp[i][j] = 0;
    }
    TofClusterHR *ptofhclbetap[4];
    for (int i=0; i<4; i++) {
      ptofhclbetap[i] = pBetaH->GetClusterHL(i);
      if (!ptofhclbetap[i]) continue;
      ntpTof.time[i] = ptofhclbetap[i]->Time;
    }
    for (int i=0; i<nTofClusterH(); i++) {
      TofClusterHR* ptofhcl = pTofClusterH(i);
      if (!ptofhcl) continue;
      int plane = ptofhcl->Layer;
      int bar = ptofhcl->Bar;
      float edep = ptofhcl->GetEdep();
      int iw = 2;                                          // default: 2:other bars in the plane
      if (ptofhclbetap[plane]) {
        if (ptofhcl==ptofhclbetap[plane]) iw=0;             // 0:Used Cluster
        else if(abs(bar-ptofhclbetap[plane]->Bar)==1) iw=1; // 1:nearby bar
      }
      nclp[plane][iw]++;
      eclp[plane][iw] += edep;
      if(iw==2 && eclp[plane][3]<edep) {eclp[plane][3] = edep;}   // 3:Maximum Cluster not in nearby bar
    }
    float eo[4] = {0};
    float rm[4] = {0};
    float rb[4] = {0};
    // int fl[34] = {0};
    // int ax = Tools::axtof(this,pBetaH,pTrTrack,fl,eo,rm,rb);
    for (int il=0;il<4;il++) eclp[il][4] = eo[il];              // 4:overlap cluster not in betah
    // BetaH track (no tracker info)
    AMSPoint tofhcoo;
    AMSDir tofhdir;
    // double time;
    // double zref = 0;
    // double path = pBetaH->TInterpolate(zref,tofhcoo,tofhdir,time,false);
    if (pBetaH->GetBeta()*tofhdir[2]>0) for(int i=0; i<3; i++) tofhdir[i]=-tofhdir[i];
    float tofhth = acos(tofhdir[2]);
    float tofhph = atan2(tofhdir[1],tofhdir[0]);
    if (tofhph<0) tofhph = tofhph+2.*M_PI;
    ntpTof.flag      = tofhflag;
    ntpTof.trk_ncl   = pBetaH->GetSumHit();
    ntpTof.beta      = isMC ? pBetaH->GetMCBeta() : pBetaH->GetBeta();
    ntpTof.beta_err  = ntpTof.beta*ntpTof.beta*pBetaH->GetEBetaV();
    ntpTof.t0        = pBetaH->GetT0();
    ntpTof.beta_ncl  = pBetaH->GetUseHit();
    ntpTof.beta_patt = tofhbetapatt;
    ntpTof.chisqcn   = pBetaH->GetNormChi2C();
    ntpTof.chisqtn   = pBetaH->GetNormChi2T();
    ntpTof.theta     = tofhth;
    ntpTof.phi       = tofhph;
    for(int i=0;i<3;i++) ntpTof.coo[i] = tofhcoo[i];
    for (int i=0;i<4;i++) {
      ntpTof.flagp[i] = tofhflagp[i];
      for (int j=0;j<3;j++) ntpTof.nclp[i][j] = nclp[i][j];
      for (int j=0;j<5;j++) ntpTof.edep[i][j] = eclp[i][j];
      ntpTof.ovresm[i] = rm[i];
      ntpTof.ovresb[i] = rb[i];
    }
    // charge
    ntpTof.q = pBetaH->GetQ(ntpTof.q_nhit,ntpTof.q_err);
    for (int il=0; il<4; il++) {
      ntpTof.q_lay       [il] = pBetaH->GetQL(il);
      ntpTof.q_lay_uncorr[il] = pBetaH->GetQL(il,2,TofRecH::kThetaCor|TofRecH::kBirkCor|TofRecH::kReAttCor|TofRecH::kQ2Q);
    }
    ntpTof.z = pBetaH->GetZ(ntpTof.z_nhit,ntpTof.z_like);
    if (ntpTof.z_like>0) ntpTof.z_like = -log(ntpTof.z_like);
    else ntpTof.z_like = -1;
    Tools::TofUnusedHits(this,pBetaH,ntpTof.clsdt,ntpTof.clsed,ntpTof.clsn);
  }
}

void Analysis::FillNtpTracker() {
  if (!pTrTrack) return;
  if (nTrRecHit()==0) return;
  // hit patterns
  ntpTracker.patty  = pTrTrack->GetBitPatternJ();
  ntpTracker.pattxy = pTrTrack->GetBitPatternXYJ();
  // sensor grid sequence
  int s0,s2,s8;
  int ks = GetTrSensorGridID(pTrTrack,s0,s2,s8);
  if (ks>=0) {
    ntpTracker.s0 = s0;
    ntpTracker.s2 = s2;
    ntpTracker.s8 = s8;
  }
  // charge
  float beta = (pBetaH) ? (isMC ? pBetaH->GetMCBeta() : pBetaH->GetBeta()) : 1;
  mean_t qold = pTrTrack->GetInnerQ_all(beta);
  ntpTracker.q_inn[0]      = qold.Mean;
  ntpTracker.q_inn_nhit[0] = qold.NPoints;
  ntpTracker.q_inn_rms[0]  = qold.RMS;
  if (!disable) {
    TrQYJTrack qyj = pTrTrack->GetQYJ_all(2,beta);
    ntpTracker.q_inn[1]      = qyj.InnerQ;
    ntpTracker.q_inn_nhit[1] = qyj.InnerQPoints;
    ntpTracker.q_inn_rms[1]  = qyj.InnerQRMS;
  }
  mean_t qhl = pTrTrack->GetInnerQH_all(2,beta);
  ntpTracker.q_inn[2]      = qhl.Mean;
  ntpTracker.q_inn_nhit[2] = qhl.NPoints;
  ntpTracker.q_inn_rms[2]  = qhl.RMS;
  for (int il = 0; il<9; il++) {
    ntpTracker.q_lay[0][il] = pTrTrack->GetLayerJQ(il+1,beta);
    if (!disable) ntpTracker.q_lay[1][il] = pTrTrack->GetLayerQYJ(il+1,2,beta);
    ntpTracker.q_lay[2][il] = pTrTrack->GetLayerJQH(il+1,2,beta);
    ntpTracker.q_lay_uncorr[0][il] = pTrTrack->GetLayerJQ(il+1);
    if (!disable) ntpTracker.q_lay_uncorr[1][il] = pTrTrack->GetLayerQYJ(il+1,2);
    ntpTracker.q_lay_uncorr[2][il] = pTrTrack->GetLayerJQH(il+1,2);
  }
  for (int il=0; il<9; il++) {
    TrRecHitR *hit = pTrTrack->GetHitLJ(il + 1);
    if (!hit)
      continue;
    if (hit->GetXCluster()) {
      ntpTracker.q_clu[0][il] = hit->GetXCluster()->GetQ(beta);
      ntpTracker.q_clu_uncorr[0][il] = hit->GetXCluster()->GetQ(1);
      ntpTracker.q_clu_status[0][il] = hit->GetXCluster()->GetQStatus();
    }
    if (hit->GetYCluster()) {
      ntpTracker.q_clu[1][il] = hit->GetYCluster()->GetQ(beta);
      ntpTracker.q_clu_uncorr[1][il] = hit->GetYCluster()->GetQ(1);
      ntpTracker.q_clu_status[1][il] = hit->GetYCluster()->GetQStatus();
    }
  }
  // all hit coordinates
  for (int il=0; il<9; il++) {
    TrRecHitR* hit = pTrTrack->GetHitLJ(il+1);
    if (!hit) continue;
    AMSPoint point = hit->GetCoord();
    for (int icoo=0; icoo<3; icoo++) ntpTracker.coo_lay[il][icoo] = point[icoo];
  }
  // interpolation
  for(int il=0; il<9; il++) {
    AMSPoint point;
    AMSDir dir;
    pTrTrack->Interpolate(tracker_layers_z[il],point,dir,id_ref_fit);
    for (int icoo=0; icoo<3; icoo++) ntpTracker.int_lay[il][icoo] = point[icoo];
  }
  // interpolation @ RICH radiator
  AMSPoint point;
  AMSDir   dir;
  pTrTrack->Interpolate(rich_radiator_z,point,dir,id_ref_fit);
  ntpTracker.int_rich_rad[0] = point[0];
  ntpTracker.int_rich_rad[1] = point[1];
  ntpTracker.theta_rich_rad  = dir.gettheta();
  ntpTracker.phi_rich_rad    = dir.getphi();
  // interpolation @ RICH PMT plane
  pTrTrack->Interpolate(rich_pmt_plane_z,point,dir,id_ref_fit);
  ntpTracker.int_rich_pmt[0] = point[0];
  ntpTracker.int_rich_pmt[1] = point[1];
  // directional cutoff
  pTrTrack->Interpolate(top_z,point,dir,id_ref_fit);
  AMSDir dirc(M_PI-dir.gettheta(),M_PI+dir.getphi());
  double rcut = 0;
  if (!GetStoermerCutoff(rcut, 1,dirc)) ntpTracker.stoermer_cutoff[0] = rcut;
  if (!GetStoermerCutoff(rcut,-1,dirc)) ntpTracker.stoermer_cutoff[1] = rcut;
  if (!GetIGRFCutoff    (rcut, 1,dirc)) ntpTracker.igrf_cutoff[0] = rcut;
  if (!GetIGRFCutoff    (rcut,-1,dirc)) ntpTracker.igrf_cutoff[1] = rcut;
  double thetaGl,phiGl;
  Loc2Gl(dirc.gettheta(),dirc.getphi(),thetaGl,phiGl);
  float radS=fHeader.RadS,latS=fHeader.ThetaS,phiS=fHeader.PhiS,thetaP=thetaGl,phiP=phiGl;
  ntpTracker.igrf_cutoff_upper[0] = rigidityCutoffPositiveUpper(radS,latS,phiS,thetaP,phiP);
  ntpTracker.igrf_cutoff_lower[0] = ntpTracker.igrf_cutoff_upper[0] -
    rigidityCutoffPositivePenumbra(latS,phiS,thetaP,phiP);
  ntpTracker.igrf_cutoff_upper[1] = rigidityCutoffNegativeUpper(radS,latS,phiS,thetaP,phiP);
  ntpTracker.igrf_cutoff_lower[1] = ntpTracker.igrf_cutoff_upper[1] -
    rigidityCutoffNegativePenumbra(latS,phiS,thetaP,phiP);
  // edep
  int   nclx[7] = {0};
  float eclx[7] = {0};
  int   ncly[7] = {0};
  float ecly[7] = {0};
  for(int il=1; il<=9; il++) {
    float maxedepx   = 0;
    float d2maxedepx = 0;
    float maxedepy   = 0;
    float d2maxedepy = 0;
    Tools::TrClusterOnLayerJ(il,this,*pTrTrack,id_ref_fit,nclx,eclx,ncly,ecly,maxedepx,d2maxedepx,maxedepy,d2maxedepy);
    for (int i=0; i<7; i++) {
      // per layers
      ntpTracker.nclu_lay[0][il-1][i] = nclx[i];
      ntpTracker.edep_lay[0][il-1][i] = eclx[i];
      ntpTracker.nclu_lay[1][il-1][i] = ncly[i];
      ntpTracker.edep_lay[1][il-1][i] = ecly[i];
    }
    ntpTracker.max_edep_lay[0][il-1]   = maxedepx;
    ntpTracker.max_edep_lay[1][il-1]   = maxedepy;
    ntpTracker.d2_max_edep_lay[0][il-1] = d2maxedepx;
    ntpTracker.d2_max_edep_lay[1][il-1] = d2maxedepy;
  }
  // magnet temperature correction
  double RigidityCorrection = 1;
  if (!isMC) {
    float fact;
    int retv(0);
    float mfcor[2] = {0};
    retv = MagnetVarp::btempcor(fact,0,1);
    if (retv==0) mfcor[0] = fact;
    retv = MagnetVarp::btempcor(fact,0,2);
    if (retv==0) mfcor[1] = fact;
    if      (mfcor[0]&&mfcor[1]) RigidityCorrection = (mfcor[0]+mfcor[1])/2;
    else if (mfcor[0])           RigidityCorrection = mfcor[0];
    else if (mfcor[1])           RigidityCorrection = mfcor[1];
  }
  if (RigidityCorrection==0) RigidityCorrection=1;
  // refit removing layer one-by-one
  int pat9inner = 0;
  for(int il=2; il<9; il++) {
    if (!pTrTrack->TestHitLayerJ(il)) continue;
    pat9inner = pat9inner + 9*pow(10,il-1);
  }
  for(int il=2; il<9; il++) {
    if (!pTrTrack->TestHitLayerJ(il)) continue;
    int patt = pat9inner - 9*pow(10,il-1);
    int fit_id = pTrTrack->iTrTrackPar(algo_ref_fit,patt,21,mass_ref_fit,char_ref_fit);
    if (fit_id<0) continue;
    const TrTrackPar &t = pTrTrack->gTrTrackPar(fit_id);
    ntpTracker.inn_rig[il-2] = t.Rigidity*(!t.BcorrFlag?RigidityCorrection:1);
    ntpTracker.inn_invrigerr[il-2] = t.ErrRinv/(!t.BcorrFlag?RigidityCorrection:1);
    AMSPoint point;
    AMSDir   dir;
    pTrTrack->Interpolate(ntpTracker.coo_lay[il-2][2],point,dir,fit_id);
    for (int icoo=0; icoo<2; icoo++) ntpTracker.inn_res[il-2][icoo] = ntpTracker.coo_lay[il-1][icoo]-point[icoo];
    ntpTracker.inn_chisqn[il-2][0] = (t.NdofX>0)? t.ChisqX/t.NdofX:0;
    ntpTracker.inn_chisqn[il-2][1] = (t.NdofY>0)? t.ChisqY/t.NdofY:0;
  }

if (Event()==3858) cout << "Analysis::FillNtpTracker" << endl;

  // other useful stuff for rare events search
  for (int j=0; j<7; j++) {
    ntpTracker.signal_ratio[j] = GetTrackerRawSignalRatio(j+2,10);
    ntpTracker.feet_dist[j] = GetTkFeetDist(j+2,pTrTrack);
  }
  //////////////////////////////////////////////////////////////////////////
  // Refits
  // 0: Max available span, Kalman fit (algo=6)
  // 1: Inner Tracker, Kalman fit (algo=6)
  // 2: Upper Inner Tracker (exclude Layer 7 and 8), Kalman fit (algo=6)
  // 3: Lower Inner Tracker (exclude Layer 2), Kalman fit (algo=6)
  // 4: Inner Tracker + Layer 1, Kalman fit (algo=6)
  // 5: Inner Tracker + Layer 9, Kalman fit (algo=6)
  // 6: Inner Tracker, Vitaly fit (algo=1)
  // 7: Inner Tracker, Vitaly fit no MS (algo=21) < IMPORTANTE
  // 8: Inner Tracker, Chickanian fit (algo=3)
  //////////////////////////////////////////////////////////////////////////
  int l26 = 0;
  int l38 = 0;
  int l34 = 0;
  int l56 = 0;
  for(int il=2; il<9; il++) {
    if (!pTrTrack->TestHitLayerJ(il)) continue;
    if (il<7) l26 += 9*pow(10,il-1);
    if (il>2) l38 += 9*pow(10,il-1);
    if ( (il!=3)&&(il!=4) ) l34 += 9*pow(10,il-1);
    if ( (il!=5)&&(il!=6) ) l56 += 9*pow(10,il-1);
  }
  int nfit = 9;
  int algo[9]       = { 6, 6,  6,  6, 6, 6, 1,21, 3};
  int patt[9]       = { 0, 3,l26,l38, 5, 6, 3, 3, 3};
  int refit_data[9] = {21,21, 21, 21,21,21,21,21,21};
  int refit_mc[9]   = {21,21, 21, 21,21,21,21,21,21};
  int fit_id[9]     = {-1,-1, -1, -1,-1,-1,-1,-1,-1};
  for (int ifit=0; ifit<nfit; ifit++) {
    int refit = (!isMC) ? refit_data[ifit] : refit_mc[ifit];
    TrFit fit;
    fit_id[ifit] = pTrTrack->iTrTrackPar(fit,algo[ifit],patt[ifit],refit,mass_ref_fit,char_ref_fit);
    if (fit_id[ifit]<0) continue;
    const TrTrackPar &t        = pTrTrack->gTrTrackPar(fit_id[ifit]);
    ntpTracker.invrigerr[ifit] = t.ErrRinv/(!t.BcorrFlag?RigidityCorrection:1);
    ntpTracker.theta[ifit]     = t.Dir.gettheta();
    ntpTracker.phi[ifit]       = t.Dir.getphi();
    AMSPoint tkcoo = t.P0;
    for (int icoo=0; icoo<3; icoo++) ntpTracker.coo[ifit][icoo] = tkcoo[icoo];
    ntpTracker.chisqn[ifit][0] = (t.NdofX>0)? t.ChisqX/t.NdofX:0;
    ntpTracker.chisqn[ifit][1] = (t.NdofY>0)? t.ChisqY/t.NdofY:0;
    AMSPoint point;
    AMSDir   dir;
    double   rigi;
    if ((algo[ifit]%10)==6) fit.InterpolateKalman(tracker_layers_z[1],point,dir,rigi);
    else pTrTrack->Interpolate(tracker_layers_z[1],point,dir,fit_id[ifit]);
    ntpTracker.int_l2[ifit][0] = point.x();
    ntpTracker.int_l2[ifit][1] = point.y();
    for (int iz=0; iz<3; iz++) ntpTracker.rig[ifit][iz] = pTrTrack->GetRigidity(fit_id[ifit],iz)*(!t.BcorrFlag?RigidityCorrection:1);
  }
  // second track stuff
  if (pTrTrack2) {
    double beta = (pBetaH) ? (isMC ? pBetaH->GetMCBeta() : pBetaH->GetBeta()) : 1;
    float q_inn = pTrTrack2->GetInnerQ(beta);
    int charge = floor(0.5+q_inn);
    if (charge<1) charge = 1;
    float mass = (charge>=2) ? TrFit::Mhelium/2*charge : TrFit::Mproton;
    ntpTracker.sec_patty = pTrTrack2->GetBitPatternJ();
    ntpTracker.sec_pattxy = pTrTrack2->GetBitPatternXYJ();
    ntpTracker.sec_inn_q = q_inn;
    int id = pTrTrack2->iTrTrackPar(algo_ref_fit,patt_ref_fit,refi_ref_fit,mass,charge);
    if (id>0) {
      const TrTrackPar &t2 = pTrTrack2->gTrTrackPar(id);
      ntpTracker.sec_inn_rig = t2.Rigidity*(!t2.BcorrFlag?RigidityCorrection:1);
      ntpTracker.sec_inn_chisqn[0] = (t2.NdofX>0)? t2.ChisqX/t2.NdofX:0;
      ntpTracker.sec_inn_chisqn[1] = (t2.NdofY>0)? t2.ChisqY/t2.NdofY:0;
    }
  }
}


void Analysis::FillNtpTrackerScattering() {

  if (!pTrTrack) return;
  if (nTrRecHit()==0) return;

  // init
  double beta[2] = {fabs(ntpTof.beta),fabs(ntpRich.beta_corrected)};
  double mass[2] = {0.938272297,1.875612928};
  double zref[2] = {27.2,-0.3};
  AMSPoint point[7];
  for (int i=0; i<7; i++) {
    TrRecHitR* hit = pTrTrack->GetHitLJ(i+2);
    if (!hit) continue;
    point[i] = hit->GetCoord();
  }

  // calculate scattering angles using beta (TOF/RICH) and points of the inner tracker
  // 1. ims=0: material between L3-L4, z=27.2, fix rigidity from beta and use the two fits 234 and 3456
  // 2. ims=1: material between L5-L6, z=-0.3, fix rigidity from beta and use the two fits 3456 and 5678
  for (int ibeta=0; ibeta<2; ibeta++) { // TOF/RICH
    if (beta[ibeta]<=0) continue;
    if (beta[ibeta]>0.999) beta[ibeta] = 0.999;
    for (int imass=0; imass<2; imass++) { // proton/deuteron
      double rigidity = mass[imass]*beta[ibeta]*sqrt(1-beta[ibeta]*beta[ibeta]);
      for (int ims=0; ims<2; ims++) {
        double angle_x[2] = {0};
        double angle_y[2] = {0};
        for (int iul=0; iul<2; iul++) { // L34/L56
          TrFit fit;
          fit.SetRigidity(rigidity);
          AMSPoint err(10e-4,10e-4,100e-4);
          int lay[2] = {0};
          if      ( (ims==0)&&(iul==0) ) { lay[0] = 2; lay[1] = 4; }
          else if ( (ims==0)&&(iul==1) ) { lay[0] = 3; lay[1] = 6; }
          else if ( (ims==1)&&(iul==0) ) { lay[0] = 3; lay[1] = 6; }
          else if ( (ims==1)&&(iul==1) ) { lay[0] = 5; lay[1] = 8; }
          int nl = 0;
          AMSPoint p1;
          AMSPoint p2;
          for (int jl=lay[0]; jl<=lay[1]; jl++) {
            if (fabs(point[jl-2].z())<=0) continue;
            fit.Add(point[jl-2],err);
            if (nl==0) p1 = point[jl-2];
            if (nl==1) p2 = point[jl-2];
            nl++;
          }
          if (nl==2) {
            AMSDir dir = p2-p1;
            AMSPoint pp = p2;
            for (int iter=0; iter<4; iter++) {
              TrProp tr(p1,dir,rigidity);
              tr.Propagate(p2.z());
              pp[0] -= tr.GetP0x()-p2.x();
              pp[1] -= tr.GetP0y()-p2.y();
              dir = pp-p1;
            }
            TrProp prop(p1,dir,rigidity);
            prop.Propagate(zref[ims]);
            angle_x[iul] = TMath::ATan(prop.GetDxDz())*TMath::RadToDeg();
            angle_y[iul] = TMath::ATan(prop.GetDyDz())*TMath::RadToDeg();
          }
          else {
            fit.AlcarazFit(1);
            fit.Propagate(zref[ims]);
            angle_x[iul] = TMath::ATan(fit.GetDxDz())*TMath::RadToDeg();
            angle_y[iul] = TMath::ATan(fit.GetDyDz())*TMath::RadToDeg();
          }
        }
        ntpTracker.scat_beta_theta[ibeta][imass][ims][0] = angle_x[1]-angle_x[0];
        ntpTracker.scat_beta_theta[ibeta][imass][ims][1] = angle_y[1]-angle_y[0];
      }
    }
  }

  // calculate scattering angles using only track information
  // 1. ims=0: material between L3-L4, z=27.2, fix rigidity with 345678 fit, then use this rigidity to fit 234 vs 345678
  // 2. ims=1: material between L5-L6, z=-0.3, fix rigidity with 23456 fit, then use this rigidity to fit 5678
  for (int ims=0; ims<2; ims++) {
    double angle_x[2] = {0};
    double angle_y[2] = {0};
    double rigidity = 0;
    // eval rigidity and direction with a fit
    {
      int lay[2][2] = { {3,8},{2,6} };
      int patt = 0;
      for (int jl=lay[ims][0]; jl<=lay[ims][1]; jl++) {
        if (!pTrTrack->TestHitLayerJ(jl)) continue;
        patt += 9*pow(10,jl-1);
      }
      TrFit fit;
      int fit_id = pTrTrack->iTrTrackPar(fit,6,patt,refi_ref_fit,mass_ref_fit,char_ref_fit); // two mass hypothesys do not matter
      if (fit_id<0) continue;
      AMSPoint point;
      AMSDir dir;
      fit.InterpolateKalman(zref[ims],point,dir,rigidity);
      angle_x[0] = TMath::ATan(fit.GetDxDz())*TMath::RadToDeg();
      angle_y[0] = TMath::ATan(fit.GetDyDz())*TMath::RadToDeg();
    }
    // eval direction with fixed rigidity
    {
      int lay[2][2] = { {2,4},{5,8} };
      TrFit fit;
      fit.SetRigidity(rigidity);
      AMSPoint err(10e-4,10e-4,100e-4);
      int nl = 0;
      AMSPoint p1;
      AMSPoint p2;
      for (int jl=lay[ims][0]; jl<=lay[ims][1]; jl++) {
        if (fabs(point[jl-2].z())<=0) continue;
        fit.Add(point[jl-2],err);
        if (nl==0) p1 = point[jl-2];
        if (nl==1) p2 = point[jl-2];
        nl++;
      }
      if (nl==2) {
        AMSDir dir = p2-p1;
        AMSPoint pp = p2;
        for (int iter=0; iter<4; iter++) {
          TrProp tr(p1,dir,rigidity);
          tr.Propagate(p2.z());
          pp[0] -= tr.GetP0x()-p2.x();
          pp[1] -= tr.GetP0y()-p2.y();
          dir = pp-p1;
        }
        TrProp prop(p1,dir,rigidity);
        prop.Propagate(zref[ims]);
        angle_x[1] = TMath::ATan(prop.GetDxDz())*TMath::RadToDeg();
        angle_y[1] = TMath::ATan(prop.GetDyDz())*TMath::RadToDeg();
      }
      else {
        fit.AlcarazFit(1);
        fit.Propagate(zref[ims]);
        angle_x[1] = TMath::ATan(fit.GetDxDz())*TMath::RadToDeg();
        angle_y[1] = TMath::ATan(fit.GetDyDz())*TMath::RadToDeg();
      }
    }
    ntpTracker.scat_rigi_theta[ims][0] = angle_x[1]-angle_x[0];
    ntpTracker.scat_rigi_theta[ims][1] = angle_y[1]-angle_y[0];
  }

}

void Analysis::FillNtpRich() {

  ///////////////////////////
  // Charge Centers
  ///////////////////////////

  // CIEMAT default charge center finding
  ntpRich.nparticle       = RichHitR::getPMTs()-RichHitR::getPMTs(false);
  ntpRich.tot_hit_uncorr  = NRichHit();
  ntpRich.tot_pmt_uncorr  = RichHitR::getPMTs();
  ntpRich.tot_p_uncorr    = RichHitR::getCollectedPhotoElectrons();
  ntpRich.max_p_uncorr[0] = RichHitR::getMaxPMTCol(true);
  ntpRich.max_p_uncorr[1] = RichHitR::getMaxPMTCol(false);

  // get the PMT calibration, if available
  RichPMTCalib* pmt_calib = (!isMC) ? RichPMTCalib::getHead() : 0;

  // occupancy
  ntpRich.occupancy_status = (pRichOcc) ? pRichOcc->GetNBad() : -1;

  // integrate on PMT-by-PMT for alternative charge center finding
  class pmt_info {
   public:
    int pmt;
    int nhit;
    float np;
    float coo[3];
    float dist;
    pmt_info() : pmt(-1), nhit(0), np(0), coo{0,0,0}, dist(0) {}
    void print() const { printf("pmt_info::print pmt:%3d nhit:%2d npe:%7.2f (x,y,z)=(%7.2f,%7.2f,%7.2f) dist=%7.2f\n",pmt,nhit,np,coo[0],coo[1],coo[2],dist); }
  };
  map<int,pmt_info> pmt_info_map;
  vector<int> crossed_pmt[2]; // primary/secondary
  for (int ihit=0; ihit<nRichHit(); ihit++) {
    RichHitR* hit = pRichHit(ihit);
    int pmt = int(hit->Channel/16);
    if ((pmt_calib)&&(find(pmt_calib->BadPMTs.begin(),pmt_calib->BadPMTs.end(),pmt)!=pmt_calib->BadPMTs.end())) continue;
    float np = hit->Npe;
    pmt_info_map[pmt].pmt = pmt;
    pmt_info_map[pmt].nhit++;
    pmt_info_map[pmt].np += np;
    for (int icoo=0; icoo<3; icoo++) pmt_info_map[pmt].coo[icoo] += np*hit->Coo[icoo];
  }
  for (map<int,pmt_info>::iterator it=pmt_info_map.begin(); it!=pmt_info_map.end(); it++) {
    for (int icoo=0; icoo<3; icoo++) it->second.coo[icoo] /= it->second.np;
    if (pTrTrack) {
      AMSPoint coo;
      AMSDir   dir;
      pTrTrack->Interpolate(it->second.coo[2],coo,dir,id_ref_fit);
      it->second.dist = sqrt(pow(it->second.coo[0]-coo.x(),2.)+pow(it->second.coo[1]-coo.y(),2.));
    }
    if (it->second.np<5) continue;
    crossed_pmt[(it->second.dist<3.5)?0:1].push_back(it->first);
  }
  typedef function<bool(pair<int,pmt_info>, pair<int,pmt_info>)> comparator;
  comparator functor = [](pair<int,pmt_info> elem1, pair<int,pmt_info> elem2) { return elem1.second.np > elem2.second.np; };
  set<pair<int,pmt_info>,comparator> pmt_info_set(pmt_info_map.begin(),pmt_info_map.end(),functor);
  int index = 0;
  for (set<pair<int,pmt_info>>::iterator it=pmt_info_set.begin(); it!=pmt_info_set.end(); it++) {
    if (index>4) continue;
    ntpRich.pmt_nhit_uncorr[index] = it->second.nhit;
    ntpRich.pmt_np_uncorr[index] = it->second.np;
    ntpRich.pmt_dist[index] = it->second.dist;
    ntpRich.pmt_pmt[index] = it->first;
    index++;
  }

  // CIEMAT ring reconstruction
  if (pRichRing) {
    // flag
    ntpRich.selection            = Tools::RichQC(pRichRing);
    // general
    ntpRich.is_naf               = pRichRing->IsNaF();
    ntpRich.status               = pRichRing->Status;
    ntpRich.correction_status    = pRichRing->PmtCorrectionsFailed();
    // ring
    ntpRich.nhit                 = pRichRing->getUsedHits();
    ntpRich.nhit_uncorr          = pRichRing->getUsedHits(false);
    ntpRich.nhit_refl            = pRichRing->getReflectedHits();
    ntpRich.npmt                 = pRichRing->NpColPMT.size();
    ntpRich.npmt_uncorr          = pRichRing->getPMTs();
    ntpRich.np                   = pRichRing->getPhotoElectrons();
    ntpRich.np_uncorr            = pRichRing->getPhotoElectrons(false);
    for(int i=0;i<10;i++)
      ntpRich.np_w[i]            = pRichRing->NpColWindow[i];
    // expected
    ntpRich.npmt_exp             = pRichRing->NpExpPMT.size();
    ntpRich.np_exp               = pRichRing->getExpectedPhotoElectrons();
    ntpRich.np_exp_elec          = RichRingR::ComputeNpExp(pRichRing->pTrTrack(),1);
    ntpRich.np_exp_uncorr        = pRichRing->getExpectedPhotoElectrons(false);
    // >>> expected number of photons in the active in an area: ACCEPTANCE EVALUATOR
    // reconstruction
    ntpRich.prob                 = pRichRing->getProb();
    ntpRich.width                = pRichRing->getWidth();
    ntpRich.udist                = pRichRing->UDist;
    ntpRich.tile_id              = pRichRing->getTileIndex();
    const float* TrackEmissionPoint = pRichRing->getTrackEmissionPoint();
    ntpRich.rad_coo[0]           = TrackEmissionPoint[0];
    ntpRich.rad_coo[1]           = TrackEmissionPoint[1];
    ntpRich.rad_theta            = pRichRing->getTrackTheta();
    ntpRich.rad_phi              = pRichRing->getTrackPhi();
    ntpRich.distance_tile_border = pRichRing->DistanceTileBorder();
    ntpRich.q                    = sqrt(pRichRing->getCharge2Estimate(true));
    ntpRich.q_consistency        = pRichRing->getPMTChargeConsistency();
    ntpRich.q_res                = pRichRing->getChargeExpectedResolution();
    ntpRich.q_rms                = pRichRing->getChargeExpectedRms();
    ntpRich.beta                 = pRichRing->getBeta();
    ntpRich.beta_raw             = pRichRing->Beta;
    ntpRich.beta_refit           = pRichRing->BetaRefit;
    ntpRich.beta_corrected       = pRichCorr->correctBeta(ntpRich.beta_raw,ntpRich.rad_coo[0],ntpRich.rad_coo[1],
                                                          ntpRich.rad_theta,ntpRich.rad_phi,
                                                          abs(ntpMCHeader.charge),(nMCEventg()>0));

    // just some check values
    if ((Run()==1381257029)&&(Event()==3752)) {
      cout << "FillNtpRich::Beta: " << ntpRich.beta << " BetaRaw:" << ntpRich.beta_raw << " BetaRefit:" << ntpRich.beta_refit << endl;
      cout << "FillNtpRich::CheckBetaCorrection Run:" << Run() << " Event:" << Event() << " " << ntpRich.beta << " "
           << ntpRich.rad_coo[0] << " " << ntpRich.rad_coo[1] << " " << ntpRich.rad_theta << " " << ntpRich.rad_phi << " "
           << abs(ntpMCHeader.charge) << " " << (nMCEventg()>0) << " ----> " << ntpRich.beta_corrected << endl;
    }

    // expected resolution (C.Delgado)
    GeomHashEnsemble* corr = pRichRing->IsNaF() ? pRichUnifNaf : pRichUnifAgl;
    float x     = pRichRing->AMSTrPars[0];
    float y     = pRichRing->AMSTrPars[1];
    float theta = pRichRing->AMSTrPars[3];
    float phi   = pRichRing->AMSTrPars[4];
    float vx    = sin(theta)*cos(phi);
    float vy    = sin(theta)*sin(phi);
    if (cos(theta)>0) { vx*=-1; vy*=-1; }
    corr->Eval(x,y,vx,vy);
    ntpRich.beta_unifcorr = pRichRing->Beta/corr->MeanPeak;
    ntpRich.beta_res = corr->MeanPeakWidth;
    ntpRich.beta_rms = corr->MeanRms;
  }

  // LIP ring reconstruction
  if (pRichRingB) {
    ntpRich.lip_q    = pRichRingB->ChargeRec;
    ntpRich.lip_beta = pRichRingB->Beta;
  }

  // Re-clustering (triggers calculation of beta-per-hit)
  if (pRichRing) {
    // clustering
    float beta_correction = (ntpRich.beta>0) ? ntpRich.beta_corrected/ntpRich.beta : 1;
    ntpRich.nclus = pRichRing->ClusterizeZ1();
    if (ntpRich.nclus>0) {
      for (int iclu=0; iclu<TMath::Min(10,ntpRich.nclus); iclu++) {
        pRichRing->GetClusters(iclu,ntpRich.clus_size[iclu],ntpRich.clus_mean[iclu],ntpRich.clus_rms[iclu]);
        ntpRich.clus_mean[iclu] *= beta_correction;
        ntpRich.clus_rms[iclu] *= beta_correction;
      }
    }
  }

  // Hit infos (including beta-per-hit)
  class hit_info {
   public:
    unsigned int status;  // Carlos status, bit 30 crossed
    unsigned int status2; // bit0: Jorge's bad PMT; bit1: occupancy table; bit2: crossed by primary (my algo.); bit3: crossed by secondary (my algo.); bit4: is not in selected ring
    int   used;           // Used in whatever ring
    int   channel;        // PMT*16+pixel
    float np;             // Number of photons (uncorrected? corrected for gain?)
    float beta[2];        // Beta hypothesys if existing
    hit_info(unsigned int s, unsigned int s2, int u, int c, float n, float b0, float b1) : status(s), status2(s2), used(u), channel(c), np(n) { beta[0] = b0; beta[1] = b1; }
    void print() const { printf("hit_info::print status:%10d status2:%3d used:%2d channel:%5d np:%7.2f beta(%7.4f,%7.4f)\n",status,status2,used,channel,np,beta[0],beta[1]); }
  };
  vector<hit_info> hit_info_list;
  vector<RichHitR*> rich_hit_on_ring_list;
  vector<RichHitR*> rich_hit_beta_hyp_list;
  map<int,RichHitR*> rich_hit_offline_map;
  if (pRichRing) {
    // selected ring
    for (int ihit=0; ihit<pRichRing->getUsedHits(false); ihit++)
      rich_hit_on_ring_list.push_back(pRichHit(pRichRing->iRichHit(ihit)));
    // beta hypothesys
    int which = 0;
    for (RichOffline::RichRawEvent* hit = new RichOffline::RichRawEvent(this); hit; hit=hit->next(),which++) {
      RichHitR* rich_hit = hit->getpointer();
      rich_hit_offline_map.insert(pair<int,RichHitR*>(which,rich_hit));
    }
    // loop on beta hypothesys
    for (int i=0; i<pRichRing->RawBetas(); i++) {
      int ihit = pRichRing->HitBeta(i);
      RichHitR* rich_hit = rich_hit_offline_map[ihit];
      if (!rich_hit) continue;
      int pmt = int(rich_hit->Channel/16);
      bool is_bad_pmt = (pmt_calib)&&(find(pmt_calib->BadPMTs.begin(),pmt_calib->BadPMTs.end(),pmt)!=pmt_calib->BadPMTs.end());
      bool is_bad_occupancy = (pRichOcc) ? (!pRichOcc->IsGood(rich_hit->Channel)) : true;
      bool is_crossed_pri = find(crossed_pmt[0].begin(),crossed_pmt[0].end(),pmt)!=crossed_pmt[0].end();
      bool is_crossed_sec = find(crossed_pmt[1].begin(),crossed_pmt[1].end(),pmt)!=crossed_pmt[1].end();
      bool is_not_in_selected_ring = (find(rich_hit_on_ring_list.begin(),rich_hit_on_ring_list.end(),rich_hit)==rich_hit_on_ring_list.end());
      int status2 = ((is_bad_pmt)?0x1:0)+((is_bad_occupancy)?0x2:0)+((is_crossed_pri)?0x4:0)+((is_crossed_sec)?0x8:0)+((is_not_in_selected_ring)?0x10:0);
      hit_info hit( rich_hit->Status,status2,pRichRing->UsedBeta(i),rich_hit->Channel,rich_hit->Npe,pRichRing->RawBeta(i,0),pRichRing->RawBeta(i,1));
      rich_hit_beta_hyp_list.push_back(rich_hit);
      hit_info_list.push_back(hit);
      // integrate beta = 1 hypothesys with only good hits
      if ((rich_hit->Status>>31)&0x1) continue;
      if ((status2&0xf)!=0) continue;
      int kind_of_tile = 0;
      if (pRichRing->IsNaF()) kind_of_tile = 1;
      float simple_resolution[2] = {1.2e-3,4e-3};
      for (int ibeta=0; ibeta<2; ibeta++) {
        if (fabs(pRichRing->RawBeta(i,ibeta)-1)>3.*simple_resolution[kind_of_tile]) continue;
        ntpRich.tot_hyp_hit_uncorr[0][ibeta]++;
        ntpRich.tot_hyp_p_uncorr[0] += rich_hit->Npe;
        if ( (pRichRing->UsedBeta(i)<0)||(pRichRing->UsedBeta(i)>1)) {
          ntpRich.tot_hyp_hit_uncorr[1][ibeta]++;
          ntpRich.tot_hyp_p_uncorr[1] += rich_hit->Npe;
        }
      }
    }
  }
  // loop on all other hits to be added to the collection
  int nrichhit = 0;
  for (int ihit=0; ihit<(int)NRichHit(); ihit++) {
    RichHitR* rich_hit = pRichHit(ihit);
    if (!rich_hit) continue;
    nrichhit++;
    if (find(rich_hit_beta_hyp_list.begin(),rich_hit_beta_hyp_list.end(),rich_hit)!=rich_hit_beta_hyp_list.end()) continue; // already included
    int pmt = int(rich_hit->Channel/16);
    bool is_bad_pmt = (pmt_calib)&&(find(pmt_calib->BadPMTs.begin(),pmt_calib->BadPMTs.end(),pmt)!=pmt_calib->BadPMTs.end());
    bool is_bad_occupancy = (pRichOcc) ? (!pRichOcc->IsGood(rich_hit->Channel)) : true;
    bool is_crossed_pri = find(crossed_pmt[0].begin(),crossed_pmt[0].end(),pmt)!=crossed_pmt[0].end();
    bool is_crossed_sec = find(crossed_pmt[1].begin(),crossed_pmt[1].end(),pmt)!=crossed_pmt[1].end();
    bool is_not_in_selected_ring = (find(rich_hit_on_ring_list.begin(),rich_hit_on_ring_list.end(),rich_hit)==rich_hit_on_ring_list.end());
    int status2 = ((is_bad_pmt)?0x1:0)+((is_bad_occupancy)?0x2:0)+((is_crossed_pri)?0x4:0)+((is_crossed_sec)?0x8:0)+((is_not_in_selected_ring)?0x10:0);
    hit_info hit(rich_hit->Status,status2,-1,rich_hit->Channel,rich_hit->Npe,0,0);
    hit_info_list.push_back(hit);
  }

  // store hits and global informations
  int nhit_stored = 0;
  for (int i1=0; i1<2; i1++) {
    for (int i2=0; i2<5; i2++) {
      ntpRich.tot_hit[i1][i2] = 0;
      ntpRich.tot_p[i1][i2] = 0;
    }
  }
  if (int(hit_info_list.size())!=(int)NRichHit()) std::cout << "Problem with hits categorization, Run:" << Run() << " Event:" << Event() << std::endl;
  for (vector<hit_info>::iterator it=hit_info_list.begin(); it!=hit_info_list.end(); it++) {
    int i1 = ((it->status2&0x10)!=0) ? 1 : 0; // 0: in-sel-ring, 1: out-of-sel-ring
    int i2 = 4;
    if      (((it->status2&0x03)!=0)||(((it->status>>28)&0x1)==0x0)) i2 = 0; // Bad (occupancy or bad PMT)
    else if ((it->status2&0x04)!=0)                                  i2 = 1; // Crossed Primary
    else if ((it->status2&0x08)!=0)                                  i2 = 2; // Crossed Secondary
    else if (((it->status2&0x0C)==0)&&(((it->status>>30)&0x1)==0x1)) i2 = 3; // Carlos-crossed only
    ntpRich.tot_hit[i1][i2]++;
    ntpRich.tot_p[i1][i2] += it->np;
    if (nhit_stored>29) continue;
    ntpRich.hit_stat[nhit_stored] = it->status;
    ntpRich.hit_stat2[nhit_stored] = it->status2;
    ntpRich.hit_used[nhit_stored] = it->used;
    ntpRich.hit_chan[nhit_stored] = it->channel;
    ntpRich.hit_np_uncorr[nhit_stored] = it->np;
    for (int ibeta=0; ibeta<2; ibeta++) ntpRich.hit_beta[nhit_stored][ibeta] = it->beta[ibeta];
    nhit_stored++;
  }

  /*
  // good hit subset
  cout << "Run:" << Run() << " Event:" << Event() << " ";
  for (int i1=0; i1<2; i1++)
    for (int i2=0; i2<5; i2++)
      std::cout << ntpRich.tot_hit[i1][i2] << " ";
  std::cout << std::endl;
  for (vector<hit_info>::iterator it=hit_info_list.begin(); it!=hit_info_list.end(); it++) {
    if (((it->status>>28)&0x1)==0x0) continue;
    if (((it->status>>30)&0x1)==0x1) continue;
    if (it->status2!=0) continue;
    if ( (it->beta[0]==0)||(it->beta[1]==0) ) continue;
    it->print();
  }
  */

  // clean-up
  hit_info_list.clear();
  rich_hit_on_ring_list.clear();
  rich_hit_beta_hyp_list.clear();
  rich_hit_offline_map.clear();

  // Veto, in case of no ring and downgoing track passing inside RICH radiator (depends on NtpTracker)
  if ( (!pRichRing)&&(pTrTrack)&&(pBetaH)&&(pBetaH->GetBeta()>0)&&(ntpTracker.IsInsideRich()>0) ) {
    ntpRich.veto_np_exp_elec = RichRingR::ComputeNpExp(pTrTrack,1,1);
    int prot_fit_id = pTrTrack->iTrTrackPar(algo_ref_fit,patt_ref_fit,21,m_p,1);
    if (prot_fit_id>=0) {
      double R = pTrTrack->GetRigidity(prot_fit_id);
      ntpRich.veto_beta_trk_prot = 1./sqrt(1.+pow(R/m_p,-2));
      ntpRich.veto_np_exp_prot = RichRingR::ComputeNpExp(pTrTrack,ntpRich.veto_beta_trk_prot,1);
    }
  }
}

void Analysis::FillNtpEcal() {
  // I ignore the associate ECAL shower.
  // We store the highest released energy showers (1st most energetic)
  int N = nEcalShower();
  int showerindex[N];
  float v[N];
  for(int ish=0; ish<nEcalShower(); ish++) {
    EcalShowerR &ecal = *(AMSEventR::pEcalShower(ish));
    v[ish] = ecal.EnergyE;
    showerindex[ish] = -1;
  }
  Tools::Mysort(N,v,showerindex);
  for(int i=0; i<min(N,2); i++) {
    int ish = showerindex[i];
    if (ish>-1) {
      EcalShowerR &ecal = *(AMSEventR::pEcalShower(ish));
      ntpEcal.energyE[i] = ecal.EnergyE;
      ntpEcal.energyD[i] = ecal.EnergyD/1000.;
      if (!disable) ntpEcal.bdt[i] = ecal.GetEcalBDT();
      for(int j=0; j<3; j++) {
        ntpEcal.entry[j][i] = ecal.Entry[j];
        ntpEcal.exit[j][i]  = ecal.Exit[j];
        ntpEcal.cog[j][i]   = ecal.CofG[j];
        ntpEcal.dir[j][i]   = ecal.Dir[j];
      }
      ntpEcal.chi2dir[i]   = ecal.Chi2Dir;
      ntpEcal.mips_tag[i]  = Tools::MipsTagCut(*this,ecal,true);
      ntpEcal.moliere[i]   = ecal.Energy3C[0];
      ntpEcal.rear_leak[i] = ecal.RearLeak;
      ntpEcal.tmax[i]      = ecal.ParProfile[1];
      ntpEcal.long_disp[i] = ecal.ShowerLongDisp;
      ntpEcal.nhit[i]      = Tools::MyEcalShowerHits(*this,ecal);
      ntpEcal.q[i]         = ecal.EcalChargeEstimator();
      for (int i2dcl=0; i2dcl<ecal.NEcal2DCluster(); i2dcl++) {
        Ecal2DClusterR *e2dcl = ecal.pEcal2DCluster(i2dcl);
        for (int icl=0; icl<e2dcl->NEcalCluster(); icl++) {
          EcalClusterR *ecl=e2dcl->pEcalCluster(icl);
          for (int ihit=0; ihit<ecl->NEcalHit(); ihit++) {
            EcalHitR *ehit = ecl->pEcalHit(ihit);
            int ilay = ehit->Plane;
            //if (ilay>3) continue;
            float edep = ehit->Edep;
            ntpEcal.nhit_lay[ilay][i]++;
            ntpEcal.edep_tot_lay[ilay][i] += edep;
            if (edep>ntpEcal.edep_max_lay[i][ilay]) ntpEcal.edep_max_lay[ilay][i] = edep;
          }
        }
      }
    }
  }
}

void Analysis::FillNtpAnti() {
  for (int ianti=0; ianti<nAntiCluster(); ianti++) {
    AntiClusterR* cluster = AMSEventR::pAntiCluster(ianti);
    int is = cluster->Sector - 1;
    ntpAnti.npair[is] = cluster->Npairs;
    ntpAnti.chisq[is] = cluster->chi;
    ntpAnti.unfz[is] = cluster->unfzeta;
    ntpAnti.edep[is] = cluster->Edep;
  }
}

double HitDistanceXY(TrRecHitR* hit, const AMSPoint &coo, int &mult) {
  mult = -1;
  if (hit->OnlyY()) return -1; // only available for XY hits
  int max_mult = 0;
  if (hit->GetXCluster()) max_mult = hit->GetXCluster()->GetMultiplicity();
  double dist = 1e+30;
  for (int imult=0; imult<max_mult; imult++) {
    double d = (coo - hit->GetCoord(imult)).norm();
    if (d<dist) {
      mult = imult;
      dist = d;
    }
  }
  if ( (mult>-1)&&(mult<=max_mult) ) return dist;
  return -2; // failed
}

void Analysis::FillNtpStandAlone() {
  // int   mult = -1;
  float dxdz[2] = {0,0};
  float dydz[2] = {0,0};
  if (pBetaH_SA) {
    ntpStandAlone.beta = isMC ? pBetaH_SA->GetMCBeta() : pBetaH_SA->GetBeta();
    ntpStandAlone.beta_err = ntpStandAlone.beta*ntpStandAlone.beta*pBetaH_SA->GetEBetaV();
    ntpStandAlone.beta_ncl = pBetaH_SA->GetUseHit();
    ntpStandAlone.beta_patt = 0;
    for (int il=0; il<4; il++) if (!(pBetaH_SA->TestExistHL(il))) ntpStandAlone.beta_patt |= (1<<il);
    ntpStandAlone.beta_chisqtn = pBetaH_SA->GetNormChi2T();
    AMSPoint tofh_coo;
    AMSDir tofh_dir;
    double tofh_time;
    // double path = pBetaH_SA->TInterpolate(0,tofh_coo,tofh_dir,tofh_time,false);
    if (ntpStandAlone.beta*tofh_dir[2]>0) for(int icoo=0; icoo<3; icoo++) tofh_dir[icoo] = -tofh_dir[icoo];
    float tofh_th = acos(tofh_dir[2]);
    float tofh_ph = atan2(tofh_dir[1],tofh_dir[0]);
    ntpStandAlone.theta = tofh_th;
    ntpStandAlone.phi = tofh_ph;
    for (int icoo=0; icoo<3; icoo++) ntpStandAlone.coo[icoo] = tofh_coo[icoo];
    ntpStandAlone.beta_q = pBetaH_SA->GetQ(ntpStandAlone.beta_nhit,ntpStandAlone.beta_q_err);
    for (int il=0; il<4; il++) ntpStandAlone.beta_q_lay[il] = pBetaH_SA->GetQL(il);
    for (int il=0; il<2; il++) {
      pBetaH_SA->TInterpolate(tracker_layers_z[(il==0)?0:8],tofh_coo,tofh_dir,tofh_time,false);
      dxdz[il] = (fabs(tofh_dir.z())>0) ? tofh_dir.x()/tofh_dir.z() : 0;
      dydz[il] = (fabs(tofh_dir.z())>0) ? tofh_dir.y()/tofh_dir.z() : 0;
      for (int icoo=0; icoo<3; icoo++) ntpStandAlone.exthit_int[il][icoo] = tofh_coo[icoo];
    }
  }
  // trd
  if (pTrdTrack_SA) {
    ntpStandAlone.trd_chisq = pTrdTrack_SA->Chi2;
    ntpStandAlone.trd_nseg = pTrdTrack_SA->NTrdSegment();
    ntpStandAlone.theta = pTrdTrack_SA->Theta;
    ntpStandAlone.phi = pTrdTrack_SA->Phi;
    for (int icoo=0; icoo<3; icoo++) ntpStandAlone.coo[icoo] = pTrdTrack_SA->Coo[icoo];
    ntpStandAlone.trd_q = pTrdTrack_SA->Q;
    AMSDir trd_dir(ntpStandAlone.theta,ntpStandAlone.phi);
    for (int il=0; il<2; il++) {
      dxdz[il] = (fabs(trd_dir.z())>0) ? trd_dir.x()/trd_dir.z() : 0;
      dydz[il] = (fabs(trd_dir.z())>0) ? trd_dir.y()/trd_dir.z() : 0;
      double z = tracker_layers_z[(il==0)?0:8];
      ntpStandAlone.exthit_int[il][0] = ntpStandAlone.coo[0] + dxdz[il]*(z-ntpStandAlone.coo[2]);
      ntpStandAlone.exthit_int[il][1] = ntpStandAlone.coo[1] + dydz[il]*(z-ntpStandAlone.coo[2]);
      ntpStandAlone.exthit_int[il][2] = z;
    }
  }
  // trdk
  TrdKCluster* trdK_SA = 0;
  if      (pTrdTrack_SA) trdK_SA = new TrdKCluster(this,pTrdTrack_SA);
  else if (pBetaH_SA)    trdK_SA = new TrdKCluster(this,pBetaH_SA);
  int    isvalid = -1;
  double lh[3] = {-1,-1,-1};
  double ll[3] = {-1,-1,-1};
  int    nh = 0;
  int    nhoff = 0;
  float  ampoff = 0;
  if (trdK_SA) {
    ntpStandAlone.trdk_is_align_ok = (trdK_SA->IsReadAlignmentOK!=0);
    ntpStandAlone.trdk_is_calib_ok = (trdK_SA->IsReadCalibOK!=0);
    for (int iopt=0; iopt<3; iopt++) {
      trdK_SA->CalculateTRDCharge(iopt);
      ntpStandAlone.trdk_q[iopt]     = trdK_SA->GetTRDCharge();
      ntpStandAlone.trdk_q_err[iopt] = trdK_SA->GetTRDChargeError();
    }
    float threshold = 15;
    float total_pathlength;
    float total_amp;
    isvalid = trdK_SA->GetLikelihoodRatio_TRDRefit(threshold,lh,ll,nh,total_pathlength,total_amp,1,1,-1,0);
    trdK_SA->GetOffTrackHit_TRDRefit(nhoff,ampoff);
    if (trdK_SA) delete trdK_SA;
  }
  ntpStandAlone.trdk_like_valid = isvalid;
  ntpStandAlone.trdk_like_e  = (ll[0]>0) ? -log(ll[0]) : -1;
  ntpStandAlone.trdk_like_p  = (ll[1]>0) ? -log(ll[1]) : -1;
  ntpStandAlone.trdk_like_he = (ll[2]>0) ? -log(ll[2]) : -1;
  if (isvalid>=0) {
    ntpStandAlone.trdk_like_nhit     = nh;
    ntpStandAlone.trdk_like_nhit_off = nhoff;
    ntpStandAlone.trdk_like_ampl_off = ampoff;
  }
  // external layers unbiased hit
  float dist_min[2] = {1e+30,1e+30};
  float charge_max[2] = {0,0};
  TrRecHitR* exthit_closest[2] = {0};
  TrRecHitR* exthit_largest[2] = {0};
  int mult_closest[2] = {-1,-1};
  int mult_largest[2] = {-1,-1};
  double beta = (pBetaH_SA) ? (isMC ? pBetaH_SA->GetMCBeta() : pBetaH_SA->GetBeta()):1;
  for (int ihit=0; ihit<nTrRecHit(); ihit++) {
    TrRecHitR* hit = pTrRecHit(ihit);
    if (!hit) continue;
    if (hit->OnlyY()) continue;
    int layer = hit->GetLayerJ();
    if ( (layer!=1)&&(layer!=9) ) continue;
    int il = (layer==1)?0:1;
    // find max charge external (let's make it independent from TRD)
    int mult = 0;
    double charge = hit->GetQ(2,beta,0,0,mult,dxdz[il],dydz[il]);
    if (charge>charge_max[il]) {
      charge_max[il] = charge;
      exthit_largest[il] = hit;
      mult_largest[il] = mult;
    }
    // find closest to TRD track
    if (!pTrdTrack_SA) continue;
    AMSPoint point(ntpStandAlone.exthit_int[il][0],ntpStandAlone.exthit_int[il][1],ntpStandAlone.exthit_int[il][2]);
    double dist = HitDistanceXY(hit,point,mult);
    if ( (dist<0)||(mult<0) ) continue; // failed
    if (dist<dist_min[il]) {
      dist_min[il] = dist;
      exthit_closest[il] = hit;
      mult_closest[il] = mult;
    }
  }
  for (int il=0; il<2; il++) {
    if (exthit_closest[il]) {
      AMSPoint point = exthit_closest[il]->GetCoord();
      for (int icoo=0; icoo<3; icoo++) ntpStandAlone.exthit_closest_coo[il][icoo] = point[icoo];
      ntpStandAlone.exthit_closest_q[0][il] = exthit_closest[il]->GetQ  (2,beta,0,0,mult_closest[il],dxdz[il],dydz[il]);
      if (!disable) ntpStandAlone.exthit_closest_q[1][il] = exthit_closest[il]->GetQYJ(2,beta,0,  mult_closest[il],dxdz[il],dydz[il]);
      ntpStandAlone.exthit_closest_q[2][il] = exthit_closest[il]->GetQH (2,beta,0,  mult_closest[il],dxdz[il],dydz[il]);
      ntpStandAlone.exthit_closest_status[il] = exthit_closest[il]->GetQStatus();
    }
    if (exthit_largest[il]) {
      AMSPoint point = exthit_largest[il]->GetCoord();
      for (int icoo=0; icoo<3; icoo++) ntpStandAlone.exthit_largest_coo[il][icoo] = point[icoo];
      ntpStandAlone.exthit_largest_q[0][il] = exthit_largest[il]->GetQ  (2,beta,0,0,mult_largest[il],dxdz[il],dydz[il]);
      if (!disable) ntpStandAlone.exthit_largest_q[1][il] = exthit_largest[il]->GetQYJ(2,beta,0,  mult_largest[il],dxdz[il],dydz[il]);
      ntpStandAlone.exthit_largest_q[2][il] = exthit_largest[il]->GetQH (2,beta,0,  mult_largest[il],dxdz[il],dydz[il]);
      ntpStandAlone.exthit_largest_status[il] = exthit_largest[il]->GetQStatus();
    }
  }
}

void Analysis::FillNtpCompact(bool use_already_calculated) {
  memset(&ntpCompact,0,sizeof(NtpCompact));
  // header
  int npart = nParticle();    if (npart>9) npart = 9;
  int nanti = nAntiCluster(); if (nanti>9) nanti = 9;
  int ntofh = nBetaH();       if (ntofh>9) ntofh = 9;
  int ntrck = nTrTrack();     if (ntrck>9) ntrck = 9;
  int ntrkh = nTrRecHit();    if (ntrkh>99) ntrkh = 99;
  int ntrdc = nTrdCluster();  if (ntrdc>99) ntrdc = 99;
  int ntofc = nTofClusterH(); if (ntofc>41) ntofc = 40;
  ntpCompact.status = npart+nanti*10+ntofh*100+ntrck*1000+ntrkh*10000+ntrdc*1000000+ntofc*100000000;
  if (nLevel1()>0) {
    Level1R* lvl1 = pLevel1(0); // anyone
    int antipatt = lvl1->AntiPatt;
    int sublvl1  = lvl1->PhysBPatt;
    int trigpatt = lvl1->JMembPatt;
    if (nMCEventg()>0) lvl1->RebuildTrigPatt(trigpatt,sublvl1,antipatt);
    ntpCompact.sublvl1  = (short)sublvl1;
    ntpCompact.trigpatt = (short)trigpatt;
  }
  if (pBetaH) {
    // beta
    ntpCompact.tof_beta = isMC ? pBetaH->GetMCBeta() : pBetaH->GetBeta();
    ntpCompact.tof_beta_ncl = pBetaH->GetUseHit();
    ntpCompact.tof_chisqtn = pBetaH->GetNormChi2T();
    ntpCompact.tof_trk_ncl = pBetaH->GetSumHit();
    ntpCompact.tof_chisqcn = pBetaH->GetNormChi2C();
    ntpCompact.tof_beta_patt = 0;
    for (int i=0; i<4; i++) if (!(pBetaH->TestExistHL(i))) ntpCompact.tof_beta_patt |= (1<<i);
    if (!pBetaH->IsGoodBeta()) ntpCompact.tof_beta_patt |= 0x10;
    if (!pBetaH->IsTkTofMatch()) ntpCompact.tof_beta_patt |= 0x20;
    // charge 
    for (int il=0; il<4; il++) ntpCompact.tof_q_lay[il] = pBetaH->GetQL(il);
    int nlay;
    ntpCompact.tof_z = pBetaH->GetZ(nlay,ntpCompact.tof_z_like);
    ntpCompact.tof_z_like = (ntpCompact.tof_z_like>0) ? -log(ntpCompact.tof_z_like) : -1;
    ntpCompact.tof_z_nhit = nlay;
    // isolation 
    float eclp[4][3] = {{0}};
    TofClusterHR *ptofhclbetap[4];
    for (int i=0; i<4; i++) ptofhclbetap[i] = pBetaH->GetClusterHL(i);
    for (int i=0; i<nTofClusterH(); i++) {
      TofClusterHR* ptofhcl = pTofClusterH(i);
      if (!ptofhcl) continue;
      int il = ptofhcl->Layer;
      int ib = ptofhcl->Bar;
      float edep = ptofhcl->GetEdep();
      int iw = 2; // default, 2: Other bars in the plane
      if (ptofhclbetap[il]) {
        if (ptofhcl==ptofhclbetap[il]) iw = 0; // 0: Used cluster
        else if (abs(ib-ptofhclbetap[il]->Bar)==1) iw = 1; // 1: Nearby bar
      }
      eclp[il][iw] += edep;
    }
    for (int il=0; il<4; il++) {
      ntpCompact.tof_edep_frac[il] = 0;
      if (eclp[il][0]+eclp[il][1]+eclp[il][2]<=0) continue;
      ntpCompact.tof_edep_frac[il] = eclp[il][0]/(eclp[il][0]+eclp[il][1]+eclp[il][2]);
    }
    // unused bars
    float clsdt[4];
    float clsed[4];
    Tools::TofUnusedHits(this,pBetaH,clsdt,clsed,ntpCompact.tof_clsn);
  }
  if (pBeta) {
    ntpCompact.tof_evgeni_beta = pBeta->Beta;
  }
  if (pTrTrack) {
    float beta = (pBetaH) ? (isMC ? pBetaH->GetMCBeta() : pBetaH->GetBeta()) : 1;
    ntpCompact.trk_patty  = pTrTrack->GetBitPatternJ();
    ntpCompact.trk_pattxy = pTrTrack->GetBitPatternXYJ();
    int fit_id = pTrTrack->iTrTrackPar(algo_ref_fit,patt_ref_fit,refi_ref_fit,mass_ref_fit,char_ref_fit);
    if (fit_id>=0) {
      int innerq_patt;
      ntpCompact.trk_q_inn = pTrTrack->GetInnerQYJ(ntpCompact.trk_q_inn_rms,innerq_patt,2,beta,fit_id); 
      for (int il = 0; il < 9; il++) {
        if (pTrTrack->GetHitLJ(il + 1))
          ntpCompact.trk_q_lay[il] = ((((pTrTrack->GetHitLJ(il+1)->GetQStatus()) & 0x10013D) == 0) ? 1 : -1) *
                                     pTrTrack->GetLayerQYJ(il+1,beta,fit_id);
      }
      // ratio
      double tredep[9][2][2] = {{{0}}};
      for (int il=0; il<9; il++) {
        TrRecHitR* hit = pTrTrack->GetHitLJ(il+1);
        TrClusterR* cluster[2] = {(!hit)?0:hit->GetXCluster(),(!hit)?0:hit->GetYCluster()};
        AMSPoint pos;
        AMSDir dir;
        pTrTrack->InterpolateLayerJ(il+1,pos,dir,fit_id); 
        for (int is=0; is<2; is++) {
          for (int iw=0; iw<2; iw++) {
            tredep[il][is][iw] = 0;
            if (iw==0) {
              if (!cluster[is]) continue;
              tredep[il][is][iw] += cluster[is]->GetEdep();
              continue;
            }
            double min_dis = 10;
            for (int icl=0; icl<nTrCluster(); icl++) {
              TrClusterR* this_cluster = pTrCluster(icl);
              if ( (is!=this_cluster->GetSide())||(il+1!=this_cluster->GetLayerJ()) ) continue;
              bool near = false;
              for (int im=0; im<this_cluster->GetMultiplicity(); im++) {
                float mul_dis = fabs(this_cluster->GetGCoord(im)-pos[is]);
               l_chisqn near |= (mul_dis<min_dis);
              }
              if (near) tredep[il][is][iw] += this_cluster->GetEdep();
            }
          }
        }
      }
      for (int il=0; il<9; il++) 
        for (int is=0; is<2; is++) 
          ntpCompact.trk_edep_frac[il][is] = (tredep[il][is][1]>0) ? tredep[il][is][0]/tredep[il][is][1] : 0; 
    }
    // Vitaly refits
    int algo[5] = {1,1,1,1,21}; // (fit 3 is important for compact definition)
    int patt[5] = {0,5,6,3,3};
    for (int irefit=0; irefit<5; irefit++) {
      int fit_id = pTrTrack->iTrTrackPar(algo[irefit],patt[irefit],refi_ref_fit,mass_ref_fit,char_ref_fit);
      if (fit_id<0) continue;
      ntpCompact.trk_fiducial[irefit] = GetPatternInsideTracker(pTrTrack,fit_id);
      const TrTrackPar &t = pTrTrack->gTrTrackPar(fit_id);
      if (!isMC) ntpCompact.trk_rig[irefit] = pTrTrack->GetCorrectedRigidity(fit_id,3); 
      else       ntpCompact.trk_rig[irefit] = pTrTrack->GetRigidity(fit_id);
      ntpCompact.trk_chisqn[irefit][0] = (t.NdofX>0)? t.ChisqX/t.NdofX:0;
      ntpCompact.trk_chisqn[irefit][1] = (t.NdofY>0)? t.ChisqY/t.NdofY:0;
      AMSPoint point;
      AMSDir   dir;
      pTrTrack->Interpolate(top_z,point,dir,fit_id);
      AMSDir dirc(M_PI-dir.gettheta(),M_PI+dir.getphi());
      double rcut = 0;
      if (!GetStoermerCutoff(rcut,1,dirc)) ntpCompact.trk_stoermer[irefit] = rcut;
      if (irefit==0) {
        AMSPoint point;
        AMSDir   dir;
        pTrTrack->Interpolate(rich_radiator_z,point,dir,fit_id);
        ntpCompact.trk_int_rich_rad[0] = point[0];
        ntpCompact.trk_int_rich_rad[1] = point[1];
      }
    }
    // ancillary fits
    fit_id = pTrTrack->iTrTrackPar(21,3,21,mass_ref_fit,char_ref_fit);
    float rigi_inner = 0;
    if (fit_id>=0) {
      if (!isMC) rigi_inner = pTrTrack->GetCorrectedRigidity(fit_id,3);
      else       rigi_inner = pTrTrack->GetRigidity(fit_id);
      fit_id = pTrTrack->iTrTrackPar(21,1,21,mass_ref_fit,char_ref_fit);
      float rigi_intop = 0;
      if (fit_id>=0) {
        if (!isMC) rigi_intop = pTrTrack->GetCorrectedRigidity(fit_id,3);
        else       rigi_intop = pTrTrack->GetRigidity(fit_id);
      }
      fit_id = pTrTrack->iTrTrackPar(21,2,21,mass_ref_fit,char_ref_fit);
      float rigi_indwn = 0;
      if (fit_id>=0) {
        if (!isMC) rigi_indwn = pTrTrack->GetCorrectedRigidity(fit_id,3);
        else       rigi_indwn = pTrTrack->GetRigidity(fit_id);
      }
      ntpCompact.trk_rho_up = rigi_intop/rigi_inner;
      ntpCompact.trk_rho_dw = rigi_indwn/rigi_inner;
    }
    // Kalman refit
    fit_id = pTrTrack->iTrTrackPar(6,patt_ref_fit,refi_ref_fit,mass_ref_fit,char_ref_fit);
    if (fit_id>=0) {
      ntpCompact.trk_kal_fiducial = GetPatternInsideTracker(pTrTrack,fit_id);
      const TrTrackPar &t = pTrTrack->gTrTrackPar(fit_id);
      for (int iz=0; iz<3; iz++) {
        if (!isMC) ntpCompact.trk_kal_rig[iz] = pTrTrack->GetCorrectedRigidity(fit_id,3,iz);
        else       ntpCompact.trk_kal_rig[iz] = pTrTrack->GetRigidity(fit_id,iz);
      }
      ntpCompact.trk_kal_chisqn[0] = (t.NdofX>0)? t.ChisqX/t.NdofX:0;
      ntpCompact.trk_kal_chisqn[1] = (t.NdofY>0)? t.ChisqY/t.NdofY:0;
      AMSPoint point;
      AMSDir   dir;
      pTrTrack->Interpolate(top_z,point,dir,fit_id);
      AMSDir dirc(M_PI-dir.gettheta(),M_PI+dir.getphi());
      double rcut = 0;
      if (!GetStoermerCutoff(rcut,1,dirc)) ntpCompact.trk_kal_stoermer = rcut;
    }
  }
  ntpCompact.rich_status = (nRichHit()<100) ? nRichHit() : 99;
  ntpCompact.rich_tot_np = (use_already_calculated) ? ntpRich.tot_p_uncorr : RichHitR::getCollectedPhotoElectrons();
  if (pRichRing) {
    ntpCompact.rich_status += (pRichRing->getPMTs()<100) ? pRichRing->getPMTs()*100 : 100*99;
    ntpCompact.rich_status += 100*100*pRichRing->getUsedHits();
    const float* TrackEmissionPoint = pRichRing->getTrackEmissionPoint();
    MCEventgR* mcpart = GetPrimaryMC();
    int charge_mc = 1;
    if (mcpart) charge_mc = mcpart->Charge;
    ntpCompact.rich_beta = (use_already_calculated) ? ntpRich.beta_corrected : pRichCorr->correctBeta(
      pRichRing->Beta,
      TrackEmissionPoint[0],TrackEmissionPoint[1],
      pRichRing->getTrackTheta(),pRichRing->getTrackPhi(),
      abs(charge_mc),(nMCEventg()>0)
    );
    ntpCompact.rich_np     = (use_already_calculated) ? ntpRich.np_uncorr     : pRichRing->getPhotoElectrons(false);
    ntpCompact.rich_np_exp = (use_already_calculated) ? ntpRich.np_exp_uncorr : pRichRing->getExpectedPhotoElectrons(false);
    ntpCompact.rich_prob   = (use_already_calculated) ? ntpRich.prob          : pRichRing->getProb();
    ntpCompact.rich_select = (use_already_calculated) ? ntpRich.selection     : Tools::RichQC(pRichRing);
    RichBDT->data->FillDataDirect(this);
    ntpCompact.rich_bdt = pRichRing->IsNaF() ? RichBDT->GetNafBDT() : RichBDT->GetAglBDT(); // GetRichBDT();
  }
  // MC Info
  if (nMCEventg()!=0) {
    MCEventgR* mcpart = GetPrimaryMC();
    ntpCompact.mc_momentum = mcpart->Momentum;
  }
}

void Analysis::FillNtpCompactStandAlone(bool use_already_calculated) {
  // BetaH
  if (pBetaH_SA) {
    ntpCompact.sa_tof_beta = (isMC) ? pBetaH_SA->GetMCBeta() : pBetaH_SA->GetBeta();
    ntpCompact.sa_tof_chisqtn = pBetaH_SA->GetNormChi2T();
    for (int il=0; il<4; il++) ntpCompact.sa_tof_q_lay[il] = pBetaH_SA->GetQL(il);
    ntpCompact.sa_tof_beta_ncl = pBetaH_SA->GetUseHit();
    ntpCompact.sa_tof_build = pBetaH_SA->GetBuildType();
    // unused bars
    float clsdt[4];
    float clsed[4];
    Tools::TofUnusedHits(this,pBetaH_SA,clsdt,clsed,ntpCompact.sa_tof_clsn);
    // isolation 
    float eclp[4][3] = {{0}};
    TofClusterHR *ptofhclbetap[4];
    for (int i=0; i<4; i++) ptofhclbetap[i] = pBetaH_SA->GetClusterHL(i);
    for (int i=0; i<nTofClusterH(); i++) {
      TofClusterHR* ptofhcl = pTofClusterH(i);
      if (!ptofhcl) continue;
      int il = ptofhcl->Layer;
      int ib = ptofhcl->Bar;
      float edep = ptofhcl->GetEdep();
      int iw = 2; // default, 2: Other bars in the plane
      if (ptofhclbetap[il]) {
        if (ptofhcl==ptofhclbetap[il]) iw = 0; // 0: Used cluster
        else if (abs(ib-ptofhclbetap[il]->Bar)==1) iw = 1; // 1: Nearby bar
      }
      eclp[il][iw] += edep;
    }
    for (int il=0; il<4; il++) {
      ntpCompact.sa_tof_edep_frac[il] = 0;
      if (eclp[il][0]+eclp[il][1]+eclp[il][2]<=0) continue;
      ntpCompact.sa_tof_edep_frac[il] = eclp[il][0]/(eclp[il][0]+eclp[il][1]+eclp[il][2]);
    }
  }
  // TRD and unbiased hit on L1
  ntpCompact.sa_trd_same = false;
  if (pTrdTrack_SA) {
    if (pTrdTrack_SA==pTrdTrack) ntpCompact.sa_trd_same = true;
    ntpCompact.sa_trd_fiducial = GetPatternInsideTracker(pTrdTrack_SA);
    ntpCompact.sa_trd_q = pTrdTrack_SA->Q;
    ntpCompact.sa_trd_chi2 = pTrdTrack_SA->Chi2;
    AMSDir trd_dir(pTrdTrack_SA->Theta,pTrdTrack_SA->Phi);
    double rcut = 0;
    if (!GetStoermerCutoff(rcut,1,trd_dir)) ntpCompact.sa_trd_stoermer = rcut;
    float dxdz = trd_dir.x()/trd_dir.z();
    float dydz = trd_dir.y()/trd_dir.z();
    AMSPoint trd_point(
      pTrdTrack_SA->Coo[0] + dxdz*(tracker_layers_z[0]-pTrdTrack_SA->Coo[2]),
      pTrdTrack_SA->Coo[1] + dydz*(tracker_layers_z[0]-pTrdTrack_SA->Coo[2]),
      tracker_layers_z[0]
    );
    float dist_min = 1e+30;
    TrRecHitR* exthit = 0;
    int mult_min = -1;
    for (int ihit=0; ihit<nTrRecHit(); ihit++) {
      TrRecHitR* hit = pTrRecHit(ihit);
      if (!hit) continue;
      if (hit->OnlyY()) continue;
      if (hit->GetLayerJ()!=1) continue;
      int mult = -1;
      double dist = HitDistanceXY(hit,trd_point,mult);
      if ( (dist<0)||(mult<0) ) continue; // failed
      if (dist<dist_min) {
        dist_min = dist;
        exthit = hit;
        mult_min = mult;
      }
    }
    if (exthit) {
      AMSPoint point = exthit->GetCoord();
      float beta = (pBetaH_SA) ? (isMC ? pBetaH_SA->GetMCBeta() : pBetaH_SA->GetBeta()) : 1;
      ntpCompact.sa_exthit_ql1 = exthit->GetQYJ(2,beta,0,mult_min,dxdz,dydz);
      ntpCompact.sa_exthit_dl1[0] = point.x()-trd_point.x();
      ntpCompact.sa_exthit_dl1[1] = point.y()-trd_point.y();
      ntpCompact.sa_exthit_status_l1 = exthit->GetQStatus();
    }
  }
  // ECAL
  for (int ish=0; ish<nEcalShower(); ish++) {
    EcalShowerR* ecal = AMSEventR::pEcalShower(ish);
    if (ecal->EnergyD/1000.<ntpCompact.sa_ecal_edepd) continue;
    ntpCompact.sa_ecal_edepd = ecal->EnergyD/1000.;
  }
}

int Analysis::Loc2Gl(double Theta, double Phi, double &ThetaGl, double &PhiGl){
  // Define transformation matrices
  ////////////
  // AMS->ISS
  double alpha = 12 * TMath::DegToRad();
  double mAMS2ISS[3][3];
  mAMS2ISS[0][0] = 0;
  mAMS2ISS[0][1] = -1;
  mAMS2ISS[0][2] = 0;
  mAMS2ISS[1][0] = -cos(alpha);
  mAMS2ISS[1][1] = 0;
  mAMS2ISS[1][2] = -sin(alpha);
  mAMS2ISS[2][0] = sin(alpha);
  mAMS2ISS[2][1] = 0;
  mAMS2ISS[2][2] = -cos(alpha);
  AMSRotMat TrAMS2ISS(mAMS2ISS);
  ////////////
  // ISS->LVLH
  float pitch,roll,yaw;
  yaw   = fHeader.Yaw;
  pitch = fHeader.Pitch;
  roll  = fHeader.Roll;
  double ca = cos(yaw);
  double sa = sin(yaw);
  double cb = cos(pitch);
  double sb = sin(pitch);
  double cg = cos(roll);
  double sg = sin(roll);
  double mISS2LVLH[3][3];
  mISS2LVLH[0][0] = ca*cb;
  mISS2LVLH[0][1] = ca*sb*sg-sa*cg;
  mISS2LVLH[0][2] = ca*sb*cg+sa*sg;
  mISS2LVLH[1][0] = sa*cb;
  mISS2LVLH[1][1] = sa*sb*sg+ca*cg;
  mISS2LVLH[1][2] = sa*sb*cg-ca*sg;
  mISS2LVLH[2][0] = -sb;
  mISS2LVLH[2][1] = cb*sg;
  mISS2LVLH[2][2] = cb*cg;
  AMSRotMat TrISS2LVLH(mISS2LVLH);
  ////////////
  // LVLH->GEO
  // float theta,phi,v,vtheta,vphi;
  float theta,phi,vtheta,vphi;
  double pi=3.1415926;
  theta  = fHeader.ThetaS;
  phi    = fHeader.PhiS;
  vtheta = fHeader.VelTheta;
  vphi   = fHeader.VelPhi;
  AMSDir amszg(pi/2+theta,phi+pi);  // z points toward the Earth's center
  AMSDir amsxg(pi/2-vtheta,vphi);   // x points along the orbital velocity
  AMSDir amsyg=amszg.cross(amsxg);
  // double prod=amsxxg.prod(amszg);
  double mLVLH2GEO[3][3];
  for(int i=0 ; i < 3 ; i++) {
    mLVLH2GEO[i][0] = amsxg[i];
    mLVLH2GEO[i][1] = amsyg[i];
    mLVLH2GEO[i][2] = amszg[i];
  }
  AMSRotMat TrLVLH2GEO(mLVLH2GEO);
  ////////////
  AMSDir _dir(Theta,Phi);
  AMSDir global = TrLVLH2GEO * (TrISS2LVLH * (TrAMS2ISS * _dir));
  ThetaGl=global.gettheta();
  PhiGl=global.getphi();
  return 0;
}
