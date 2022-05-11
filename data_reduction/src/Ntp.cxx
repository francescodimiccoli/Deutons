#include "Ntp.h"
#include "RichBDT.h"
#include "TrackerBDT.h"
#include "Classifier.h"
#include "RichOccupancy.h"

#include "TMath.h"

ClassImp(RTIInfo);
ClassImp(FileInfo);
ClassImp(FileMCInfo);
ClassImp(ProcInfo);
ClassImp(NtpSHeader);
ClassImp(NtpHeader);
ClassImp(NtpMCHeader);
ClassImp(NtpTrd);
ClassImp(NtpTof);
ClassImp(NtpTracker);
ClassImp(NtpRich);
ClassImp(NtpEcal);
ClassImp(NtpAnti);
ClassImp(NtpStandAlone);
ClassImp(Event);

bool RTIInfo::Select() {
  bool cut0 = (ntrig/nev>0.98);
  bool cut1 = (npart/ntrig>0.07/1600*ntrig)&&(npart/ntrig<0.25);
  bool cut2 = (lf>0.5);
  bool cut3 = (zenith<25);
  bool cut4 = (nerr>=0)&&(nerr/nev<0.1);
  bool cut5 = (npart>0)&&(nev<1800);
  bool cut6 = (fabs(dl1l9[0][1])<35)&&(fabs(dl1l9[1][1])<45);
  bool cut7 = (!isinsaa);
  bool tcut = (cut0&&cut1&&cut2&&cut3&&cut4&&cut5&&cut6&&cut7);
  return tcut;
}

bool RTIInfo::IsBadRun() {
  if ( (run==1306219312)||
       (run==1306219522)||
       (run==1306233745)||
       ( (run>=1307125541)&&(run<=1307218054) )||
       (run==1321198167) ) return true;
  return false;
}

int FileMCInfo::GetNGen() {
  int ngen1 = event[1];
  if (event[0]>1000) ngen1 -= event[0];
  int ngen2 = ngen_datacard;
  int ngen = 0;
  // if file and datacard agree at percent, use datacard
  if (fabs(ngen2-ngen1)<0.01*ngen2) ngen = ngen2;
  // otherwise use last-first
  else ngen = ngen1;
  // if number of generated is too low return 0
  if (ngen<1000) return 0;
  return ngen;
}

double FileMCInfo::GetRMin() {
  if (abs(charge)<=0) return 0;
  // by default take it from datacard
  double rmin = fabs(mom_datacard[0]/charge);
  // otherwise from file name
  if (rmin<=0) rmin = fabs(mom_filename[0]/charge);
  return rmin;
}

double FileMCInfo::GetRMax() {
  if (abs(charge)<=0) return 0;
  // by default take it from datacard
  double rmax = fabs(mom_datacard[1]/charge);
  // otherwise from file name
  if (rmax<=0) rmax = fabs(mom_filename[1]/charge);
  return rmax;
}

double NtpTof::GetQ(int &n, double& rms, int pattern, bool good_flag, bool uncorr) {
  n = 0;
  rms = 0;
  double mea = 0;
  for (int i=0; i<4; i++) {
    if ((pattern&(1<<i)) == 0) continue;
    if (q_lay[i]<=0) continue;
    if ( (good_flag)&&(flagp[i]!=0) ) continue;
    mea += (uncorr)?q_lay_uncorr[i]:q_lay[i];
    rms += pow((uncorr)?q_lay_uncorr[i]:q_lay[i],2);
    n++;
  }
  if (n<=0) return 0;
  mea /= n;
  rms /= n;
  rms = sqrt(rms-mea*mea);
  return mea;
}

void NtpTrd::Interpolate(float z, float xy[2]) {
  xy[0] = coo[0] + tan(theta)*cos(phi)*(z-coo[2]);
  xy[1] = coo[1] + tan(theta)*sin(phi)*(z-coo[2]);
}

static float cut_aerogel_external_border = 3500.;   // Aerogel external border (r**2)
static float cut_aerogel_naf_border[2] = {17.,19.}; // NaF/Aerogel border

int NtpTracker::IsInsideRich() {
  float x = int_rich_rad[0];
  float y = int_rich_rad[1];
  if ((x*x+y*y)>cut_aerogel_external_border) return 0;
  if      (max(fabs(x),fabs(y))<cut_aerogel_naf_border[0]) return 1;
  else if (max(fabs(x),fabs(y))>cut_aerogel_naf_border[1]) return 2;
  return 0;
}

static float rich_radiator_z = -71.87;

int NtpTrd::IsInsideRich() {
  float coo[2] = {0};
  Interpolate(rich_radiator_z,coo);
  float x = coo[0];
  float y = coo[1];
  if ((x*x+y*y)>cut_aerogel_external_border) return 0;
  if      (max(fabs(x),fabs(y))<cut_aerogel_naf_border[0]) return 1;
  else if (max(fabs(x),fabs(y))>cut_aerogel_naf_border[1]) return 2;
  return 0;
}

static float tracker_planes_edges[9][4] = { // tracker edges xmin,ymin,xmax,ymax (xmax is used also as disk radius for L1, ..., L8)
  {-62.14,  -47.40,   62.14,   47.40},
  {-62.14,  -40.10,   62.14,   40.10},
  {-49.70,  -43.75,   49.70,   43.75},
  {-49.72,  -43.75,   49.72,   43.75},
  {-49.71,  -36.45,   49.70,   36.45},
  {-49.72,  -36.45,   49.72,   36.45},
  {-49.72,  -43.75,   49.71,   43.75},
  {-49.72,  -43.75,   49.71,   43.75},
  {-45.62,  -29.48,   45.55,   29.53}
};

int NtpTracker::GetPatternInsideTracker() {
  int pattern = 0;
  for (int ilayer=0; ilayer<9; ilayer++) {
    float x = int_lay[ilayer][0];
    float y = int_lay[ilayer][1];
    bool isinlayer = false;
    if ( (x>tracker_planes_edges[ilayer][0])&&(x<tracker_planes_edges[ilayer][2])&&
         (y>tracker_planes_edges[ilayer][1])&&(y<tracker_planes_edges[ilayer][3]) ) {
      if ((ilayer+1)==9) isinlayer = true;
      else {
        if ( (sqrt(x*x+y*y)<tracker_planes_edges[ilayer][2]) ) isinlayer = true;
      }
    }
    if (isinlayer) pattern |= (1<<ilayer);
  }
  return pattern;
}

static float tracker_layers_z[9] = {158.920,53.060,29.228,25.212,1.698,-2.318,-25.212,-29.228,-135.882};

int NtpTrd::GetPatternInsideTracker() {
  int pattern = 0;
  for (int ilayer=0; ilayer<9; ilayer++) {
    float coo[2];
    Interpolate(tracker_layers_z[ilayer],coo);
    float x = coo[0];
    float y = coo[1];
    bool isinlayer = false;
    if ( (x>tracker_planes_edges[ilayer][0])&&(x<tracker_planes_edges[ilayer][2])&&
         (y>tracker_planes_edges[ilayer][1])&&(y<tracker_planes_edges[ilayer][3]) ) {
      if ((ilayer+1)==9) isinlayer = true;
      else {
        if ( (sqrt(x*x+y*y)<tracker_planes_edges[ilayer][2]) ) isinlayer = true;
      }
    }
    if (isinlayer) pattern |= (1<<ilayer);
  }
  return pattern;
}

static int   ntrd_center = 9;
static float trd_accep_centerx[9] = {-80.0, -47.0, 47.0, 80.0,  80.0,  47.0, -47.0, -80.0, -80.0};
static float trd_accep_centery[9] = { 43.5,  75.5, 75.5, 43.5, -43.5, -75.5, -75.5, -43.5,  43.5};
static int   ntrd_top = 37;
static float trd_accep_topx[37] = {
  -99.0,-89.0,-89.0,-78.7,-78.7,-67.8,-67.8,-57.7,-57.7, 57.7, 57.7, 67.8, 67.8, 78.7, 78.7, 89.0, 89.0, 99.0,
   99.0, 89.0, 89.0, 78.7, 78.7, 67.8, 67.8, 57.7, 57.7,-57.7,-57.7,-67.8,-67.8,-78.7,-78.7,-89.0,-89.0,-99.0,-99.0
};
static float trd_accep_centerz = 100;
static float trd_accep_topy[37] = {
   54.5, 54.5, 62.5, 62.5, 74.0, 74.0, 84.0, 84.0, 95.3, 95.3, 84.0, 84.0, 74.0, 74.0, 62.5, 62.5, 54.5, 54.5,
  -51.7,-51.7,-62.2,-62.2,-72.0,-72.0,-82.5,-82.5,-92.5,-92.5,-82.5,-82.5,-72.0,-72.0,-62.2,-62.2,-51.7,-51.7, 54.5
};
static float trd_accep_topz = 143.925;

bool NtpTracker::IsInsideTRD() {
  double x0 = int_lay[0][0]; // l1
  double y0 = int_lay[0][1]; // l1
  double z0 = int_lay[0][2]; // l1
  double dxdz = (int_lay[1][0]-int_lay[0][0])/(int_lay[1][2]-int_lay[0][2]); // l1-l2
  double dydz = (int_lay[1][1]-int_lay[0][1])/(int_lay[1][2]-int_lay[0][2]); // l1-l2
  float trd_coo[2][2] = {
    {static_cast<float>( x0 + dxdz*(trd_accep_centerz-z0) ),
     static_cast<float>( y0 + dydz*(trd_accep_centerz-z0) )},
    {static_cast<float>( x0 + dxdz*(trd_accep_topz-z0) ),
     static_cast<float>( y0 + dydz*(trd_accep_topz-z0) )}
  };
  bool passTrdCenter = TMath::IsInside(trd_coo[0][0],trd_coo[0][1],ntrd_center,trd_accep_centerx,trd_accep_centery);
  bool passTrdTop    = TMath::IsInside(trd_coo[1][0],trd_coo[1][1],ntrd_top,   trd_accep_topx,   trd_accep_topy);
  return passTrdCenter&&passTrdTop;
}

bool NtpTrd::IsInsideTRD() {
  float trd_coo[2][2] = {{0}};
  Interpolate(trd_accep_centerz,trd_coo[0]);
  Interpolate(trd_accep_topz,trd_coo[1]);
  bool passTrdCenter = TMath::IsInside(trd_coo[0][0],trd_coo[0][1],ntrd_center,trd_accep_centerx,trd_accep_centery);
  bool passTrdTop    = TMath::IsInside(trd_coo[1][0],trd_coo[1][1],ntrd_top,   trd_accep_topx,   trd_accep_topy);
  return passTrdCenter&&passTrdTop;
}

static float ecal_topz = -142.732;
static float ecal_botz = -159.382;
static float ecal_xmin =  -32.270;
static float ecal_xmax =   32.530;
static float ecal_ymin =  -32.470;
static float ecal_ymax =   32.330;

bool NtpTracker::IsInsideECAL() {
  float ecal_coo[2][2] = {{0}};
  float coo_l8[3] = {static_cast<float>(int_lay[7][0]),static_cast<float>(int_lay[7][1]),tracker_layers_z[7]};
  float coo_l9[3] = {static_cast<float>(int_lay[8][0]),static_cast<float>(int_lay[8][1]),tracker_layers_z[8]};
  for (int i=0; i<2; i++) {
    ecal_coo[0][i] = coo_l8[i]+(coo_l9[i]-coo_l8[i])/(coo_l9[2]-coo_l8[2])*(ecal_topz-coo_l8[2]);
    ecal_coo[1][i] = coo_l8[i]+(coo_l9[i]-coo_l8[i])/(coo_l9[2]-coo_l8[2])*(ecal_botz-coo_l8[2]);
  }
  bool isinside = true;
  for (int i=0; i<2; i++) {
    isinside &= (ecal_coo[i][0]>ecal_xmin)&&(ecal_coo[i][0]<ecal_xmax);
    isinside &= (ecal_coo[i][1]>ecal_ymin)&&(ecal_coo[i][1]<ecal_ymax);
  }
  return isinside;
}

bool NtpTrd::IsInsideECAL() {
  float ecal_coo[2][2] = {{0}};
  Interpolate(ecal_topz,ecal_coo[0]);
  Interpolate(ecal_botz,ecal_coo[1]);
  bool isinside = true;
  for (int i=0; i<2; i++) {
    isinside &= (ecal_coo[i][0]>ecal_xmin)&&(ecal_coo[i][0]<ecal_xmax);
    isinside &= (ecal_coo[i][1]>ecal_ymin)&&(ecal_coo[i][1]<ecal_ymax);
  }
  return isinside;
}

double NtpRich::GetNpBetaHyp(int& nhit, double beta, int used, double* res) {
  vector<int> cross;
  for (int ipmt=0; ipmt<5; ipmt++) {
    if (pmt_np_uncorr[ipmt]<5.) continue;
    cross.push_back(pmt_pmt[ipmt]);
  }
  double window[2] = {4e-3,1.2e-3};
  if (res!=0) for (int i=0; i<2; i++) window[i] = res[i];
  nhit = 0;
  double np = 0;
  for (int ihit=0; ihit<30; ihit++) {
    if (hit_np_uncorr[ihit]<=0) continue;
    if ( ((used%10)==1)&&(hit_used[ihit]>=0)&&(hit_used[ihit]<=1) ) continue;
    if ( (((int(used/10)%10)==1)||((int(used/10)%10)==3))&&((hit_stat[ihit]>>31)&0x1) ) continue;
    if ( (((int(used/10)%10)==2)||((int(used/10)%10)==3))&&
         (find(cross.begin(),cross.end(),int(hit_chan[ihit]/16))!=cross.end()) ) continue;
    for (int ibeta=0; ibeta<2; ibeta++) {
      if (fabs(hit_beta[ihit][ibeta]-beta)>3*window[(is_naf)?0:1]) continue;
      np += hit_np_uncorr[ihit];
      nhit++;
    }
  }
  cross.clear();
  return np;
}

double NtpTrd::trdk_like_phe(int i) {
  if ( (i<0)||(i>2) ) return -1;
  if (trdk_like_he[i]<0) return -2;
  if (trdk_like_p[i]<0) return -3;
  double he = exp(-trdk_like_he[i]);
  double p = exp(-trdk_like_p[i]);
  return -log(p/(he+p));
}

double NtpTrd::trdk_like_ep(int i) {
  if ( (i<0)||(i>2) ) return -1;
  if (trdk_like_e[i]<0) return -2;
  if (trdk_like_p[i]<0) return -3;
  double e = exp(-trdk_like_e[i]);
  double p = exp(-trdk_like_p[i]);
  return -log(e/(e+p));
}

double NtpStandAlone::GetSAQTof(int &n, double& rms, int pattern) {
  n = 0;
  rms = 0;
  double mea = 0;
  for (int i=0; i<4; i++) {
    if ((pattern&(1<<i)) == 0) continue;
    if (beta_q_lay[i]<=0) continue;
    mea += beta_q_lay[i];
    rms += rms*rms;
    n++;
  }
  if (n<=0) return 0;
  mea /= n;
  rms /= n;
  rms = sqrt(mea*mea-rms);
  return mea;
}

Event::Event() {
  SHeader  = nullptr;
  Header   = nullptr;
  Trd      = nullptr;
  Tof      = nullptr;
  Tracker  = nullptr;
  Rich     = nullptr;
  Ecal     = nullptr;
  Anti     = nullptr;
  SA       = nullptr;
  MCHeader = nullptr;
  RTI      = nullptr;

  RichBDT = new RichBDTMgr();
  TrackerBDT = new TrackerBDTMgr();
  Classifier = new ClassifierManager();
}

Event::~Event() {
  delete SHeader;
  delete Header;
  delete Trd;
  delete Tof;
  delete Tracker;
  delete Rich;
  delete Ecal;
  delete Anti;
  delete SA;
  delete MCHeader;
}

Double_t Event::GetRichBDT(bool debug){
  // cout << "Called Event::GetRichBDT(), using mgr at " << RichBDT << endl;
  static unsigned int currRun = 0;
  static int          currEv  = 0;
  static double       lastBDT = 0;
  //Only if we are changing event!
  if( SHeader->run == currRun && SHeader->event == currEv ) return lastBDT;
  RichBDT->data->FillData(this);
  if(debug) RichBDT->data->Dump();
  //Set current event
  currRun = SHeader->run;
  currEv  = SHeader->event;
  lastBDT = Rich->is_naf ? RichBDT->GetNafBDT() : RichBDT->GetAglBDT();
  if(debug) std::cout << "Rich BDT result: " << lastBDT << std::endl;
  return lastBDT;
}

Double_t Event::GetTrackerBDT(bool debug){
  // cout << "Called Event::GetTrackerBDT(), using mgr at " << TrackerBDT << endl;
  static unsigned int currRun = 0;
  static int          currEv  = 0;
  static double       lastBDT = 0;
  //Only if we are changing event!
  if( SHeader->run == currRun && SHeader->event == currEv ) return lastBDT;
  TrackerBDT->data->FillData(this);
  if(debug) TrackerBDT->data->Dump();
  //Set current event
  currRun = SHeader->run;
  currEv  = SHeader->event;
  lastBDT = TrackerBDT->GetBDT();
  if(debug) std::cout << "Tracker BDT result: " << lastBDT << std::endl;
  return lastBDT;
}

Double_t Event::GetClassifier(int classifier_type) {
  return Classifier->GetClassifier(this,classifier_type);
}

Double_t Event::GetMinusLogLikelihood(int ll_type) { 
  return Classifier->GetMinusLogLikelihood(this,ll_type); 
}

bool NtpCompact::IsStandalone() {
  bool good_beta = (sa_tof_beta_ncl==4)&&(sa_tof_chisqtn<10);
  bool good_trd = (sa_trd_chi2<10)&&(sa_tof_build==2);
  bool trd_is_inside = (sa_trd_fiducial&0xff)==0xff;
  return good_beta&&good_trd&&trd_is_inside;
}

bool NtpCompact::IsAnalysis() {
  bool good_trk_inn = ((trk_patty&0x2)!=0)&&((trk_patty&0xc)!=0)&&((trk_patty&0x30)!=0)&&((trk_patty&0xc0)!=0);
  bool trk_inn_is_inside = (trk_fiducial[3]&0xff)==0xff;
  return good_trk_inn&&trk_inn_is_inside;
}

int NtpCompact::IsInsideRich() {
  float x = trk_int_rich_rad[0];
  float y = trk_int_rich_rad[1];
  if ((x*x+y*y)>cut_aerogel_external_border) return 0;
  if      (max(fabs(x),fabs(y))<cut_aerogel_naf_border[0]) return 1;
  else if (max(fabs(x),fabs(y))>cut_aerogel_naf_border[1]) return 2;
  return 0;
}

