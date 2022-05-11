#include "Tools.h"

float Tools::GetValueDatacard(TFile* file, const char* par_name, int par_numb) {
  float value = 0;
  if (!file) return 0;
  TObjString* DataCards = (TObjString*) file->Get("DataCards");
  TObjArray*  Lines = (TObjArray*) DataCards->GetString().Tokenize("\n");
  bool done = false;
  for (int i = 0; i < Lines->GetEntries(); i++) {
    TObjString* Line = (TObjString*) Lines->At(i);
    // 1st possibility: PART=67
    TObjArray* Equals = (TObjArray*) Line->GetString().Tokenize("=");
    if ( (Equals->GetEntries() > 1) && (((TObjString*)Equals->At(0))->String().CompareTo(par_name) == 0) && (!done) ) {
      // ((TObjString*)Equals->At(0))->Print(); ((TObjString*)Equals->At(1))->Print();
      value = ((TObjString*)Equals->At(1))->String().Atof();
      done = true;
    }
    if (Equals) delete Equals;
    // 2nd possibility: HADR 1
    TObjArray* Spaces = (TObjArray*) Line->GetString().Tokenize(" ");
    if ( (Spaces->GetEntries() > 1) && (((TObjString*)Spaces->At(0))->String().CompareTo(par_name) == 0) &&
         (((TObjString*)Spaces->At(1))->String().CompareTo("=") == 0) && (!done) ) {
      value = ((TObjString*)Spaces->At(1))->String().Atof();
      done = true;
    }
    // 3rd possibility:
    if ( (Spaces->GetEntries() > 1) && (((TObjString*)Spaces->At(0))->String().CompareTo(par_name) == 0) && (!done) ) {
      for (int j = 0; j < Spaces->GetEntries(); j++) {
        TObjString* Space = (TObjString*) Spaces->At(j);
        TObjArray* Equals = (TObjArray*) Space->GetString().Tokenize("=");
        if ( (Equals->GetEntries() == 2) && (((TObjString*)Equals->At(0))->String().Atoi() == par_numb) && (!done) ) {
          value = ((TObjString*)Equals->At(1))->String().Atof();
          done = true;
        }
        if (Equals) delete Equals;
      }
    }
    if (Spaces) delete Spaces;
  }
  if (Lines) delete Lines;
  if (DataCards) delete DataCards; // !
  //if (par_numb>=0) printf("%s %d=%f\n",par_name,par_numb,value);
  //else             printf("%s=%f\n",par_name,value);
  return value;
}

int Tools::RichQC(RichRingR *prich) {
  RichRingR &rich = *prich;
  float cut_prob = 0.01;                      //  Kolmogorov test probability
  float cut_pmt = 3;                          //  number of pmts
  float cut_collovertotal = 0.4;              //  ring photoelctrons / total photoelectrons in the event
  float cut_chargeconsistency = 10;           //  hit/PMT charge consistency test
  float cut_betaconsistency[2] = {0.01, 0.005}; //  beta_lip vs beta_ciemat consistency ([0]=NaF, [1]=aerogel)
  float cut_expphe[2] = {1, 2};               //  expected number of photoelectrons   ([0]=NaF, [1]=aerogel)
  float cut_aerogelexternalborder = 3500.;    //  aerogel external border (r**2)
  float cut_aerogel_nafborder[2] = {17., 19.}; //  aerogel/NaF border                  ([0]=NaF, [1]=aerogel)

  int nbadtiles = 5;
  int kbadtile[nbadtiles];
  kbadtile[0] = 3;
  kbadtile[1] = 7;
  kbadtile[2] = 87;
  kbadtile[3] = 100;
  kbadtile[4] = 108; //  tiles with bad beta reconstruction

  int mask = 0;

  if (!rich.IsGood() || !rich.IsClean()) mask |= (1 << 0);
  if (rich.getProb() < cut_prob) mask |= (1 << 1);
  if (rich.getPMTs() < cut_pmt) mask |= (1 << 2);
  if (rich.getPhotoElectrons() / RichHitR::getCollectedPhotoElectrons() < cut_collovertotal) mask |= (1 << 3);
  if (rich.getPMTChargeConsistency() > cut_chargeconsistency) mask |= (1 << 4);

  float x = rich.getTrackEmissionPoint()[0];
  float y = rich.getTrackEmissionPoint()[1];

  if (rich.IsNaF()) {
    if (rich.getExpectedPhotoelectrons() < cut_expphe[0]) mask |= (1 << 5);
    if (rich.getBetaConsistency() > cut_betaconsistency[0]) mask |= (1 << 6);
    if (max(abs(x), abs(y)) > cut_aerogel_nafborder[0]) mask |= (1 << 7);
  }
  else {
    if (rich.getExpectedPhotoelectrons() < cut_expphe[1]) mask |= (1 << 5);
    if (rich.getBetaConsistency() > cut_betaconsistency[1]) mask |= (1 << 6);
    if (x * x + y * y            > cut_aerogelexternalborder) mask |= (1 << 7);
    if (max(abs(x), abs(y)) < cut_aerogel_nafborder[1]) mask |= (1 << 8);
    for (int kbad = 0; kbad < nbadtiles; kbad++) {
      if (rich.getTileIndex() == kbadtile[kbad]) mask |= (1 << 9);
    }
  }

  if (mask != 0) return -mask;
  if (rich.IsNaF()) return 1;
  else return 2;
}

void Tools::Mysort(int N, float *v, int *sort) {
  bool mask[N]; for (int i = 0; i < N; i++) mask[i] = false;
  for (int i = 0; i < N; i++) {
    double maxv = -10000000.;
    for (int j = 0; j < N; j++) {
      if (mask[j]) continue;
      if (v[j] > maxv) {maxv = v[j]; sort[i] = j;}
    }
    mask[sort[i]] = true;
  }
}

bool Tools::GetInnerNHits(TrTrackR* track, int &nxy, int &ny) {
  nxy = 0;
  ny = 0;
  if (!track) return false;
  int patt_xy = track->GetBitPatternXYJ();
  int patt_y  = track->GetBitPatternJ();
  for (int ilay = 0 + 1; ilay < 9 - 1; ilay++) {
    if ((patt_xy & (1 << ilay)) > 0) nxy++;
    if ((patt_y & (1 << ilay)) > 0) ny++;
  }
  return true;
}

int Tools::axtof(AMSEventR* pev, BetaHR* pbetah, TrTrackR* ptrtk, int *flagoverlap, float *edepoverlap, float *resmeasoverlap, float *resborderoverlap) {
  // init output
  for (int il = 0; il < 4; il++) {edepoverlap[il] = 0; resmeasoverlap[il] = 99; resborderoverlap[il] = 99;}
  for (int i = 0; i < 34; i++) {flagoverlap[i] = -1;}
  // TofHBeta
  if (pev->nTofClusterH() == 0) return -1;
  if (!pbetah) return -1;
  BetaHR &betah = *pbetah;
  // TofHBeta OK
  if (!betah.IsGoodBeta() || !betah.IsTkTofMatch()) return -2;
  // Inner Tracker
  if (pbetah->pTrTrack()) ptrtk = pbetah->pTrTrack();
  int idinner = -1;
  if (ptrtk) idinner = ptrtk->iTrTrackPar(1, 3, 1);
  if (idinner < 0) return -3;
  int iaxtof = 0;
  static float width[4][10] = { // Bar's width
    {22.5, 12., 12., 12., 12., 12., 12., 22.5, 0.,   0.},
    {25.5, 12., 12., 12., 12., 12., 12., 25.5, 0.,   0.},
    {18.5, 12.,  12., 12., 12., 12., 12., 12., 12.,  18.5},
    {26.,  12., 12., 12., 12., 12., 12., 26.,  0.,   0.}
  };
  // TofHPlanes
  for (int il = 0; il < 4; il++) { // loop on layers
    if (betah.TestExistHL(il) && !betah.IsIsolationHL(il)) { // check layers with non isolated clusters
      TofClusterHR* ptofbhcl = betah.GetClusterHL(il);
      AMSPoint fitcoo;
      AMSDir   fitdir;
      double   time;
      if (ptrtk)
	ptrtk->Interpolate(ptofbhcl->Coo[2], fitcoo, fitdir, idinner);
      else
	pbetah->TInterpolate(ptofbhcl->Coo[2], fitcoo, fitdir, time, false);
 
     if (betah.IsInOverlap(il, fitcoo[0], fitcoo[1], 2)) { // check only clusters in the overlap region
        iaxtof++;
        //      printf(" iaxtof %i --- Plane-Bar %i-%i: Tk Coo: %6.3f %6.3f %6.3f \n",iaxtof,il,ptofbhcl->Bar,fitcoo[0],fitcoo[1],fitcoo[2]);
        int n = 0;
        float rm[2] = {99, 99};
        float rb[2] = {99, 99};
        float eb[2] = {0};
        int   cl[2] = {0, 0};
        int km = (il == 0 || il == 3) ? 0 : 1;
        int kb = (il == 0 || il == 3) ? 1 : 0;
        for (int icl = 0; icl < pev->nTofClusterH(); icl++) { // look for neighbor clusters in the overlap region
          TofClusterHR* pcl = pev->pTofClusterH(icl);
          if (!pcl) continue;
          if (pcl->Layer != il) continue;
          if (abs(pcl->Bar - ptofbhcl->Bar) == 1) {         // ... neighbor
            if (ptrtk) 
	      ptrtk->Interpolate(pcl->Coo[2], fitcoo, fitdir, idinner);
	    else 
	      pbetah->TInterpolate(pcl->Coo[2], fitcoo, fitdir, time, false);
            cl[n] = icl;
            rm[n] = abs(pcl->Coo[km] - fitcoo[km]);
            rb[n] = abs(pcl->Coo[kb] - fitcoo[kb]) - width[pcl->Layer][pcl->Bar] / 2.;
            eb[n] = pcl->GetEdep();
            //      printf(" Cluster %i (L%iB%i): Coo: %6.3f %6.3f %6.3f --- Edep %6.3f Res: %6.3f %6.3f EdepOverlap %6.3f\n"
            //             ,icl,pcl->Layer,pcl->Bar,pcl->Coo[0],pcl->Coo[1],pcl->Coo[2],pcl->GetEdep(),rm[n],rb[n],eb[n]);
            n++;
          }
        }
        int side = (abs(rb[0]) < abs(rb[1])) ? 0 : 1;
        resmeasoverlap[il]  = rm[side];
        resborderoverlap[il] = rb[side];
        edepoverlap[il]     = eb[side];
        flagoverlap[cl[side]] = 1;
      }
    }
  }
  return iaxtof;
}

int Tools::TkSelect(AMSEventR* pEvent, TrTrackR* pTrTrack, BetaHR* pBetaH) {
  if (!pTrTrack) return 0;
  float tof_betah = (pBetaH) ? pBetaH->GetBeta() : 1;
  int   tk_charge = TMath::Max(int(pTrTrack->GetInnerQ(tof_betah)), 1);
  int tk_hitb[2];
  tk_hitb[0] = pTrTrack->GetBitPatternXYJ();
  tk_hitb[1] = pTrTrack->GetBitPatternJ();
  int nhit[2] = {0};
  for (int ilay = 0 + 1; ilay < 9 - 1; ilay++) {
    for (int ixy = 0; ixy < 2; ixy++) {
      if ((tk_hitb[ixy] & (1 << ilay)) > 0)nhit[ixy]++;
    }
  }
  int haslay1 = 0;
  if ((tk_hitb[1] & (1 << 1)) > 0) haslay1 = 1;
  int algo  = 1; // 1:Choutko  2:Alcaraz  (+10 no multiple scattering)
  int patt  = 3; // 3:Inner only
  int refit = 1; //
  int mfit;
  mfit = pTrTrack->iTrTrackPar(algo, patt, refit); //Choutko+Proton Mass MulScat Refit
  if (mfit < 0) return -4;
  float tk_rigidity = pTrTrack->GetRigidity(mfit);
  //--Vitaly+Alcaraz+ChikanF Calulate Mass to do MulScat Refit
  float aovz = 1;
  if (tk_charge > 1) aovz = 2.;
  float mass = (fabs(tof_betah) > 0.9) ? tk_charge * aovz * 0.938272297 : 0; //can replace TofBeta with RichBeta // or directly use Study Particle Mass
  float tk_rigidityI[3];
  float tk_rigidityIUL[3][2];
  float tk_chisI[3][3];//inner chis
  for (int ialg = 1; ialg <= 3; ialg++) {
    patt = 3; //inter
    mfit = pTrTrack->iTrTrackPar(ialg, patt, refit, mass, tk_charge, tof_betah);
    if (mfit < 0) {
      tk_rigidityI[ialg - 1] = 9999;
      tk_chisI[ialg - 1][0] = tk_chisI[ialg - 1][1] = tk_chisI[ialg - 1][2] = 9999;
    }
    else {
      tk_rigidityI[ialg - 1] = pTrTrack->GetRigidity(mfit);
      tk_chisI[ialg - 1][0]  = pTrTrack->GetNormChisqX(mfit);
      tk_chisI[ialg - 1][1]  = pTrTrack->GetNormChisqY(mfit);
      tk_chisI[ialg - 1][2]  = pTrTrack->GetChisq(mfit);
    }
    patt = 1;//inner up
    mfit = pTrTrack->iTrTrackPar(ialg, patt, refit, mass, tk_charge, tof_betah);
    if (mfit < 0) tk_rigidityIUL[ialg - 1][0] = 9999;
    else        tk_rigidityIUL[ialg - 1][0] = pTrTrack->GetRigidity(mfit);
    patt = 2;//inn/r down
    mfit = pTrTrack->iTrTrackPar(ialg, patt, refit, mass, tk_charge, tof_betah);
    if (mfit < 0) tk_rigidityIUL[ialg - 1][1] = 9999;
    else        tk_rigidityIUL[ialg - 1][1] = pTrTrack->GetRigidity(mfit);
  }
  //--cut All Fitting Algorithem exist
  int rigf = 1;
  float cutval = 999;
  for (int ialg = 0; ialg < 3; ialg++) if (tk_rigidityI[ialg] == 0 || fabs(tk_rigidityI[ialg]) > cutval) rigf = 0;
  int righf = 1;
  for (int ialg = 0; ialg < 3; ialg++) {
    for (int iud = 0; iud < 2; iud++) {
      if (tk_rigidityIUL[ialg][iud] == 0) righf = 0;
      if (fabs(tk_rigidityIUL[ialg][iud]) > cutval) righf = 0;
    }
  }
  if (rigf == 0 || righf == 0) return -5;
  //-----Cut Value
  double disru      = tk_rigidityIUL[0][0] / tk_rigidityI[0]; //UpRig/AllRig
  double disrd      = tk_rigidityIUL[0][1] / tk_rigidityI[0]; //DoRig/AllRig
  double disrignovc = (tk_rigidityI[0] - tk_rigidity) / (tk_rigidityI[0] + tk_rigidity); // Proton MassFit-Nucleus MassFit
  double disrigvc   = (tk_rigidityI[1] - tk_rigidityI[0]) / (tk_rigidityI[1] + tk_rigidityI[0]); //Alcaraz-Choutko
  double disrigch   = (tk_rigidityI[2] - tk_rigidityI[0]) / (tk_rigidityI[2] + tk_rigidityI[0]); //Chikan-Choutko
  //--Track Edep
  const int NTKL = 9;
  const int NTKS = 2;
  const int NTKDIS = 4;
  float tk_dedx[9][2] = {{0}}; //track dedx(X+Y)
  float tk_dedx_ns[9][2][NTKDIS] = {{{0}}}; //near 1cm +2cm+4cm+8cm sum all
  const float TKDIS[NTKDIS] = {1, 2, 4, 6}; //dedx dis(cm) near tk
  for (int ilay = 0; ilay < NTKL; ilay++) {
    AMSPoint postr;
    AMSDir dirtr;
    algo = 1;
    patt = 3; //inner
    mfit = pTrTrack->iTrTrackPar(algo, patt, refit, mass, tk_charge, tof_betah);
    pTrTrack->InterpolateLayerJ(ilay + 1, postr, dirtr, mfit);
    TrRecHitR *tkhit = pTrTrack->GetHitLJ(ilay + 1);
    TrClusterR *tkcl[2] = {0};
    if (tkhit) {
      tkcl[0] = tkhit->GetXCluster(); tkcl[1] = tkhit->GetYCluster();
      for (int ixy = 0; ixy < NTKS; ixy++) {
        if (tkcl[ixy]) tk_dedx[ilay][ixy] = tkcl[ixy]->GetEdep(); //each layer dedx
      }
    }
    for (int icl = 0; icl < pEvent->nTrCluster(); icl++) {
      TrClusterR *cltr = pEvent->pTrCluster(icl);
      if ((cltr == 0) || (cltr == tkcl[0]) || (cltr == tkcl[1]) || cltr->GetLayerJ() != ilay + 1) {continue;} //Non used
      float cldis = 999999;
      for (int m = 0; m < cltr->GetMultiplicity(); m++) { //reso multi
        float muldis = cltr->GetGCoord(m) - postr[cltr->GetSide()];
        if (fabs(muldis) < cldis) cldis = fabs(muldis);
      }
      for (int idis = 0; idis < NTKDIS; idis++) {
        if (cldis < TKDIS[idis]) tk_dedx_ns[ilay][cltr->GetSide()][idis] += cltr->GetEdep(); //near 1+2+4+8cm No-Track dedx
      }
    }
  }
  //-Cut Value
  double mindep = 9999999;
  double maxdep = 0;
  double sumdep = 0;
  int uhity = 0;
  for (int ilay = 0 + 1; ilay < 9 - 1; ilay++) { //Inner Edep
    if (tk_dedx[ilay][1] == 0)continue;
    if (tk_dedx[ilay][1] < mindep) mindep = tk_dedx[ilay][1];
    if (tk_dedx[ilay][1] > maxdep) maxdep = tk_dedx[ilay][1];
    sumdep += tk_dedx[ilay][1];
    uhity++;
  }
  double disdep = (maxdep - mindep) / (maxdep + mindep);
  //---Near max
  double maxdep_n[NTKDIS] = {0};
  double sumdep_n[NTKDIS] = {0};
  for (int idis = 0; idis < NTKDIS; idis++) {
    for (int ilay = 0 + 1; ilay < 9 - 1; ilay++) {
      if (tk_dedx_ns[ilay][1][idis] == 0) continue;
      if (tk_dedx_ns[ilay][1][idis] > maxdep_n[idis]) maxdep_n[idis] = tk_dedx_ns[ilay][1][idis];
      sumdep_n[idis] += tk_dedx_ns[ilay][1][idis];
    }
  }
  //////--Alcut
  bool cut[100];
  cut[20] = (tk_chisI[0][0] < 5); //X Chis cut
  cut[21] = (tk_chisI[0][1] < 4); //Y Chis cut
  cut[22] = (disru > 0.8 && disru < 1.2); //URig/ARig cut
  cut[23] = (disrd > 0.93 && disrd < 1.1); //DRig/ARig
  cut[24] = (disdep < 0.5); //TkEdep  (Max-Min)/(Max+Min)
  cut[25] = (maxdep_n[1] < 0.2); //near 2cm Edep<0.2MeV isolate
  cut[26] = (sumdep_n[1] / (sumdep + sumdep_n[1]) < 0.1); // 2cm Edep/(TkEdep+2cm Edep)<0.1 near isolate
  cut[27] = (fabs(disrignovc) < 0.005); //
  cut[28] = (fabs(disrigvc) < 0.01); //Alcarz-VC
  cut[29] = (nhit[0] >= 4); //Inner Hit cut
  cut[30] = (nhit[1] >= 6);
  cut[31] = (fabs(disrigch) < 0.05); //Chis-VC
  cut[32] = (haslay1 == 1); //Layer 2 Hit Cut
  cut[33] = (mindep / (sumdep / uhity) > 0.6); //TkEdp Min/(Average Edep)>0.6
  //---rig cut
  bool cuthit  = (cut[29] && cut[30] && cut[32]);
  bool cutchis = (cut[20] && cut[21]);
  bool cutalg  = (cut[27] && cut[28] && cut[31]);
  bool cutud   = (cut[22] && cut[23]);
  bool cutdep  = (cut[24] && cut[25] && cut[26] && cut[33]);
  bool tkcut1  = (cutdep && cutchis && cutud && cuthit && cutalg);
  if (!cuthit)  return -11;
  if (!cutchis) return -12;
  if (!cutalg)  return -13;
  if (!cutud)   return -14;
  if (!cutdep)  return -15;
  if (!tkcut1)  return -10;
  return 1;
}

void Tools::TrClusterOnLayerJ(int layerJ, AMSEventR* pev, TrTrackR &trtk, int idfit, int *nclx, float *eclx, int *ncly, float* ecly, float &maxedepx, float &d2maxedepx, float &maxedepy, float &d2maxedepy) {
  float window[5] = {0.1, 1., 2., 5., 10.};
  for (int i = 0; i < 7; i++) {
    nclx[i] = 0;
    eclx[i] = 0;
    ncly[i] = 0;
    ecly[i] = 0;
  }
  maxedepx   = 0;
  maxedepy   = 0;
  d2maxedepx = 0;
  d2maxedepy = 0;
  // identify the clusters used in trtrack
  TrClusterR *ptrtkclx = trtk.GetHitLJ(layerJ) != NULL ? trtk.GetHitLJ(layerJ)->GetXCluster() : 0;
  TrClusterR *ptrtkcly = trtk.GetHitLJ(layerJ) != NULL ? trtk.GetHitLJ(layerJ)->GetYCluster() : 0;
  // Get track interpolation to layerJ
  AMSPoint fitcoo;
  AMSDir   fitdir;
  trtk.InterpolateLayerJ(layerJ, fitcoo, fitdir, idfit);
  // Loop on clusters
  //  float thrnoise=50;
  float thrnoise = 0;
  for (int icl = 0; icl < pev->nTrCluster(); icl++) {
    TrClusterR* ptrcl = (TrClusterR *) pev->pTrCluster(icl);
    if (!ptrcl) continue;
    if (ptrcl->GetLayerJ() != layerJ) continue; // only clusters in the required layer
    float edep = 1000 * ptrcl->GetEdep();
    if (ptrcl == ptrtkclx) { nclx[0]++; eclx[0] += edep; continue; } //[0]: x-cluster in track
    if (ptrcl == ptrtkcly) { ncly[0]++; ecly[0] += edep; continue; } //[0]: y-cluster in tracl
    if (edep < thrnoise) continue; // apply a threshold to account for noise
    int side = ptrcl->GetSide();
    // Deal with multiplicities: Choose closest to track
    float dcl2tk = 99999.;
    for (int im = 0; im < ptrcl->GetMultiplicity(); im++) {
      float coo = ptrcl->GetGCoord(im);
      if (abs(fitcoo[side] - coo) < dcl2tk) {
        dcl2tk = abs(fitcoo[side] - coo);
      }
    }
    // Count clusters, Edep
    if (side == 0) {
      nclx[6]++;
      eclx[6] += edep;
      if (edep > maxedepx) {
        maxedepx = edep;
        d2maxedepx = dcl2tk;
      }
      for (int iw = 0; iw < 5; iw++) {
        if (abs(dcl2tk) < window[iw]) {
          nclx[iw + 1]++;
          eclx[iw + 1] += edep;
        }
      }
    }
    else {
      ncly[6]++;
      ecly[6] += edep;
      if (edep > maxedepy) {
        maxedepy = edep;
        d2maxedepy = dcl2tk;
      }
      for (int iw = 0; iw < 5; iw++) {
        if (abs(dcl2tk) < window[iw]) {
          ncly[iw + 1]++;
          ecly[iw + 1] += edep;
        }
      }
    }
  }
}

int Tools::MyEcalShowerHits(AMSEventR &event, EcalShowerR &ecal) {
  int Nhits = 0;
  for (int j = 0; j < ecal.NEcal2DCluster(); j++) {
    if (!ecal.pEcal2DCluster(j)) continue;
    Ecal2DClusterR &cluster2D = *(ecal.pEcal2DCluster(j));
    if (!(cluster2D.Status & (1 << 5))) continue;
    for (int k = 0; k < cluster2D.NEcalCluster(); k++) {
      if (!cluster2D.pEcalCluster(k)) continue;
      EcalClusterR &cluster = *(cluster2D.pEcalCluster(k));
      if (!(cluster.Status & (1 << 5))) continue;
      for (int l = 0; l < cluster.NEcalHit(); l++) {
        if (!cluster.pEcalHit(l)) continue;
        EcalHitR &hit = *(cluster.pEcalHit(l));
        if (!(hit.Status & (1 << 5))) continue;
        Nhits++;
      }
    }
  }
  return Nhits;
}

bool Tools::MipsTagCut(AMSEventR &event, EcalShowerR &ecal, bool looseCut) {
  double alpha = 1;
  bool mipcut[18];
  bool NOMIPS[18];
  if (looseCut) alpha = 0.020;
  else          alpha = 0.020;
  int laymaxfrac = 0;
  double maxfrac = 0;
  for (int i = 1; i < 18; i++) {
    if (ecal.EnergyFractionLayer[i] > maxfrac) {
      maxfrac = ecal.EnergyFractionLayer[i];
      laymaxfrac = i;
    }
  }
  for (int i = 0; i < 18; i++) mipcut[i] = laymaxfrac > i ? ecal.EnergyD * ecal.EnergyFractionLayer[i] / 1000 > alpha : 1; //0.020 GeV
  NOMIPS[0] = true;
  for (int i = 1; i < 18; i++) NOMIPS[i] = NOMIPS[i - 1] && mipcut[i];
  return NOMIPS[17];
}

bool Tools::IsInsideTRD(float TRDCoo[2][2]) {
  static int nTrdCenter = 9;
  static float AccepCenterX[9] = { -80.0, -47.0, 47.0, 80.0,  80.0,  47.0, -47.0, -80.0, -80.0};
  static float AccepCenterY[9] = { 43.5,  75.5, 75.5, 43.5, -43.5, -75.5, -75.5, -43.5,  43.5};
  static int nTrdTop = 37;
  static float AccepTopX[37] = { -99.0, -89.0, -89.0, -78.7, -78.7, -67.8, -67.8, -57.7, -57.7, 57.7, 57.7, 67.8, 67.8, 78.7, 78.7, 89.0, 89.0, 99.0,
    99.0, 89.0, 89.0, 78.7, 78.7, 67.8, 67.8, 57.7, 57.7,-57.7,-57.7,-67.8,-67.8,-78.7,-78.7,-89.0,-89.0,-99.0,-99.0};
  float AccepTopY[37] = { 54.5, 54.5, 62.5, 62.5, 74.0, 74.0, 84.0, 84.0, 95.3, 95.3, 84.0, 84.0, 74.0, 74.0, 62.5, 62.5, 54.5, 54.5,
    -51.7,-51.7,-62.2,-62.2,-72.0,-72.0,-82.5,-82.5,-92.5,-92.5,-82.5,-82.5,-72.0,-72.0,-62.2,-62.2,-51.7,-51.7, 54.5};
  bool passTrdCenter = TMath::IsInside(TRDCoo[0][0], TRDCoo[0][1], nTrdCenter, AccepCenterX, AccepCenterY);
  bool passTrdTop    = TMath::IsInside(TRDCoo[1][0], TRDCoo[1][1], nTrdTop,   AccepTopX,   AccepTopY);
  return (passTrdCenter && passTrdTop);
}

bool Tools::IsInsideTkL1(float  xin, float yin ){
  double yp=7.3;
  double xp=4.14;
  float xx=fabs(xin);
  float yy=fabs(yin);
  if      ( yy >6.5*yp         ||       xx > 15  *xp) return false;
  else  if( yy<=2.5*yp &&               xx<= 15  *xp) return true;
  else  if( yy >2.5*yp && yy<=3.5*yp && xx<= 14.5*xp) return true;
  else 	if( yy >3.5*yp && yy<=4.5*yp && xx<= 14. *xp) return true;
  else  if( yy >4.5*yp && yy<=5.5*yp && xx<= 13  *xp) return true;
  else  if( yy >5.5*yp && yy<=6.5*yp && xx<= 12.5*xp) return true;
  else  return false ;
}

void Tools::TofUnusedHits(AMSEventR* pev, BetaHR* pBetaH, 
  float clsdt[4], //[in_time/of_time(2top+2bot)],aver.times wrt beta-used hit
  float clsed[4], //[in_time/of_time(2top+2bot)],tot.edep
  short int clsn[4] //[in_time/of_time(2top+2bot)],clusters number
) {
  for (int i=0; i<4; i++) {
    clsdt[i] = 0;
    clsed[i] = 0;
    clsn[i] = 0;
  }
  //-->look around(time) used hits:
  for (int il=0; il<4; il++) { //layer loop
    int itb = (il<2) ? 0 : 1;
    if (pBetaH->TestExistHL(il)) { //hit exists
      float ltime = pBetaH->GetTime(il); //ns
      int ntfcls = pev->nTofClusterH();
      for (int icl=0; icl<ntfcls; icl++) {
        TofClusterHR &tfcl = pev->TofClusterH(icl);
        if (tfcl.Layer!=il) continue;
        if (tfcl.NBetaHUsed()>0) continue; //used by BetaHtfcl.Time();
        if (!tfcl.IsGoodTime()) continue;
        float cltime = tfcl.Time; //ns
        float ed = tfcl.GetEdep(); //MeV
        float dt = cltime-ltime; //later cluster has positive dt
        float tcut(0);
        if (itb==0) tcut = 10; //ns,top
        else tcut = 4; //ns, bot
        int itm(-1);
        if (fabs(dt)<=tcut) itm = 0; //around-time hits ("in time")
        if (dt>tcut) itm = 1; //later hits ("off time")
        if (itm>=0) {
          int indx = itm+2*itb;
          if(indx>=0 && indx<4){
            clsdt[indx] += dt;
            clsed[indx] += ed;
            clsn[indx] += 1;
          }
        }     
      }//-->endof "secondary clust-loop"
    }//-->endof "beta-hit exists"
  }//-->endof "Tof layer loop"
  for (int i=0; i<4; i++) { //make averages
    if (clsn[i]>0) {
      clsdt[i] /= clsn[i];
      clsed[i] /= clsn[i];
    }
  }
}

