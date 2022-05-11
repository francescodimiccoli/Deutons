#include "TrackerBDT.h"

TrackerBDTMgr::TrackerBDTMgr() : initialized(false) {
  data = new TrackerBDTData();
}

TrackerBDTMgr::~TrackerBDTMgr() { delete data; }

void TrackerBDTMgr::Init() {
  TMVA::Tools::Instance();
  // TString opt = "Silent";
  TString opt = "V";
  Reader = new TMVA::Reader(opt);

  SetReaderVariables(Reader);

  TString weightPath = "/cvmfs/ams.cern.ch/Offline/dbar/public/release_v4/"
                       "e1_vdev_180213/aux_data/TrackerBDTWeights/";
  Reader->BookMVA("BDT", weightPath + "HighMass.xml");

  initialized = true;
}

void TrackerBDTMgr::SetReaderVariables(TMVA::Reader *reader) {
  reader->AddVariable("TrEdepL2On", &(data->TrEdepL2On));
  reader->AddVariable("TrEdepL2Off", &(data->TrEdepL2Off));
  reader->AddVariable("TrEdepInnOff0", &(data->TrEdepInnOff0));
  reader->AddVariable("TrEdepInnOff1", &(data->TrEdepInnOff1));
  reader->AddVariable("TrChiSqChY", &(data->TrChiSqChY));
  reader->AddVariable("TrChiSqChX", &(data->TrChiSqChX));
  reader->AddVariable("TrChiSqKaY", &(data->TrChiSqKaY));
  reader->AddVariable("TrChiSqKaX", &(data->TrChiSqKaX));
  reader->AddVariable("TrHalfRUpCh", &(data->TrHalfRUpCh));
  reader->AddVariable("TrHalfRDwCh", &(data->TrHalfRDwCh));
  reader->AddVariable("TrHalfRUpKa", &(data->TrHalfRUpKa));
  reader->AddVariable("TrHalfRDwKa", &(data->TrHalfRDwKa));
  reader->AddVariable("L2ChRes_x", &(data->L2ChRes_x));
  reader->AddVariable("L2ChRes_y", &(data->L2ChRes_y));
  reader->AddVariable("L2KaRes_x", &(data->L2KaRes_x));
  reader->AddVariable("L2KaRes_y", &(data->L2KaRes_y));
  reader->AddVariable("L34Scat_x", &(data->L34Scat_x));
  reader->AddVariable("L34Scat_y", &(data->L34Scat_y));
  reader->AddVariable("L56Scat_x", &(data->L56Scat_x));
  reader->AddVariable("L56FeetDist", &(data->L56FeetDist));
  reader->AddVariable("L2FeetDist", &(data->L2FeetDist));
  reader->AddVariable("TrQMin", &(data->TrQMin));
  reader->AddVariable("TrQAsymm", &(data->TrQAsymm));
  reader->AddVariable("L56Scat_y*(Beta*Rigidity)", &(data->L56ScatBR_y));

  reader->AddSpectator("Category", &(data->Category));
  reader->AddSpectator("Mass", &(data->Mass));
  reader->AddSpectator("Beta", &(data->Beta));
  reader->AddSpectator("Rigidity", &(data->Rigidity));
  reader->AddSpectator("RichBDT", &(data->RichBDT));
}

double TrackerBDTMgr::GetBDT() {
  if (!initialized)
    Init();

  return Reader->EvaluateMVA("BDT");
};

void TrackerBDTData::FillData(Event *event) {
  TrEdepL2On = event->Tracker->edep_lay[1][1][0] / 1000;
  TrEdepL2Off = event->Tracker->edep_lay[1][1][6] / 1000;

  for (int ilay = 1; ilay < 8; ilay++) {
    if (event->Tracker->q_lay[1][ilay] <= 0) {
      TrEdepInnOff1 += event->Tracker->edep_lay[1][ilay][6] / 1000;
    } else {
      TrEdepInnOff0 += event->Tracker->edep_lay[1][ilay][6] / 1000;
    }
  }

  float TrQMax = 0;
  for (int ilay = 1; ilay < 8; ilay++) {
    if (event->Tracker->q_lay[1][ilay] <= 0)
      continue;
    TrQMax = std::max(TrQMax, event->Tracker->q_lay[1][ilay]);
    TrQMin = std::min(TrQMin, event->Tracker->q_lay[1][ilay]);
  }
  TrQAsymm = (TrQMax - TrQMin) / event->Tracker->q_inn[0];

  TrChiSqChY = event->Tracker->chisqn[6][1];
  TrChiSqChX = event->Tracker->chisqn[6][0];
  TrChiSqKaY = event->Tracker->chisqn[1][1];
  TrChiSqKaX = event->Tracker->chisqn[1][0];
  L34Scat_x = event->Tracker->scat_rigi_theta[0][0];
  L34Scat_y = event->Tracker->scat_rigi_theta[0][1];
  L56Scat_x = event->Tracker->scat_rigi_theta[1][0];
  float L56Scat_y = event->Tracker->scat_rigi_theta[1][1];
  L56FeetDist =
      std::min(event->Tracker->feet_dist[3], event->Tracker->feet_dist[4]);
  L2FeetDist = event->Tracker->feet_dist[1];

  L2ChRes_x =
      1e4 * (event->Tracker->coo_lay[1][0] - event->Tracker->int_l2[6][0]);
  L2ChRes_y =  
      1e4 * (event->Tracker->coo_lay[1][1] - event->Tracker->int_l2[6][1]);
  L2KaRes_x =
      1e4 * (event->Tracker->coo_lay[1][0] - event->Tracker->int_l2[1][0]);
  L2KaRes_y =
      1e4 * (event->Tracker->coo_lay[1][1] - event->Tracker->int_l2[1][1]);

  if (fabs(event->Tracker->rig[2][1]) > 0) {
    TrHalfRUpCh =
        (fabs(event->Tracker->rig[1][1]) > 0)
            ? fabs(1 - event->Tracker->rig[1][1] / event->Tracker->rig[2][1])
            : 0;
  }
  if (fabs(event->Tracker->rig[3][1]) > 0) {
    TrHalfRDwCh =
        (fabs(event->Tracker->rig[1][1]) > 0)
            ? fabs(1 - event->Tracker->rig[1][1] / event->Tracker->rig[3][1])
            : 0;
  }
  if (fabs(event->Tracker->rig[2][1]) > 0) {
    TrHalfRUpKa =
        (fabs(event->Tracker->rig[1][1]) > 0)
            ? fabs(1 - event->Tracker->rig[1][1] / event->Tracker->rig[2][1])
            : 0;
  }
  if (fabs(event->Tracker->rig[3][1]) > 0) {
    TrHalfRDwKa =
        (fabs(event->Tracker->rig[1][1]) > 0)
            ? fabs(1 - event->Tracker->rig[1][1] / event->Tracker->rig[3][1])
            : 0;
  }
  Rigidity = event->Tracker->rig[1][1];
  Beta = event->Rich->beta_corrected;
  L56ScatBR_y = L56Scat_y * Beta * Rigidity;

  // dummy spectators
  Category = 1;
  Mass = 1;
  RichBDT = 1;
}

void TrackerBDTData::Dump() {}
