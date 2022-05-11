#include <bitset>

#include "RichBDT.h"

RichBDTMgr::RichBDTMgr() : initialized(false) { data = new RichBDTData(); }

RichBDTMgr::~RichBDTMgr() { delete data; }

void RichBDTMgr::Init() {
  TMVA::Tools::Instance();
  TString opt = "Silent";
  nafReader = new TMVA::Reader(opt);
  aglReader = new TMVA::Reader(opt);

  SetReaderVariables(nafReader);
  SetReaderVariables(aglReader);

  TString weightPath = "/cvmfs/ams.cern.ch/Offline/dbar/public/release_v5/e1_vdev_181025/aux_data/RichBDTWeights/";
  nafReader->BookMVA("BDT", weightPath + "QualityNaF_BDT.weights.xml");
  aglReader->BookMVA("BDT", weightPath + "QualityAgl_BDT.weights.xml");

  initialized = true;
}

void RichBDTMgr::SetReaderVariables( TMVA::Reader* reader ){
  reader->AddVariable("Richtotused"             , &(data->Richtotused));
  reader->AddVariable("RichPhEl"                , &(data->RichPhEl));
  reader->AddVariable("RICHprob"                , &(data->RICHprob));
  reader->AddVariable("RICHcollovertotal"       , &(data->RICHcollovertotal));
  reader->AddVariable("RICHLipBetaConsistency"  , &(data->RICHLipBetaConsistency));
  reader->AddVariable("RICHTOFBetaConsistency"  , &(data->RICHTOFBetaConsistency));
  reader->AddVariable("RICHChargeConsistency"   , &(data->RICHChargeConsistency));
  reader->AddVariable("RICHPmts"                , &(data->RICHPmts));
  reader->AddVariable("RICHgetExpected"         , &(data->RICHgetExpected));
  reader->AddVariable("tot_hyp_p_uncorr"        , &(data->tot_hyp_p_uncorr));
  reader->AddVariable("Bad_ClusteringRICH"      , &(data->Bad_ClusteringRICH));
  reader->AddVariable("NSecondariesRICHrich"    , &(data->NSecondariesRICHrich));
  reader->AddVariable("HitHVoutdir"             , &(data->HitHVoutdir));
  reader->AddVariable("HitHVoutrefl"            , &(data->HitHVoutrefl));
  // dummy vars. ROOT wants all the spectators to be present...
  reader->AddSpectator("R", &(data->Spectator1));
  reader->AddSpectator("BetaRICH_new", &(data->Spectator2));
}

double RichBDTMgr::GetNafBDT(){
  if (!initialized) Init();
  return nafReader->EvaluateMVA("BDT");
};

double RichBDTMgr::GetAglBDT(){
  if (!initialized) Init();
  return aglReader->EvaluateMVA("BDT");
};

double RichBDTMgr::GetNafMinusLogRarity(){
  if (!initialized) Init();
  double rarity = nafReader->GetRarity("BDT");
  return (rarity<=0) ? 300 : -log(rarity);
};

double RichBDTMgr::GetAglMinusLogRarity(){
  if (!initialized) Init();
  double rarity = aglReader->GetRarity("BDT");
  return (rarity<=0) ? 300 : -log(rarity);
};

void RichBDTData::FillData( Event* event ){
  RICHprob = event->Rich->prob;
  RICHChargeConsistency = event->Rich->q_consistency;
  RICHPmts = event->Rich->npmt;
  RICHgetExpected = event->Rich->np_exp_uncorr;
  tot_hyp_p_uncorr = event->Rich->tot_hyp_p_uncorr[1];
  Richtotused = event->Header->nrichhit - event->Rich->nhit;
  if (event->Rich->np_uncorr > 0)
    RichPhEl = event->Rich->np_exp_uncorr / event->Rich->np_uncorr;
  else
    RichPhEl = 0;
  if (event->Rich->tot_p_uncorr > 0)
    RICHcollovertotal = event->Rich->np_uncorr / event->Rich->tot_p_uncorr;
  else
    RICHcollovertotal = 0;
  RICHLipBetaConsistency =
      fabs(event->Rich->lip_beta - event->Rich->beta_corrected);
  if (event->Rich->beta_corrected > 0)
    RICHTOFBetaConsistency =
        fabs(event->Rich->beta_corrected - event->Tof->beta) /
        event->Rich->beta_corrected;
  else
    RICHTOFBetaConsistency = 0;

  Bad_ClusteringRICH = 0;
  for (int is = 0; is < 10; is++)
    if ((event->Rich->clus_mean[is] - event->Rich->beta) > 0.01)
      Bad_ClusteringRICH++;

  NSecondariesRICHrich = 0;
  for (int is = 0; is < 5; is++)
    if ((event->Rich->pmt_np_uncorr[is] > 5) &&
        (event->Rich->pmt_dist[is] > 3.5))
      NSecondariesRICHrich++;
  HitHVoutdir = event->Rich->tot_hyp_hit_uncorr[1][0];
  HitHVoutrefl = event->Rich->tot_hyp_hit_uncorr[1][1];
}

void RichBDTData::Dump(){
  // std::cout << "Chisquare: " << Chisquare << std::endl;
  std::cout << "RICHprob: " << RICHprob << std::endl;
  std::cout << "RICHPmts: " << RICHPmts << std::endl;
  std::cout << "RICHgetExpected: " << RICHgetExpected << std::endl;
  std::cout << "tot_hyp_p_uncorr: " << tot_hyp_p_uncorr << std::endl;
  std::cout << "Richtotused: " << Richtotused << std::endl;
  std::cout << "RichPhEl: " << RichPhEl << std::endl;
  std::cout << "RICHcollovertotal: " << RICHcollovertotal << std::endl;
  std::cout << "RICHLipBetaConsistency: " << RICHLipBetaConsistency
            << std::endl;
  std::cout << "RICHTOFBetaConsistency: " << RICHTOFBetaConsistency
            << std::endl;
  std::cout << "RICHprob: " << RICHprob << std::endl;
  std::cout << "RICHChargeConsistency: " << RICHChargeConsistency << std::endl;
  std::cout << "Bad_ClusteringRICH: " << Bad_ClusteringRICH << std::endl;
  std::cout << "NSecondariesRICHrich: " << NSecondariesRICHrich << std::endl;
  std::cout << "HitHVoutdir: " << HitHVoutdir << std::endl;
  std::cout << "HitHVoutrefl: " << HitHVoutrefl << std::endl;
}
