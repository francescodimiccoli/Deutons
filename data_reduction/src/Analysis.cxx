#include <functional>

#include "Analysis.h"

// some globals
static string        current_filename;
static int           current_entry = 0;
static unsigned int  current_rti_time = 0;
static clock_t       clock_start;
static clock_t       clock_end;

int   Analysis::char_ref_fit = 1;
float Analysis::mass_ref_fit = TrFit::Mproton;
int   Analysis::algo_ref_fit = 1; // Choutko
int   Analysis::patt_ref_fit = 0; // Max span
int   Analysis::refi_ref_fit = 21; // Refit if does not exist
int   Analysis::id_ref_fit = -1;

int   Analysis::pres_strategy      = 1;

static bool disable_fill = false;

GHBManager*       Analysis::pRichCorr = 0;
GeomHashEnsemble* Analysis::pRichUnifAgl = 0;
GeomHashEnsemble* Analysis::pRichUnifNaf = 0;

void Analysis::InitPrescaler() {
  //                        name                method           false  true
  vPrescale.AddCondition("Acc",          &Analysis::Accepted      , 0   , 1);
  vPrescale.AddCondition("HasBeta",      &Analysis::HasBeta       , 100,  1);
  vPrescale.AddCondition("DownG",        &Analysis::IsDowngoing   , 100 , 1);
  vPrescale.AddCondition("InL1",         &Analysis::IsInsideL1    , 2  ,  1);

  vPrescale.AddCondition("HasTk",        &Analysis::HasTrack      , 100 , 1);

  vPrescale.AddCondition("Nega",         &Analysis::IsNegative    , 1   , 1);

  vPrescale.AddCondition("HasRich",      &Analysis::Is_rich       , 10 ,  1);
  vPrescale.AddCondition("IsNaF",        &Analysis::IsNaF         , 1 ,   1);

  vPrescale.AddCondition("IsZgt1",       &Analysis::IsZgt1        , 10 ,  1);
  vPrescale.AddCondition("IsZgt2",       &Analysis::IsZgt2        , 2   , 1);

  vPrescale.PrintCond();

  // Bad
  vPrescale.AddCategory("0000000000", "1000000000");
  // Good Nobeta notrack
  vPrescale.AddCategory("1000000000", "1100100000",0);
  // Good Nobeta track
  vPrescale.AddCategory("1000100000", "1100100000");

  //=============================================================
  // accepted beta down inside
  // notrack
  vPrescale.AddCategory("1111000000", "1111100000");
  // Track, associated negativ
  vPrescale.AddCategory("1111110000", "1111110000");

  // Track,  positiv  !RICH !HZ  !HZZ
  vPrescale.AddCategory("1111100000", "1111111010");
  // Track,  positiv  !RICH  HZ  !HZZ
  vPrescale.AddCategory("1111100010", "1111111011");
  // Track,  positiv  !RICH  HZ  HZZ
  vPrescale.AddCategory("1111100011", "1111111011");

  // Track,  positiv RICH NAF !HZ  !HZZ
  vPrescale.AddCategory("1111101100", "1111111110");
  // Track,  positiv RICH NAF HZ   !HZZ
  vPrescale.AddCategory("1111101110", "1111111111");
  // Track,  positiv RICH NAF HZ   HZZ
  vPrescale.AddCategory("1111101111", "1111111111");

  // Track,  positiv RICH !NAF !HZ  !HZZ
  vPrescale.AddCategory("1111101000", "1111111110");
  // Track,  positiv RICH !NAF HZ   !HZZ
  vPrescale.AddCategory("1111101010", "1111111111");
  // Track,  positiv RICH !NAF HZ   HZZ
  vPrescale.AddCategory("1111101011", "1111111111");

  //=============================================================
  // accepted beta down outside
  // notrac
  vPrescale.AddCategory("1110000000", "1111100000");
  // Track,  negativ
  vPrescale.AddCategory("1110110000", "1111110000",1);

  // Track,  positiv !RICH !HZ  !HZZ
  vPrescale.AddCategory("1110100000", "11111111010");
  // Track,  positiv  !RICH  HZ  !HZZ
  vPrescale.AddCategory("1110100010", "11111111011");
  // Track,  positiv L1 !RICH  HZ  HZZ
  vPrescale.AddCategory("1110100011", "11111111011");

  // Track,  positiv  RICH  NAF !HZ  !HZZ
  vPrescale.AddCategory("1110101100", "1111111110");
  // Track,  positiv  RICH  NAF HZ   !HZZ
  vPrescale.AddCategory("1110101110", "1111111111");
  // Track,  positiv  RICH  NAF HZ   HZZ
  vPrescale.AddCategory("1110101111", "1111111111");

  // Track,  positiv  RICH  !NAF !HZ  !HZZ
  vPrescale.AddCategory("1110101000", "1111111110");
  // Track,  positiv  RICH  !NAF HZ   !HZZ
  vPrescale.AddCategory("1110101010", "1111111111");
  // Track,  positiv  RICH  !NAF HZ   HZZ
  vPrescale.AddCategory("1110101011", "1111111111");

  //===================================================================================
  // accepted beta UP inside
  // notrac
  vPrescale.AddCategory("1100000000", "1110100000",1000);
  // Track,  negativ
  vPrescale.AddCategory("1100110000", "1110110000");

  // Track,  positiv  RICH NAF
  vPrescale.AddCategory("1100100000", "1110110000");

  vPrescale.PrintCateg();
  vPrescale.CheckConsistency();
  vPrescale.BuildTree();

  if (pres_strategy==0)
    printf("PRESCALING IN DUMMY MODE (pres_strategy==0)  NO PRESCALING ACTUALLY APPLIED!!\n");
}

Analysis::Analysis(bool is_mc) {
  pFileOutput = nullptr;
  pTreeRTI = nullptr;
  pTreeFile = nullptr;
  pTreeEvent = nullptr;
  pTreeProc = nullptr;
  pTreeComp = nullptr;
  isMC = is_mc;
  pBeta = nullptr;
  pBetaH = nullptr;
  pTrTrack = nullptr;
  pTrTrack2 = nullptr;
  pTrdTrack = nullptr;
  pTrdHTrack = nullptr;
  pRichRing = nullptr;
  pRichRingB = nullptr;
  pEcalShower = nullptr;
  pBetaH_SA = nullptr;
  pTrdTrack_SA = nullptr;
  current_filename.clear();
  current_entry = 0;
  pRichOcc = nullptr;
  RichBDT = std::unique_ptr<RichBDTMgr>(new RichBDTMgr());
}

Analysis::~Analysis() {
  ClearData();
  ClearSelection();
  pFileOutput = 0;
  pTreeRTI = 0;
  pTreeFile = 0;
  pTreeEvent = 0;
  pTreeProc = 0;
  pTreeComp = 0;
}

void Analysis::ClearData() {
  memset(&rtiInfo, 0, sizeof(RTIInfo));
  memset(&fileInfo, 0, sizeof(FileInfo));
  memset(&fileMCInfo, 0, sizeof(FileMCInfo));
  memset(&procInfo, 0, sizeof(ProcInfo));
  memset(&ntpHeader, 0, sizeof(NtpHeader));
  memset(&ntpSHeader, 0, sizeof(NtpSHeader));
  memset(&ntpMCHeader, 0, sizeof(NtpMCHeader));
  memset(&ntpTrd, 0, sizeof(NtpTrd));
  memset(&ntpTof, 0, sizeof(NtpTof));
  memset(&ntpTracker, 0, sizeof(NtpTracker));
  memset(&ntpRich, 0, sizeof(NtpRich));
  memset(&ntpEcal, 0, sizeof(NtpEcal));
  memset(&ntpAnti, 0, sizeof(NtpAnti));
  memset(&ntpStandAlone, 0, sizeof(NtpStandAlone));
  memset(&ntpCompact, 0, sizeof(NtpCompact));
}

void Analysis::ClearSelection() {
  pBeta = nullptr;
  pBetaH = nullptr;
  pTrTrack = nullptr;
  pTrTrack2 = nullptr;
  pTrdTrack = nullptr;
  pTrdHTrack = nullptr;
  pRichRing = nullptr;
  pRichRingB = nullptr;
  pEcalShower = nullptr;
  pBetaH_SA = nullptr;
  pTrdTrack_SA = nullptr;
}

void Analysis::UBegin() {

  printf("\n\n Analysis::UBegin =============================================\n");

  clock_start = clock();

  // prepare output
  pFileOutput = TFile::Open(GetOption(), "recreate");
  pFileOutput->SetCompressionLevel(7);

  if (!isMC) {
    pTreeRTI = new TTree("RTI", "RTI");
    pTreeRTI->SetAutoFlush(0);
    pTreeRTI->SetAutoSave(1000000000);
    pTreeRTI->Branch("SHeader", "SHeader", &ntpSHeader);
    pTreeRTI->Branch("RTIInfo", "RTIInfo", &rtiInfo);
  }

  pTreeFile = new TTree("File", "File");
  pTreeFile->SetAutoFlush(0);
  pTreeFile->SetAutoSave(1000000000);
  pTreeFile->Branch("FileInfo", "FileInfo", &fileInfo);
  if (isMC) pTreeFile->Branch("FileMCInfo", "FileMCInfo", &fileMCInfo);

  if (compact!=2) {
    pTreeEvent = new TTree("Event", "Event");
    pTreeEvent->SetAutoFlush(0);
    pTreeEvent->SetAutoSave(1000000000);

    pTreeEvent->Branch("SHeader", "SHeader", &ntpSHeader);
    pTreeEvent->Branch("Header", "Header", &ntpHeader);
    if (isMC) pTreeEvent->Branch("MCHeader", "MCHeader", &ntpMCHeader);
    pTreeEvent->Branch("Trd", "Trd", &ntpTrd);
    pTreeEvent->Branch("Tof", "Tof", &ntpTof);
    pTreeEvent->Branch("Tracker", "Tracker", &ntpTracker);
    pTreeEvent->Branch("Rich", "Rich", &ntpRich);
    pTreeEvent->Branch("Ecal", "Ecal", &ntpEcal);
    pTreeEvent->Branch("Anti", "Anti", &ntpAnti);
    pTreeEvent->Branch("SA", "SA", &ntpStandAlone);
  }

  pTreeProc = new TTree("Processing", "Processing");
  pTreeProc->SetAutoFlush(0);
  pTreeProc->SetAutoSave(1000000000);
  pTreeProc->Branch("Processing", "Processing", &procInfo);

  if (compact!=0) {
    pTreeComp = new TTree("Compact", "Compact");
    pTreeComp->SetAutoFlush(0);
    pTreeComp->SetAutoSave(1000000000);
    pTreeComp->Branch("SHeader", "SHeader", &ntpSHeader);
    pTreeComp->Branch("Compact", "Compact", &ntpCompact);
  }

  ClearData();

  procInfo.utctime[0] = time(0);
  procInfo.is_mc = isMC;
  procInfo.pres_strategy = pres_strategy;
  sprintf(procInfo.ams_ver, "vdev_181025");
  procInfo.lib_ver[0] = 7; // to be incremented for each production release (if ntuple structure has been modified)
  procInfo.lib_ver[1] = 0; // to be incremented for each commit in a production branch (if ntuple structure has been modified)
  procInfo.exe_ver[0] = 7; // to be incremented for each production release
  procInfo.exe_ver[1] = 0; // to be incremented for each commit in a production branch
  printf(" === procInfo.utctime[0]:    %d\n", procInfo.utctime[0]);
  printf(" === procInfo.is_mc:         %d\n", procInfo.is_mc);
  printf(" === procInfo.pres_strategy: %d\n", procInfo.pres_strategy);
  printf(" === procInfo.ams_ver:       %s\n", procInfo.ams_ver);
  printf(" === procInfo.lib_ver:       %d.%d\n", procInfo.lib_ver[0], procInfo.lib_ver[1]);
  printf(" === procInfo.exe_ver:       %d.%d\n", procInfo.exe_ver[0], procInfo.exe_ver[1]);

  InitPrescaler();

  clock_end = clock();
  procInfo.time[0] += float(clock_end - clock_start) / CLOCKS_PER_SEC;

  clock_start = clock();
}

void Analysis::ProcessFile() {
  TFile* file = Tree()->GetCurrentFile();
  string name(file->GetName());
  if (current_filename != name) {
    // fill
    if (!current_filename.empty()) {
      cout << " Analysis::ProcessFile fill information for file " << name << endl;
      pTreeFile->Fill();
      pTreeFile->Show(pTreeFile->GetEntries() - 1);
    }
    // init
    current_filename = name;
    fileInfo.run = Run();
    fileInfo.event[0] = Event();
    fileInfo.event[1] = Event();
    fileInfo.nentries = 0;
    fileInfo.utime[0] = fHeader.Time[0];
    fileInfo.utime[1] = fHeader.Time[0];
    if (isMC) {
      fileMCInfo.event[0] = Event();
      fileMCInfo.event[1] = Event();
      MCEventgR* primary = (nMCEventg() != 0) ? pMCEventg(0) : 0;
      fileMCInfo.charge = (primary) ? primary->Charge : 0;
      fileMCInfo.mass = (primary) ? primary->Mass : 0;
      fileMCInfo.momentum[0] = FLT_MAX;
      fileMCInfo.momentum[1] = 0;
      fileMCInfo.pid_datacard = Tools::GetValueDatacard(file, "PART");
      fileMCInfo.mom_datacard[0] = Tools::GetValueDatacard(file, "PMIN");
      fileMCInfo.mom_datacard[1] = Tools::GetValueDatacard(file, "PMAX");
      fileMCInfo.ngen_datacard = int(Tools::GetValueDatacard(file, "TRIG"));
      if      (current_filename.find(".208000") != std::string::npos) { fileMCInfo.mom_filename[0] =   20; fileMCInfo.mom_filename[1] =  8000; }
      else if (current_filename.find(".2800")   != std::string::npos) { fileMCInfo.mom_filename[0] =    2; fileMCInfo.mom_filename[1] =   800; }
      else if (current_filename.find(".8032000") != std::string::npos) { fileMCInfo.mom_filename[0] =   80; fileMCInfo.mom_filename[1] = 32000; }
      else if (current_filename.find(".1800")   != std::string::npos) { fileMCInfo.mom_filename[0] =    1; fileMCInfo.mom_filename[1] =   800; }
      else if (current_filename.find(".808000") != std::string::npos) { fileMCInfo.mom_filename[0] =   80; fileMCInfo.mom_filename[1] =  8000; }
      else if (current_filename.find(".1010000") != std::string::npos) { fileMCInfo.mom_filename[0] =   10; fileMCInfo.mom_filename[1] = 10000; }
      else if (current_filename.find(".1212000") != std::string::npos) { fileMCInfo.mom_filename[0] =   12; fileMCInfo.mom_filename[1] = 12000; }
      else if (current_filename.find(".1616000") != std::string::npos) { fileMCInfo.mom_filename[0] =   16; fileMCInfo.mom_filename[1] = 16000; }
      else if (current_filename.find(".1664000") != std::string::npos) { fileMCInfo.mom_filename[0] =   16; fileMCInfo.mom_filename[1] = 64000; }
      else if (current_filename.find(".24000")  != std::string::npos) { fileMCInfo.mom_filename[0] =    2; fileMCInfo.mom_filename[1] =  4000; }
      else if (current_filename.find(".216000") != std::string::npos) { fileMCInfo.mom_filename[0] =    2; fileMCInfo.mom_filename[1] = 16000; }
      else if (current_filename.find(".8032000") != std::string::npos) { fileMCInfo.mom_filename[0] =   80; fileMCInfo.mom_filename[1] = 32000; }
      else if (current_filename.find(".36000")  != std::string::npos) { fileMCInfo.mom_filename[0] =    3; fileMCInfo.mom_filename[1] =  6000; }
      else if (current_filename.find(".66000")  != std::string::npos) { fileMCInfo.mom_filename[0] =    6; fileMCInfo.mom_filename[1] =  6000; }
      else if (current_filename.find(".624000") != std::string::npos) { fileMCInfo.mom_filename[0] =    6; fileMCInfo.mom_filename[1] = 24000; }
      else if (current_filename.find(".48000")  != std::string::npos) { fileMCInfo.mom_filename[0] =    4; fileMCInfo.mom_filename[1] =  8000; }
      else if (current_filename.find(".832000") != std::string::npos) { fileMCInfo.mom_filename[0] =    8; fileMCInfo.mom_filename[1] = 32000; }
      else if (current_filename.find(".510000") != std::string::npos) { fileMCInfo.mom_filename[0] =    5; fileMCInfo.mom_filename[1] = 10000; }
      else if (current_filename.find(".1040000") != std::string::npos) { fileMCInfo.mom_filename[0] =   10; fileMCInfo.mom_filename[1] = 40000; }
      else if (current_filename.find(".612000") != std::string::npos) { fileMCInfo.mom_filename[0] =    6; fileMCInfo.mom_filename[1] = 12000; }
      else if (current_filename.find(".1248000") != std::string::npos) { fileMCInfo.mom_filename[0] =   12; fileMCInfo.mom_filename[1] = 48000; }
      else if (current_filename.find(".714000") != std::string::npos) { fileMCInfo.mom_filename[0] =    7; fileMCInfo.mom_filename[1] = 14000; }
      else if (current_filename.find(".1456000") != std::string::npos) { fileMCInfo.mom_filename[0] =   14; fileMCInfo.mom_filename[1] = 56000; }
      else if (current_filename.find(".816000") != std::string::npos) { fileMCInfo.mom_filename[0] =    8; fileMCInfo.mom_filename[1] = 16000; }
      else if (current_filename.find(".1664000") != std::string::npos) { fileMCInfo.mom_filename[0] =   16; fileMCInfo.mom_filename[1] = 64000; }
      else if (current_filename.find(".1200")   != std::string::npos) { fileMCInfo.mom_filename[0] =    1; fileMCInfo.mom_filename[1] =   200; }
      else if (current_filename.find(".0_5200") != std::string::npos) { fileMCInfo.mom_filename[0] =  0.5; fileMCInfo.mom_filename[1] =   200; }
      else if (current_filename.find(".400")    != std::string::npos) { fileMCInfo.mom_filename[0] =  400; fileMCInfo.mom_filename[1] =   400; }
      else if (current_filename.find(".0_510")  != std::string::npos) { fileMCInfo.mom_filename[0] =  0.5; fileMCInfo.mom_filename[1] =    10; }
      else if (current_filename.find(".10200")  != std::string::npos) { fileMCInfo.mom_filename[0] =   10; fileMCInfo.mom_filename[1] =   200; }
      else if (current_filename.find(".2016000") != std::string::npos) { fileMCInfo.mom_filename[0] =   20; fileMCInfo.mom_filename[1] = 16000; }
      else if (current_filename.find(".204000") != std::string::npos) { fileMCInfo.mom_filename[0] =   20; fileMCInfo.mom_filename[1] =  4000; }
      else if (current_filename.find(".0_252")  != std::string::npos) { fileMCInfo.mom_filename[0] = 0.25; fileMCInfo.mom_filename[1] =     2; }
      fileMCInfo.isl1_filename = (current_filename.find(".l1." ) != std::string::npos);
      fileMCInfo.isl19_filename = false;
      if ( (current_filename.find(".l19." ) != std::string::npos) ||
           (current_filename.find(".l1a9.") != std::string::npos) ) fileMCInfo.isl19_filename = true;
      fileMCInfo.istb_filename = (current_filename.find(".tb.") != std::string::npos);
    }
  }
  // process current entry
  fileInfo.event[1] = Event();
  fileInfo.utime[1] = fHeader.Time[0];
  fileInfo.nentries++;
  if (isMC) {
    if (Event() < fileMCInfo.event[0]) fileMCInfo.event[0] = Event();
    if (Event() > fileMCInfo.event[1]) fileMCInfo.event[1] = Event();
    MCEventgR* primary = (nMCEventg() != 0) ? pMCEventg(0) : 0;
    float p = (primary) ? primary->Momentum : 0;
    if (p < fileMCInfo.momentum[0]) fileMCInfo.momentum[0] = p;
    if (p > fileMCInfo.momentum[1]) fileMCInfo.momentum[1] = p;
  }
}

void Analysis::ProcessRTI() {
  AMSSetupR::RTI rti;
  GetRTI(rti);
  if (current_rti_time != rti.utime) {
    FillRTI();
    FillNtpSHeader();
    pTreeRTI->Fill();
    current_rti_time = rti.utime;
  }
}


bool Analysis::IsCompactStandalone() {
  // events with TOF and TRD in fiducial volume
  if ( (pBetaH_SA)&&(pTrdTrack_SA) ) {
    bool good_beta = (pBetaH_SA->GetUseHit()==4)&&(pBetaH_SA->GetNormChi2T()<10);
    bool good_trd = (pTrdTrack_SA->Chi2<10)&&(pBetaH_SA->GetBuildType()==2);
    int trd_fiducial_pattern = GetPatternInsideTracker(pTrdTrack_SA);
    bool trd_is_inside = (trd_fiducial_pattern&0xff)==0xff;
    return good_beta&&good_trd&&trd_is_inside;
  }
  return false;
}

bool Analysis::IsCompactAnalysis() {
  // events with a L1+Inner good pattern track in fiducial volume
  if (pTrTrack) {
    int trk_patty = pTrTrack->GetBitPatternJ();
    bool good_trk_inn = ((trk_patty&0x2)!=0)&&((trk_patty&0xc)!=0)&&((trk_patty&0x30)!=0)&&((trk_patty&0xc0)!=0);
    int fit_id = pTrTrack->iTrTrackPar(1,3,refi_ref_fit,mass_ref_fit,char_ref_fit); // Choutko Inner
    int trk_inn_fiducial = GetPatternInsideTracker(pTrTrack,fit_id);
    bool trk_inn_is_inside = (trk_inn_fiducial&0xff)==0xff;
    return good_trk_inn&&trk_inn_is_inside;
  }
  return false;
}

void Analysis::DefineStandAlone() {
  // 0) Do a rebuild allowing recon. without TrTrack
  // 0-normal Build
  // 1-Exclude Track-Association Build
  // 2-Exclude Track+Trd Association
  // 3-Build Continue After Track Association Finding
  // 11110 Trk-Trd-Ecal-TOF-Other
  // 31110 Trk(Track Find Continue)-Trd-Ecal-TOF-Other
  // -1 Calibration Build.
  TofRecH::BuildOpt = 31110; 
  TofRecH::ReBuild();
  // 1) Best BetaH with TrdTrack
  double min_chisqtn = 1e+30;
  for (int ibetah = 0; ibetah<nBetaH(); ibetah++) {
    BetaHR* betah = AMSEventR::pBetaH(ibetah);
    if ( (betah->GetBuildType()!=2)||(!betah->pTrdTrack()) ) continue;
    if (betah->NTofClusterH()<3) continue;
    if (betah->GetNormChi2T()<min_chisqtn) {
      pBetaH_SA = betah;
      pTrdTrack_SA = pBetaH_SA->pTrdTrack();
      min_chisqtn = betah->GetNormChi2T();
    }
  }
  // 2) if no TrdTrack, take best BetaH without TrTrack
  if (!pBetaH_SA) {
    min_chisqtn = 1e+30;
    for (int ibetah = 0; ibetah<nBetaH(); ibetah++) {
      BetaHR* betah = AMSEventR::pBetaH(ibetah);
      if (betah->GetBuildType()==1) continue;
      if (betah->NTofClusterH()<3) continue;
      if (betah->GetNormChi2T()<min_chisqtn) {
        pBetaH_SA = betah;
        min_chisqtn = betah->GetNormChi2T();
      }
    }
  }
  // 3) Tuning
  if (pBetaH_SA) if (isMC) pBetaH_SA->DoMCtune();
}

static bool got_it = false;

bool Analysis::UProcessCut() {

  if (singleevent>=0) {
    if (Event()!=singleevent) return kFALSE;
    else cout << "Analysis::UProcessCut-Found-Event: " << Event() << endl;
  }

  if (Entry()<firstentry) return kFALSE;
  if (!got_it) cout << "Analysis::UProcessCut-First-Event-Processed Entry:" << Entry() << " Run:" << Run() << " Event:" << Event() << endl;
  got_it = true;

  ClearSelection();

  accepted       = false;
  a_track        = false;
  a_betah        = false;
  is_associated  = false;
  is_downgoing   = false;
  is_positive    = false;
  is_tofh        = false;
  is_rich        = false;
  is_in_lay1     = false;
  is_zgt1        = false;
  is_zgt2        = false;
  has_tkl1       = false;
  is_naf         = false;

  double trck_rig   = 0;
  double tofh_beta  = 0;
  double rich_beta  = 0;
  double trck_inn_q = 0;

  clock_end = clock();
  procInfo.time[1] += float(clock_end - clock_start) / CLOCKS_PER_SEC;
  clock_start = clock();
  procInfo.nevents[0]++;

  // report
  current_entry++;
  if ((current_entry % 10000) == 0) {
    cout << "Processed entry " << current_entry << ", run " << Run() << " event " << Event() << endl;
  }

  // run
  ProcessFile();

  // RTI
  if (!isMC) ProcessRTI();

  // cut really bad stuff
  if ( (nDaqEvent() < 1) || (nLevel1() < 1) ) {
    clock_end = clock();
    procInfo.time[2] += float(clock_end - clock_start) / CLOCKS_PER_SEC;
    clock_start = clock();
    accepted = false;
  } else {
    accepted = true;

    procInfo.nevents[1]++;

    ////////////////////////////////////////////////////
    // Selection
    ////////////////////////////////////////////////////

    // 1) select the particle with a track
    // - drop particles without a beta with less than 3 points
    // - drop particles in which pointed track in BetaH is different from track pointed in Particle
    // - if more than one track select the one with the maximum number of inner tracker XY hits
    // - if the tracks have the same XY number of hits then look for the most Y clusters
    int ipart_select = -1;
    int ibetah_select = -1;
    int max_hits = 0;
    for (int ipart = 0; ipart < NParticle(); ipart++) {
      ParticleR* particle = pParticle(ipart);
      if (!particle->pBetaH()) continue;
      if (particle->pBetaH()->NTofClusterH()<3) continue;
      if (!particle->pTrTrack()) continue;
      TrTrackR* track = particle->pTrTrack();
      if (particle->pBetaH()->pTrTrack() != track) continue;
      int nxy, ny;
      Tools::GetInnerNHits(track,nxy,ny);
      if ((ny*10+nxy)>max_hits) {
        ipart_select = ipart;
        max_hits = ny*10+nxy;
      }
    }
    if (ipart_select >= 0) {
      pBeta       = pParticle(ipart_select)->pBeta();
      pBetaH      = pParticle(ipart_select)->pBetaH();
      pTrTrack    = pParticle(ipart_select)->pTrTrack();
      pTrdTrack   = pParticle(ipart_select)->pTrdTrack();
      pTrdHTrack  = pParticle(ipart_select)->pTrdHTrack();
      pRichRing   = pParticle(ipart_select)->pRichRing();
      pRichRingB  = pParticle(ipart_select)->pRichRingB();
      pEcalShower = pParticle(ipart_select)->pEcalShower();
      if ((pBetaH) && (pBetaH->pTrTrack() != pTrTrack)) cout << "TEST KO" << endl;
    }
    else {
      // 2) if there is no particle, then look for the track using BetaH container
      int max_hits = 0;
      for (int ibetah = 0; ibetah < nBetaH(); ibetah++) {
        BetaHR* betah = AMSEventR::pBetaH(ibetah);
        if (betah->NTofClusterH()<3) continue;
        if (!betah->pTrTrack()) continue;
        TrTrackR* track = betah->pTrTrack();
        int nxy, ny;
        Tools::GetInnerNHits(track, nxy, ny);
        if ((ny*10+nxy)>max_hits) {
          ibetah_select = ibetah;
          max_hits = ny*10+nxy;
        }
      }
      if (ibetah_select >= 0) {
        pBetaH      = AMSEventR::pBetaH(ibetah_select);
        pTrTrack    = pBetaH->pTrTrack();
        pTrdTrack   = pBetaH->pTrdTrack();
        for (int ibeta = 0; ibeta < nBeta();      ibeta++) if (AMSEventR::pBeta     (ibeta)->pTrTrack() == pTrTrack) pBeta      = AMSEventR::pBeta     (ibeta);
        for (int irich = 0; irich < nRichRing();  irich++) {
          if (AMSEventR::pRichRing(irich)->pTrTrack() == pTrTrack) {
            if (!pRichRing) pRichRing = AMSEventR::pRichRing(irich);
            else if (!(pRichRing->IsClean())) pRichRing = AMSEventR::pRichRing(irich);
          }
        }
        for (int irich = 0; irich < nRichRingB(); irich++) if (AMSEventR::pRichRingB(irich)->pTrTrack() == pTrTrack) pRichRingB = AMSEventR::pRichRingB(irich);
        pEcalShower = pBetaH->pEcalShower();
      }
    }
    // 3) if there is still no track, loop on TrTrack container (this case is needed for beta efficiency)
    if (!pTrTrack) {
      int max_hits = 0;
      for (int itrtrack = 0; itrtrack < nTrTrack(); itrtrack++) {
        TrTrackR* track = AMSEventR::pTrTrack(itrtrack);
        int nxy, ny;
        Tools::GetInnerNHits(track, nxy, ny);
        if ((ny*10+nxy) > max_hits) {
          pTrTrack = track;
          max_hits = ny*10+nxy;
        }
      }
    }
    // 4) if there is still no beta loop on BetaH container
    // - presence of track is required, otherwise we will end up with wrong associations
    // - look for the best time chi2
    if ( (!pBetaH) && (!pTrTrack) ) {
      float min_chisqtn = 1e+30;
      for (int ibetah = 0; ibetah < nBetaH(); ibetah++) {
        BetaHR* betah = AMSEventR::pBetaH(ibetah);
        if (betah->NTofClusterH() < 3) continue;
        if (betah->GetNormChi2T() < min_chisqtn) {
          min_chisqtn = betah->GetNormChi2T();
          pBetaH = betah;
        }
      }
    }
    procInfo.nevents[2]++;

    // Tracker settings
    if (pTrTrack) {
      // spatial resolution corrections/improvements
      if (isMC) {
        SetDefaultMCTuningParameters();
        TrExtAlignDB::SmearExtAlign(); // Outer layer smearing
        TRCLFFKEY.UseSensorAlign = 0;  // Cancel inner sensors disalignment
        TRCLFFKEY.ClusterCofGOpt = 1;  // Mask of negative signal strip
        TRFITFFKEY.Zshift = -1;
      }
      else {
        TrLinearEtaDB::SetLinearCluster(); // New Eta uniformity (only for ISS data)
        TRFITFFKEY.Zshift = 2;
      }
      TRFITFFKEY.ErcHeY = 0; // Use realistic charge-based errors
    }

    // BetaH settings
    if ( (isMC)&&(pBetaH) ) {
      short particle = fabs(GetPrimaryMC()->Particle);
      if ( (particle==47)||(particle==49)||(particle==147)||(particle==149) ) {
        TofMCPar::MCtuneDT = -87; // 120; // MC Smear Par(ps) : 120ps, increase beta resolution from 0.03 to 0.04
        TofMCPar::MCtuneST =  10; // -40; // MC Shift Par(ps) : 40ps shift, shifting 1/beta ~ -0.01
        pBetaH->DoMCtune();
      }
    }

    ////////////////////////////////////////////////////
    // Prescaling condition build
    ////////////////////////////////////////////////////

    // beta and rigidity
    a_track        = (pTrTrack != 0);
    a_betah        = (pBetaH != 0);
    is_associated  = ((pBetaH) && (pTrTrack)) ? (pBetaH->pTrTrack() == pTrTrack) : false;
    tofh_beta      = (pBetaH) ? (isMC ? pBetaH->GetMCBeta() : pBetaH->GetBeta()) : 1;
    trck_inn_q     = (pTrTrack) ? pTrTrack->GetInnerQ_all(tofh_beta,0).Mean : 0;
    if ((pTrTrack)&&(trck_inn_q>1.5)) {  
      TrQYJTrack trck_q_yj = pTrTrack->GetQYJ_all(2,tofh_beta,0);
      if (trck_q_yj.InnerQ>8.5) trck_inn_q = trck_q_yj.InnerQ;
    }
    char_ref_fit   = floor(0.5 + trck_inn_q);
    if (char_ref_fit<1) char_ref_fit = 1;
    mass_ref_fit   = (char_ref_fit>=2) ? TrFit::Mhelium/2*char_ref_fit : TrFit::Mproton;
    refi_ref_fit   = 23;
    if ( (!isMC)&&(pTrTrack)&&(char_ref_fit==int(pTrTrack->GetAdvancedFitCharge()+0.5)) ) refi_ref_fit = 21;
    id_ref_fit     = (pTrTrack) ? pTrTrack->iTrTrackPar(algo_ref_fit,patt_ref_fit,refi_ref_fit,mass_ref_fit,char_ref_fit) : -1;
    trck_rig       = (id_ref_fit >= 0) ? pTrTrack->GetRigidity(id_ref_fit) : 0;
    rich_beta      = (pRichRing) ? pRichRing->getBeta() : 0;
    is_downgoing   = (tofh_beta > 0);
    is_positive    = ((trck_rig > 0) && (tofh_beta>=0))|| ((trck_rig < 0) && (tofh_beta<0));
    is_tofh        = (fabs(tofh_beta) > 0);
    is_rich        = (fabs(rich_beta) > 0);
    is_naf         = (pRichRing) ? pRichRing->IsNaF() : 0;
    is_zgt1        = (trck_inn_q > 1.8);
    is_zgt2        = (trck_inn_q > 2.7);
    has_tkl1       = (pTrTrack) ? (pTrTrack->GetHitLJ(1) !=0)   : false;

    // L1 geometry
    AMSPoint lay1_point(0, 0, 0);
    AMSDir lay1_dir;
    double z_lay1 = tracker_layers_z[0];
    double lay1_time;
    is_in_lay1 = false;
    if (pTrTrack) {
      pTrTrack->Interpolate(z_lay1, lay1_point, lay1_dir);
      is_in_lay1 = Tools::IsInsideTkL1(lay1_point[0],lay1_point[1]);
    }
    else if (pTrdTrack) {
      pTrdTrack->Interpolate(z_lay1, lay1_point, lay1_dir);
      is_in_lay1 = Tools::IsInsideTkL1(lay1_point[0],lay1_point[1]);
    }
    else if (pBetaH) {
      pBetaH->TInterpolate(z_lay1, lay1_point, lay1_dir, lay1_time, false);
      is_in_lay1=Tools::IsInsideTkL1(lay1_point[0],lay1_point[1]);
    } else {
      is_in_lay1=false;
    }
  }

  // pre-scaling
  long int ret = vPrescale.PrescaleEvent(*this);

  // compact
  if ( ( (ret<0)&&(pres_strategy>0) )||(compact==2) ) {
    if ( (pTreeComp)&&(!disable_fill) ) {
      ClearCompact();
      FillNtpSHeader();
      bool is_analysis = IsCompactAnalysis();
      if (is_analysis) FillNtpCompact(false);
      DefineStandAlone();  
      bool is_standalone = IsCompactStandalone();     
      FillNtpCompactStandAlone(false);
      if ( (is_analysis||is_standalone) ) pTreeComp->Fill();
    }
    clock_end = clock();
    procInfo.time[2] += float(clock_end - clock_start) / CLOCKS_PER_SEC;
    clock_start = clock();
    return false;
  }

  long int ret2 = (ret<0) ? ret*-1 : ret;
  ntpHeader.pres_weight     = vPrescale.FindCat(ret2)->GetPrf();
  ntpHeader.pres_trck_rig   = trck_rig;
  ntpHeader.pres_tofh_beta  = tofh_beta;
  ntpHeader.pres_rich_beta  = rich_beta;
  ntpHeader.pres_trck_inn_q = trck_inn_q;
  ntpHeader.pres_patt       = ret2;
  procInfo.nevents[3]++;

  // look for the track with the largest momentum that is not the choosen track
  double max_mom = 0;
  for (int itrtrack = 0; itrtrack < nTrTrack(); itrtrack++) {
    TrTrackR* track = AMSEventR::pTrTrack(itrtrack);
    if (pTrTrack == track) continue;
    float mom = fabs(track->GetInnerQ()*track->GetRigidity());
    if (mom > max_mom) {
      pTrTrack2 = track;
      max_mom = mom;
    }
  }
  procInfo.nevents[4]++;
  clock_end = clock();
  procInfo.time[2] += float(clock_end - clock_start)/CLOCKS_PER_SEC;
  return true;
}

void Analysis::UProcessFill() {

  procInfo.nevents[5]++;

  if (disable_fill) return;

  ////////////////////////////////////////////////////
  // Event
  ////////////////////////////////////////////////////

  clock_start = clock();
  ClearNtp();
  clock_end = clock();
  procInfo.time[3] += float(clock_end - clock_start) / CLOCKS_PER_SEC;

  clock_start = clock();
  FillNtpSHeader();
  FillNtpHeader();
  clock_end = clock();
  procInfo.time[4] += float(clock_end - clock_start) / CLOCKS_PER_SEC;

  clock_start = clock();
  if (isMC) FillNtpMCHeader();
  clock_end = clock();
  procInfo.time[5] += float(clock_end - clock_start) / CLOCKS_PER_SEC;

  clock_start = clock();
  FillNtpTrd();
  clock_end = clock();
  procInfo.time[6] += float(clock_end - clock_start) / CLOCKS_PER_SEC;

  clock_start = clock();
  FillNtpTof();
  clock_end = clock();
  procInfo.time[7] += float(clock_end - clock_start) / CLOCKS_PER_SEC;

  clock_start = clock();
  FillNtpTracker();
  clock_end = clock();
  procInfo.time[8] += float(clock_end - clock_start) / CLOCKS_PER_SEC;

  clock_start = clock();
  FillNtpRich();
  clock_end = clock();
  procInfo.time[9] += float(clock_end - clock_start) / CLOCKS_PER_SEC;

  clock_start = clock();
  FillNtpTrackerScattering();
  clock_end = clock();
  procInfo.time[8] += float(clock_end - clock_start) / CLOCKS_PER_SEC;

  clock_start = clock();
  FillNtpEcal();
  clock_end = clock();
  procInfo.time[10] += float(clock_end - clock_start) / CLOCKS_PER_SEC;

  clock_start = clock();
  FillNtpAnti();
  clock_end = clock();
  procInfo.time[11] += float(clock_end - clock_start) / CLOCKS_PER_SEC;

  clock_start = clock();
  ClearCompact();
  FillNtpSHeader();
  bool is_analysis = IsCompactAnalysis();
  if (is_analysis) FillNtpCompact(true);
  clock_end = clock();
  procInfo.time[13] += float(clock_end - clock_start) / CLOCKS_PER_SEC;

  ////////////////////////////////////////////////////
  // Stand-alone
  ////////////////////////////////////////////////////

  clock_start = clock();
  DefineStandAlone();
  FillNtpStandAlone();
  clock_end = clock();
  procInfo.time[12] += float(clock_end - clock_start) / CLOCKS_PER_SEC;

  clock_start = clock();
  bool is_standalone = IsCompactStandalone();
  FillNtpCompactStandAlone(true);
  clock_end = clock();
  procInfo.time[13] += float(clock_end - clock_start) / CLOCKS_PER_SEC;

  ////////////////////////////////////////////////////
  // Fill
  ////////////////////////////////////////////////////

  clock_start = clock();

  if ( (is_analysis||is_standalone)&&(pTreeComp) ) pTreeComp->Fill(); 
  if (pTreeEvent) pTreeEvent->Fill(); 
  clock_end = clock();
  procInfo.time[14] += float(clock_end - clock_start) / CLOCKS_PER_SEC;

  clock_start = clock();
}

void Analysis::UTerminate() {
  clock_end = clock();
  procInfo.time[1] += float(clock_end - clock_start) / CLOCKS_PER_SEC;
  printf("\n\n Analysis::UTerminate =============================================\n");
  procInfo.utctime[1] = time(0);
  printf(" === stop: %d\n", procInfo.utctime[1]);
  printf(" === pTreeFile->Show \n");
  pTreeFile->Fill();
  pTreeFile->Show(pTreeFile->GetEntries() - 1);
  printf(" === pTreeProc->Show \n");
  pTreeProc->Fill();
  pTreeProc->Show(pTreeProc->GetEntries() - 1);
  if (pTreeEvent) printf(" === pTreeEvent->GetEntries: %d\n",(int)pTreeEvent->GetEntries());
  if (pTreeComp) printf(" === pTreeComp->GetEntries: %d\n",(int)pTreeComp->GetEntries());
  // save
  if (pFileOutput) {
    pFileOutput->cd();
    if (pTreeRTI) {
      int ret = pTreeRTI->BuildIndex("SHeader.utime");
      printf("Producing Index for rti size %d\n", ret);
      pTreeRTI->Write();
    }
    if (pTreeFile)  pTreeFile ->Write();
    if (pTreeEvent) pTreeEvent->Write();
    if (pTreeProc)  pTreeProc->Write();
    if (pTreeComp)  pTreeComp->Write();
    vPrescale.WriteArchive(pFileOutput);
    if (pRichOcc)   pRichOcc->Write();
    pFileOutput->Close();
  }
  if (pres_strategy==0)
    printf("PRESCALING IN DUMMY MODE (pres_strategy==0)  NO PRESCALING ACTUALLY APPLIED!!\n");
  vPrescale.Report();
  if (pres_strategy==0)
    printf("\n\nPRESCALING IN DUMMY MODE (pres_strategy==0)  NO PRESCALING ACTUALLY APPLIED!!\n\n");
}

#include "RichCharge.h"
#include "richidOff.h"
#include "richradidOff.h"
#include "richrecOff.h"
void RichBDTData::FillDataDirect(Analysis *event) {

  Richtotused = event->nRichHit() - event->pRichRing->getUsedHits();
  if (event->pRichRing->getPhotoElectrons(false) > 0)
    RichPhEl = event->pRichRing->getExpectedPhotoElectrons(false) / event->pRichRing->getPhotoElectrons(false);
  else
    RichPhEl = 0;
  RICHprob = event->pRichRing->getProb();
  if (RichHitR::getCollectedPhotoElectrons() > 0)
    RICHcollovertotal = event->pRichRing->getPhotoElectrons(false) / RichHitR::getCollectedPhotoElectrons();
  else
    RICHcollovertotal = 0;

  // getting beta_corrected is a pain, better to take it from the compact tree
  // for now...
  RICHLipBetaConsistency = event->pRichRingB ? fabs(event->pRichRingB->Beta - event->ntpCompact.rich_beta) : 0;
  if (event->ntpCompact.rich_beta > 0)
    RICHTOFBetaConsistency =
        fabs(event->ntpCompact.rich_beta - (event->isMC ? event->pBetaH->GetMCBeta() : event->pBetaH->GetBeta())) /
        event->ntpCompact.rich_beta;
  else
    RICHTOFBetaConsistency = 0;
  RICHChargeConsistency = event->pRichRing->getPMTChargeConsistency();
  RICHPmts = event->pRichRing->NpColPMT.size();
  RICHgetExpected = event->pRichRing->getExpectedPhotoElectrons(false);

  // float Bad_ClusteringRICH;
  float beta_correction = (event->pRichRing->Beta > 0) ? event->ntpCompact.rich_beta / event->pRichRing->Beta : 1;
  auto nclus = event->pRichRing->ClusterizeZ1();
  int clus_size[10];
  float clus_mean[10];
  float clus_rms[10];
  Bad_ClusteringRICH = 0;
  if (nclus > 0) {
    for (int iclu = 0; iclu < TMath::Min(10, nclus); iclu++) {
      event->pRichRing->GetClusters(iclu, clus_size[iclu], clus_mean[iclu], clus_rms[iclu]);
      clus_mean[iclu] *= beta_correction;
      clus_rms[iclu] *= beta_correction;
      if ((clus_mean[iclu] - event->ntpCompact.rich_beta) > 0.01)
        Bad_ClusteringRICH++;
    }
  }

  // it seemed so easy until this point...
  // =========================
  // ==== HERE BE DRAGONS ====
  // =========================
  float tot_hyp_hit_uncorr[2] = {0, 0};
  tot_hyp_p_uncorr = 0;
  RichPMTCalib *pmt_calib = RichPMTCalib::getHead();
  vector<RichHitR *> rich_hit_on_ring_list;
  map<int, RichHitR *> rich_hit_offline_map;
  for (int ihit = 0; ihit < event->pRichRing->getUsedHits(false); ihit++)
    rich_hit_on_ring_list.push_back(event->pRichHit(event->pRichRing->iRichHit(ihit)));
  int which = 0;
  for (RichOffline::RichRawEvent *hit = new RichOffline::RichRawEvent(event); hit; hit = hit->next(), which++) {
    RichHitR *rich_hit = hit->getpointer();
    rich_hit_offline_map.insert(pair<int, RichHitR *>(which, rich_hit));
  }

  class pmt_info {
  public:
    int pmt;
    int nhit;
    float np;
    float coo[3];
    float dist;
    pmt_info() : pmt(-1), nhit(0), np(0), coo{0, 0, 0}, dist(0) {}
    void print() const {
      printf("pmt_info::print pmt:%3d nhit:%2d npe:%7.2f "
             "(x,y,z)=(%7.2f,%7.2f,%7.2f) dist=%7.2f\n",
             pmt, nhit, np, coo[0], coo[1], coo[2], dist);
    }
  };
  map<int, pmt_info> pmt_info_map;
  vector<int> crossed_pmt[2]; // primary/secondary
  for (int ihit = 0; ihit < event->nRichHit(); ihit++) {
    RichHitR *hit = event->pRichHit(ihit);
    int pmt = int(hit->Channel / 16);
    if ((pmt_calib) && (find(pmt_calib->BadPMTs.begin(), pmt_calib->BadPMTs.end(), pmt) != pmt_calib->BadPMTs.end()))
      continue;
    float np = hit->Npe;
    pmt_info_map[pmt].pmt = pmt;
    pmt_info_map[pmt].nhit++;
    pmt_info_map[pmt].np += np;
    for (int icoo = 0; icoo < 3; icoo++)
      pmt_info_map[pmt].coo[icoo] += np * hit->Coo[icoo];
  }
  for (auto it : pmt_info_map) {
    for (int icoo = 0; icoo < 3; icoo++)
      it.second.coo[icoo] /= it.second.np;
    if (event->pTrTrack) {
      AMSPoint coo;
      AMSDir dir;
      event->pTrTrack->Interpolate(it.second.coo[2], coo, dir, event->id_ref_fit);
      it.second.dist = sqrt(pow(it.second.coo[0] - coo.x(), 2.) + pow(it.second.coo[1] - coo.y(), 2.));
    }
    if (it.second.np < 5)
      continue;
    crossed_pmt[(it.second.dist < 3.5) ? 0 : 1].push_back(it.first);
  }

  float pmt_np_uncorr[5];
  float pmt_dist[5];
  typedef std::function<bool(pair<int, pmt_info>, pair<int, pmt_info>)> comparator;
  comparator functor = [](pair<int, pmt_info> elem1, pair<int, pmt_info> elem2) {
    return elem1.second.np > elem2.second.np;
  };
  set<pair<int, pmt_info>, comparator> pmt_info_set(pmt_info_map.begin(), pmt_info_map.end(), functor);
  int index = 0;
  for (auto it : pmt_info_set) {
    if (index > 4)
      continue;
    pmt_np_uncorr[index] = it.second.np;
    pmt_dist[index] = it.second.dist;
    index++;
  }
  NSecondariesRICHrich = 0;
  for (int is = 0; is < 5; is++)
    if ((pmt_np_uncorr[is] > 5) && (pmt_dist[is] > 3.5))
      NSecondariesRICHrich++;

  // Hit infos (including beta-per-hit)
  class hit_info {
  public:
    unsigned int status;  // Carlos status, bit 30 crossed
    unsigned int status2; // bit0: Jorge's bad PMT; bit1: occupancy table; bit2:
                          // crossed by primary (my algo.); bit3: crossed by
                          // secondary (my algo.); bit4: is not in selected ring
    int used;             // Used in whatever ring
    int channel;          // PMT*16+pixel
    float np;             // Number of photons (uncorrected? corrected for gain?)
    float beta[2];        // Beta hypothesys if existing
    hit_info(unsigned int s, unsigned int s2, int u, int c, float n, float b0, float b1)
        : status(s), status2(s2), used(u), channel(c), np(n) {
      beta[0] = b0;
      beta[1] = b1;
    }
    void print() const {
      printf("hit_info::print status:%10d status2:%3d used:%2d channel:%5d "
             "np:%7.2f beta(%7.4f,%7.4f)\n",
             status, status2, used, channel, np, beta[0], beta[1]);
    }
  };
  for (int i = 0; i < event->pRichRing->RawBetas(); i++) {
    int ihit = event->pRichRing->HitBeta(i);
    RichHitR *rich_hit = rich_hit_offline_map[ihit];
    if (!rich_hit)
      continue;
    int pmt = int(rich_hit->Channel / 16);
    bool is_bad_pmt =
        (pmt_calib) && (find(pmt_calib->BadPMTs.begin(), pmt_calib->BadPMTs.end(), pmt) != pmt_calib->BadPMTs.end());
    bool is_bad_occupancy = (event->pRichOcc) ? (!event->pRichOcc->IsGood(rich_hit->Channel)) : true;
    bool is_crossed_pri = find(crossed_pmt[0].begin(), crossed_pmt[0].end(), pmt) != crossed_pmt[0].end();
    bool is_crossed_sec = find(crossed_pmt[1].begin(), crossed_pmt[1].end(), pmt) != crossed_pmt[1].end();
    bool is_not_in_selected_ring =
        (find(rich_hit_on_ring_list.begin(), rich_hit_on_ring_list.end(), rich_hit) == rich_hit_on_ring_list.end());
    int status2 = ((is_bad_pmt) ? 0x1 : 0) + ((is_bad_occupancy) ? 0x2 : 0) + ((is_crossed_pri) ? 0x4 : 0) +
                  ((is_crossed_sec) ? 0x8 : 0) + ((is_not_in_selected_ring) ? 0x10 : 0);
    hit_info hit(rich_hit->Status, status2, event->pRichRing->UsedBeta(i), rich_hit->Channel, rich_hit->Npe,
                 event->pRichRing->RawBeta(i, 0), event->pRichRing->RawBeta(i, 1));
    // rich_hit_beta_hyp_list.push_back(rich_hit);
    // hit_info_list.push_back(hit);
    // integrate beta = 1 hypothesys with only good hits
    if ((rich_hit->Status >> 31) & 0x1)
      continue;
    if ((status2 & 0xf) != 0)
      continue;
    int kind_of_tile = 0;
    if (event->pRichRing->IsNaF())
      kind_of_tile = 1;
    float simple_resolution[2] = {1.2e-3, 4e-3};
    for (int ibeta = 0; ibeta < 2; ibeta++) {
      if (fabs(event->pRichRing->RawBeta(i, ibeta) - 1) > 3. * simple_resolution[kind_of_tile])
        continue;
      if ((event->pRichRing->UsedBeta(i) < 0) || (event->pRichRing->UsedBeta(i) > 1)) {
        tot_hyp_hit_uncorr[ibeta]++;
        tot_hyp_p_uncorr += rich_hit->Npe;
      }
    }
  }

  HitHVoutdir = tot_hyp_hit_uncorr[0];
  HitHVoutrefl = tot_hyp_hit_uncorr[0];
}
