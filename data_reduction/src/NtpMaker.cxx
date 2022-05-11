#include "Analysis.h"
#include "RichOccupancy.h"
#include "RichOccupancyAnalysis.h"

#include "Tofrec02_ihep.h"
#include "TkDBc.h"
#include "TrExtAlignDB.h"
#include "TrCharge.h"
#include "TrRecon.h"
#include "RichCharge.h"
#include "RichConfig.h"
#include "TrdKCluster.h"
#include "AntiPG.h"
#include "MagField.h"
#include "GM_SubLibrary.h"
#include "bcorr.h"
#include "TrGainDB.h"
#include "HistoMan.h"
#include "RichBeta.h"
#include "RichTools.h"

#include "point.h"
#include "root.h"
#include "amschain.h"

#include "TTree.h"
#include "TFile.h"
#include "TString.h"

#include <cstdio>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>
#include <ctime>

using namespace std;

void usage() {
    std::cout
    << "\n\nUsage: NtpMaker -d <InputDir> -r <RunNumber> [-f <start>] [-n <end>] [-o <OutFile>] [-p <PrescalingStrategy>] [-e <Event>] [-R <RichOccupancy>] [-c <Compact>]\n"
    << "\n\t-d: Input data folder.                      Required."
    << "\n\t-r: Run number.                             Required.\n"
    << "\n\t-n: Number of events to process.            Optional, default 1000000000000."
    << "\n\t-f: Starting event number.                  Optional, default 0."
    << "\n\t-o: Output root file.                       Optional, default is <RunNumer>.root"
    << "\n\t-p: Prescaling.                             Optional: 1: 0 No; !=0 Yes. Default 1."
    << "\n\t-e: Select a single event to be processed.  Optional, default process all events."
    << "\n\t-R: Perform the RICH occupancy calculation. Optional: 0 No; !=0 Yes. Default 1."
    << "\n\t-c: Store the compact tree.                 Optional: 0 No; 1 Yes; 2 yes compact, no event tree. Default 1."
    << std::endl;
    exit(-1);
}

void ExpandFilelist(TString infilename, TChain* chain, TString basedir);


int main(int argc, char ** argv){

  ////////////////////
  // Process Arguments
  ////////////////////

  //vars
  char * dirname = (char*) "";
  char * ntupleRootFile = (char*) "";
  char * runlist = (char*) "";
  unsigned int run = 0;
  Long64_t maxentries = 1000000000000LL;
  Long64_t firstentry = 0;
  int pres_strategy = 1;
  Long64_t singleevent = -1;
  int occupancy = 1;
  int compact = 1;

  // Parsing with getopt
  char c;
  while ((c = getopt(argc, argv, "d:r:n:f:o:p:e:R:c:")) != -1) {
    if      (c == 'd') dirname = strdup(optarg);
    else if (c == 'r') runlist = strdup(optarg);
    else if (c == 'n') maxentries = atol(optarg);
    else if (c == 'f') firstentry = atol(optarg);
    else if (c == 'o') ntupleRootFile = strdup(optarg);
    else if (c == 'p') pres_strategy = atoi(optarg);
    else if (c == 'e') singleevent = atol(optarg);
    else if (c == 'R') occupancy = atoi(optarg);
    else if (c == 'c') compact = atoi(optarg);
  }

  // Checking required arguments
  if (strlen(dirname) == 0)  usage();
  if (strlen(runlist) == 0)  usage();

  Bool_t flag_filelist = kFALSE, flag_runnumber = kFALSE;
  if( ((TString) runlist).Contains(".txt") ){
    flag_filelist = kTRUE;
  }
  else if( ((TString) runlist).Contains(".root") ){
    flag_filelist = kFALSE;
  }
  else if( atoi(runlist) ){
    flag_filelist = kFALSE;
    flag_runnumber = kTRUE;
    run = atoi(runlist);
  } else {
    cout << "ERROR: not a valid inputfile" << endl;
    usage();
  }

  if (strlen(ntupleRootFile) == 0) { ntupleRootFile = strdup(Form("%d.root", run)); }
  bool is_test_beam = (strstr(dirname,"BT.")!=0);

  ////////////////////////////////////////////////////////////////////
  // settings
  ////////////////////////////////////////////////////////////////////

  if (!is_test_beam) TkDBc::UseFinal();

  TRMCFFKEY.init();
  TRMCFFKEY_DEF::ReadFromFile = 0;
  TRFITFFKEY_DEF::ReadFromFile = 0;
  TKGEOMFFKEY.MaxAlignedRun= 1527491433-1; // Updated on Oct 11, 2018

  TofRecH::RebuildBetaHInReadHeader = false;

  RichRingR::useEffectiveTemperatureCorrection = true; // Additional (effective) temperature corrections
  RichRingR::reloadRunTag = true;   // force load Config & Status from ext. files (needed for pass6 run>=1407139304 && run<=1411991495)

  // AntiRecoPG* Acci = AntiRecoPG::gethead();

  // chain
  AMSChain chain;

  TString basedir = dirname;
  if( basedir.EndsWith("/") ) basedir = basedir(0, basedir.Length()-1);

  if( flag_filelist ){
    ExpandFilelist( runlist, &chain, basedir );
  }
  else if( !flag_runnumber ){
    chain.Add( Form( "%s/%s", basedir.Data(), runlist ) );
  }
  else if( flag_runnumber ){
    chain.Add( Form( "%s/%d*root", basedir.Data(), run ) );
  }
  // ams.AddFromFile(runlist,first,last+1); // ,false,10);
  AMSEventR* event = chain.GetEvent(0);
  bool is_mc = (event->nMCEventg()>0);
  printf("Entries in the chain: %d\n",(int)chain.GetEntries());

  // RTI
  if (!is_mc) AMSSetupR::RTI::UseLatest(7);

  // TrdK stuff
  if (is_mc) {
    TrdKCluster::ForceReadAlignment=0;
    TrdKCluster::ForceReadCalibration=0;
    TrdKCluster::ForceReadXePressure=0;
    TrdKCluster::SetDefaultMCXePressure(900);
  }

  // run RICH occupancy creation
  RichOccupancyAnalysis rich_occupancy;
  if ((!is_mc)&&(occupancy!=0)) chain.Process(&rich_occupancy,ntupleRootFile); 

  // data reduct analysis
  Analysis analysis(is_mc);
  analysis.firstentry = firstentry;
  analysis.singleevent = singleevent;
  analysis.pres_strategy = pres_strategy;
  analysis.compact = compact;

  // set more RICH stuff
  TString richCorrFilePath = "root://eosams.cern.ch//eos/ams/group/dbar/data/"; // eos/ams/user/j/jorgec/betacorr/";
  analysis.pRichCorr = new GHBManager(richCorrFilePath + "dataHashBeta.root", richCorrFilePath + "mcHashBeta.root");
  TFile* hashFile = TFile::Open(Form("%s/v5.00/RichBetaUniformityCorrection.root",getenv("AMSDataDir")));
  analysis.pRichUnifAgl = (GeomHashEnsemble*) hashFile->Get("BetaAgl");
  analysis.pRichUnifNaf = (GeomHashEnsemble*) hashFile->Get("BetaNaF");
  analysis.pRichOcc = rich_occupancy.pOccupancy;
  hashFile->Close();

  // run data reduction
  chain.Process(&analysis,ntupleRootFile,maxentries);

  return 0;
}


void ExpandFilelist(TString infilename, TChain* chain, TString basedir){
  ifstream infilelist( infilename.Data() );
  TString bufname;
  if( !basedir.EndsWith("/") ) basedir += "/";
  while( !infilelist.eof() ){
    infilelist >> bufname;
    if( !bufname.CompareTo("") ) continue;
    if(bufname.Contains(".root")){
      cout << "Adding file to chain: " << bufname << endl;
      chain->Add( basedir + bufname );
    }
    else if(bufname.Contains(".txt")){
      ExpandFilelist(bufname, chain, basedir);
    }
  }
}
