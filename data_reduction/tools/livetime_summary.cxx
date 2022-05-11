#include "Ntp.h"

#include "GM_SubLibrary.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h" 

#include <string>
#include <iostream>
#include <fstream>
#include <cstdio>

using namespace std; 

static int type = 1; // use the so called IGRF-(RTI) 

string basename(const string& str) {
  size_t found = str.find_last_of("/\\");
  return str.substr(found+1);
} 

void add(TChain* chain, const char* runlist, int start, int end) {
  ifstream file(runlist);
  if (!file.is_open()) return; 
  int index = 0;
  while (!file.eof()) {
    string filename;
    getline(file,filename); 
    if ( (!filename.empty())&&(index>=start)&&(index<=end) ) { 
      chain->Add(filename.c_str()); 
      printf("Add %s\n",filename.c_str()); 
    } 
    index++; 
  }
  file.close(); 
}

int main(int argc, char ** argv) {
  if (argc<4) {
    printf("%s <Runlist> <Start> <End>\n",argv[0]);
    return 1; 
  }
  // open
  TChain* chain = new TChain("RTI");
  add(chain,argv[1],atoi(argv[2]),atoi(argv[3])); 
  RTIInfo* rti = new RTIInfo();
  if (chain->GetBranch("RTIInfo")) chain->SetBranchAddress("RTIInfo",&(rti));
  int nentries = chain->GetEntries();
  printf("%s::nentries=%d\n",argv[0],nentries);
  // Theta_{CGM} vs LiveTime 
  TH2F* thetam_vs_livetime[2][4][3] = {{{0}}};
  for (int ifov=0; ifov<4; ifov++) {
    for (int isign=0; isign<3; isign++) {
      for (int isel=0; isel<2; isel++) {
        thetam_vs_livetime[isel][ifov][isign] = new TH2F(
          Form("thetam_vs_livetime_isel%1d_ifov%1d_isign%1d",isel,ifov,isign),
          "; log_{10}(R_{c}/GV); #lambda_{CGM} [degrees]; T [s]",300,-1,2,30,-1.5,1.5
        );
      }
    }
  }
  for (int i=0; i<nentries; i++) {
    chain->GetEntry(i);
    if ((i%10000)==0) cout << i << "/" << nentries << endl;
    double deg2rad  = 1.74532925199432955e-02; 
    double Re       = 6371.2; // km Earth radius
    double Altitude = rti->r/1.e5-Re;
    double ThetaISS = rti->theta/deg2rad;
    double PhiISS   = rti->phi/deg2rad;
    time_t Utime    = rti->utime;
    float thetaM    = GM_GetThetaM(Utime,Altitude,ThetaISS,PhiISS);
    thetaM = thetaM*(-1)*deg2rad;
    double dt       = rti->lf*rti->nev/(rti->nev+rti->nerr); // corrected time
    for (int ifov=0; ifov<4; ifov++) {
      for (int isign=0; isign<3; isign++) {
        double Rc = (isign==2) ? 
          fabs(std::max(fabs(rti->cf[type][ifov][0]),fabs(rti->cf[type][ifov][1]))) :
          fabs(rti->cf[type][ifov][isign]);
        thetam_vs_livetime[0][ifov][isign]->Fill(log10(Rc),thetaM,dt);
        if (!rti->Select()) continue;
        if (rti->IsBadRun()) continue;
        thetam_vs_livetime[1][ifov][isign]->Fill(log10(Rc),thetaM,dt);
      }
    }
  }
  // save
  string runlist(argv[1]);      
  TFile* file = TFile::Open(Form("livetime_summary_%s_%06d_%06d.root",basename(runlist).c_str(),atoi(argv[2]),atoi(argv[3])),"recreate");
  for (int ifov=0; ifov<4; ifov++) {
    for (int isign=0; isign<3; isign++) {
      for (int isel=0; isel<2; isel++) {
        if (thetam_vs_livetime[isel][ifov][isign]) thetam_vs_livetime[isel][ifov][isign]->Write();
      } 
    }
  }
  file->Close();
  return 0;
}

