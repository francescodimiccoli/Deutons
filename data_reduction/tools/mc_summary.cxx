#include "Ntp.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include <string>
#include <iostream>
#include <fstream>
#include <cstdio>

using namespace std; 

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
  TChain* chain = new TChain("File");
  add(chain,argv[1],atoi(argv[2]),atoi(argv[3]));
  FileInfo* fileInfo = new FileInfo();
  FileMCInfo* fileMCInfo = new FileMCInfo();
  if (chain->GetBranch("FileInfo"))   chain->SetBranchAddress("FileInfo",  &(fileInfo));
  if (chain->GetBranch("FileMCInfo")) chain->SetBranchAddress("FileMCInfo",&(fileMCInfo));
  int nentries = chain->GetEntries();
  printf("%s::nentries=%d\n",argv[0],nentries);
  string runlist(argv[1]);
  FILE* txtfile = fopen(Form("mc_summary_%s_%06d_%06d.txt",basename(runlist).c_str(),atoi(argv[2]),atoi(argv[3])),"w");
  for (int ifile=0; ifile<nentries; ifile++) {
    chain->GetEntry(ifile);
    if ((ifile%100)==0) printf("%d/%d\n",ifile,nentries);
    if ((ifile%1000)==0) cout << ifile << "/" << nentries << endl; 
    fprintf(txtfile,"%10d %10d %10d %10d %10d %10d %10d %2d %2d %10.4f %10.4f %10.4f %10.4f\n",
      fileInfo->run,fileInfo->nentries,
      fileInfo->event[0],fileInfo->event[1],
      fileMCInfo->event[0],fileMCInfo->event[1],
      fileMCInfo->ngen_datacard,
      fileMCInfo->pid_datacard,
      fileMCInfo->charge,
      fileMCInfo->mom_datacard[0],
      fileMCInfo->mom_datacard[1],
      fileMCInfo->momentum[0],
      fileMCInfo->momentum[1]
    );
  }
  fclose(txtfile);
  return 0;
}

