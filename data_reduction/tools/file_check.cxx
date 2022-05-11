#include "TFile.h"
#include "TTree.h"

#include <iostream>

using namespace std; 

int main(int argc, char ** argv) {
  if (argc<2) return 1;  
  cout << argv[1];
  TFile* file = TFile::Open(argv[1],"read");
  if (!file) return 2;
  if (file->IsZombie()) return 3; 
  if (!file->Get("Processing")) { cout << endl; return 4; }
  cout << " " << ((TTree*)file->Get("Processing"))->GetEntries();
  if (!file->Get("File")) { cout << endl; return 5; } 
  cout << " " << ((TTree*)file->Get("File"))->GetEntries();
  if (!file->Get("Event")) { cout << endl; return 6; }
  cout << " " << ((TTree*)file->Get("Event"))->GetEntries();
  if (((TTree*)file->Get("Event"))->GetEntries()<=0) { cout << endl; return 7; }
  if (((TTree*)file->Get("Event"))->GetEntry(0)==0) { cout << endl; return 8; }
  if (!file->Get("Compact")) { cout << endl; return 9; }
  cout << " " << ((TTree*)file->Get("Compact"))->GetEntries();
  if (!file->Get("RTI")) { cout << endl; return 10; } 
  cout << " " << ((TTree*)file->Get("RTI"))->GetEntries(); 
  if (!file->Get("VPSArchive")) { cout << endl; return 11; } 
  if (!file->Get("RichOccupancy")) { cout << endl; return 12; }
  cout << endl; 
  return 0; // ok 
}
