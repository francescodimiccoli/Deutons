#include <string>
#include <vector>
#include <iostream>
#include <getopt.h>

void UsageAndExit() {
    std::cout << "Usage: ProduceHistograms -d inputData.root -m imputMC.root -o outputFile.root -c calibrations.root"
              << "\n\t-d: root file with DATA ntuples\n"
              << "\n\t-m: root file with MC ntuples\n"
              << "\n\t-o: output file with produced histograms\n"
              << "\n\t-c: file with calibrations\n"
              << std::endl;
    exit(-1);
}

//#include "TH1F.h"


void Loop(NTuple * ntuple, Distributions * dist) {
    TupleVars tup(ntuple);
    
    for(int i=0; i<ntuple->GetEntries(); i++) {
        ntuple->GetEvent(i);
        dist->Fill(tup);
    }
}


int main(int argc, char **argv){ 

    std::string outFname;
    std::string calFname;
    std::string inDataFname;
    std::string inMCFname;

    int c;
    while((c = getopt(argc, argv, "m:o:c:d:")) != -1) {
             if(c == 'o') outFname    = std::string(optarg);
        else if(c == 'c') calFname    = std::string(optarg);
        else if(c == 'd') inDataFname = std::string(optarg);
        else if(c == 'm') inMCFname   = std::string(optarg);
    }
    if(    outFname == "" ) UsageAndExit();
    if(    calFname == "" ) UsageAndExit();
    if( inDataFname == "" ) UsageAndExit();
    if(   inMCFname == "" ) UsageAndExit();

    std::cout << "Input DATA file: \""             << inDataFname << "\"\n";
    std::cout << "Input MC file: \""               << inMCFname   << "\"\n";
    std::cout << "Input Calibrations file: \""     << calFname    << "\"\n";
    std::cout << "Output file with histograms: \"" << outFname    << "\"\n";

    Calibrations calibs(calFname);
    AnalysisBins bins(20,20,20,20);

    Distributions * distTrigMC = Distributions::TriggerMC(bins, calibs);
    Distributions * distQualMC = Distributions::QualMC   (bins, calibs);

    TFile * fileTirgMC = TFile::Open(inMCFname.c_str(), "READ");
        TNtuple * ntupTrigMC = static_cast<TNtuple *>(fileTrigMC->Get("trig"));
        TNtuple * ntupQualMC = static_cast<TNtuple *>(fileTrigMC->Get("grandezzesepd"));
        Loop(ntupTrigMC, distTrigMC);
        Loop(ntupQualMC, distQualMC);
    fileTrigMC->Close();

    fileOut = TFile::Open(outFname.c_str(), "RECREATE");
        distTrigMC->Write(fileOut);
        distQualMC->Write(fileOut);
    fileOut->Close();

    return 0;
}
