#include "particle.h"
#include "binning.h"
#include "histogram.h"
#include "spectrum.h"

// Preferred compilation command : g++ -g -std=c++11 $(root-config --libs --cflags)  -o GCRtest GCRtest.cpp GCR_data.cpp -fmax-errors=10 -Wno-pointer-arith

int main() {
   Particle proton(0.9382720813); // proton mass 938 MeV
   Particle deuton(0.18756129, 1, 2);   // deuterium mass 1876 MeV
   deuton.FillFromEk(5);
   Binning bins(proton);
   bins.setBinsFromEk(15, 2, 16);
   Histogram histo;
   histo.fillWithGalpropFile("CRDB_ProtonsAMS_R.galprop");
   histo.printContent();
   std::cout << " ---------------------------------------" << std::endl;
   histo.normalize();
   histo.printContent();
   std::cout << " ---------------------------------------" << std::endl;
   Spectrum spectrum(bins);
   spectrum.rebinHistoInRig(histo);
   std::vector<float> flux=spectrum.getRebinnedHisto();
   for (int i=0; i<bins.size()-1; i++)
      std::cout << i << " " << bins.RigBinCent(i) << " " << flux[i] << std::endl;

   return 0;
}


void printEdgesContent(std::vector<float> bins, std::vector<float> content) {
   

   
}
