#include "particle.h"
#include "binning.h"
#include "histogram.h"
#include "spectrum.h"

// Preferred compilation command : g++ -g -std=c++11 $(root-config --libs --cflags)  -o GCRtest GCRtest.cpp GCR_data.cpp -fmax-errors=10 -Wno-pointer-arith


int main() {
   // Inits
   Particle proton(0.9382720813); // proton mass 938 MeV
   Particle deuton(0.18756129, 1, 2);   // deuterium mass 1876 MeV
   deuton.FillFromEk(5);
   Binning bins(proton);
   bins.setBinsFromEk(15, 2, 16);

   // Get normalized protons histogram
   Histogram histo;
   histo.fillWithGalpropFile("CRDB_ProtonsAMS_R.galprop");
   std::cout << " #####  AMS-02 proton file, raw histogram  ##### " << std::endl;
   histo.printContent();
   std::cout << " #####  AMS-02 proton file, normalized histogram  ##### " << std::endl;
   histo.normalize();
   histo.printContent();
   Spectrum spectrum(bins);
   spectrum.rebinHistoInRig(histo);
   std::cout << " #####  Migration matrix (transposed)  ##### " << std::endl;
   spectrum.printTMatrix();
   std::cout << " #####  AMS-02 proton file, rebinned  ##### " << std::endl;
   spectrum.printContentInRig();

   // Mock MC lognorm histo
   std::cout << " #####  AMS-02 generated lognorm spectrum ##### " << std::endl;
   Histogram MCgen;
   MCgen.genMCLogNormFlux();
   MCgen.normalize();
   MCgen.printContent();
   Spectrum mcspectrum(bins);
   mcspectrum.rebinHistoInRig(MCgen);
   mcspectrum.printContentInRig();
   

   // Flux ratio
   int size=bins.size();
   std::vector<float> MCOverP(size);
   for (int i=0; i<size; i++)
      MCOverP[i] = mcspectrum.getRebinnedHisto()[i]/spectrum.getRebinnedHisto()[i];
   Histogram hRatio;
   bool fillok=hRatio.fillWithVectors(bins.RigBins(), MCOverP);
   if (fillok) {
      std::cout << " #####  Normalized to unity ##### " << std::endl;
      hRatio.normalize();
      hRatio.printContent();
      std::cout << " #####  Normalized to the number of bins (average=1) ##### " << std::endl;
      hRatio.normalizeToAllBins();
      hRatio.printContent();
      std::cout << " #####  Normalized to the number of non-empty bins ##### " << std::endl;
      hRatio.normalizeToNonEmptyBins();
      hRatio.printContent();
   } else {
      std::cout << "Error fill" << std::endl;
   }
   

   return 0;
}
