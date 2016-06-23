#ifndef REWEIGHT_H
#define REWEIGHT_H

#include "spectrum.h"

class Reweight {
   public:
      Reweight();
      Reweight(Binning binning) : bins(binning) { init(); };
      float MCToDataForRig(float rig);
      float DataToMCForRig(float rig);
      std::vector<float>  getMCOverData() {return MCOverData; }
      std::vector<float>  getDataOverMC() {return DataOverMC; }
      std::vector<float>  getRigBins()    {return bins.RigBinsCent(); }
      
   private:
      Binning bins;
      void init();
      void normalizeVectors();
      void getVectorsFromSpectra(Spectrum data, Spectrum mc);
      Spectrum loadDataSpectrum();
      Spectrum loadMCSpectrum();
      std::vector<float> MCOverData;
      std::vector<float> DataOverMC;
};


Reweight::Reweight() : bins(0.9382720813)
{
   bins.setBinsFromEk(43, 0.5, 100);
   init();
}



void Reweight::init()
{
   Spectrum dataspectrum=loadDataSpectrum();
   Spectrum mcspectrum = loadMCSpectrum();
   getVectorsFromSpectra(dataspectrum, mcspectrum);
   normalizeVectors();
   return;
}


Spectrum Reweight::loadDataSpectrum()
{
   Histogram hdata;
   hdata.fillWithGalpropFile("/storage/gpfs_ams/ams/users/fdimicco/Deutons/include/CRDB_ProtonsAMS_R.galprop");
   Spectrum dataspectrum(bins);
   dataspectrum.rebinHistoInRig(hdata);
   return dataspectrum;
}

Spectrum Reweight::loadMCSpectrum()
{
   Histogram MCgen;
   MCgen.genMCLogNormFlux(bins);
   Spectrum mcspectrum(bins);
   mcspectrum.rebinHistoInRig(MCgen);
   return mcspectrum;
}

void Reweight::getVectorsFromSpectra(Spectrum data, Spectrum mc) {
   uint size=bins.size();
   MCOverData.assign(size, 0);
   DataOverMC.assign(size, 0);

   // If division by 0 I chose to keep 0 in the histogram instead of "std::numeric_limits<double>::infinity();"
   // Because it might correspond to uncomplete spectrum, etc.
   // A change should be propagated to the histogram "normalize" methods.

   for (int i=0; i<bins.size(); i++) {
      float binData=data.getRebinnedHisto()[i];
      float binMC  =mc  .getRebinnedHisto()[i];
      if (binMC!=0)   DataOverMC[i]=binData/binMC;
      if (binData!=0) MCOverData[i]=binMC/binData;
   }
   return;
   
}

void Reweight::normalizeVectors() {
   Histogram hDataOMC;
   hDataOMC.fillWithVectors(bins.RigBins(), DataOverMC);
   hDataOMC.normalizeToNonEmptyBins();
   DataOverMC=hDataOMC.getContent();

   Histogram hMCOData;
   hMCOData.fillWithVectors(bins.RigBins(), MCOverData);
   hMCOData.normalizeToNonEmptyBins();
   MCOverData=hMCOData.getContent();

   return;
}


float Reweight::MCToDataForRig(float rig) {
   int bin=bins.GetRBin(rig);
   if (bin==-1) return 0;
   return DataOverMC[bin];
   // Yes ! 1 MC would correspond (in data evts) to the inverse of this ratio
}

float Reweight::DataToMCForRig(float rig) {
   int bin=bins.GetRBin(rig);
   if (bin==-1) return 0;
   return MCOverData[bin];
}


#endif
