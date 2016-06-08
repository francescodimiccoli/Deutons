#ifndef SPECTRUM_H
#define SPECTRUM_H

#include "binning.h"
#include "histogram.h"

class Spectrum
{
   public:
      Spectrum (Binning thosebins) : bins (thosebins) {};
      bool rebinHistoInRig (Histogram htorebin);
      bool rebinHistoInEtotPerMass (Histogram htorebin);
      std::vector<float> getRebinnedHisto() {return rebinnedVector;}
      Binning getBins() {return bins;}
         
   private:
         bool rebinHisto(Histogram htorebin, std::vector<float> binGoodUnit );
         void resetMatrix() { matrix.clear(); } ///< Reset the matrix transition
         Binning bins;
         void setMatrix (std::vector<float> vinput, std::vector<float> vref);
         void rebinVector (std::vector<float> vinput);
         std::vector< std::vector<float> > matrix; ///< Migration matrix for fluxes recorded and our binning.
         std::vector<float> rebinnedVector;

};





bool Spectrum::rebinHistoInRig (Histogram htorebin) 
{
   return rebinHisto(htorebin, bins.RigBins());
}


bool Spectrum::rebinHistoInEtotPerMass (Histogram htorebin) // Format of Galprop :/
{
   return rebinHisto(htorebin, bins.EtotPerMassBins());
}

bool Spectrum::rebinHisto(Histogram htorebin, std::vector<float> binGoodUnit) {
   if (! htorebin.checkIntegrity()) return false;
   resetMatrix();
   setMatrix( htorebin.getEdges(), binGoodUnit );
   rebinVector(htorebin.getContent());
}


void Spectrum::rebinVector (std::vector<float> vinput) {
   // What we really do is multiply a matrix by a vector
   int nmatlines=matrix[0].size();
   rebinnedVector.resize(nmatlines);
   for (uint ibin=0; ibin<vinput.size(); ibin++)
      for (uint obin=0; obin<nmatlines; obin++)
         rebinnedVector[obin] += matrix[ibin][obin] * vinput[ibin];
}



void Spectrum::setMatrix (std::vector<float> vinput, std::vector<float> vref)
{

   //if (matrix.size()!=0) return;     // Already done

   for (uint ib=0; ib<vinput.size()-1; ib++) { // prefix / suffix b for incoming binning
      float bmin=vinput[ib], bmax=vinput[ib+1];
      float bwidth=bmax-bmin;
      std::vector<float> column;

      for (uint it=0; it<vref.size()-1; it++) { // prefix / suffix b for This binning
         float tmin=vref[it], tmax=vref[it+1];
         float weight=0;
         // Which fraction of the bbin is in tbin? 4 possibilities:
         if      (tmin > bmax || tmax < bmin ) weight=0;                      // either tbin is outside bbin
         else if (tmin > bmin && tmax < bmax ) weight = (tmax-tmin) / bwidth; // or fully included into bbin
         else if (tmin <=bmin && tmax < bmax ) weight = (tmax-bmin) / bwidth; // or partly included to the right
         else if (tmin > bmin && tmax >=bmax ) weight = (bmax-tmin) / bwidth; // ...or to the left
         column.push_back (weight);
      }

      matrix.push_back (column);
   }
   return;
}



/*Class SpectrumWeithing
 Binning B
 LoadDataFromCSV(string filename)
 vector<float> generated_binsClass DataSpectrumWeighting: public SpectrumWeighting
 LoadDataFromGalprop(string filename)
 LoadDataFromCRDB(some options, maybe a date)Class MCSpectrumWeighting: public SpectrumWeighting
 vector <float> impingingbinsClass SpectrumWeightingMaker :
  Binning B
  MCSpectrumWeighting mcweight
  MCSpectrumWeighting Data
  float GetRatioForRigidity (R)*/


#endif /* SPECTRUM_H */ 




