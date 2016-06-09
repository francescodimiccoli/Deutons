#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include "GCR_data.h"
#include "printmatrix.h"
#include <string>
#include <algorithm>

/** A histogram is a combination of two vectors, edges and contents
 * Edges is sorted, contents has one element less (both are checked by CheckIntegrity)
 *
 * */

class Histogram {
   public:
      std::vector<float> getEdges()   {return edges  ;}
      std::vector<float> getContent() {return content;}

      bool checkIntegrity();
      bool fillWithVectors (std::vector<float> edge, std::vector<float> cont );
      void fillWithGalpropFile (std::string filename);
      void genMCLogNormFlux (int nbins=100, float rmin=0.5, float rmax=10);
      void normalize();
      void printContent();

   private:
      std::vector<float> fluxtoHisto (std::vector<float> inContent);
      bool checkIntegrity ( std::vector<float> edge, std::vector<float> cont );
      std::vector<float> edges ;
      std::vector<float> content ;
};





bool Histogram::checkIntegrity (std::vector<float> edge, std::vector<float> cont  ) {
   if ( !std::is_sorted (edge.begin(),edge.end() ) ) return false;
   if ( edge.size() != cont.size() +1 )              return false;
   return true;
}

bool Histogram::checkIntegrity() {
   return checkIntegrity (edges, content);
}

bool Histogram::fillWithVectors (std::vector<float> edge, std::vector<float> cont) {
   bool goodEntryVectors=checkIntegrity (edge, cont);
   if (goodEntryVectors) {
      edges=edge;
      content=cont;
   }
   return goodEntryVectors;
}


void Histogram::fillWithGalpropFile (std::string filename)
{
   GCR_data galprop;
   galprop.read ( (char *) filename.data(), ( char*) "m2", ( char*) "GeV");
   int nentries=galprop.n;
   std::vector<float> fluxContent (nentries);
   edges.assign(nentries, 0);

   for (int i=0; i<nentries; i++) {
      edges      [i] = galprop.E_low_input[i];
      fluxContent[i] = galprop.value_input[i];
   }

   edges.push_back (galprop.E_high_input[nentries-1]); // last bin
   content=fluxtoHisto (fluxContent);
   return;
}


/** In our jargon, a histo is an (integrated) number of events in a bin, as opposed to a flux (derivative, independent of binning)
 *  This function turns a flux and histo
 * */
std::vector<float> Histogram::fluxtoHisto (std::vector<float> inContent) {
   int nbins=inContent.size();
   std::vector<float> outContent (nbins);
   for (int i=0; i<nbins; i++)
      outContent[i]=inContent[i] * (edges[i+1] - edges[i]);
   return outContent;
}

/** MC AMS LogNorm flux has the same number of entries in a normalised binning in momentum.
 *  This function creates such a basic histogram
 * */
void Histogram::genMCLogNormFlux (int nbins, float rmin, float rmax) {
   Binning bins (1); // Mass not important, we only work in rigidity
   bins.setBinsFromRigidity (nbins, rmin, rmax);
   edges=bins.RigBins();
   content.assign (bins.size(), 1);
   return;
}

void Histogram::normalize() {
   float integral=0;
   for (int i=0; i<content.size(); i++)  integral += content[i];
   for (int i=0; i<content.size(); i++)  content[i] /= integral;
   return;
}

void Histogram::printContent() {
   printMatrix::print( {edges, content} );
}

#endif
