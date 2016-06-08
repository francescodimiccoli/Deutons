#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include "GCR_data.h"
#include <string>
#include <algorithm>

/** A histogram is a combination of two vectors, edges and contents
 * Edges is sorted, contents has one element less (both are checked by CheckIntegrity)
 *
 * */

class Histogram {
   public:
      bool checkIntegrity();
      bool fillWithVectors (std::vector<float> edge, std::vector<float> cont );
      void fillWithGalpropFile(std::string filename);
      std::vector<float> getEdges()   {return edges  ;}
      std::vector<float> getContent() {return content;}
      
   private:
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


void Histogram::fillWithGalpropFile(std::string filename)
{
   GCR_data galprop;
   galprop.read((char *)filename.data(), ( char*) "m2", ( char*) "GeV");
   int nentries=galprop.n;
   for (int i=0; i<nentries; i++) {
      edges.  push_back (galprop.E_low_input[i]);
      content.push_back (galprop.value_input[i]);
   }
   edges.push_back (galprop.E_high_input[nentries-1]); // last bin

   return;
}

#endif 
