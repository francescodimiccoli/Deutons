#ifndef HISTUTILS_H
#define HISTUTILS_H

#include "histogram.h"
#include "binning.h"
#include "GCR_data.h"

inline Histogram loadGalpropFile (std::string filename)
{
    GCR_data galprop;
    galprop.read ( (char *) filename.data(), ( char*) "m2", ( char*) "GeV");
    int nentries=galprop.n;

    std::vector<float> content(nentries);
    std::vector<float>   edges(nentries);
    edges.assign(nentries, 0);

    for (int i = 0; i < nentries; i++) {
        edges  [i] = galprop.E_low_input[i];
        content[i] = galprop.value_input[i];
    }

    edges.push_back (galprop.E_high_input[nentries-1]); // last bin
    return Histogram(edges, content);
}

inline Histogram makeLogUniform(int nbins, float min, float max) {
    Binning binning;
    binning.setBinsFromRigidity (nbins, min, max,ResponseAgl,4.28781e-05,67.8521);
    Histogram h( binning.RigBins() );
    std::vector<float> counts(nbins, 1);
    h.fillWithCounts(counts);
    return h;
}
#endif
