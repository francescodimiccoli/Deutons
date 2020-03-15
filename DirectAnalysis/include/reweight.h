#ifndef REWEIGHT_H
#define REWEIGHT_H

#include "histogram.h"
int FindTimeIndex(std::string input);

class Reweighter {
public:
    Reweighter(){ }	
    Reweighter(Histogram from, Histogram to);
    Histogram getFrom() const { return from; }
    Histogram getTo  () const { return to;   }
    float getWeight(float value, float def = 0.0) const;
private:
    Histogram from, to;
};

#endif
