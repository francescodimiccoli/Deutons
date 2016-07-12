#ifndef REWEIGHT_H
#define REWEIGHT_H

#include "histogram.h"

class Reweighter {
public:
    Reweighter(Histogram from, Histogram to);
    Histogram getFrom() const { return from; }
    Histogram getTo  () const { return to;   }
    inline float getWeight(float value, float def = 0.0) const;
private:
    Histogram from, to;
};


Reweighter::Reweighter(Histogram f, Histogram t): from(f), to(t) {}

float Reweighter::getWeight(float value, float def) const {
    float fromV = from.at(value);
    float   toV =   to.at(value);
    if(fromV == 0.0) return def;
    return toV / fromV;
}

#endif
