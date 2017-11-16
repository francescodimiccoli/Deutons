#include "reweight.h"

Reweighter::Reweighter(Histogram f, Histogram t): from(f), to(t) {}

float Reweighter::getWeight(float value, float def) const {
    float fromV = from.at(value);
    float   toV =   to.at(value);
    if(fromV == 0.0) return def;
    return toV / fromV;
}


