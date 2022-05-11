#ifndef __Fiducial_h__ 
#define __Fiducial_h__

#include "root.h"

int GetPatternInsideTracker(TrTrackR* track, int fit_id = 0);
int GetPatternInsideTracker(TrdTrackR* track);

#endif
