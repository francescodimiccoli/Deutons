#include "Fiducial.h"

static float tracker_layers_z[9] = {158.920,53.060,29.228,25.212,1.698,-2.318,-25.212,-29.228,-135.882};

static float tracker_planes_edges[9][4] = { // tracker edges xmin,ymin,xmax,ymax (xmax is used also as disk radius for L1, ..., L8)
  {-62.14,  -47.40,   62.14,   47.40},
  {-62.14,  -40.10,   62.14,   40.10},
  {-49.70,  -43.75,   49.70,   43.75},
  {-49.72,  -43.75,   49.72,   43.75},
  {-49.71,  -36.45,   49.70,   36.45},
  {-49.72,  -36.45,   49.72,   36.45},
  {-49.72,  -43.75,   49.71,   43.75},
  {-49.72,  -43.75,   49.71,   43.75},
  {-45.62,  -29.48,   45.55,   29.53}
};

int GetPatternInsideTracker(TrTrackR* track, int fit_id) {
  if ((!track)||(fit_id<0)) return 0;
  int pattern = 0;
  for (int ilayer=0; ilayer<9; ilayer++) {
    AMSPoint point;
    AMSDir dir;
    track->Interpolate(tracker_layers_z[ilayer],point,dir,fit_id);
    float x = point.x();
    float y = point.y();
    bool isinlayer = false;
    if ( (x>tracker_planes_edges[ilayer][0])&&(x<tracker_planes_edges[ilayer][2])&&
         (y>tracker_planes_edges[ilayer][1])&&(y<tracker_planes_edges[ilayer][3]) ) {
      if ((ilayer+1)==9) isinlayer = true;
      else {
        if ( (sqrt(x*x+y*y)<tracker_planes_edges[ilayer][2]) ) isinlayer = true;
      }
    }
    if (isinlayer) pattern |= (1<<ilayer);
  }
  return pattern;
}

int GetPatternInsideTracker(TrdTrackR* track) {
  if (!track) return 0; 
  int pattern = 0;
  for (int ilayer=0; ilayer<9; ilayer++) {
    float coo[3] = {0};
    for (int i=0; i<3; i++) coo[i] = track->Coo[i];
    float theta = track->Theta;
    float phi = track->Phi;
    float z = tracker_layers_z[ilayer];
    float x = coo[0] + tan(theta)*cos(phi)*(z-coo[2]);
    float y = coo[1] + tan(theta)*sin(phi)*(z-coo[2]);
    bool isinlayer = false;
    if ( (x>tracker_planes_edges[ilayer][0])&&(x<tracker_planes_edges[ilayer][2])&&
         (y>tracker_planes_edges[ilayer][1])&&(y<tracker_planes_edges[ilayer][3]) ) {
      if ((ilayer+1)==9) isinlayer = true;
      else {
        if ( (sqrt(x*x+y*y)<tracker_planes_edges[ilayer][2]) ) isinlayer = true;
      }
    }
    if (isinlayer) pattern |= (1<<ilayer);
  }
  return pattern;
}


