#ifndef _PGTRACK_
#define _PGTRACK_
#include "TrTrack.h"
#endif

#include <root.h>

typedef struct mchit{
 int pid;
 float pos[3];
 float mom[3];
 float edep;
 mchit():pid(-1),edep(-1){
 	pos[0]=-1000;pos[1]=-1000;;pos[2]=-1000;
 	mom[0]=-1000;mom[1]=-1000;;mom[2]=-1000;
 }

} mchit;
/* Finds closest MC cluster to a given AMSPoint 
 * the clusters beyond match_raidus are discarded
 * */
TrMCClusterR * findClosestMCCluster(AMSEventR * event, const AMSPoint & point, double match_radius=0.1);


/* Gets the Geant PID for a closest MC Cluster at a give JLayer.
 * Returns 0 if nothing is found. 
 * */
mchit getJLayerGeantID(AMSEventR * event, int ntrack, int jlayer, double match_raidus=0.1);


/* Gets layer 1-to-4 Geant PIDs of the closest MC clusters and packs them 
 * into a 32-bit integer:
 *
 * | PID Layer 4 |  PID Layer 3 | PID Layer 2 | PID Layer 1 |
 * 32            24             16            8             0
 *
 * */
UInt_t getPackedLayers_1to4(AMSEventR * event, int ntrack, double match_radius=0.1);
