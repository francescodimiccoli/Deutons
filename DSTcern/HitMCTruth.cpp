#include "HitMCTruth.h"

TrMCClusterR * findClosestMCCluster(AMSEventR * event, const AMSPoint & point, double match_radius) {
    if( ! event ) return NULL;
    if( event->nMCEventgC() == 0 ) return NULL;

    double rmin = 10000;
    TrMCClusterR * closest = NULL;

    // Finding cluster with minimal distance
    for(int c = 0; c < event->NTrMCCluster(); c++) {
        TrMCClusterR * cluster = event->pTrMCCluster(c);
        double r = point.dist(cluster->GetXgl());
        if(r < rmin) {
            rmin = r;
            closest = cluster;
        }
    }

    // Fail if minimal distance is not in match_radius 
    if(rmin > match_radius) return NULL;

    return closest; 
}


int getJLayerGeantID(AMSEventR * event, int ntrack, int jlayer, double match_radius) {
    TrTrackR * track = event->pTrTrack(ntrack); 
    if(!track) return 0;
    TrRecHitR * hit = track->GetHitLJ(jlayer);
    if(!hit) return 0;

    TrMCClusterR * cluster = findClosestMCCluster(event, hit->GetCoord(), match_radius);
    if(!cluster) return 0;

    return cluster->GetPart();
}

UInt_t getPackedLayers_1to4(AMSEventR * event, int ntrack, double match_radius) {
    TrTrackR * track = event->pTrTrack(ntrack); 
    if(!track) return 0;

    UInt_t id1 = getJLayerGeantID(event, ntrack, 1,  match_radius);
    UInt_t id2 = getJLayerGeantID(event, ntrack, 2,  match_radius);
    UInt_t id3 = getJLayerGeantID(event, ntrack, 3,  match_radius);
    UInt_t id4 = getJLayerGeantID(event, ntrack, 4,  match_radius);

    id1 = id1 < 127 ? id1 : 0;
    id2 = id2 < 127 ? id2 : 0;
    id3 = id3 < 127 ? id3 : 0;
    id4 = id4 < 127 ? id4 : 0;

    return id1 + (id2 << 8) + (id3 << 16) + (id4 << 24);
}
