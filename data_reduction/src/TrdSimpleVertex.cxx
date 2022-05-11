#include "TrdSimpleVertex.h" 

void TrdSimpleVertexRecon::Clear() {
  for (int ixy=0; ixy<2; ixy++) {
    for (vector<TrdSimpleVertex2D*>::iterator it=vertices2D[ixy].begin(); it!=vertices2D[ixy].end(); it++) delete (*it);
    vertices2D[ixy].clear();
  }
  for (vector<TrdSimpleVertex3D*>::iterator it=vertices3D.begin(); it!=vertices3D.end(); it++) delete (*it);
  vertices3D.clear();
}

void TrdSimpleVertexRecon::Print() {
  for (int ixy=0; ixy<2; ixy++) 
    for (vector<TrdSimpleVertex2D*>::iterator it=vertices2D[ixy].begin(); it!=vertices2D[ixy].end(); it++) 
      cout << "  TrdSimpleVertexRecon::Print2D ixy:" << ixy << " NSeg:" << (*it)->NSeg << " Coo(xy,z):(" << (*it)->Coo[0] << "," << (*it)->Coo[1] << ") D2:" << " " << (*it)->Dist2 << endl;             
  for (vector<TrdSimpleVertex3D*>::iterator it=vertices3D.begin(); it!=vertices3D.end(); it++) 
    cout << "  TrdSimpleVertexRecon::Print3D NSeg(X,Y):(" << (*it)->NSegX << "," << (*it)->NSegY << ") Coo(x,y,z):(" << (*it)->Coo[0] << "," << (*it)->Coo[1] << "," << (*it)->Coo[2] << ") D2:" << " " << (*it)->Dist2 << endl;
} 

bool TrdSimpleVertexRecon::Build(AMSEventR* event, double max_dist, double max_d2, int debug) {
  Clear();
  if (!event) return false;
  if (debug>0) cout << "TrdSimpleVertexRecon::Build-Processing nTrdSegment:" << event->nTrdSegment() << endl;

  /////////////////////////
  // 2D    
  /////////////////////////

  // init
  int nseg = event->nTrdSegment();
  int* used = new int[nseg]; // 0: not used, 1: used, 2: bad chi2, 3: tmp used (to be confirmed)
  for (int iseg=0; iseg<nseg; iseg++) used[iseg] = 0;
  // start point
  for (int iseg1=0; iseg1<nseg-1; iseg1++) {
    double z_ref = 0;
    double x_ref = 0;
    double y_ref = 0;
    double d2_ref = 0;
    TrdSegmentR* seg1 = event->pTrdSegment(iseg1);
    if (seg1->Chi2<=0) {
      used[iseg1] = 2;
      continue;
    }
    vector<TrdSegmentR*> list;
    list.push_back(seg1);
    used[iseg1] = 3;
    int ixy = seg1->Orientation;
    if (debug>1) cout << "TrdSimpleVertexRecon::Build-Seed-Segment index:" << iseg1 << " ixy:" << ixy << endl; 
    // try to add one
    for (int iseg2=iseg1+1; iseg2<nseg; iseg2++) {
      TrdSegmentR* seg2 = event->pTrdSegment(iseg2);
      if (ixy!=seg2->Orientation) continue;
      if (seg2->Chi2<=0) {
        used[iseg2] = 2;
        continue;
      }
      if (used[iseg2]) continue;
      // calculate new vertex position
      double x = 0;
      double y = 0;
      double z = 0;               
      double d2 = 0;
      list.push_back(seg2);
      used[iseg2] = 3;
      Fit3D(x,y,z,d2,list);
      // check that is not too far 
      double dist = sqrt(pow(x-x_ref,2)+pow(y-y_ref,2)+pow(z-z_ref,2));
      if ( (debug>1)&&(list.size()>2) ) cout << "TrdSimpleVertexRecon::Build-Try-to-Add index:" << iseg2 << " (x,y,z):(" << x << "," << y << "," << z << ") dist(cm):" << dist << endl;
      if ( (list.size()>2)&&(dist>max_dist) ) {
        used[iseg2] = 0;
        list.pop_back();
        continue;
      }
      // update
      used[iseg2] = 3;
      Fit3D(x_ref,y_ref,z_ref,d2_ref,list);
      if (debug>1) cout << "TrdSimpleVertexRecon::Build-Update-Cand-Vertex N:" << (int)list.size() << " (x,y,z):(" << x_ref << "," << y_ref << "," << z_ref << ") sqrt(d2)(cm):" << sqrt(d2_ref) << endl;
    }
    // drop only 2 
    if (list.size()<=2) {
      for (int iseg=0; iseg<nseg; iseg++) if (used[iseg]==3) used[iseg] = 0;
      list.clear();
      continue;
    }
    // take it
    for (int iseg=0; iseg<nseg; iseg++) if (used[iseg]==3) used[iseg] = 1;        
    vertices2D[ixy].push_back(new TrdSimpleVertex2D((int)list.size(),ixy,z_ref,(ixy==0)?x_ref:y_ref,d2_ref,list));
    list.clear();
  }
  // best on top
  for (int ixy=0; ixy<2; ixy++) sort(vertices2D[ixy].begin(),vertices2D[ixy].end(),SortVertex2D);

  /////////////////////////
  // 3D    
  /////////////////////////

  // all combinations, excluding too bad ones 
  for (int ix=0; ix<(int)vertices2D[0].size(); ix++) {
    for (int iy=0; iy<(int)vertices2D[1].size(); iy++) {
      vector<TrdSegmentR*> list;
      int nsegx = (int) vertices2D[0].at(ix)->List.size();
      int nsegy = (int) vertices2D[1].at(iy)->List.size();
      for (int isx=0; isx<nsegx; isx++) list.push_back(vertices2D[0].at(ix)->List.at(isx));
      for (int isy=0; isy<nsegy; isy++) list.push_back(vertices2D[1].at(iy)->List.at(isy));
      double coo[3];
      double dist2;
      Fit3D(coo[0],coo[1],coo[2],dist2,list);
      if (dist2>max_d2) continue; 
      vertices3D.push_back(new TrdSimpleVertex3D(nsegx,nsegy,coo,dist2));
    }
  }
  // best on top 
  sort(vertices3D.begin(),vertices3D.end(),SortVertex3D);
  if (debug>0) Print();

  // clean-up and return
  delete [] used;
  return true;
}

bool TrdSimpleVertexRecon::Intersection2D(double &z, double& xy, TrdSegmentR* seg1, TrdSegmentR* seg2) {
  z = 0;
  xy = 0;
  if ( (!seg1)||(!seg2) ) return false;
  double a = seg1->FitPar[0];
  double b = seg2->FitPar[0];
  double c = seg1->FitPar[1];
  double d = seg2->FitPar[1];
  if (fabs(a-b)<=0) return false;
  z = (d-c)/(a-b);
  xy = (a*d-b*c)/(a-b);  
  return true;
}

bool TrdSimpleVertexRecon::Fit3D(double& x, double &y, double& z, double& d2, vector<TrdSegmentR*> list) {
  x = 0;
  y = 0;
  z = 0;
  d2 = 0;
  if (list.empty()) return false;
  double sw[2] = {0};
  double sk[2] = {0};
  double sm[2] = {0};
  double smk[2] = {0};
  double sm2[2] = {0};
  for (vector<TrdSegmentR*>::iterator it=list.begin(); it!=list.end(); it++) {
    double m = (*it)->FitPar[0];
    double k = (*it)->FitPar[1];
    double w = 1./(1.+m*m);
    int ixy = (*it)->Orientation; 
    sw[ixy] += w;
    sk[ixy] += k*w;
    sm[ixy] += m*w;
    smk[ixy] += k*m*w;
    sm2[ixy] += m*m*w;
  }
  if (sw[0]>0) sw[0] = 1/sw[0];
  if (sw[1]>0) sw[1] = 1/sw[1];  
  z = (sk[0]*sm[0]*sw[0]+sk[1]*sm[1]*sw[1]-smk[0]-smk[1])/(sm2[0]+sm2[1]-sm[0]*sm[0]*sw[0]-sm[1]*sm[1]*sw[1]); 
  x = (sk[0]+sm[0]*z)*sw[0];
  y = (sk[1]+sm[1]*z)*sw[1];
  d2 = 0;
  for (vector<TrdSegmentR*>::iterator it=list.begin(); it!=list.end(); it++) {
    double m = (*it)->FitPar[0];
    double k = (*it)->FitPar[1];
    int ixy = (*it)->Orientation;
    d2 += (ixy==0) ?
      (m*z+k-x)*(m*z+k-x)/(1+m*m) :
      (m*z+k-y)*(m*z+k-y)/(1+m*m);
  }
  return true;
}

bool TrdSimpleVertexRecon::SortVertex2D(TrdSimpleVertex2D* a, TrdSimpleVertex2D* b) { 
  if ( (!a)||(!b) ) return false;
  return (a->Dist2<b->Dist2);
}

bool TrdSimpleVertexRecon::SortVertex3D(TrdSimpleVertex3D* a, TrdSimpleVertex3D* b) {
  if ( (!a)||(!b) ) return false;
  return (a->Dist2<b->Dist2);
}


