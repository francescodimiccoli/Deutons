#ifndef __TrdSimpleVertex_h__
#define __TrdSimpleVertex_h__

#include <root.h>

#include <algorithm>
#include <iostream>
#include <vector>

using namespace std;

class TrdSimpleVertex2D {

public:
  vector<TrdSegmentR *> List;

  int NSeg;
  int Orien;
  double Coo[2];
  double Dist2;

  TrdSimpleVertex2D(int nseg, int orien, double z, double xy, double d2)
      : NSeg(nseg), Orien(orien), Coo{xy, z}, Dist2(d2) {}
  TrdSimpleVertex2D(int nseg, int orien, double z, double xy, double d2,
                    vector<TrdSegmentR *> list)
      : List(list), NSeg(nseg), Orien(orien), Coo{xy, z}, Dist2(d2) {}
};

class TrdSimpleVertex3D {

public:
  int NSegX;
  int NSegY;
  double Coo[3];
  double Dist2;

  TrdSimpleVertex3D(int nsegx, int nsegy, double coo[3], double dist2)
      : NSegX(nsegx), NSegY(nsegy), Coo{coo[0], coo[1], coo[2]}, Dist2(dist2) {}
};

class TrdSimpleVertexRecon {

public:
  vector<TrdSimpleVertex2D *> vertices2D[2];
  vector<TrdSimpleVertex3D *> vertices3D;

  TrdSimpleVertexRecon() {
    vertices2D[0].clear();
    vertices2D[1].clear();
    vertices3D.clear();
  }

  ~TrdSimpleVertexRecon() { Clear(); }
  void Clear();
  void Print();

  bool Build(AMSEventR *event, double max_dist = 2, double max_d2 = 1000,
             int debug = 0);

  static bool Intersection2D(double &z, double &xy, TrdSegmentR *seg1,
                             TrdSegmentR *seg2);
  static bool Fit3D(double &x, double &y, double &z, double &d2,
                    vector<TrdSegmentR *> list);
  static bool SortVertex2D(TrdSimpleVertex2D *a, TrdSimpleVertex2D *b);
  static bool SortVertex3D(TrdSimpleVertex3D *a, TrdSimpleVertex3D *b);
};

#endif
