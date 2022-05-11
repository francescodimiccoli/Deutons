#ifndef _RICHANALYSISMORE_
#define _RICHANALYSISMORE_

#include <vector>
#include <fstream>
#include <assert.h>
#include "TObject.h"

using namespace std;

///// A geometrical hash useful for multidimensional studies
#include<stdlib.h>
#include<math.h>

namespace GeomHashBetaConstants{
  const int nBuffers=4;
};

using namespace GeomHashBetaConstants;

// Faster (hopefully much faster) implementation of the tool
class GeomHashBeta: public TObject{
 public:
  static bool storeEntries;
  static bool storeTemplates;
  static bool storeTemplatesRms;
  static bool storeMean;
  static bool storeRms;
  static bool storePeak;
  static bool storePeakWidth;
  static bool computePeakAsMean;
  static float peakFinderFraction; 

  static void setDefaults(){
    GeomHashBeta::storeEntries=true;
    GeomHashBeta::storeTemplates=true;
    GeomHashBeta::storeTemplatesRms=true;
    GeomHashBeta::storeMean=true;
    GeomHashBeta::storeRms=true;
    GeomHashBeta::storePeak=true;
    GeomHashBeta::storePeakWidth=true;
    GeomHashBeta::computePeakAsMean=false;
    GeomHashBeta::peakFinderFraction=0.68;
  }

  int dimension;
  int numNodes;
  int minSize;              // min size allowed to bin
  vector<unsigned char> points;       // Chosen direction
  vector<float> limit;     // The limit to separate among them
  vector<int>   nodes[2];   // If the node j is terminal, node[0][j]=-1, node[1][j]=entries in bin, point[0][j*dimension]=mean value, point[1][j*dimension]=rms                                                   
  vector<float> Means;
  vector<float> Rms;
  vector<int>   Entries;
  vector<float> Template;
  vector<float> TemplateRms;
  vector<float> Peak;
  vector<float> PeakWidth;

  GeomHashBeta(int d=1);

  float *getTemplate(int h);
  float *getTemplateRms(int h);

  inline int offset(int which){return which*dimension;}
  //Evaluating
  int hash(float *point);
  int get(double x,...);

  // Virtual function
  void store(int *pointer,int size,int parent);                                           // Storing information in the tree

  // Growing
  vector<float> samples;  //! Vector storing all the samples with the format x0,x1,x2...xn,y,weight
  vector<float> values;   //! Vector storing the guiding values
  void push(float value,float *x);
  void fill(double x,...);


  void grow(int min_size=0);

  void grow(int *pointers,                          // Buffer pointing to the points indexes
	    int *scratch,
	    int min_size=0);        // Scratch region to store the distances to the points


  // grow a tree trying tries time for optimization, and with minSize minimum elements per node. size is the number of samples to be used
  // stored in buffer (as indexes to samples) and scratch is a temporary scratch space
  void  grow_internal(int *pointers,
                      int  size,
		      int *scratch,
                      int parent=0);



  float getMean(int node);
  float getRms(int node);
  float getPeak(int node);
  float getPeakWidth(int node);

  int getEntries(int node);


  // Some internal buffers
  static float *buffer[GeomHashBetaConstants::nBuffers]; //!
  static int    bufferSize;       //!
  void checkBuffers();

#ifdef VERSION6
#else
#pragma omp threadprivate(fgIsA)
#endif
  ClassDef(GeomHashBeta,1);
};


class GeomHashBetaEnsemble: public GeomHashBeta{
 public:
  GeomHashBetaEnsemble(int d=1):GeomHashBeta(d){}

  vector<GeomHashBeta> hashes;
  void growOne(int minSize=10,bool bootstrap=false);
  void Eval(float *x);
  void Eval(double x,...);

  // Result from last evaluation
  double MeanValue;
  double ValueRms;
  double MeanRms;
  double RmsRms;
  double MeanPeak;
  double PeakRms;
  double MeanPeakWidth;
  double PeakWidthRms;
  double MeanEntries;
  double EntriesRms;
  int    Hashes;

  int numHashes();
  GeomHashBeta &getHash(int i);

  ClassDef(GeomHashBetaEnsemble,1);
#ifdef VERSION6
#else
#pragma omp threadprivate(fgIsA)
#endif
};

// A small tool
void merge(GeomHashBetaEnsemble &receiver,GeomHashBetaEnsemble &small);


// Manager to deal with the default GeomHashBeta
#include "TFile.h"
#include "TRandom.h"
class GHBManager{
 public:
  TFile *dataF;
  TFile *mcF;
  GeomHashBetaEnsemble *data;
  GeomHashBetaEnsemble *mc;
  
  GHBManager(TString fnameData="dataHashBeta.root",TString fnameMC="mcHashBeta.root");
  ~GHBManager(){dataF->Close();mcF->Close();}


  void getParameters(float x,float y,float theta,float phi, bool isMC,float &betaCorrection,float &width);
  float correctBeta(float beta,float x,float y,float theta,float phi, int charge=1,bool isMC=false);
};






#endif





#ifndef _CUTOFF_TOOLS_
#define _CUTOFF_TOOLS_

#include <vector>
#include <fstream>
#include <assert.h>
#include "TNamed.h"

using namespace std;

///// A geometrical hash useful for multidimensional studies
#include<stdlib.h>
#include<math.h>

namespace GeomHashNamedConstants{
  const int nBuffers=4;
};

using namespace GeomHashNamedConstants;

// Faster (hopefully much faster) implementation of the tool
class GeomHashNamed: public TNamed{
 public:
  static bool storeEntries;
  static bool storeTemplates;
  static bool storeTemplatesRms;
  static bool storeMean;
  static bool storeRms;
  static bool storePeak;
  static bool storePeakWidth;
  static bool computePeakAsMean;
  static float peakFinderFraction; 
  static bool sampleAllDirections;

  static void setDefaults(){
    GeomHashNamed::storeEntries=true;
    GeomHashNamed::storeTemplates=true;
    GeomHashNamed::storeTemplatesRms=true;
    GeomHashNamed::storeMean=true;
    GeomHashNamed::storeRms=true;
    GeomHashNamed::storePeak=true;
    GeomHashNamed::storePeakWidth=true;
    GeomHashNamed::computePeakAsMean=false;
    GeomHashNamed::peakFinderFraction=0.68;
    GeomHashNamed::sampleAllDirections=false;
  }

  int dimension;
  int numNodes;
  int minSize;              // min size allowed to bin
  vector<unsigned char> points;       // Chosen direction
  vector<float> limit;     // The limit to separate among them
  vector<int>   nodes[2];   // If the node j is terminal, node[0][j]=-1, node[1][j]=entries in bin, point[0][j*dimension]=mean value, point[1][j*dimension]=rms                                                   
  vector<float> Means;
  vector<float> Rms;
  vector<int>   Entries;
  vector<float> Template;
  vector<float> TemplateRms;
  vector<float> Peak;
  vector<float> PeakWidth;

  GeomHashNamed(int d=1);

  float *getTemplate(int h);
  float *getTemplateRms(int h);

  inline int offset(int which){return which*dimension;}
  //Evaluating
  int hash(float *point);
  int get(double x,...);

  // Virtual function
  void store(int *pointer,int size,int parent);                                           // Storing information in the tree

  // Growing
  vector<float> samples;  //! Vector storing all the samples with the format x0,x1,x2...xn,y,weight
  vector<float> values;   //! Vector storing the guiding values
  void push(float value,float *x);
  void fill(double x,...);


  void grow(int min_size=0);

  void grow(int *pointers,                          // Buffer pointing to the points indexes
	    int *scratch,
	    int min_size=0);        // Scratch region to store the distances to the points


  // grow a tree trying tries time for optimization, and with minSize minimum elements per node. size is the number of samples to be used
  // stored in buffer (as indexes to samples) and scratch is a temporary scratch space
  void  grow_internal(int *pointers,
                      int  size,
		      int *scratch,
                      int parent=0);



  float getMean(int node);
  float getRms(int node);
  float getPeak(int node);
  float getPeakWidth(int node);

  int getEntries(int node);


  // Some internal buffers
  static float *buffer[GeomHashNamedConstants::nBuffers]; //!
  static int    bufferSize;       //!
  void checkBuffers();

#ifdef VERSION6
#else
#pragma omp threadprivate(fgIsA)
#endif
  ClassDef(GeomHashNamed,1);
};


class GeomHashNamedEnsemble: public GeomHashNamed{
 public:
  GeomHashNamedEnsemble(int d=1):GeomHashNamed(d){}

  vector<GeomHashNamed> hashes;
  void growOne(int minSize=10,bool bootstrap=false);
  void Eval(float *x);
  void Eval(double x,...);

  // Result from last evaluation
  double MeanValue;
  double ValueRms;
  double MeanRms;
  double RmsRms;
  double MeanPeak;
  double PeakRms;
  double MeanPeakWidth;
  double PeakWidthRms;
  double MeanEntries;
  double EntriesRms;
  int    Hashes;

  int numHashes();
  GeomHashNamed &getHash(int i);

  ClassDef(GeomHashNamedEnsemble,1);
#ifdef VERSION6
#else
#pragma omp threadprivate(fgIsA)
#endif
};



#endif





