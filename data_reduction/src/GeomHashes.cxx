#include "GeomHashes.h"
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <iostream>
#include <algorithm>
#undef _DEBUG_
#define ASSERT(x) 

using namespace std;

ClassImp(GeomHashBeta);

using namespace GeomHashBetaConstants;

float *GeomHashBeta::buffer[GeomHashBetaConstants::nBuffers]={0,0,0,0};
int GeomHashBeta::bufferSize=0;

bool GeomHashBeta::storeEntries=true;
bool GeomHashBeta::storeTemplates=true;
bool GeomHashBeta::storeTemplatesRms=true;
bool GeomHashBeta::storeMean=true;
bool GeomHashBeta::storeRms=true;
bool GeomHashBeta::storePeak=true;
bool GeomHashBeta::storePeakWidth=true;
bool GeomHashBeta::computePeakAsMean=false;
float GeomHashBeta::peakFinderFraction=0.68;

void GeomHashBeta::checkBuffers(){
  if(dimension<=bufferSize) return;
  bufferSize=dimension;
  for(int i=0;i<GeomHashBetaConstants::nBuffers;i++){
    if(buffer[i]) delete[] buffer[i];
    buffer[i]=new float[dimension+1];
  }
}


GeomHashBeta::GeomHashBeta(int d){
  if(d>255){cout<<"Initializing to 255 variables."<<endl;d=255;}
  if(d<0) d=0;
  dimension=d;
  numNodes=0;
  points.reserve(1024);
  nodes[0].reserve(1024);
  nodes[1].reserve(1024);
  limit.reserve(1024);

  // Grow the buffers if needed
  //  checkBuffers();
}


int GeomHashBeta::hash(float *point){
  checkBuffers();
  if(numNodes==0) return 0;
  int current=0;
  for(;;){
    int nodeNumber=point[points[current]]<limit[current]?0:1;
    if(nodes[nodeNumber][current]<0) return -nodes[nodeNumber][current]-1;
    current=nodes[nodeNumber][current];
  }
}

void GeomHashBeta::push(float value,float *x){
  // Set the minimum capacity
  if(int(samples.capacity())<1024*dimension) samples.reserve(1024*dimension);

  // Store the point
  for(int i=0;i<int(dimension);i++)  samples.push_back(x[i]);

  values.push_back(value);
}

#include <stdarg.h>
void GeomHashBeta::fill(double x,...){
  checkBuffers();
  va_list ap;
  va_start(ap, x);
  float value=x;
  for(int i=0;i<dimension;i++) buffer[2][i]=va_arg(ap, double);
  va_end(ap);
  push(value,buffer[2]);
}

int GeomHashBeta::get(double x,...){
  checkBuffers();
  va_list ap;
  va_start(ap, x);
  buffer[2][0]=x;
  for(int i=1;i<dimension;i++) buffer[2][i]=va_arg(ap, double);
  va_end(ap);
  return hash(buffer[2]);
}

int *_WallB_;  // Just in case

void GeomHashBeta::grow(int min_size){
  points.clear();
  limit.clear();
  nodes[0].clear(); nodes[1].clear();
  int *pointers=new int[samples.size()/dimension];
  int *scratch=new  int[samples.size()/dimension];

  int size=samples.size()/dimension;
  _WallB_=pointers+size;
  grow(pointers,scratch,min_size);

  delete[] pointers;
  delete[] scratch;
}


void GeomHashBeta::grow(int *pointers,int *scratch,int min_size){
  int entries=samples.size()/dimension;

  // Init the buffer
  for(int i=0;i<entries;i++) pointers[i]=i;

  // Start the real growing process
  minSize=2*min_size;
  grow_internal(pointers,entries,scratch); // Adjust the min_size to obtain a better approximation

  // Liberates the space used by the samples and values
  samples.resize(0);
  values.resize(0);
  vector<float> empty; 
  samples=empty;
  values=empty;
}

void GeomHashBeta::store(int *pointers,int size,int parent){
  checkBuffers();
  int currentHashNumber=numNodes++;					
  int trueParent=fabs(parent)-1;					
  int parentNode=parent<0?0:1;						
  nodes[parentNode].at(trueParent)=-currentHashNumber-1;
  
  // Here we can store the information concerning the data, for example mean and rms
  if(storeEntries) Entries.push_back(size);

  double sum=0,sum2=0;
  for(int i=0;i<size;i++){
    sum+=values[pointers[i]];
    sum2+=values[pointers[i]]*values[pointers[i]];
  }
  sum/=size;
  sum2/=size;
  sum2-=sum*sum;
  sum2=sqrt(fabs(sum2));

  if(storeMean) Means.push_back(sum);
  if(storeRms) Rms.push_back(sum2);

  // Fill the templates
  for(int i=0;i<dimension;i++) buffer[0][i]=buffer[1][i]=0;
  for(int j=0;j<dimension;j++){
    double sum=0;
    double sum2=0;
    for(int i=0;i<size;i++){
      int &element=pointers[i];
      float &v=samples[offset(element)+j];
      sum+=v;
      sum2+=v*v;
    }
    sum/=size;
    sum2/=size;
    sum2-=sum*sum;
    sum2=sqrt(fabs(sum2));
    if(storeTemplates) Template.push_back(sum);
    if(storeTemplatesRms) TemplateRms.push_back(sum2);
  }

  // Peak search and width determination
  vector<float> ordered;
  ordered.reserve(size);
  for(int i=0;i<size;i++) ordered.push_back(values[pointers[i]]);
  sort(ordered.begin(),ordered.end());
  int window=int(ceil(size*peakFinderFraction))-1; 
  if(window==0) window=1;

#ifdef _DEBUG_
  cout<<"PEAK FINDER FOR"<<endl;
  for(int i=0;i<size;i++) cout<<ordered[i]<<" ";
  cout<<endl;
  cout<<"Window size "<<window<<endl;
  cout<<"Sample size "<<size<<endl;
  cout<<"Mean "<<Means[Means.size()-1];
#endif


  float bestWidth=INFINITY;
  float peakPos=0;
  for(int i=0;i<size;i++){
    if(i+window>=size) break;
    if(i && ordered[i]==ordered[i-1]) continue; // skip identical events
    int finalIndex=i+window;
    for(int j=finalIndex+1;j<size;j++){
      if(ordered[j]!=ordered[finalIndex]) break;  
      finalIndex=j;
    }

    float width=ordered[finalIndex]-ordered[i];
    if(width<bestWidth){
      bestWidth=width;
      peakPos=0.5*(ordered[finalIndex]+ordered[i]);
#ifdef _DEBUG_
      cout<<"CURRENT BEST "<<peakPos<<" from "<<i<<" to "<<finalIndex<<" "<<ordered[i]<<" "<<ordered[finalIndex]<<endl;
#endif
    }
  }

  bestWidth/=2; // width meaning should be equivalent to sigma 

  if(computePeakAsMean){
    double sum=0;
    double sum2=0;
    int total=0; 
    for(int i=0;i<size;i++){
      if(ordered[i]>peakPos+bestWidth) continue;
      if(ordered[i]<peakPos-bestWidth) continue;
      sum+=ordered[i];
      sum2+=ordered[i]*ordered[i];
      total++;
    }
    sum/=total;
    sum2/=total;
    sum2-=sum*sum;
    sum2=sqrt(fabs(sum2));
    peakPos=sum;
    bestWidth=sum2;
  }

  if(storePeak) Peak.push_back(peakPos);
  if(storePeakWidth) PeakWidth.push_back(bestWidth);
}

float *GeomHashBeta::getTemplate(int h){
  if(!Template.size()) return 0;
  if(h<0 || h>=numNodes) return 0;
  checkBuffers();
  for(int i=0;i<dimension;i++) buffer[0][i]=Template[offset(h)+i];
  return buffer[0];
}

float *GeomHashBeta::getTemplateRms(int h){
  if(!TemplateRms.size()) return 0;
  if(h<0 || h>=numNodes) return 0;
  checkBuffers();
  for(int i=0;i<dimension;i++) buffer[1][i]=TemplateRms[offset(h)+i];
  return buffer[1];
}


float GeomHashBeta::getMean(int node){
  if(!Means.size()) return 0;
  if(node>=numNodes) return 0;
  return Means.at(node);
}

float GeomHashBeta::getRms(int node){
  if(!Rms.size()) return 0;
  if(node>=numNodes) return 0;
  return Rms.at(node);
}


float GeomHashBeta::getPeak(int node){
  if(!Peak.size()) return 0;
  if(node>=numNodes) return 0;
  return Peak.at(node);
}

float GeomHashBeta::getPeakWidth(int node){
  if(!PeakWidth.size()) return 0;
  if(node>=numNodes) return 0;
  return PeakWidth.at(node);
}

int GeomHashBeta::getEntries(int node){
  if(!Entries.size()) return 0;
  if(node>=numNodes) return 0;
  return Entries.at(node);
}


void GeomHashBeta::grow_internal(int *pointers,int size,int *scratch,int parent){
  const float epsilon=5e-6;
  // Check if we should proceed
  if(size==1 || size<minSize) {store(pointers,size,parent);return;}

  // Fast search binning
  //  const int numBins=1000;  // Number of bins to be used during spliting
  const int numBins=10;  // Number of bins to be used during spliting
  int entries[numBins];
  double mean[numBins];
  double mean2[numBins]; 
  for(int i=0;i<numBins;i++) entries[i]=0;
  for(int i=0;i<numBins;i++) mean[i]=0;
  for(int i=0;i<numBins;i++) mean2[i]=0;

  int me=nodes[0].size();

  // Pick randomly a direction
  int direction=int(dimension*(rand() / (RAND_MAX + 1.0)));

  // Select the sorting range
  double minValue=INFINITY;
  double maxValue=-INFINITY;

  // DEBUG
  for(int index=0;index<size;index++){
    int element=pointers[index];
    float x=samples[offset(element)+direction];
    if(x<minValue) minValue=x;
    if(x>maxValue) maxValue=x;
  }


  // No binning
  if(fabs((maxValue-minValue)/((maxValue+minValue)/2))<epsilon) {store(pointers,size,parent);return;}

  double bw=(maxValue-minValue)/numBins;
  maxValue+=bw/2;
  minValue-=bw/2;
  bw=(maxValue-minValue)/numBins;

  // No binning
  if(bw<=0) {store(pointers,size,parent);return;}

  // Fill the arrays
  double totalMean=0;
  double totalMean2=0;
  bool equals=true;

  double prevValue=values[0];
  for(int index=0;index<size;index++){
    int element=pointers[index];

    float x=samples[offset(element)+direction];
    float value=values[element];
    float value2=value*value;
    int bin=int(floor((x-minValue)/bw));

    if(value!=prevValue) equals=false;

    entries[bin]++;
    mean[bin]+=value;
    mean2[bin]+=value2;
    totalMean+=value;
    totalMean2+=value2;
  }


  // Search the best split position
  const double rms01=totalMean2/size-totalMean/size*totalMean/size;
  const int n01=size;

  int accEntries=0;
  double acc=0;
  double acc2=0;
  float bestValue=FP_INFINITE;
  float splitPoint;
  int splitCounter=0;
  bool fail=true;

  if(!equals)
  for(int i=0;i<numBins-1;i++){
    accEntries+=entries[i];
    acc+=mean[i];
    acc2+=mean2[i];
    if(size-accEntries<minSize/2) break;
    if(accEntries<minSize/2) continue;

    int n0=accEntries;
    double rms0=acc2/n0-acc/n0*acc/n0;
    int n1=size-n0; 
    double rms1=(totalMean2-acc2)/n1-(totalMean-acc)/n1*(totalMean-acc)/n1; 
    double value=n0*rms0+n1*rms1-n01*rms01;
    
    // Currently this is disabled on purpose by adding a negative sign into the first comparison
    // The reason for this is to avoid too unbalanced trees which give rise to problems
    // Somehow this has to be corrected in some moment
    if(0)
    if(value<-0.05*n01*rms01 && n0>0.5*minSize && n1>0.5*minSize)  // Large enough correction and enough entries on it
      if(value<bestValue){
	bestValue=value;
	splitCounter=n0;
	splitPoint=minValue+bw*(i+1);
	fail=false;
      }
  }


  if(fail){
    // Try the standard split
    int accEntries=entries[0];
    for(int i=1;i<numBins-1;i++){
      if(accEntries+entries[i]>size/2){
	splitPoint=minValue+bw*(i+1);
	splitCounter=accEntries;
	fail=false;
	break;
      }
      accEntries+=entries[i];
    }
    // No way, store it a go ahead
    if(fail) {store(pointers,size,parent);return;}
  }

  if(splitCounter==size) {store(pointers,size,parent);return;}


  // Copy the pointer to the scratch area to sort them 
  for(int index=0;index<size;index++) scratch[index]=pointers[index];

  // Sort
  int lefties=0,righties=size-1;
  splitCounter=0;
  for(int index=0;index<size;index++){
    int element=scratch[index];
    float x=samples[offset(element)+direction];
    if(x<splitPoint){
      pointers[lefties++]=scratch[index];
      splitCounter++;
    }else{
      pointers[righties--]=scratch[index];
    }
  }


  // Store the current information
  limit.push_back(splitPoint);
  points.push_back(direction);

  nodes[0].push_back(0);  nodes[1].push_back(0);
  nodes[0][me]=nodes[0].size();
  grow_internal(pointers,splitCounter,scratch,-(me+1));
  nodes[1][me]=nodes[1].size();
  grow_internal(pointers+splitCounter,size-splitCounter,scratch,me+1);
}


ClassImp(GeomHashBetaEnsemble);


void GeomHashBetaEnsemble::growOne(int minSize,bool bootstrap){
  hashes.push_back(GeomHashBeta(dimension));
  GeomHashBeta &hash=hashes.back();

  if(!bootstrap){
    hash.samples=samples;
    hash.values=values;
  }else{
    checkBuffers();
    for(unsigned int i=0;i<values.size();i++){
      int pointer=int(values.size()*(rand() / (RAND_MAX + 1.0)));
      for(int j=0;j<dimension;j++) buffer[0][j]=samples[offset(pointer)+j];
      hash.push(values[pointer],buffer[0]);
    }
  }
  hash.grow(minSize);
}


void GeomHashBetaEnsemble::Eval(float *x){
  MeanValue=0;
  ValueRms=0;
  MeanRms=0;
  RmsRms=0;
  MeanPeak=0;
  PeakRms=0;
  MeanPeakWidth=0;
  PeakWidthRms=0;
  Hashes=0;
  
  for(unsigned int i=0;i<hashes.size();i++){
    Hashes++;
    GeomHashBeta &hash=hashes[i];
    int h=hash.hash(x);
    double xx=hash.getMean(h);
    MeanValue+=xx;
    ValueRms+=xx*xx;

    xx=hash.getRms(h);
    MeanRms+=xx;
    RmsRms+=xx*xx;

    xx=hash.getPeak(h);
    MeanPeak+=xx;
    PeakRms+=xx*xx;

    xx=hash.getPeakWidth(h);
    MeanPeakWidth+=xx;
    PeakWidthRms+=xx*xx;

    xx=hash.getEntries(h);
    MeanEntries+=xx;
    EntriesRms+=xx*xx;
  }

#define Do(_x) {Mean##_x/=Hashes; _x##Rms /=Hashes; _x##Rms-=Mean##_x*Mean##_x; _x##Rms=sqrt(fabs(_x##Rms));}
  Do(Value);
  Do(Rms);
  Do(Peak);
  Do(PeakWidth);
  Do(Entries);
#undef Do
}
/*
void GeomHashBetaEnsemble::Eval(float *x){
  MeanValue=0;
  ValueRms=0;
  MeanRms=0;
  RmsRms=0;
  MeanPeak=0;
  double MeanPeakW=0;
  PeakRms=0;
  MeanPeakWidth=0;
  PeakWidthRms=0;
  Hashes=0;
  
  for(unsigned int i=0;i<hashes.size();i++){
    Hashes++;
    GeomHashBeta &hash=hashes[i];
    int h=hash.hash(x);
    double x=hash.getMean(h);
    MeanValue+=x;
    ValueRms+=x*x;

    x=hash.getRms(h);
    MeanRms+=x;
    RmsRms+=x*x;

    x=hash.getPeak(h);
    //    MeanPeak+=x;
    double w=hash.getPeakWidth(h);
    MeanPeak+=x/w/w;
    MeanPeakW+=1/w/w;
    PeakRms+=x*x;

    x=hash.getPeakWidth(h);
    MeanPeakWidth+=x;
    PeakWidthRms+=x*x;

    x=hash.getEntries(h);
    MeanEntries+=x;
    EntriesRms+=x*x;
  }

  double p=MeanPeak/MeanPeakW;

#define Do(_x) {Mean##_x/=Hashes; _x##Rms /=Hashes; _x##Rms-=Mean##_x*Mean##_x; _x##Rms=sqrt(fabs(_x##Rms));}
  Do(Value);
  Do(Rms);
  Do(Peak);
  Do(PeakWidth);
  Do(Entries);
#undef Do

  MeanPeak=p;

}
*/



void GeomHashBetaEnsemble::Eval(double x,...){
  checkBuffers();
  va_list ap;
  va_start(ap, x);
  buffer[2][0]=x;
  for(int i=1;i<dimension;i++) buffer[2][i]=va_arg(ap, double);
  va_end(ap);
  Eval(buffer[2]);
}

int GeomHashBetaEnsemble::numHashes(){
  return hashes.size();
}

GeomHashBeta& GeomHashBetaEnsemble::getHash(int i){
  return hashes.at(i);
}


//////////////// TOOL

void merge(GeomHashBetaEnsemble &receiver,GeomHashBetaEnsemble &small){
  if(receiver.dimension!=small.dimension){
    cout<<"DIMENSIONS DOES NOT MATCH"<<endl;
    return;
  }

  for(int i=0;i<small.numHashes();i++){
    receiver.hashes.push_back(small.getHash(i));
  }
}

////////////// Manager to use it on data and MC

GHBManager::GHBManager(TString fnameData,TString fnameMC){
  dataF= TFile::Open(fnameData);
  mcF= TFile::Open(fnameMC);
  // Open files
  if(!dataF || ! mcF){
    cerr<<"Problem with input files: "<<(!dataF?fnameData:"")<<" "<<(!mcF?fnameMC:"")<<endl;
    exit(1);
  }
  // Retrieve the hashes
  cout<<"GHBManager::READING DATA ..."<<endl;
  data=(GeomHashBetaEnsemble*)dataF->Get("GeomHashBetaEnsemble");  
  cout<<"GHBManager::READING MC ..."<<endl;
  mc=(GeomHashBetaEnsemble*)mcF->Get("GeomHashBetaEnsemble"); 
  if(!data || ! mc){
    cerr<<"Problem with input files: "<<(!data?fnameData:"")<<" "<<(!mc?fnameMC:"")<<endl;
    exit(1);
  }
}


void GHBManager::getParameters(float x,float y,float theta,float phi, bool isMC,float &betaCorrection,float &width){
  // Retrieve geom hash
  GeomHashBetaEnsemble *hash=isMC?mc:data;

  // Prepare input
  float vx=sin(theta)*cos(phi);
  float vy=sin(theta)*sin(phi);
  if(cos(theta)>0) {vx*=-1;vy*=-1;} // Keep a coherent definition

  // Evaluate
  hash->Eval(x,y,vx,vy);

  // Retrieve results
  betaCorrection=hash->MeanPeak;
  width=hash->MeanPeakWidth;
}


float GHBManager::correctBeta(float beta,float x,float y,float theta,float phi, int charge,bool isMC){
  // Retrieve the correction
  float peak,width;
  getParameters(x,y,theta,phi,isMC,peak,width);

  float result=beta/peak;
  
  // If it is MC, apply a gaussian smearing to improve the agreement of the width
  if(isMC){
    float c=charge*0.5;
    float peakD,widthD;
    getParameters(x,y,theta,phi,false,peakD,widthD);
    float sigma=widthD*widthD-width*width;
    //    sigma*=1.5; // Correction factor
    sigma/=c*c;
    if(sigma>0) result+=gRandom->Gaus()*sqrt(sigma);
  }

  return result;
}

#include "GeomHashes.h"
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <iostream>
#include <algorithm>
#undef _DEBUG_
#define ASSERT(x) 

using namespace std;

ClassImp(GeomHashNamed);

using namespace GeomHashNamedConstants;

float *GeomHashNamed::buffer[GeomHashNamedConstants::nBuffers]={0,0,0,0};
int GeomHashNamed::bufferSize=0;

bool GeomHashNamed::storeEntries=true;
bool GeomHashNamed::storeTemplates=true;
bool GeomHashNamed::storeTemplatesRms=true;
bool GeomHashNamed::storeMean=true;
bool GeomHashNamed::storeRms=true;
bool GeomHashNamed::storePeak=true;
bool GeomHashNamed::storePeakWidth=true;
bool GeomHashNamed::computePeakAsMean=false;
float GeomHashNamed::peakFinderFraction=0.68;
bool GeomHashNamed::sampleAllDirections=false;

void GeomHashNamed::checkBuffers(){
  if(dimension<=bufferSize) return;
  bufferSize=dimension;
  for(int i=0;i<GeomHashNamedConstants::nBuffers;i++){
    if(buffer[i]) delete[] buffer[i];
    buffer[i]=new float[dimension+1];
  }
}


GeomHashNamed::GeomHashNamed(int d){
  if(d>255){cout<<"Initializing to 255 variables."<<endl;d=255;}
  if(d<0) d=0;
  dimension=d;
  numNodes=0;
  points.reserve(1024);
  nodes[0].reserve(1024);
  nodes[1].reserve(1024);
  limit.reserve(1024);

  // Grow the buffers if needed
  //  checkBuffers();
}


int GeomHashNamed::hash(float *point){
  checkBuffers();
  if(numNodes==0) return 0;
  int current=0;
  for(;;){
    int nodeNumber=point[points[current]]<limit[current]?0:1;
    if(nodes[nodeNumber][current]<0) return -nodes[nodeNumber][current]-1;
    current=nodes[nodeNumber][current];
  }
}

void GeomHashNamed::push(float value,float *x){
  // Set the minimum capacity
  if(samples.capacity()<1024*dimension) samples.reserve(1024*dimension);

  // Store the point
  for(int i=0;i<dimension;i++)  samples.push_back(x[i]);

  values.push_back(value);
}

#include <stdarg.h>
void GeomHashNamed::fill(double x,...){
  checkBuffers();
  va_list ap;
  va_start(ap, x);
  float value=x;
  for(int i=0;i<dimension;i++) buffer[2][i]=va_arg(ap, double);
  va_end(ap);
  push(value,buffer[2]);
}

int GeomHashNamed::get(double x,...){
  checkBuffers();
  va_list ap;
  va_start(ap, x);
  buffer[2][0]=x;
  for(int i=1;i<dimension;i++) buffer[2][i]=va_arg(ap, double);
  va_end(ap);
  return hash(buffer[2]);
}

int *_WallN_;  // Just in case

void GeomHashNamed::grow(int min_size){
  points.clear();
  limit.clear();
  nodes[0].clear(); nodes[1].clear();
  int *pointers=new int[samples.size()/dimension];
  int *scratch=new  int[samples.size()/dimension];

  int size=samples.size()/dimension;
  _WallN_=pointers+size;
  grow(pointers,scratch,min_size);

  delete[] pointers;
  delete[] scratch;
}


void GeomHashNamed::grow(int *pointers,int *scratch,int min_size){
  int entries=samples.size()/dimension;

  // Init the buffer
  for(int i=0;i<entries;i++) pointers[i]=i;

  // Start the real growing process
  minSize=2*min_size;
  grow_internal(pointers,entries,scratch); // Adjust the min_size to obtain a better approximation

  // Liberates the space used by the samples and values
  samples.resize(0);
  values.resize(0);
  vector<float> empty; 
  samples=empty;
  values=empty;
}

void GeomHashNamed::store(int *pointers,int size,int parent){
  checkBuffers();
  int currentHashNumber=numNodes++;					
  int trueParent=abs(parent)-1;					
  int parentNode=parent<0?0:1;						
  nodes[parentNode].at(trueParent)=-currentHashNumber-1;
  
  // Here we can store the information concerning the data, for example mean and rms
  if(storeEntries) Entries.push_back(size);

  double sum=0,sum2=0;
  for(int i=0;i<size;i++){
    sum+=values[pointers[i]];
    sum2+=values[pointers[i]]*values[pointers[i]];
  }
  sum/=size;
  sum2/=size;
  sum2-=sum*sum;
  sum2=sqrt(fabs(sum2));

  if(storeMean) Means.push_back(sum);
  if(storeRms) Rms.push_back(sum2);

  // Fill the templates
  for(int i=0;i<dimension;i++) buffer[0][i]=buffer[1][i]=0;
  for(int j=0;j<dimension;j++){
    double sum=0;
    double sum2=0;
    for(int i=0;i<size;i++){
      int &element=pointers[i];
      float &v=samples[offset(element)+j];
      sum+=v;
      sum2+=v*v;
    }
    sum/=size;
    sum2/=size;
    sum2-=sum*sum;
    sum2=sqrt(fabs(sum2));
    if(storeTemplates) Template.push_back(sum);
    if(storeTemplatesRms) TemplateRms.push_back(sum2);
  }

  // Peak search and width determination
  vector<float> ordered;
  ordered.reserve(size);
  for(int i=0;i<size;i++) ordered.push_back(values[pointers[i]]);
  sort(ordered.begin(),ordered.end());
  int window=int(ceil(size*peakFinderFraction))-1; 
  if(window==0) window=1;

#ifdef _DEBUG_
  cout<<"PEAK FINDER FOR"<<endl;
  for(int i=0;i<size;i++) cout<<ordered[i]<<" ";
  cout<<endl;
  cout<<"Window size "<<window<<endl;
  cout<<"Sample size "<<size<<endl;
  cout<<"Mean "<<Means[Means.size()-1];
#endif


  float bestWidth=INFINITY;
  float peakPos=0;
  for(int i=0;i<size;i++){
    if(i+window>=size) break;
    if(i && ordered[i]==ordered[i-1]) continue; // skip identical events
    int finalIndex=i+window;
    for(int j=finalIndex+1;j<size;j++){
      if(ordered[j]!=ordered[finalIndex]) break;  
      finalIndex=j;
    }

    float width=ordered[finalIndex]-ordered[i];
    if(width<bestWidth){
      bestWidth=width;
      peakPos=0.5*(ordered[finalIndex]+ordered[i]);
#ifdef _DEBUG_
      cout<<"CURRENT BEST "<<peakPos<<" from "<<i<<" to "<<finalIndex<<" "<<ordered[i]<<" "<<ordered[finalIndex]<<endl;
#endif
    }
  }

  bestWidth/=2; // width meaning should be equivalent to sigma 

  if(computePeakAsMean){
    double sum=0;
    double sum2=0;
    int total=0; 
    for(int i=0;i<size;i++){
      if(ordered[i]>peakPos+bestWidth) continue;
      if(ordered[i]<peakPos-bestWidth) continue;
      sum+=ordered[i];
      sum2+=ordered[i]*ordered[i];
      total++;
    }
    sum/=total;
    sum2/=total;
    sum2-=sum*sum;
    sum2=sqrt(fabs(sum2));
    peakPos=sum;
    bestWidth=sum2;
  }

  if(storePeak) Peak.push_back(peakPos);
  if(storePeakWidth) PeakWidth.push_back(bestWidth);
}

float *GeomHashNamed::getTemplate(int h){
  if(!Template.size()) return 0;
  if(h<0 || h>=numNodes) return 0;
  checkBuffers();
  for(int i=0;i<dimension;i++) buffer[0][i]=Template[offset(h)+i];
  return buffer[0];
}

float *GeomHashNamed::getTemplateRms(int h){
  if(!TemplateRms.size()) return 0;
  if(h<0 || h>=numNodes) return 0;
  checkBuffers();
  for(int i=0;i<dimension;i++) buffer[1][i]=TemplateRms[offset(h)+i];
  return buffer[1];
}


float GeomHashNamed::getMean(int node){
  if(!Means.size()) return 0;
  if(node>=numNodes) return 0;
  return Means.at(node);
}

float GeomHashNamed::getRms(int node){
  if(!Rms.size()) return 0;
  if(node>=numNodes) return 0;
  return Rms.at(node);
}


float GeomHashNamed::getPeak(int node){
  if(!Peak.size()) return 0;
  if(node>=numNodes) return 0;
  return Peak.at(node);
}

float GeomHashNamed::getPeakWidth(int node){
  if(!PeakWidth.size()) return 0;
  if(node>=numNodes) return 0;
  return PeakWidth.at(node);
}

int GeomHashNamed::getEntries(int node){
  if(!Entries.size()) return 0;
  if(node>=numNodes) return 0;
  return Entries.at(node);
}


void GeomHashNamed::grow_internal(int *pointers,int size,int *scratch,int parent){
  const float epsilon=5e-6;
  // Check if we should proceed
  if(size==1 || size<minSize) {store(pointers,size,parent);return;}

  // Fast search binning
  const int numBins=1000;  // Number of bins to be used during spliting
  //const int numBins=10;  // Number of bins to be used during spliting /* CJC */
  int entries[numBins];
  double mean[numBins];
  double mean2[numBins]; 
  for(int i=0;i<numBins;i++) entries[i]=0;
  for(int i=0;i<numBins;i++) mean[i]=0;
  for(int i=0;i<numBins;i++) mean2[i]=0;

  int me=nodes[0].size();

  // Pick randomly a direction
  int direction=int(dimension*(rand() / (RAND_MAX + 1.0)));

  int tries=0;
 loop:
  tries++;

  // Select the sorting range
  double minValue=INFINITY;
  double maxValue=-INFINITY;

  // DEBUG
  int minP=-1;
  int maxP=-1;

  for(int index=0;index<size;index++){
    int element=pointers[index];
    float x=samples[offset(element)+direction];
    if(x<minValue) {minValue=x;minP=index;}
    if(x>maxValue) {maxValue=x;maxP=index;}
  }

  // No binning
  if(fabs((maxValue-minValue)/((maxValue+minValue)/2))<epsilon) {
    if (sampleAllDirections && tries<dimension) goto loop;
    else {store(pointers,size,parent);return;}
  }

  double bw=(maxValue-minValue)/numBins;
  maxValue+=bw/2;
  minValue-=bw/2;
  bw=(maxValue-minValue)/numBins;

  // No binning
  if(bw<=0) {
    if (sampleAllDirections && tries<dimension) goto loop;
    else {store(pointers,size,parent);return;}
  }

  // Fill the arrays
  double totalMean=0;
  double totalMean2=0;
  bool equals=true;

  double prevValue=values[0];
  for(int index=0;index<size;index++){
    int element=pointers[index];

    float x=samples[offset(element)+direction];
    float value=values[element];
    float value2=value*value;
    int bin=int(floor((x-minValue)/bw));

    if(value!=prevValue) equals=false;

    entries[bin]++;
    mean[bin]+=value;
    mean2[bin]+=value2;
    totalMean+=value;
    totalMean2+=value2;
  }


  // Search the best split position
  const double rms01=totalMean2/size-totalMean/size*totalMean/size;
  const int n01=size;

  int accEntries=0;
  double acc=0;
  double acc2=0;
  float bestValue=FP_INFINITE;
  float splitPoint;
  int splitCounter=0;
  bool fail=true;

  if(!equals)
  for(int i=0;i<numBins-1;i++){
    accEntries+=entries[i];
    acc+=mean[i];
    acc2+=mean2[i];
    if(size-accEntries<minSize/2) break;
    if(accEntries<minSize/2) continue;

    int n0=accEntries;
    double rms0=acc2/n0-acc/n0*acc/n0;
    int n1=size-n0; 
    double rms1=(totalMean2-acc2)/n1-(totalMean-acc)/n1*(totalMean-acc)/n1; 
    double value=n0*rms0+n1*rms1-n01*rms01;

    if(fabs(value/n01*rms01)<-epsilon)  // Minimum improvement to split /*CJC BUG??*/
      if(value<bestValue){
	bestValue=value;
	splitCounter=n0;
	splitPoint=minValue+bw*(i+1);
	fail=false;
      }
  }

  if(fail){
    // Try the standard split
    int accEntries=entries[0];
    for(int i=1;i<numBins-1;i++){
      //if(accEntries+entries[i]>size/2){ /* CJC */
      if(accEntries+entries[i]>=size/2){
	splitPoint=minValue+bw*(i+1);
	splitCounter=accEntries;
	fail=false;
	break;
      }
      accEntries+=entries[i];
    }
    // No way, store it a go ahead
    if(fail) {store(pointers,size,parent);return;}
  }

  if(splitCounter==size) {store(pointers,size,parent);return;}


  // Copy the pointer to the scratch area to sort them 
  for(int index=0;index<size;index++) scratch[index]=pointers[index];

  // Sort
  int lefties=0,righties=size-1;
  splitCounter=0;
  for(int index=0;index<size;index++){
    int element=scratch[index];
    float x=samples[offset(element)+direction];
    if(x<splitPoint){
      pointers[lefties++]=scratch[index];
      splitCounter++;
    }else{
      pointers[righties--]=scratch[index];
    }
  }


  // Store the current information
  limit.push_back(splitPoint);
  points.push_back(direction);

  nodes[0].push_back(0);  nodes[1].push_back(0);
  nodes[0][me]=nodes[0].size();
  grow_internal(pointers,splitCounter,scratch,-(me+1));
  nodes[1][me]=nodes[1].size();
  grow_internal(pointers+splitCounter,size-splitCounter,scratch,me+1);
}


ClassImp(GeomHashNamedEnsemble);


void GeomHashNamedEnsemble::growOne(int minSize,bool bootstrap){
  hashes.push_back(GeomHashNamed(dimension));
  GeomHashNamed &hash=hashes.back();

  if(!bootstrap){
    hash.samples=samples;
    hash.values=values;
  }else{
    checkBuffers();
    for(int i=0;i<values.size();i++){
      int pointer=int(values.size()*(rand() / (RAND_MAX + 1.0)));
      for(int j=0;j<dimension;j++) buffer[0][j]=samples[offset(pointer)+j];
      hash.push(values[pointer],buffer[0]);
    }
  }
  hash.grow(minSize);
}


void GeomHashNamedEnsemble::Eval(float *x){
  MeanValue=0;
  ValueRms=0;
  MeanRms=0;
  RmsRms=0;
  MeanPeak=0;
  PeakRms=0;
  MeanPeakWidth=0;
  PeakWidthRms=0;
  Hashes=0;
  
  for(int i=0;i<hashes.size();i++){
    Hashes++;
    GeomHashNamed &hash=hashes[i];
    int h=hash.hash(x);
    double xx=hash.getMean(h);
    MeanValue+=xx;
    ValueRms+=xx*xx;

    xx=hash.getRms(h);
    MeanRms+=xx;
    RmsRms+=xx*xx;

    xx=hash.getPeak(h);
    MeanPeak+=xx;
    PeakRms+=xx*xx;

    xx=hash.getPeakWidth(h);
    MeanPeakWidth+=xx;
    PeakWidthRms+=xx*xx;

    xx=hash.getEntries(h);
    MeanEntries+=xx;
    EntriesRms+=xx*xx;
  }

#define Do(_x) {Mean##_x/=Hashes; _x##Rms /=Hashes; _x##Rms-=Mean##_x*Mean##_x; _x##Rms=sqrt(fabs(_x##Rms));}
  Do(Value);
  Do(Rms);
  Do(Peak);
  Do(PeakWidth);
  Do(Entries);
#undef Do
}

void GeomHashNamedEnsemble::Eval(double x,...){
  checkBuffers();
  va_list ap;
  va_start(ap, x);
  buffer[2][0]=x;
  for(int i=1;i<dimension;i++) buffer[2][i]=va_arg(ap, double);
  va_end(ap);
  Eval(buffer[2]);
}

int GeomHashNamedEnsemble::numHashes(){
  return hashes.size();
}

GeomHashNamed& GeomHashNamedEnsemble::getHash(int i){
  return hashes.at(i);
}


