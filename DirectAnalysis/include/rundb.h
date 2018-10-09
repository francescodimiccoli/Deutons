#ifndef rundb_h
#define rundb_h
#include <map>
#include <cstdlib>
#include  <cstdio>
#include <cmath>
#include <iostream>

class rundb_el{
public:

  int run;
  float Rpmin;
  float Rpmax;
  int RTrig;
  float  pmin;
  float pmax;
  int Trig;
  int Events;

  rundb_el(){
    run=0;
    Rpmin=0;
    Rpmax=0;
    RTrig=0;
    pmin=0;
    pmax=0;
    Trig=0;
    Events=0;
  }  
  virtual ~rundb_el(){}
  float get_evnorm()  const{
    float ll,uu;
    if(pmin!=0) ll=pmin; else ll=Rpmin; 
    if(pmax!=0) uu=pmax; else uu=Rpmax; 
    return Trig/log(uu/ll);
  }
  int GetTrig() {return Trig;}
  int GetEv() {return Events;} 
  void Print(){
    printf("%11d %6.1f %6.1f  %11d %9.3f %9.3f %7d %9d  %7.3f\n",
	   run,Rpmin,Rpmax,RTrig, pmin,pmax,Trig,Events,get_evnorm());
    
  }
};


class rundb{
public:
  std::map<int,rundb_el> db;

  void Add(const rundb_el&  el){
    db[el.run]=el;
  }

  rundb_el* find(int run){
    std::map<int,rundb_el>::iterator it;
    it=db.find(run);
    if(it==db.end())
      return 0;
    else return &(it->second);
  }
  int size(){return db.size();}
  int readdb(const char* fname);
  void Print(){
    std::map<int,rundb_el>::iterator it;
    it=db.begin();
    while(it!=db.end())
      {it->second.Print();it++;}
  }

  void Summary(){
    std::map<int,rundb_el>::iterator it;
    it=db.begin();
    double aa=0,aa2=0;
    int nn=0;
    while(it!=db.end())
      {
	aa+=it->second.get_evnorm();
	aa2+=pow(it->second.get_evnorm(),2);
	nn++;
	it++;}
    double av=aa/nn;
    double rms=sqrt(aa2/nn-pow(av,2));
    double err=rms/sqrt(nn);
    double rrms=rms/av;
    double rerr=err/av;
    printf("NN: %d Average: %7.1f [ev/ln(R)]  RMS: %7.1f (%4.1f\%)   Error: %7.2f (%5.3f\%)\n",nn,av,rms,rrms*100,err,rerr*100);

  }
  float GetTrigRate(){
    std::map<int,rundb_el>::iterator it;
    it=db.begin();
    double trigrate=0;     
    float nn=0;
    while(it!=db.end())
    {
           trigrate+= it->second.GetEv()/(float)it->second.GetTrig();
		nn++;
		it++;
	}
    return trigrate/nn;
  }    

  float GetTotEvents(){
    std::map<int,rundb_el>::iterator it;
    it=db.begin();
    long int totev = 0;     
    while(it!=db.end())
    {
		    totev+= it->second.GetEv();
		    it++;
	}
    return totev;
  }    

  float GetTotTrigg(){
    std::map<int,rundb_el>::iterator it;
    it=db.begin();
    long int tottrig = 0;     
    while(it!=db.end())
    {
       	   
		    tottrig+= it->second.GetTrig();
		    it++;
	     
    }
    return tottrig;
  }    



  rundb(){}
  virtual ~rundb(){db.clear();}
};


#endif
