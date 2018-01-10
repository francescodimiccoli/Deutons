#ifndef FILESAVER_H
#define FILESAVER_H

#include <iostream>
#include <string>
#include "TFile.h"
#include "TObject.h"
#include "TObjArray.h"

using namespace std;

class FileSaver {
   public:
      FileSaver(bool Isfinal=false) : filename("") {fArr=new TObjArray(); IsFinal = Isfinal;}
      FileSaver (string fname,bool Isfinal=false) : filename (fname) {fArr=new TObjArray(); IsFinal = Isfinal;}
      void writeObjsInFolder(string folder, bool recreate = false);
      void writeObjs();
      void setName(string fname) {filename=fname;}
      void Add(TObject* obj) {fArr->Add(obj);}
      string setName() {return filename;}
      bool CheckFile();	      
      TObject * Get(std::string name);	
      TFile * GetFile();
	
   private:
      string filename;
      TObjArray* fArr;
      bool IsFinal;
};

#endif
