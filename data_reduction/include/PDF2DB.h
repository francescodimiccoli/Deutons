#ifndef __PDF2DB_h__
#define __PDF2DB_h__

#include "PDF2.h"

#include "TNamed.h"

#include <cstdio>
#include <cmath>
#include <map>
#include <string>

//! Class for PDF2 handling.
class PDF2DB : public TNamed {

 public:

  //! PDF2 map
  std::map<std::string,PDF2*> Pdf2Map;

  //! c-tor
  PDF2DB() { Pdf2Map.clear(); }
  //! d-tor
  virtual ~PDF2DB() { Clear(); }
  //! clear
  void Clear(Option_t *opt="");
  //! print
  void Print(Option_t *opt="") const;
  //! Draw
  void Draw(Option_t *opt="");
  //! add
  void Add(std::string name, PDF2* pdf);
  //! get
  PDF2* Get(std::string name);

  ClassDef(PDF2DB,1);
};

#endif
