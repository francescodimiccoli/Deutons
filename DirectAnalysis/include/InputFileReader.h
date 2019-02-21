#ifndef INPUTFILEREADER_H
#define INPUTFILEREADER_H

#include <TChain.h>
#include <iostream>
#include <fstream>
#include <sstream>

TChain * InputFileReader(std::string listfilename,std::string treename);


#endif
