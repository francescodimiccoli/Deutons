#ifndef FILEREADER_H
#define FILEREADER_H

#include "TChain.h"
#include <fstream>
#include <sstream>

TChain * InputFileReader(std::string listfilename,std::string treename){

	TChain * chain=new TChain(treename.c_str());
	std::ifstream infile(listfilename.c_str());
	std::string line;
	while (std::getline(infile, line))
		{
    		std::istringstream iss(line);
    		std::cout<<line.c_str()<<std::endl;
		chain->Add(line.c_str());
	}
	std::cout<<chain->GetEntries()<<std::endl;
	return chain;
}


#endif
