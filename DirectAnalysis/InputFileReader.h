
#include <TChain.h>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

TChain * InputFileReader(std::string listfilename,std::string treename)
{

	TChain * chain=new TChain(treename.c_str());
	std::ifstream infile(listfilename.c_str());
	std::string line;
	while (std::getline(infile, line))
		{
    		std::istringstream iss(line);
    		std::cout<<line.c_str()<<std::endl;
		chain->Add(line.c_str());
	}
	cout<<"Chain: "<<chain<<endl;
	cout<<chain->GetEntries()<<endl;
	if(chain->GetEntries()==0) return 0; 
	else return chain;
}


