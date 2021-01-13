#include "InputFileReader.h"

using namespace std;

TChain * InputFileReader(std::string listfilename,std::string treename)
{

	TChain * chain=new TChain(treename.c_str());
	std::ifstream infile(listfilename.c_str());
	std::string line;
	int max =0;
	while (std::getline(infile, line))
		{
    		max++;
		std::istringstream iss(line);
    		std::cout<<line.c_str()<<std::endl;
		if(max<=100)chain->Add(line.c_str());
		else{cout<<"Warning: File List too long!"<<endl;}	
	}
	cout<<"Chain: "<<chain<<endl;
	cout<<chain->GetEntries()<<endl;
	if(chain->GetEntries()==0) return 0; 
	else return chain;
}


