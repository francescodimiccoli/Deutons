

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <iterator>


using namespace std;


bool IsProtonMC    (Variables * vars){ return (vars->Massa_gen<1&&vars->Massa_gen>0);}
bool IsDeutonMC    (Variables * vars){ return (vars->Massa_gen<2&&vars->Massa_gen>1);}
bool IsData        (Variables * vars){ return (vars->Massa_gen==0);}
bool IsPreselected (Variables * vars){ return (vars->joinCutmask&187)==187;}
bool IsOnlyFromToF (Variables * vars){ return ((vars->joinCutmask>>11)!=512 && (vars->joinCutmask>>11)!=0 );}
bool IsFromNaF     (Variables * vars){ return ((vars->joinCutmask>>11)==512);}
bool IsFromAgl     (Variables * vars){ return ((vars->joinCutmask>>11)==0);}
bool L1LooseCharge1(Variables * vars){ return (vars->qL1>0 && vars->qL1<1.6);} 

template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}


bool ApplyCuts(std::string cut, Variables * Vars){

	std::vector<std::string> spl = split(cut,'&');
	
	bool IsPassed = true;
	for(int i=0;i<spl.size();i++){
		if(spl[i]=="IsProtonMC"	   ) IsPassed=IsPassed && IsProtonMC    (Vars);
		if(spl[i]=="IsDeutonMC"	   ) IsPassed=IsPassed && IsDeutonMC    (Vars);
		if(spl[i]=="IsData"	   ) IsPassed=IsPassed && IsData        (Vars);
		if(spl[i]=="IsPreselected" ) IsPassed=IsPassed && IsPreselected (Vars); 
		if(spl[i]=="IsOnlyFromToF" ) IsPassed=IsPassed && IsOnlyFromToF (Vars);
		if(spl[i]=="IsFromNaF"	   ) IsPassed=IsPassed && IsFromNaF     (Vars);
		if(spl[i]=="IsFromAgl"	   ) IsPassed=IsPassed && IsFromAgl     (Vars);
		if(spl[i]=="L1LooseCharge1") IsPassed=IsPassed && L1LooseCharge1(Vars);

	}

	return IsPassed;
}



