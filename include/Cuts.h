

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <iterator>


using namespace std;



bool IsProtonMC    (Variables * vars){ return (vars->Massa_gen<1&&vars->Massa_gen>0);}
bool IsDeutonMC    (Variables * vars){ return (vars->Massa_gen<2&&vars->Massa_gen>1);}
bool IsPrimary	   (Variables * vars){ return (vars->R>1.0*vars->Rcutoff); }
bool IsMC          (Variables * vars){ return (vars->Massa_gen>0);} 
bool IsData        (Variables * vars){ return (vars->Massa_gen==0);}
bool IsPreselected (Variables * vars){ return (((int)vars->joinCutmask&187)==187&&(vars->R_L1>0||vars->EdepL1>0));}
bool IsFromNaF     (Variables * vars){ return (((int)vars->joinCutmask>>11)==512&&vars->BetaRICH_new>0);}
bool IsFromAgl     (Variables * vars){ return (((int)vars->joinCutmask>>11)==0&&vars->BetaRICH_new>0);}
bool IsOnlyFromToF (Variables * vars){ return !((IsFromNaF(vars))||(IsFromAgl(vars)));}
bool L1LooseCharge1(Variables * vars){ return (vars->qL1>0 && vars->qL1<1.6);} 
bool ProtonsMassCut(Variables * vars){ return GetRecMassTOF(vars)>0.5&&GetRecMassTOF(vars)<1.5;}


bool Qualitycut(Variables * vars, float cutvariable, float cutTOF, float cutNaF, float cutAgl){

        bool IsQual=false;
        if(IsOnlyFromToF(vars) && cutvariable<cutTOF)  IsQual=true;
        if(IsFromNaF(vars)     && cutvariable<cutNaF)  IsQual=true;
        if(IsFromAgl(vars)     && cutvariable<cutAgl)  IsQual=true;

        return IsQual;
}

bool DistanceCut   (Variables * vars){ return (Qualitycut(vars,vars->DistP,3,4,4)||Qualitycut(vars,vars->DistD,3,4,4));}
bool LikelihoodCut (Variables * vars){ return  Qualitycut(vars,log(1-vars->Likelihood),-1.2,-2.6,-3.2);}


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
		if(spl[i]=="IsPrimary"	   ) IsPassed=IsPassed && IsPrimary     (Vars);
		if(spl[i]=="IsMC"          ) IsPassed=IsPassed && IsMC          (Vars);
		if(spl[i]=="IsData"	   ) IsPassed=IsPassed && IsData        (Vars);
		if(spl[i]=="IsPreselected" ) IsPassed=IsPassed && IsPreselected (Vars); 
		if(spl[i]=="IsOnlyFromToF" ) IsPassed=IsPassed && IsOnlyFromToF (Vars);
		if(spl[i]=="IsFromNaF"	   ) IsPassed=IsPassed && IsFromNaF     (Vars);
		if(spl[i]=="IsFromAgl"	   ) IsPassed=IsPassed && IsFromAgl     (Vars);
		if(spl[i]=="L1LooseCharge1") IsPassed=IsPassed && L1LooseCharge1(Vars);
		if(spl[i]=="DistanceCut")    IsPassed=IsPassed && DistanceCut(Vars);
		if(spl[i]=="LikelihoodCut")  IsPassed=IsPassed && LikelihoodCut(Vars);
		if(spl[i]=="ProtonsMassCut") IsPassed=IsPassed && ProtonsMassCut(Vars);
	}

	return IsPassed;
}



