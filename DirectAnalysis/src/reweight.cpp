#include "reweight.h"

int FindTimeIndex(std::string input){

	input.erase(0,input.find("1"));	
	input.erase(input.find("-"),input.length());	

	std::cout<<"TIME STRING: "<<input.c_str()<<std::endl;
	
	return std::stoi(input);
}	

Reweighter::Reweighter(Histogram f, Histogram t): from(f), to(t) {}

float Reweighter::getWeight(float value, float def) const {
    float fromV = from.at(value);
    float   toV =   to.at(value);
    if(fromV == 0.0) return def;
    return toV / fromV;
}


