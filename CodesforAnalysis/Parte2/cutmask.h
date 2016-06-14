#ifndef CUTMASK_H
#define CUTMASK_H

#include <iostream>
using namespace std;

class Cutmask{
	public:
		Cutmask() {cmask=0;}
		Cutmask(int cutmask)   : cmask(cutmask) {}
		Cutmask(float cutmask) : cmask((int)(cutmask+0.5)) {}
		void setMask (int cutmask)   {cmask=cutmask;}
		void setMask (float cutmask) {cmask=(int)(cutmask+0.5);}
		int getMask () {return cmask;}
		
		bool isRichMeasureFromNaF() {return cmask>>11 == 512;}
		bool isRichMeasureFromAgl() {return cmask>>11 == 0;}
		bool notPassed(uint s);
		bool passed (uint s);
		void print();

	private:
		int cmask;
	
};

void Cutmask::print() {
	cout << "Cutmask description:" << endl;
	if (isRichMeasureFromNaF()) cout << "RICH measure taken from NaF" <<endl;
}

bool Cutmask::notPassed(uint i) {
	if (i>2) return 0;
	int notpassed[3]= {155,139,11};
	return (cmask&notpassed[i])==notpassed[i];
}

bool Cutmask::passed(uint i) {
	if (i>2) return 0;
	int passed[3]= {187,155,139};
	return (cmask&passed[i])==passed[i];
}


#endif /* CUTMASK_H */ 
