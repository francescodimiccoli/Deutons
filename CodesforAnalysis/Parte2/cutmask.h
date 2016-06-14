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
		void print();

	private:
		int cmask;
	
};

void Cutmask::print() {
	cout << "Cutmask description:" << endl;
	if (isRichMeasureFromNaF()) cout << "RICH measure taken from NaF" <<endl;
}



#endif /* CUTMASK_H */ 
