#ifndef CUTMASK_H
#define CUTMASK_H

#include <iostream>
using namespace std;

class Cutmask{
	public:
		Cutmask(int cutmask)   : cmask(cutmask) {}
		Cutmask(float cutmask) : cmask((int)(cutmask+0.5)) {}
		bool isRichMeasureFromNaF() {return (cmask>>11==512);}
		void print();

	private:
		int cmask;
	
};

void Cutmask::print() {
	cout << "Cutmask description:" << endl;
	if (isRichMeasureFromNaF()) cout << "RICH measure taken from NaF" <<endl;
	
	
}



#endif /* CUTMASK_H */ 
