#ifndef CUTMASK_H
#define CUTMASK_H

#include <iostream>
#include <iomanip>
using namespace std;

class Cutmask {
   public:
      Cutmask() {cmask=0;}
      Cutmask (int cutmask)   : cmask (cutmask) {}
      Cutmask (float cutmask) : cmask ( (int) (cutmask+0.5) ) {}
      void setMask (int cutmask)   {cmask=cutmask;}
      void setMask (float cutmask) {cmask= (int) (cutmask+0.5);}
      int getMask () {return cmask;}


      bool isMinimumBiasTrigger()       { return cmask& (1<<0);}
      bool isMinimumBiasToF3or4Layers() { return cmask& (1<<1);}
      bool isMinimumBiasTRD()           { return cmask& (1<<2);}
      bool isMinimumBiasTracker()       { return cmask& (1<<3);}
      bool isGoldenTracker()            { return cmask& (1<<4);}
      bool isGoldenToF3or4Layers()      { return cmask& (1<<5);}
      bool isGoldenTRD()                { return cmask& (1<<6);}
      bool hasSingleTrTrack()           { return cmask& (1<<7);}
      bool isMinimumBiasToF4Layers()    { return cmask& (1<<8);}
      bool isGoldenToF4Layers()         { return cmask& (1<<9);}
      // Bit 10 is always set
      bool isFromNaF() {return true;/*(cmask>>11) == 512;*/}
      bool isFromAgl() {return true;/*(cmask>>11) == 0;*/}
      bool isOnlyFromToF() {return ( !isFromAgl() && !isFromNaF() );}

      bool isPreselected(); ///< This 187 thing

      bool notPassed (uint s);
      bool passed (uint s);
      void print();

   private:
      int cmask;
      void printANSIYesNo (string legend, bool condition, int width=25);

};


bool Cutmask::isPreselected() {
	return isMinimumBiasTracker()
	    && isMinimumBiasToF3or4Layers()
	    && isMinimumBiasTracker()
	    && isGoldenTracker()
	    && isGoldenToF3or4Layers()
	    && hasSingleTrTrack()
	    //&& isMinimumBiasTRD()
            //&& isGoldenTRD()
            ;		
}

	

void Cutmask::print()
{
   const int nbits=10;
   string legend[nbits]= {
      "Min. bias trigger",
      "Min. bias ToF >=3 lay. ",
      "Min. bias TRD",
      "Min. bias Tracker",
      "Golden Tracker",
      "Golden ToF >=3 lay.",
      "Golden TRD",
      "Single Tracker track",
      "Min. bias ToF 4 lay.",
      "Golden ToF 4 layers"
   };

	cout << "Cutmask description:" << endl;
   for (int ibit=0; ibit<nbits; ibit++) {
      printANSIYesNo (legend[ibit], (cmask& (1<<ibit)) );
   }
   printANSIYesNo ("RICH meas. from NaF", isFromNaF());
   printANSIYesNo ("RICH meas. from Agl", isFromAgl());
   printANSIYesNo ("Preselected (code 187)", isPreselected());
   
	return;
}


void Cutmask::printANSIYesNo (string legend, bool condition, int width) {
	std::cout << std::right << std::setw (width) << std::setfill (' ') << legend;
   if (condition) cout << "\033[32m YES \033[0m" << endl;
   else        	cout << "\033[31m NO  \033[0m" << endl;
   return;
}

bool Cutmask::notPassed (uint i)
{
   if (i>2) return 0;
   int notpassed[3]= {139,171,11};
   return (cmask&notpassed[i]) ==notpassed[i];
}

bool Cutmask::passed (uint i)
{
   if (i>2) return 0;
   int passed[3]= {171,187,139};
   return (cmask&passed[i]) ==passed[i];
}


#endif /* CUTMASK_H */
