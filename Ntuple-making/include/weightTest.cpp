#include "reweight.h"


// Preferred compilation command :
// g++ -g -std=c++11 $(root-config --libs --cflags)  -o wTest weightTest.cpp GCR_data.cpp -Wno-pointer-arith  -fmax-errors=5


int main() {
   // Inits
   Reweight weight;
   std::vector<float> MoD=weight.getMCOverData();
   std::vector<float> DoM=weight.getDataOverMC();
   printMatrix::print(
            {weight.getRigBins(), MoD, DoM},
            {"Bin cent", "MC/data", "data/MC"}
             );
   
   return 0;
}
