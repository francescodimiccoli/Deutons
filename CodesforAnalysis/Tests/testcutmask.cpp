#include "../Parte2/cutmask.h"


int main(int argc, char * argv[])
{
   Cutmask cmask;
   cmask.print();
   for (int i=0; i<200; i++)
      for (uint j=0; j<3; j++) {
         cmask.setMask(i);
         bool notpassed=cmask.notPassed(j);
         bool passed=cmask.passed(j);
         if (passed )    cout << i << " " << j << " passed" << endl;
         if (notpassed)  cout << i << " " << j << " notpassed" << endl;
      }
      

   return 0;
}
