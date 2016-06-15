#include "../Parte2/cutmask.h"


int main(int argc, char * argv[])
{
   Cutmask cmask;
   //cmask.print();
   for (int i : {1, 3, 11})
      {
         cout << "####### " << i << endl;
         cmask.setMask(i);
         cmask.print();
      }
      
   return 0;
}
