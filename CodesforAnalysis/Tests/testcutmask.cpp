#include "../Parte2/cutmask.h"


int main(int argc, char * argv[])
{
   Cutmask cmask;
   cmask.print();
   for (int i=0; i<200; i++)
      {
         cmask.setMask(i);
         cmask.print();
      }
      
   return 0;
}
