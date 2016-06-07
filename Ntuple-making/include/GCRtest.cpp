#include "GCR_data.h"

int main() {
   GCR_data galprop;
   galprop.read("galprop.dat", "m2", "GeV");
   return 0;
}
