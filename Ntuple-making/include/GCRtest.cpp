#include "Reweighting.cpp"

int main() {
   Particle proton(0.9382720813); // proton mass 938 MeV
   Particle deuton(0.18756129);   // deuterium mass 1876 MeV
   Binning bins(proton);
   GCR_data galprop;
   galprop.read("galprop.dat", "m2", "GeV");
   return 0;
}
