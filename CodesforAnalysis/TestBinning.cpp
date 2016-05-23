#include <iostream>
#include "Parte2/Binning.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TGraph.h"

int main(int argc, char **argv)
{
   // General creation and info dumping
   DBinning mrbin;
   mrbin.Setbins(10, 0.2, 11, 2);
   mrbin.Print();
   cout << endl;

   // Checking taht the inverse of the sum of bins is around 50
   std::vector<float> mcweight=mrbin.Rebin(mrbin.LoadData("MCweight.data"));
   float sum=0;
   for (int i=0; i<mcweight.size(); i++) {
      cout << mcweight[i]	<< " ";
      sum += mcweight[i]*(mrbin.RigBin(i+1)-mrbin.RigBin(i))/mcweight.size();
   }
   cout << endl;
   cout << "Sum vector : " << sum << " inverse " << 1/sum << endl;

   // Plotting
   TCanvas c1;
   TH1F h("h", "h", mcweight.size(), mrbin.RigBins().data());
   for (int i=0; i<mcweight.size(); i++)
      h.SetBinContent(i, mcweight[i]);

   h.Draw();
   c1.SaveAs("TestBin.png");

   // Proton flux in 4 different quantities.

   c1.Clear();
   c1.SetLogx();
   PBinning pb;
   pb.Setbins(100, 0.5, 100);
   std::vector<float> protflux=pb.Rebin(mrbin.LoadCRDB("AMSProtons.dat"));


   TGraph g[4];
   g[0]=TGraph(pb.size(), pb.  EkBins().data(), protflux.data());
   g[1]=TGraph(pb.size(), pb. MomBins().data(), protflux.data());
   g[2]=TGraph(pb.size(), pb. RigBins().data(), protflux.data());
   g[3]=TGraph(pb.size(), pb.BetaBins().data(), protflux.data());

   g[0].SetLineColor(4);
   g[1].SetLineColor(1);
   g[2].SetLineColor(2);
   g[3].SetLineColor(3);

   g[0].Draw("apl");
   g[1].Draw("pl");
   g[2].Draw("pl");
   g[3].Draw("pl");

   c1.SaveAs("TestCanvas.png");
   

   return 0;
}

