#include <iostream>
#include "Parte2/Binning.h"
#include "TH1F.h"
#include "TCanvas.h"

int main(int argc, char **argv)
{
	DBinning mrbin;
	mrbin.Setbins(10, 0.2, 11, 2);
	mrbin.Print();
	cout << endl;
	std::vector<float> vec=mrbin.Rebin(mrbin.LoadCRDB("AMSProtons.dat"));
	std::vector<float> mcweight=mrbin.Rebin(mrbin.LoadData("MCweight.data"));
	float sum=0;
	for (int i=0; i<mcweight.size(); i++) {
		cout << mcweight[i]	<< " ";
		sum += mcweight[i]*(mrbin.RigBins()[i+1]-mrbin.RigBins()[i])/mcweight.size();
	}
	cout << endl;
	cout << "Sum vector : " << sum << " inverse " << 1/sum << endl;

	TCanvas c1;
	TH1F h("h", "h", mcweight.size(), mrbin.RigBins().data());
	for (int i=0; i<mcweight.size(); i++)
		h.SetBinContent(i, mcweight[i]);

	h.Draw();
	c1.SaveAs("TestBin.png");
	
	
	return 0;
}

