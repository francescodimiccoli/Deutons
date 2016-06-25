#include "PlottingFunctions/DATApreSeleff_Plot.h"

using namespace std;


LATcorr *LATpreSelDATA = new LATcorr("LATpreSelDATA",3);


void DATApreSeleff_Fill(int zona)
{

	if(Tup.Unbias!=0||Tup.R_pre<=0||Tup.Beta_pre>protons->Eval(Tup.R_pre)+0.1||Tup.Beta_pre<protons->Eval(Tup.R_pre)-0.1) return;
	if(!(Tup.EdepL1>0&&Tup.EdepL1<EdepL1beta->Eval(Tup.Beta_pre)+0.1&&Tup.EdepL1>EdepL1beta->Eval(Tup.Beta_pre)-0.1)) return;
	if(!Herejcut) return;  

	for(int S=0; S<3; S++) {
		int Kbin=PRB.GetRBin(fabs(Tup.R_pre));
		if(cmask.notPassed(S)) ((TH3 *)LATpreSelDATA->beforeR)->Fill(Kbin,zona,S);
		if(cmask.passed(S))    ((TH3 *)LATpreSelDATA->afterR )->Fill(Kbin,zona,S);
	}
	return;
}


void DATApreSeleff_Write()
{
	LATpreSelDATA->Write();
	return;
}


void DATApreSeleff(string filename)
{
	cout<<"********************** DATA PRESELECTIONS EFFICIENCIES ****************************************"<<endl;

	cout<<"*** Reading  P1 file ****"<<endl;
        TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");

	LATcorr * LATpreSelDATA = new LATcorr(inputHistoFile,"LATpreSelDATA",3);


	cout<<"********************** DATA PRESELECTIONS EFFICIENCIES ****************************************"<<endl;

	LATpreSelDATA -> Eval_Efficiency();
	TH3F *LATpreSelDATA_R = (TH3F *) LATpreSelDATA     -> effR -> Clone();
	cout<<"********************** LAT. Eff. CORRECTION **************************************************"<<endl;

	LATpreSelDATA -> Eval_LATcorr(3);

	// (LATpreSelDATA   -> LATcorrR))
	TH2F * preSelLATcorr (static_cast<TH2F *>(LATpreSelDATA   -> LATcorrR));
	TH2F * preSelLATcorr_fit(static_cast<TH2F *>(LATpreSelDATA   -> LATcorrR_fit));

	cout<<"*** Updating P1 file ****"<<endl;

	cout << LATpreSelDATA_R ->GetName() << " : " <<  LATpreSelDATA_R ->GetEntries() << endl;
	cout << preSelLATcorr ->GetName() << " : " <<  preSelLATcorr ->GetEntries() << endl;
	cout << preSelLATcorr_fit ->GetName() << " : " <<  preSelLATcorr_fit ->GetEntries() << endl;
	cout << ((TH2F *) LATpreSelDATA   -> LATcorrR_fit)->GetName() << " : " <<
		((TH2F *) LATpreSelDATA   -> LATcorrR_fit)->GetEntries() << endl;
	cout << ((TH2F *) LATpreSelDATA   -> LATcorrR)->GetName() << " : " <<
		((TH2F *) LATpreSelDATA   -> LATcorrR)->GetEntries() << endl;



	cout << LATpreSelDATA_R ->GetName() << " : " <<  LATpreSelDATA_R ->GetEntries() << endl;
	cout << preSelLATcorr ->GetName() << " : " <<  preSelLATcorr ->GetEntries() << endl;
	cout << preSelLATcorr_fit ->GetName() << " : " <<  preSelLATcorr_fit ->GetEntries() << endl;
	cout << ((TH2F *) LATpreSelDATA   -> LATcorrR_fit)->GetName() << " : " <<
		((TH2F *) LATpreSelDATA   -> LATcorrR_fit)->GetEntries() << endl;
	cout << ((TH2F *) LATpreSelDATA   -> LATcorrR)->GetName() << " : " <<
		((TH2F *) LATpreSelDATA   -> LATcorrR)->GetEntries() << endl;

	LATpreSelDATA_R  -> Write();
	preSelLATcorr    -> Write();
	preSelLATcorr_fit-> Write();
	inputHistoFile->Close();


	finalHistos.Add(LATpreSelDATA_R  );
        finalHistos.Add(preSelLATcorr    );
        finalHistos.Add(preSelLATcorr_fit);

	finalHistos.writeObjsInFolder("Results");

        cout<<"*** Plotting ...  ****"<<endl;

	DATApreSeleff_Plot(LATpreSelDATA_R  ,
                           preSelLATcorr    ,
                           preSelLATcorr_fit
	);


	return;

}


