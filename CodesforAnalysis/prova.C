using namespace std;

#include <iostream>

#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TH3.h"
#include "TF2.h"
#include "TVector3.h"
#include "TMath.h"
#include "TKey.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSpline.h"
#include "TFractionFitter.h"
#include "THStack.h"
#include "TNtuple.h"
#include "TObject.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TObjArray.h"

int size=3;
TH1F * p1[10][3];
TH2F * Histo = new TH2F("","",18,-1,1,100,0,1e-3);

void Scan(int i,int size,int comb[]){
	TH1F * Molt;
	for(int c=0;c<10;c++){
		comb[i]=c;

		for(int z=0;z<size;z++)
			cout<<comb[z]<<" ";
			cout<<endl;

		for(int z=0;z<size;z++)
                        cout<<p1[comb[z]][z]<<" ";
                        cout<<endl;

		Molt =(TH1F*) p1[comb[0]][0]->Clone();
		Molt->Multiply(p1[comb[1]][1]);
		Molt->Multiply(p1[comb[2]][2]);
		for(int l=0;l<18;l++)
		Histo->Fill(l,Molt->GetBinContent(l+1));

		if(i+1<size) Scan(i+1,size,comb);
		}	
	return;
}


int prova(){
	for(int i=0;i<10;i++) 
		for(int n=0;n<size;n++) {p1[i][n]=new TH1F("","",18,-1,1); p1[i][n]->FillRandom("gaus");p1[i][n]->Scale(1/p1[i][n]->GetEntries());}

	int comb[3]={0,0,0};
	Scan(0,size,comb);

	Histo->Draw("col");
	return 0;
}
