#include <iostream>

#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TH3.h"
#include "TF2.h"
#include <vector>
#include <string>
#include <sstream>

#include "TVector3.h"
#include "TMath.h"
#include "TKey.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TNtuple.h"
#include "TObject.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TObjArray.h"

#include "../include/filesaver.h"

using namespace std;

TH2F * PopulateHisto(TNtuple * Q,std::string variable1,std::string variable2,std::string binning1,std::string binning2,std::string cut,std::string name){
	cout<<(variable2+":"+variable1+">>"+name+"("+binning1+","+binning2+")").c_str()<<", "<<cut.c_str()<<endl;	
	Q->Draw((variable2+":"+variable1+">>"+name+"("+binning1+","+binning2+")").c_str(),cut.c_str());
	TH2F *Histo= (TH2F*)gDirectory->Get(name.c_str());
	Histo -> GetXaxis() ->SetTitle(variable1.c_str());
	Histo -> GetYaxis() ->SetTitle(variable2.c_str());
	return Histo;
};

class CorrelationCalculator{
	private:
	TNtuple * Q;
	std::vector<std::string> Variables;
	std::vector<std::string> Binnings;
	std::string MinimumCut;
	std::string Basename;
	std::vector<std::vector<TH2F*>> Matrix;
	TH2F * CorrMatrix;
	public:
	CorrelationCalculator(TNtuple * q,std::string basename,std::vector<std::string> variables,std::vector<std::string> binnings,std::string Cut){
		Basename=basename;
		Variables=variables;
		Binnings=binnings;
		MinimumCut = Cut;
		Q = q;
		CorrMatrix = new TH2F((basename+"_Correl.").c_str(),(basename+"_Correl.").c_str(),variables.size(),0,variables.size(),variables.size(),0,variables.size());
		for(int i=0;i<Variables.size();i++){
			CorrMatrix->GetXaxis()->SetBinLabel(i+1,Variables[i].c_str());
			CorrMatrix->GetYaxis()->SetBinLabel(i+1,Variables[i].c_str());
		}

		for(int i=0;i<variables.size();i++){
			Matrix.push_back(std::vector<TH2F *>());
			for(int j=0;j<variables.size();j++){
				 Matrix[i].push_back(new TH2F());
				}	
		}
	}
	void PopulateHistos();
	void Eval_CorrelMatrix();
	void SaveHistos(FileSaver finalhistos);	
};


void CorrelationCalculator::PopulateHistos(){
	for(int i=0;i<Variables.size();i++)
		for(int j=0;j<Variables.size();j++)
			if(i>=j){
				Matrix[i][j]=PopulateHisto(Q,Variables[i],Variables[j],Binnings[i],Binnings[j],MinimumCut.c_str(),to_string(i)+"vs"+to_string(j).c_str());
				Matrix[i][j]->SetName((Variables[i]+"vs"+Variables[j]).c_str());	
				Matrix[i][j]->SetTitle((Variables[i]+"vs"+Variables[j]).c_str());
			}

	return;
}

void CorrelationCalculator::Eval_CorrelMatrix(){
	for(int i=0;i<Variables.size();i++)
                for(int j=0;j<Variables.size();j++)
                        if(i>=j)
			       CorrMatrix->SetBinContent(i+1,j+1,fabs(Matrix[i][j]->GetCorrelationFactor()));
			else   CorrMatrix->SetBinContent(i+1,j+1,-10);
	return;			
}

void CorrelationCalculator::SaveHistos(FileSaver finalhistos){
	for(int i=0;i<Variables.size();i++)
                for(int j=0;j<Variables.size();j++)
			if(i>=j)
				finalhistos.Add(Matrix[i][j]);
	finalhistos.Add(CorrMatrix);
	finalhistos.writeObjsInFolder(Basename.c_str());
	return;
}	

int Correlations(){
	cout<<"************************ READING DATA ***************************"<<endl;

	string inputfile = "../Ntuple-making/Ntuples/MC/NtupleMC24.root";
	TFile * input = TFile::Open(inputfile.c_str());
	TNtuple * Q = (TNtuple *)input->Get("Q");
	cout<<input<<endl;

	FileSaver finalHistos;
        finalHistos.setName("SelectionCorrelations.root");

	cout<<"************************ CUTS & VARIABLES DEFINITIONS ***************************"<<endl;
	std::string IsPreselected = "(Cutmask&187)==187&&R>0";
	std::string L1Hit         = "qL1>0";
	std::string IsOnlyFromTOF = "(((Cutmask>>11)!=0&&(Cutmask>>11)!=512)||BetaRICH_new<0)";
	std::string IsFromNaF     = "(Cutmask>>11)==512&&BetaRICH_new>0";
	std::string IsFromAgl     = "(Cutmask>>11)==0&&BetaRICH_new>0";	
	std::string IsProtonMC    = "Massa_gen<1";
	std::string IsDeutonMC    = "Massa_gen>1&&Massa_gen<2";



	std::vector<std::string> QualityVariables ={"Dist5D_P","LDiscriminant","qL1","qUtof","qLtof","qInner"};
	std::vector<std::string> QualityBinnings  ={"200,0,18","100,0,1","100,0,3","100,0,3","100,0,3","100,0,3"};
	
	std::vector<std::string> RICHVariables ={"qL1","qUtof","qLtof","qInner",IsFromNaF,IsFromAgl};
	std::vector<std::string> RICHBinnings  ={"100,0,3","100,0,3","100,0,3","100,0,3","2,0,2","2,0,2"};

	std::vector<std::string> PresVariables ={"qL1","qUtof","qLtof","qInner",IsPreselected,L1Hit,IsFromAgl};
	std::vector<std::string> PresBinnings  ={"100,0,3","100,0,3","100,0,3","100,0,3","2,0,2","2,0,2","2,0,2"};


	

	cout<<"************************ CALCULATION of CORRELATION MATRIXES ***************************"<<endl;
	
	CorrelationCalculator * QualityTOF = new CorrelationCalculator(Q,"QualityTOF",QualityVariables,QualityBinnings,(IsPreselected+"&&"+L1Hit+"&&"+IsOnlyFromTOF+"&&"+IsProtonMC).c_str());
	
	QualityTOF->PopulateHistos();
	QualityTOF->Eval_CorrelMatrix();
	QualityTOF->SaveHistos(finalHistos);

	CorrelationCalculator * QualityAgl = new CorrelationCalculator(Q,"QualityAgl",QualityVariables,QualityBinnings,(IsPreselected+"&&"+L1Hit+"&&"+IsFromAgl+"&&"+IsProtonMC).c_str());
	
	QualityAgl->PopulateHistos();
	QualityAgl->Eval_CorrelMatrix();
	QualityAgl->SaveHistos(finalHistos);
	
	CorrelationCalculator * RICH = new CorrelationCalculator(Q,"RICH",RICHVariables,RICHBinnings,(IsPreselected+"&&"+L1Hit+"&&"+IsProtonMC).c_str());
	
	RICH->PopulateHistos();
	RICH->Eval_CorrelMatrix();
	RICH->SaveHistos(finalHistos);

	CorrelationCalculator * PreSel = new CorrelationCalculator(Q,"PreSel",PresVariables,PresBinnings,IsProtonMC.c_str());
	
	PreSel->PopulateHistos();
	PreSel->Eval_CorrelMatrix();
	PreSel->SaveHistos(finalHistos);



	return 0;
}
