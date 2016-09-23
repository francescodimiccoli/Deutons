#include "PlottingFunctions/MCQcheck_Plot.h"



using namespace std;

class QDist_check{

	public:	
	TH3 * effQ;
	TH3 * bkgndQ;	        

	TH1*  effDist; 
	TH1*  bkgndDist;

	//creation constructor	
	QDist_check(std::string basename){
		effQ    = new TH3F((basename+"_effQ" ).c_str(),(basename+"_effQ" ).c_str(),400,0,4,400,0,4,400,0,4);
		bkgndQ  = new TH3F((basename+"_bkndQ").c_str(),(basename+"_bkndQ").c_str(),400,0,4,400,0,4,400,0,4);
		
		effDist  = new TH1F((basename+"_effDist"  ).c_str(),(basename+"_effDist"  ).c_str(),500,0,100);
		bkgndDist = new TH1F((basename+"_bkndDist" ).c_str(),(basename+"_bkndDist" ).c_str(),500,0,100);
	}
	//reading constructor
	QDist_check(TFile * file, std::string basename){	
		effQ   = (TH3 *)file->Get((basename + "_effQ"     ).c_str());
		bkgndQ = (TH3 *)file->Get((basename + "_bkndQ"   ).c_str());

		effDist   = (TH1 *)file->Get((basename + "_effDist"     ).c_str());
		bkgndDist = (TH1 *)file->Get((basename + "_bkndDist"   ).c_str());

	}

	void Write(){
		effQ     ->Write(); 	
                bkgndQ   ->Write(); 
                effDist  ->Write(); 
                bkgndDist ->Write(); 
	}


	TGraph * PlotROC_Q();
	TGraph * PlotROC_Dist();

};


TGraph * QDist_check::PlotROC_Q(){

	TGraph * ROC_Q = new TGraph();
	int a=0;
	int b=2;
	int bin_a=0;
	int bin_b= effQ->GetXaxis()->FindBin(2);
	int nbins=bin_b - bin_a;

	float eff,bk=0;

	for(int i=0;i<nbins/2;i++){
		
		eff=effQ->Integral(bin_a+i,bin_b-i,bin_a+i,bin_b-i,bin_a+i,bin_b-i);
		//for(int bin=bin_a+i;bin<bin_b-i;bin++) eff+=effQ->GetBinContent(bin+1,bin+1,bin+1);
		eff/=effQ->Integral();
		
		bk=bkgndQ->Integral(bin_a+i,bin_b-i,bin_a+i,bin_b-i,bin_a+i,bin_b-i);
        //        for(int bin=bin_a+i;bin<bin_b-i;bin++) bk+=bkgndQ->GetBinContent(bin+1,bin+1,bin+1);
		bk/=bkgndQ->Integral();		
		cout<<bin_a+i<<" "<<bin_b-i<<eff<<" "<<1-bk<<endl;
		ROC_Q->SetPoint(i,eff,1-bk);
	}
	return ROC_Q;

}

	

TGraph * QDist_check::PlotROC_Dist(){

        TGraph * ROC_Dist = new TGraph();
        
	float eff,bk=0;

        for(int i=effDist->GetNbinsX();i>0;i--){

                eff = effDist-> Integral(0,i);
		eff/= effDist-> Integral();
		bk = bkgndDist-> Integral(0,i);
		bk/= bkgndDist-> Integral();
		ROC_Dist->SetPoint(i,eff,1-bk);
        }
        return ROC_Dist;

}




QDist_check * PrejcheckTOF =   new QDist_check("PrejcheckTOF");
QDist_check * PrejcheckNaF =   new QDist_check("PrejcheckNaF");	
QDist_check * PrejcheckAgl =   new QDist_check("PrejcheckAgl");

QDist_check * HeFragcheckTOF = new QDist_check("HeFragcheckTOF");
QDist_check * HeFragcheckNaF = new QDist_check("HeFragcheckNaF");
QDist_check * HeFragcheckAgl = new QDist_check("HeFragcheckAgl");


void MCQcheck_Fill() {

	float mass = 0;
	//cuts
	if(Tup.Beta<=0||Tup.R<=0) return;
	if(!(Likcut)) return;
	if(!(Betastrongcut)) return;

	if(cmask.isOnlyFromToF()){	
		mass = ((Tup.R/Tup.Beta)*pow((1-pow(Tup.Beta,2)),0.5));	
		if(mass>1.9){
			if(Massa_gen>1&&Massa_gen<2) PrejcheckTOF->effQ  ->Fill(Tup.qUtof,Tup.qInner,Tup.qLtof);
			if(Massa_gen<1) 	     PrejcheckTOF->bkgndQ->Fill(Tup.qUtof,Tup.qInner,Tup.qLtof);
			if(Massa_gen>1&&Massa_gen<2) PrejcheckTOF->effDist  ->Fill(Tup.Dist5D);
			if(Massa_gen<1) 	     PrejcheckTOF->bkgndDist->Fill(Tup.Dist5D);
		}
			if(Massa_gen<1) 	     HeFragcheckTOF->effQ  ->Fill(Tup.qUtof,Tup.qInner,Tup.qLtof);
			if(Massa_gen>2) 	     HeFragcheckTOF->bkgndQ->Fill(Tup.qUtof,Tup.qInner,Tup.qLtof);
			if(Massa_gen<1) 	     HeFragcheckTOF->effDist  ->Fill(Tup.Dist5D);
			if(Massa_gen>2) 	     HeFragcheckTOF->bkgndDist->Fill(Tup.Dist5D);

	}

	if(cmask.isFromNaF()){
                mass = ((Tup.R/Tup.BetaRICH)*pow((1-pow(Tup.BetaRICH,2)),0.5));
                if(mass>1.9){
                        if(Massa_gen>1&&Massa_gen<2) PrejcheckNaF->effQ  ->Fill(Tup.qUtof,Tup.qInner,Tup.qLtof);
                        if(Massa_gen<1)              PrejcheckNaF->bkgndQ->Fill(Tup.qUtof,Tup.qInner,Tup.qLtof);
			if(Massa_gen>1&&Massa_gen<2) PrejcheckNaF->effDist  ->Fill(Tup.Dist5D);
			if(Massa_gen<1) 	     PrejcheckNaF->bkgndDist->Fill(Tup.Dist5D);
		
		}
			
			if(Massa_gen<1) 	     HeFragcheckNaF->effQ  ->Fill(Tup.qUtof,Tup.qInner,Tup.qLtof);
			if(Massa_gen>2) 	     HeFragcheckNaF->bkgndQ->Fill(Tup.qUtof,Tup.qInner,Tup.qLtof);
			if(Massa_gen<1) 	     HeFragcheckNaF->effDist  ->Fill(Tup.Dist5D);
			if(Massa_gen>2) 	     HeFragcheckNaF->bkgndDist->Fill(Tup.Dist5D);

			
        }

	if(cmask.isFromAgl()){
                mass = ((Tup.R/Tup.BetaRICH)*pow((1-pow(Tup.BetaRICH,2)),0.5));
                if(mass>1.9){
                        if(Massa_gen>1&&Massa_gen<2) PrejcheckAgl->effQ  ->Fill(Tup.qUtof,Tup.qInner,Tup.qLtof);
                        if(Massa_gen<1)              PrejcheckAgl->bkgndQ->Fill(Tup.qUtof,Tup.qInner,Tup.qLtof);
			if(Massa_gen>1&&Massa_gen<2) PrejcheckAgl->effDist  ->Fill(Tup.Dist5D);
			if(Massa_gen<1) 	     PrejcheckAgl->bkgndDist->Fill(Tup.Dist5D);
			
	  	}
			if(Massa_gen<1) 	     HeFragcheckAgl->effQ  ->Fill(Tup.qUtof,Tup.qInner,Tup.qLtof);
			if(Massa_gen>2) 	     HeFragcheckAgl->bkgndQ->Fill(Tup.qUtof,Tup.qInner,Tup.qLtof);
			if(Massa_gen<1) 	     HeFragcheckAgl->effDist  ->Fill(Tup.Dist5D);
			if(Massa_gen>2) 	     HeFragcheckAgl->bkgndDist->Fill(Tup.Dist5D);

        }

	return;

}


void MCQcheck_Write() {

	PrejcheckTOF ->Write(); 
	PrejcheckNaF ->Write(); 
	PrejcheckAgl ->Write();
	
	HeFragcheckTOF ->Write(); 
	HeFragcheckNaF ->Write(); 
	HeFragcheckAgl ->Write(); 
 
}


void MCQcheck(std::string filename){


	cout<<"******* MC Q DISTANCE Check ********"<<endl;

	cout<<"*** Reading  P1 file ****"<<endl;
	TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");


	QDist_check * PrejcheckTOF =   new QDist_check(inputHistoFile, "PrejcheckTOF");
	QDist_check * PrejcheckNaF =   new QDist_check(inputHistoFile, "PrejcheckNaF");	
	QDist_check * PrejcheckAgl =   new QDist_check(inputHistoFile, "PrejcheckAgl");

	QDist_check * HeFragcheckTOF = new QDist_check(inputHistoFile, "HeFragcheckTOF");
	QDist_check * HeFragcheckNaF = new QDist_check(inputHistoFile, "HeFragcheckNaF");
	QDist_check * HeFragcheckAgl = new QDist_check(inputHistoFile, "HeFragcheckAgl");

	cout<<"******* MC Q DISTANCE Check ********"<<endl;

	cout<<"TOF"<<endl;
	TGraph * Pbckgnd_ROC_QTOF = PrejcheckTOF ->  PlotROC_Q();
	TGraph * Pbckgnd_ROC_DistTOF = PrejcheckTOF ->  PlotROC_Dist();

	TGraph * Hebckgnd_ROC_QTOF = HeFragcheckTOF ->  PlotROC_Q();
	TGraph * Hebckgnd_ROC_DistTOF = HeFragcheckTOF ->  PlotROC_Dist();
	
	cout<<"NaF"<<endl;
	TGraph * Pbckgnd_ROC_QNaF = PrejcheckNaF ->  PlotROC_Q();
	TGraph * Pbckgnd_ROC_DistNaF = PrejcheckNaF ->  PlotROC_Dist();

	TGraph * Hebckgnd_ROC_QNaF = HeFragcheckNaF ->  PlotROC_Q();
	TGraph * Hebckgnd_ROC_DistNaF = HeFragcheckNaF ->  PlotROC_Dist();
	
	cout<<"Agl"<<endl;
	TGraph * Pbckgnd_ROC_QAgl = PrejcheckAgl ->  PlotROC_Q();
	TGraph * Pbckgnd_ROC_DistAgl = PrejcheckAgl ->  PlotROC_Dist();

	TGraph * Hebckgnd_ROC_QAgl = HeFragcheckAgl ->  PlotROC_Q();
	TGraph * Hebckgnd_ROC_DistAgl = HeFragcheckAgl ->  PlotROC_Dist();
	
	cout<<"**** Plotting ********"<<endl;
	
	MCQcheck_Plot(
		 	Pbckgnd_ROC_QTOF,
			Pbckgnd_ROC_QNaF,   
			Pbckgnd_ROC_QAgl,
			Pbckgnd_ROC_DistTOF,
                        Pbckgnd_ROC_DistNaF,
                        Pbckgnd_ROC_DistAgl,
			Hebckgnd_ROC_QTOF,
			Hebckgnd_ROC_QNaF,   
			Hebckgnd_ROC_QAgl,
			Hebckgnd_ROC_DistTOF,
                        Hebckgnd_ROC_DistNaF,
                        Hebckgnd_ROC_DistAgl 
		      
		     );	


}




