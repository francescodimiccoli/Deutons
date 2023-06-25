#include "RangeMerger.h"


void RangeMerger::UseRTOIEdges(){
	Global_D.UseRTOIEdges();
	Global_P.UseRTOIEdges();
	ToFP.UseRTOIEdges();
	NaFP.UseRTOIEdges();
	AglP.UseRTOIEdges();
	ToFD.UseRTOIEdges();
	NaFD.UseRTOIEdges();
	AglD.UseRTOIEdges();
}

void RangeMerger::UseBetaTOIEdges(){
	Global_D.UseBetaTOIEdges();
	Global_P.UseBetaTOIEdges();
	ToFP.UseBetaTOIEdges();
	NaFP.UseBetaTOIEdges();
	AglP.UseBetaTOIEdges();
	ToFD.UseBetaTOIEdges();
	NaFD.UseBetaTOIEdges();
	AglD.UseBetaTOIEdges();
}

void RangeMerger::UseREdges(){
	Global_D.UseREdges();
	Global_P.UseREdges();
	ToFP.UseREdges();
	NaFP.UseREdges();
	AglP.UseREdges();
	ToFD.UseREdges();
	NaFD.UseREdges();
	AglD.UseREdges();
}

void RangeMerger::UseBetaEdges(){
	Global_D.UseBetaEdges();
	Global_P.UseBetaEdges();
	ToFP.UseBetaEdges();
	NaFP.UseBetaEdges();
	AglP.UseBetaEdges();
	ToFD.UseBetaEdges();
	NaFD.UseBetaEdges();
	AglD.UseBetaEdges();
}

void RangeMerger::Reset(){
		Global_D.Reset();
		Global_P.Reset();
		ToFD.Reset();
		ToFP.Reset();
		NaFD.Reset();
		NaFP.Reset();
		AglD.Reset();
		AglP.Reset();

	}

void RangeMerger::setBinsFromRDatacard(std::string datacard, TF1 * ResponseTOF, TF1 * ResponseNaF, TF1 * ResponseAgl,float a, float b){
		float betamin=0.555; float betamax=0.87;
		ToFD.setBinsFromRDatacard (datacard.c_str(), betamin, betamax ,ResponseTOF,a,b);
		ToFP.setBinsFromRDatacard (datacard.c_str(), betamin, betamax ,ResponseTOF,a,b);
		betamin=0.86, betamax=0.977;
		NaFD.setBinsFromRDatacard (datacard.c_str(), betamin, betamax ,ResponseNaF,a,b);
		NaFP.setBinsFromRDatacard (datacard.c_str(), betamin, betamax ,ResponseNaF,a,b);
		betamin=0.97, betamax=0.995;
		AglD.setBinsFromRDatacard (datacard.c_str(), betamin, betamax ,ResponseAgl,a,b);
		AglP.setBinsFromRDatacard (datacard.c_str(), betamin, betamax ,ResponseAgl,a,b);
		
		Global_D.setBinsFromRDatacard (datacard.c_str(), 0.55, 0.995 ,ResponseTOF,a,b);
		Global_P.setBinsFromRDatacard (datacard.c_str(), 0.55, 0.995 ,ResponseTOF,a,b);
		
	}



int RangeMerger::GetGlobalBinRTOICenter(int n)    {return n;}
int RangeMerger::GetGlobalBinEkinTOICenter(int n) {return n;}
int RangeMerger::GetGlobalBinBetaTOICenter(int n) {return n;}

int RangeMerger::GetGlobalBinRCenter(int n)	  {return n;}
int RangeMerger::GetGlobalBinEkinCenter(int n) {return n;}
int RangeMerger::GetGlobalBinBetaCenter(int n) {return n;}

int RangeMerger::GetToFBinP(int globalbin){ return ToFP.GetRTOIBin( Global_P.rigbincent_TOI[globalbin]);}
int RangeMerger::GetNaFBinP(int globalbin){ return NaFP.GetRTOIBin( Global_P.rigbincent_TOI[globalbin]);}
int RangeMerger::GetAglBinP(int globalbin){ return AglP.GetRTOIBin( Global_P.rigbincent_TOI[globalbin]);}
int RangeMerger::GetToFBinD(int globalbin){ return ToFD.GetRTOIBin( Global_D.rigbincent_TOI[globalbin]);}
int RangeMerger::GetNaFBinD(int globalbin){ return NaFD.GetRTOIBin( Global_D.rigbincent_TOI[globalbin]);}
int RangeMerger::GetAglBinD(int globalbin){ return AglD.GetRTOIBin( Global_D.rigbincent_TOI[globalbin]);}





TH1F * RangeMerger::MergeSubDResult_P(TH1F * ResultTOF, TH1F * ResultNaF, TH1F * ResultAgl,bool nafpriority) {
		std::cout<<ResultTOF->GetName()<<" "<<ResultNaF<<" "<<ResultAgl<<std::endl;
			
		std::cout<<ResultTOF<<" "<<ResultNaF<<" "<<ResultAgl<<std::endl;
		TH1F * Merged = new TH1F("MergedRange_RTOI_P","MergedRange_RTOI_P",Global_P.size(),0,Global_P.size());
		std::cout<<"GLOBAL P SIZE: "<<Global_P.size()<<std::endl;	
		for(int j =0; j<Global_P.size(); j++){
			int i=0;
			if(GetToFBinP(j)>=0) {
				i=GetToFBinP(j);
				if(i>=0&&ResultTOF->GetBinContent(i+1)>0){
					Merged->SetBinContent(j+1,ResultTOF->GetBinContent(i+1));
					Merged->SetBinError(j+1,ResultTOF->GetBinError(i+1));
				}
				continue;
			}
			//NaF gives the priority to other subd.
			if(GetNaFBinP(j)>=0) {
				i=GetNaFBinP(j);
				if(i>=0&&ResultNaF->GetBinContent(i+1)>0){
					Merged->SetBinContent(j+1,ResultNaF->GetBinContent(i+1));
					Merged->SetBinError(j+1,ResultNaF->GetBinError(i+1));
				}
			/*	if(GetToFBinP(j)>=0) {
					i=GetToFBinP(j);
					if(i>=0&&ResultTOF->GetBinContent(i+1)>0){
						Merged->SetBinContent(j+1,ResultTOF->GetBinContent(i+1));
						Merged->SetBinError(j+1,ResultTOF->GetBinError(i+1));
					}
				}
				if(GetAglBinP(j)>=0&&!nafpriority) {
					i=GetAglBinP(j);
					if(i>=0&&ResultAgl->GetBinContent(i+1)>0){
						Merged->SetBinContent(j+1,ResultAgl->GetBinContent(i+1));
						Merged->SetBinError(j+1,ResultAgl->GetBinError(i+1));
					}
				}*/
				continue;
			}
			if(GetAglBinP(j)>=0) {
				i=GetAglBinP(j);
				if(i>=0&&ResultAgl->GetBinContent(i+1)>0){
					Merged->SetBinContent(j+1,ResultAgl->GetBinContent(i+1));
					Merged->SetBinError(j+1,ResultAgl->GetBinError(i+1));
				}
			}
		}
		std::cout<<"end merging"<<std::endl;
	
		return Merged;	

}


TH1F * RangeMerger::MergeSubDResult_D(TH1F * ResultTOF, TH1F * ResultNaF, TH1F * ResultAgl,bool nafpriority) {
	
		std::cout<<ResultTOF->GetName()<<" "<<ResultNaF<<" "<<ResultAgl<<std::endl;
		
		TH1F * Merged = new TH1F("MergedRange_RTOI_D","MergedRange_RTOI_D",Global_D.size(),0,Global_D.size());
		std::cout<<"GLOBAL D SIZE: "<<Global_D.size()<<std::endl;	
		for(int j =0; j<Global_D.size(); j++){
			int i=0;
			if(GetToFBinD(j)>=0) {
				i=GetToFBinD(j);
				if(i>=0&&ResultTOF->GetBinContent(i+1)>0){
					Merged->SetBinContent(j+1,ResultTOF->GetBinContent(i+1));
					Merged->SetBinError(j+1,ResultTOF->GetBinError(i+1));
				}
				continue;
			}
			//NaF gives the priority to other subd.
			if(GetNaFBinD(j)>=0) {
				i=GetNaFBinD(j);
				if(i>=0&&ResultNaF->GetBinContent(i+1)>0){
					Merged->SetBinContent(j+1,ResultNaF->GetBinContent(i+1));
					Merged->SetBinError(j+1,ResultNaF->GetBinError(i+1));
				}
				/*if(GetToFBinD(j)>=0) {
					i=GetToFBinD(j);
					if(i>=0&&ResultTOF->GetBinContent(i+1)>0){
						Merged->SetBinContent(j+1,ResultTOF->GetBinContent(i+1));
						Merged->SetBinError(j+1,ResultTOF->GetBinError(i+1));
					}
				}
			*/	if(GetAglBinD(j)>=0&&!nafpriority) {
					i=GetAglBinD(j);
					if(i>=0&&ResultAgl->GetBinContent(i+1)>0){
						Merged->SetBinContent(j+1,ResultAgl->GetBinContent(i+1));
						Merged->SetBinError(j+1,ResultAgl->GetBinError(i+1));
					}
				}

				continue;
			}
			if(GetAglBinD(j)>=0) {
				i=GetAglBinD(j);
				if(i>=0&&ResultAgl->GetBinContent(i+1)>0){
					Merged->SetBinContent(j+1,ResultAgl->GetBinContent(i+1));
					Merged->SetBinError(j+1,ResultAgl->GetBinError(i+1));
				}
			}
		}
		std::cout<<"end merging"<<std::endl;
		return Merged;	

}


TH1F * RangeMerger::MergedRatio(TH1F * Result_D, TH1F * Result_P){
	int nbins = Global_D.size();
	
	TH1F * Ratio = new TH1F("Ratio_R","Ratio_R",nbins,0,nbins);
	std::cout<<"*************** MERGING RATIO in Ekin TOI ********************"<<std::endl;
	TF1 * protonmodel = new TF1("protonModel","[0]*x^[1]");
	protonmodel->SetParameter(0,14788.883);
	protonmodel->SetParameter(1,-2.77);

	for(int i=0;i<nbins;i++){
		float Rtoicenter = Global_D.RigTOIBinsCent()[i];
		float bin_P = Global_P.GetRTOIBin(Rtoicenter);
		float bin_D = Global_D.GetRTOIBin(Rtoicenter);	
		std::cout<<"merging: "<<bin_D<<" "<<bin_P<<std::endl;
		if(bin_P>=0 && bin_D>=0){
			Ratio->SetBinContent(i+1,Result_D->GetBinContent(bin_D+1)/Result_P->GetBinContent(bin_P+1));
			Ratio->SetBinError(i+1,pow(pow(Result_D->GetBinError(bin_D+1)/Result_D->GetBinContent(bin_D+1),2)+
				           pow(Result_P->GetBinError(bin_P+1)/Result_P->GetBinContent(bin_P+1),2),0.5)*Ratio->GetBinContent(i+1));
		}
	
		if(bin_P<0 && bin_D>=0)	{
			Ratio->SetBinContent(i+1,Result_D->GetBinContent(bin_D+1)/protonmodel->Eval(Global_D.RigTOIBinsCent()[bin_D]));
			Ratio->SetBinError(i+1,2*pow(pow(Result_D->GetBinError(bin_D+1)/Result_D->GetBinContent(bin_D+1),2)+
				           pow(Result_P->GetBinError(Global_P.size())/Result_P->GetBinContent(Global_P.size()),2),0.5)*Ratio->GetBinContent(i+1));
			
		}	
	}
        return Ratio;	

}

TH1F * RangeMerger::MergedRatio_Ekin(TH1F * Result_D, TH1F * Result_P){
	int nbins = Global_D.size();
	
	TH1F * Ratio = new TH1F("Ratio_Ekin","Ratio_Ekin",nbins,0,nbins);
	std::cout<<"*************** MERGING RATIO in Ekin TOI ********************"<<std::endl;
	for(int i=0;i<nbins;i++){
		float Betatoicenter = Global_D.BetaTOIBinsCent()[i];
		float bin_P = Global_P.GetBetaTOIBin(Betatoicenter);
		float bin_D = Global_D.GetBetaTOIBin(Betatoicenter);	
		std::cout<<"merging: "<<bin_D<<" "<<bin_P<<std::endl;
		if(bin_P>=0 && bin_D>=0){
			if(Result_P->GetBinContent(bin_P+1)>0){
			Ratio->SetBinContent(i+1,Result_D->GetBinContent(bin_D+1)/Result_P->GetBinContent(bin_P+1));
			Ratio->SetBinError(i+1,pow(pow(Result_D->GetBinError(bin_D+1)/Result_D->GetBinContent(bin_D+1),2)+
				           pow(Result_P->GetBinError(bin_P+1)/Result_P->GetBinContent(bin_P+1),2),0.5)*Ratio->GetBinContent(i+1));
			}
		}
	}
        return Ratio;	

}
