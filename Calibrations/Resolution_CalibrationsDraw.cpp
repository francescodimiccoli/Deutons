#include <bitset>
#include "TROOT.h"
#include "TNtuple.h"
#include <TSpline.h>
#include "../include/binning.h"
#include "TFile.h"
#include "TH2.h"
#include "TF2.h"
#include <TVector3.h>
#include "TMath.h"
#include <TFile.h>
#include "TFile.h"
#include "TH2.h"
#include "TF2.h"
#include <TVector3.h>
#include "TMath.h"
#include "TGraphErrors.h"

#include "../include/GlobalBinning.h"

#include "../Ntuple-making/Commonglobals.cpp"
#include "../Ntuple-making/Variables.hpp"
#include "../include/Cuts.h"

#include "../include/filesaver.h"

#include "../include/FitError.h"
#include "../include/Resolution.h"



bool Calculate_Resolution( Resolution * Reso, bool checkfile, TTree * treeMC,  FileSaver finalHistos , std::vector<float> ExpValues={-1}, bool spline = false){
	
	if(!checkfile){
                Reso->Fill(treeMC);
		Reso->Normalize();
        }

        else Reso = new Resolution(finalHistos,Reso->GetName(),Reso->GetBinning());

        if(Reso->CheckHistos()){
                Reso->Eval_Resolution(ExpValues);
		Reso->Save(finalHistos);
        }
   	finalHistos.Add(Reso->Get_Means()	);
   	finalHistos.Add(Reso->Get_Sigmas()	);
	finalHistos.Add(Reso->Get_Resolutions()	);

	if(spline)   finalHistos.Add(Reso->ModelSigmasWithSpline());
	else 	     finalHistos.Add(Reso->ModelSigmasWithPoly());

	if(spline)   finalHistos.Add(Reso->ModelMeansWithSpline());
	else 	     finalHistos.Add(Reso->ModelMeansWithPoly());

	
   	finalHistos.writeObjsInFolder((Reso->GetName()+"/Fit Results").c_str());				
	
	return Reso->CheckHistos();
}


int main(int argc, char * argv[])
{


	cout<<"****************************** FILES OPENING ***************************************"<<endl;

	string INPUT1(argv[1]);
	string INPUT2(argv[2]);
        string OUTPUT(argv[3]);

	FileSaver finalHistos;
        finalHistos.setName(OUTPUT.c_str());
	bool checkfile = finalHistos.CheckFile();	

	TFile *fileDT =TFile::Open(INPUT1.c_str());
        TFile *fileMC =TFile::Open(INPUT2.c_str());
        
	TTree *treeMC = (TTree *)fileMC->Get("parametri_geo");
	TTree *treeDT = (TTree *)fileDT->Get("parametri_geo");


	cout<<"****************************** BINS ***************************************"<<endl;

	SetBins();	

	PResB.Print();

	cout<<"**TOF**"<<endl;
	ToFResB.Print();

	cout<<"**NaF**"<<endl;
	NaFResB.Print();

	cout<<"**Agl**"<<endl;
	AglResB.Print();

	ToFResB.UseBetaEdges();
	NaFResB.UseBetaEdges();
	AglResB.UseBetaEdges();

	PResB.UseREdges();


	cout<<endl;

	cout<<"****************************** VARIABLES ***************************************"<<endl;
	
/*	Variables * varsMC = new Variables;
        varsMC->ReadBranches(treeMC);
*/
	std::string IsProtonMC    = "GenMass<1&&GenMass>0"  ;
	std::string IsDeutonMC    = "GenMass<2&&GenMass>1"  ;
	std::string IsData 	  = "GenMass==0"  	    ;
	std::string IsPreselected = "(CUTMASK&187)==187"    ;
	std::string IsOnlyFromToF = "(RICHmask&1023)!=512&&(RICHmask&1023)>0";
	std::string IsFromNaF 	  = "(RICHmask&1023)==512"  ;
	std::string IsFromAgl     = "(RICHmask&1023)==0"    ;
	std::string Beta_gen	  = "((GenMomentum/GenMass)^2/((GenMomentum/GenMass)^2+1))^0.5"    ;
	std::string InverseBeta_gen = "(((GenMomentum/GenMass)^2 + 1)/(GenMomentum/GenMass)^2)^0.5";
	std::string EdepTrack     = "(trtrack_edep[1]+trtrack_edep[2]+trtrack_edep[3]+trtrack_edep[4]+trtrack_edep[5]+trtrack_edep[6]+trtrack_edep[7])";

	cout<<"****************************** ANALYSIS ***************************************"<<endl;

	//rigidity resolution vs rigidity	
	Resolution * RigidityResolution_P = new Resolution("RvsR Resolution (P)",PResB,(IsPreselected + "&&" + IsProtonMC).c_str(),1000,-0.5,1.5,"1/R - 1/GenMomentum","GenMomentum");
	Calculate_Resolution( RigidityResolution_P, checkfile, treeMC, finalHistos,PResB.RigBinsCent(),true);

	/*Resolution * RigidityResolution_D = new Resolution("RvsR Resolution (D)",PResB,(IsPreselected + "&&" + IsDeutonMC).c_str(),"1/R - 1/GenMomentum","GenMomentum",1000,-0.5,1.5);
        Calculate_Resolution( RigidityResolution_D, checkfile, treeMC, finalHistos,PResB.RigBinsCent(),true);


	//rigidity resolution vs beta
	Resolution * RigidityTOFResolution_P = new Resolution("RvsBetaTOF Resolution (P)",ToFResB,(IsPreselected+ "&&" + IsProtonMC).c_str(),"1/R - 1/GenMomentum",Beta_gen.c_str(),1000,-0.5,1.5);
	Calculate_Resolution( RigidityTOFResolution_P, checkfile, treeMC, finalHistos,ToFResB.RigBinsCent());

	Resolution * RigidityNaFResolution_P = new Resolution("RvsBetaNaF Resolution (P)",NaFResB,(IsPreselected + "&&" + IsProtonMC + "&&" + IsFromNaF ).c_str(),"1/R - 1/GenMomentum",Beta_gen.c_str(),1000,-0.5,1.5);
	Calculate_Resolution( RigidityNaFResolution_P, checkfile, treeMC, finalHistos,NaFResB.RigBinsCent());

	Resolution * RigidityAglResolution_P = new Resolution("RvsBetaAgl Resolution (P)",AglResB,(IsPreselected + "&&" + IsProtonMC + "&&" + IsFromAgl ).c_str(),"1/R - 1/GenMomentum",Beta_gen.c_str(),1000,-0.5,1.5);
	Calculate_Resolution( RigidityAglResolution_P, checkfile, treeMC, finalHistos,AglResB.RigBinsCent());

	Resolution * RigidityTOFResolution_D = new Resolution("RvsBetaTOF Resolution (D)",ToFResB,(IsPreselected+ "&&" + IsDeutonMC).c_str(),"1/R - 1/GenMomentum",Beta_gen.c_str(),1000,-0.5,1.5);
	Calculate_Resolution( RigidityTOFResolution_D, checkfile, treeMC, finalHistos,ToFResB.RigBinsCent());

	Resolution * RigidityNaFResolution_D = new Resolution("RvsBetaNaF Resolution (D)",NaFResB,(IsPreselected + "&&" + IsDeutonMC + "&&" + IsFromNaF ).c_str(),"1/R - 1/GenMomentum",Beta_gen.c_str(),1000,-0.5,1.5);
	Calculate_Resolution( RigidityNaFResolution_D, checkfile, treeMC, finalHistos,NaFResB.RigBinsCent());

	Resolution * RigidityAglResolution_D = new Resolution("RvsBetaAgl Resolution (D)",AglResB,(IsPreselected + "&&" + IsDeutonMC + "&&" + IsFromAgl ).c_str(),"1/R - 1/GenMomentum",Beta_gen.c_str(),1000,-0.5,1.5);
	Calculate_Resolution( RigidityAglResolution_D, checkfile, treeMC, finalHistos,AglResB.RigBinsCent());


	//Beta resolution vs rigidity	
	Resolution * BetaTOF_RResolution_P = new Resolution("BetaTOFvsR Resolution (P)",PResB,(IsPreselected+ "&&" + IsProtonMC + "&&" + IsOnlyFromToF).c_str(),("1/BetaHR -" + InverseBeta_gen).c_str(),"GenMomentum",1000,-0.5,1.5);
	Calculate_Resolution( BetaTOF_RResolution_P, checkfile, treeMC, finalHistos,PResB.BetaBinsCent(),true);

	Resolution * BetaNaF_RResolution_P = new Resolution("BetaNaFvsR Resolution (P)",PResB,(IsPreselected + "&&" + IsProtonMC + "&&" + IsFromNaF   ).c_str(),("1/BetaRICH -" + InverseBeta_gen).c_str(),"GenMomentum",250,-0.05,0.15);
	Calculate_Resolution( BetaNaF_RResolution_P, checkfile, treeMC, finalHistos,PResB.BetaBinsCent(),true);

	Resolution * BetaAgl_RResolution_P = new Resolution("BetaAglvsR Resolution (P)",PResB,(IsPreselected + "&&" + IsProtonMC + "&&" + IsFromAgl   ).c_str(),("1/BetaRICH -" + InverseBeta_gen).c_str(),"GenMomentum",250,-0.075,0.075);
	Calculate_Resolution( BetaAgl_RResolution_P, checkfile, treeMC, finalHistos,PResB.BetaBinsCent(),true);

	Resolution * BetaTOF_RResolution_D = new Resolution("BetaTOFvsR Resolution (D)",PResB,(IsPreselected+ "&&" + IsDeutonMC + "&&" + IsOnlyFromToF).c_str(),("1/BetaHR -" + InverseBeta_gen).c_str(),"GenMomentum",1000,-0.5,1.5);
	Calculate_Resolution( BetaTOF_RResolution_D, checkfile, treeMC, finalHistos,PResB.BetaBinsCent(),true);

	Resolution * BetaNaF_RResolution_D = new Resolution("BetaNaFvsR Resolution (D)",PResB,(IsPreselected + "&&" + IsDeutonMC + "&&" + IsFromNaF   ).c_str(),("1/BetaRICH -" + InverseBeta_gen).c_str(),"GenMomentum",250,-0.05,0.15);
	Calculate_Resolution( BetaNaF_RResolution_D, checkfile, treeMC, finalHistos,PResB.BetaBinsCent(),true);

	Resolution * BetaAgl_RResolution_D = new Resolution("BetaAglvsR Resolution (D)",PResB,(IsPreselected + "&&" + IsDeutonMC + "&&" + IsFromAgl   ).c_str(),("1/BetaRICH -" + InverseBeta_gen).c_str(),"GenMomentum",250,-0.075,0.075);
	Calculate_Resolution( BetaAgl_RResolution_D, checkfile, treeMC, finalHistos,PResB.BetaBinsCent(),true);
		


	//beta resolution vs beta (ToF, NaF, Agl)
	Resolution * BetaTOFResolution_P = new Resolution("BetaTOFvsBeta Resolution (P)",ToFResB,(IsPreselected+ "&&" + IsProtonMC).c_str(),("1/BetaHR -" + InverseBeta_gen).c_str(),Beta_gen.c_str(),1000,-0.5,1.5);
	Calculate_Resolution( BetaTOFResolution_P, checkfile, treeMC, finalHistos);

	Resolution * BetaNaFResolution_P = new Resolution("BetaNaFvsBeta Resolution (P)",NaFResB,(IsPreselected + "&&" + IsProtonMC + "&&" + IsFromNaF ).c_str(),("1/BetaRICH -" + InverseBeta_gen).c_str(),Beta_gen.c_str(),250,-0.05,0.15);
	Calculate_Resolution( BetaNaFResolution_P, checkfile, treeMC, finalHistos);

	Resolution * BetaAglResolution_P = new Resolution("BetaAglvsBeta Resolution (P)",AglResB,(IsPreselected + "&&" + IsProtonMC + "&&" + IsFromAgl ).c_str(),("1/BetaRICH -" + InverseBeta_gen).c_str(),Beta_gen.c_str(),250,-0.075,0.075);
	Calculate_Resolution( BetaAglResolution_P, checkfile, treeMC, finalHistos);

	Resolution * BetaTOFResolution_D = new Resolution("BetaTOFvsBeta Resolution (D)",ToFResB,(IsPreselected+ "&&" + IsDeutonMC).c_str(),("1/BetaHR -" + InverseBeta_gen).c_str(),Beta_gen.c_str(),1000,-0.5,1.5);
	Calculate_Resolution( BetaTOFResolution_D, checkfile, treeMC, finalHistos);

	Resolution * BetaNaFResolution_D = new Resolution("BetaNaFvsBeta Resolution (D)",NaFResB,(IsPreselected + "&&" + IsDeutonMC + "&&" + IsFromNaF ).c_str(),("1/BetaRICH -" + InverseBeta_gen).c_str(),Beta_gen.c_str(),250,-0.05,0.15);
	Calculate_Resolution( BetaNaFResolution_D, checkfile, treeMC, finalHistos);

	Resolution * BetaAglResolution_D = new Resolution("BetaAglvsBeta Resolution (D)",AglResB,(IsPreselected + "&&" + IsDeutonMC + "&&" + IsFromAgl ).c_str(),("1/BetaRICH -" + InverseBeta_gen).c_str(),Beta_gen.c_str(),250,-0.075,0.075);
	Calculate_Resolution( BetaAglResolution_D, checkfile, treeMC, finalHistos);

	

	//E. dep. U. TOF vs beta (ToF)
	
	Resolution * EdepUTOFResolution_P = new Resolution("EdepUTOFvsBeta Resolution (P)",ToFResB,(IsPreselected+ "&&" + IsProtonMC).c_str(),"2/(TOFEndep[0]+TOFEndep[1])",Beta_gen.c_str(),1000,-1.5,1.5);
	Calculate_Resolution( EdepUTOFResolution_P, checkfile, treeMC, finalHistos);

	Resolution * EdepUTOFResolution_D = new Resolution("EdepUTOFvsBeta Resolution (D)",ToFResB,(IsPreselected+ "&&" + IsDeutonMC).c_str(),"2/(TOFEndep[0]+TOFEndep[1])",Beta_gen.c_str(),1000,-1.5,1.5);
	Calculate_Resolution( EdepUTOFResolution_D, checkfile, treeMC, finalHistos);
	
	
	//E. dep. L. TOF vs beta (ToF)
	
	Resolution * EdepLTOFResolution_P = new Resolution("EdepLTOFvsBeta Resolution (P)",ToFResB,(IsPreselected+ "&&" + IsProtonMC).c_str(),"2/(TOFEndep[2]+TOFEndep[3])",Beta_gen.c_str(),1000,-1.5,1.5);
	Calculate_Resolution( EdepLTOFResolution_P, checkfile, treeMC, finalHistos);

	Resolution * EdepLTOFResolution_D = new Resolution("EdepLTOFvsBeta Resolution (D)",ToFResB,(IsPreselected+ "&&" + IsDeutonMC).c_str(),"2/(TOFEndep[2]+TOFEndep[3])",Beta_gen.c_str(),1000,-1.5,1.5);
	Calculate_Resolution( EdepLTOFResolution_D, checkfile, treeMC, finalHistos);

	//E. dep. Inner Tracker vs beta (ToF)
	
	Resolution * EdepTrack_P = new Resolution("EdepTrackvsBeta Resolution (P)",ToFResB,(IsPreselected+ "&&" + IsProtonMC).c_str(),("7/"+EdepTrack).c_str(),Beta_gen.c_str(),500,-1.5,15);
	Calculate_Resolution( EdepTrack_P, checkfile, treeMC, finalHistos);

	Resolution * EdepTrack_D = new Resolution("EdepTrackvsBeta Resolution (D)",ToFResB,(IsPreselected+ "&&" + IsDeutonMC).c_str(),("7/"+EdepTrack).c_str(),Beta_gen.c_str(),500,-1.5,15);
	Calculate_Resolution( EdepTrack_D, checkfile, treeMC, finalHistos);


	// E.dep. calibration

	//E. dep. U. TOF vs beta (ToF)
	
	Resolution * EdepUTOFMC_P = new Resolution("EdepUTOFvsBeta Measured MC",ToFResB,(IsPreselected+ "&&" + IsProtonMC).c_str(),"2/(TOFEndep[0]+TOFEndep[1])","BetaHR",1000,-1.5,1.5);
	Calculate_Resolution( EdepUTOFMC_P, checkfile, treeMC, finalHistos);

	Resolution * EdepUTOFDT_P = new Resolution("EdepUTOFvsBeta Measured DT",ToFResB,(IsPreselected).c_str(),"2/(TOFEndep[0]+TOFEndep[1])","BetaHR",1000,-1.5,1.5);
	Calculate_Resolution( EdepUTOFDT_P, checkfile, treeDT, finalHistos);

	//E. dep. L. TOF vs beta (ToF)
	
	Resolution * EdepLTOFMC_P = new Resolution("EdepLTOFvsBeta Measured MC",ToFResB,(IsPreselected+ "&&" + IsProtonMC).c_str(),"2/(TOFEndep[2]+TOFEndep[3])","BetaHR",1000,-1.5,1.5);
	Calculate_Resolution( EdepLTOFMC_P, checkfile, treeMC, finalHistos);

	Resolution * EdepLTOFDT_P = new Resolution("EdepLTOFvsBeta Measured DT",ToFResB,(IsPreselected).c_str(),"2/(TOFEndep[2]+TOFEndep[3])","BetaHR",1000,-1.5,1.5);
	Calculate_Resolution( EdepLTOFDT_P, checkfile, treeDT, finalHistos);


	//E. dep. Inner Tracker vs beta (ToF)
	
	Resolution * EdepTrackMC_P = new Resolution("EdepTrackvsBeta Measured MC",ToFResB,(IsPreselected+ "&&" + IsProtonMC).c_str(),("7/"+EdepTrack).c_str(),"BetaHR",500,-1.5,15);
	Calculate_Resolution( EdepTrackMC_P, checkfile, treeMC, finalHistos);

	Resolution * EdepTrackDT_P = new Resolution("EdepTrackvsBeta Measured DT",ToFResB,(IsPreselected).c_str(),("7/"+EdepTrack).c_str(),"BetaHR",500,-1.5,15);
	Calculate_Resolution( EdepTrackDT_P, checkfile, treeDT, finalHistos);

	*/





	return 0;
}


