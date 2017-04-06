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
	

	if(spline)  finalHistos.Add(Reso->ModelSigmaWithSpline());
	else 	    finalHistos.Add(Reso->ModelSigmaWithPoly());
	
   	finalHistos.writeObjsInFolder((Reso->GetName()+"/Fit Results").c_str());				
	
	return Reso->CheckHistos();
}


int main(int argc, char * argv[])
{


	cout<<"****************************** FILES OPENING ***************************************"<<endl;

	string INPUT(argv[1]);
        string OUTPUT(argv[2]);

	FileSaver finalHistos;
        finalHistos.setName(OUTPUT.c_str());
	bool checkfile = finalHistos.CheckFile();	


        TFile *fileMC =TFile::Open(INPUT.c_str());
        TTree *treeMC = (TTree *)fileMC->Get("parametri_geo");


	cout<<"****************************** BINS ***************************************"<<endl;

	SetBins();	

	PRB.Print();
	DRB.Print();

	cout<<"**TOF**"<<endl;
	ToFDB.Print();
	ToFPB.Print();

	cout<<"**NaF**"<<endl;
	NaFDB.Print();
	NaFPB.Print();

	cout<<"**Agl**"<<endl;
	AglDB.Print();
	AglPB.Print();

	ToFDB.UseBetaEdges();
	ToFPB.UseBetaEdges();
	NaFDB.UseBetaEdges();
	NaFPB.UseBetaEdges();
	AglDB.UseBetaEdges();
	AglPB.UseBetaEdges();

	DRB.UseREdges();
	PRB.UseREdges();


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
	Resolution * RigidityResolution_P = new Resolution("RvsR Resolution (P)",PRB,(IsPreselected + "&&" + IsProtonMC).c_str(),"1/R - 1/GenMomentum","GenMomentum",1000,-0.5,1.5);
	Calculate_Resolution( RigidityResolution_P, checkfile, treeMC, finalHistos,PRB.RigBinsCent(),true);

	Resolution * RigidityResolution_D = new Resolution("RvsR Resolution (D)",PRB,(IsPreselected + "&&" + IsDeutonMC).c_str(),"1/R - 1/GenMomentum","GenMomentum",1000,-0.5,1.5);
        Calculate_Resolution( RigidityResolution_D, checkfile, treeMC, finalHistos,PRB.RigBinsCent(),true);


	//rigidity resolution vs beta
	Resolution * RigidityTOFResolution_P = new Resolution("RvsBetaTOF Resolution (P)",ToFPB,(IsPreselected+ "&&" + IsProtonMC).c_str(),"1/R - 1/GenMomentum",Beta_gen.c_str(),1000,-0.5,1.5);
	Calculate_Resolution( RigidityTOFResolution_P, checkfile, treeMC, finalHistos,ToFPB.RigBinsCent());

	Resolution * RigidityNaFResolution_P = new Resolution("RvsBetaNaF Resolution (P)",NaFPB,(IsPreselected + "&&" + IsProtonMC + "&&" + IsFromNaF ).c_str(),"1/R - 1/GenMomentum",Beta_gen.c_str(),1000,-0.5,1.5);
	Calculate_Resolution( RigidityNaFResolution_P, checkfile, treeMC, finalHistos,NaFPB.RigBinsCent());

	Resolution * RigidityAglResolution_P = new Resolution("RvsBetaAgl Resolution (P)",AglPB,(IsPreselected + "&&" + IsProtonMC + "&&" + IsFromAgl ).c_str(),"1/R - 1/GenMomentum",Beta_gen.c_str(),1000,-0.5,1.5);
	Calculate_Resolution( RigidityAglResolution_P, checkfile, treeMC, finalHistos,AglPB.RigBinsCent());

	Resolution * RigidityTOFResolution_D = new Resolution("RvsBetaTOF Resolution (D)",ToFDB,(IsPreselected+ "&&" + IsDeutonMC).c_str(),"1/R - 1/GenMomentum",Beta_gen.c_str(),1000,-0.5,1.5);
	Calculate_Resolution( RigidityTOFResolution_D, checkfile, treeMC, finalHistos,ToFDB.RigBinsCent());

	Resolution * RigidityNaFResolution_D = new Resolution("RvsBetaNaF Resolution (D)",NaFDB,(IsPreselected + "&&" + IsDeutonMC + "&&" + IsFromNaF ).c_str(),"1/R - 1/GenMomentum",Beta_gen.c_str(),1000,-0.5,1.5);
	Calculate_Resolution( RigidityNaFResolution_D, checkfile, treeMC, finalHistos,NaFDB.RigBinsCent());

	Resolution * RigidityAglResolution_D = new Resolution("RvsBetaAgl Resolution (D)",AglDB,(IsPreselected + "&&" + IsDeutonMC + "&&" + IsFromAgl ).c_str(),"1/R - 1/GenMomentum",Beta_gen.c_str(),1000,-0.5,1.5);
	Calculate_Resolution( RigidityAglResolution_D, checkfile, treeMC, finalHistos,AglDB.RigBinsCent());


	//Beta resolution vs rigidity	
	Resolution * BetaTOF_RResolution_P = new Resolution("BetaTOFvsR Resolution (P)",PRB,(IsPreselected+ "&&" + IsProtonMC + "&&" + IsOnlyFromToF).c_str(),("1/BetaHR -" + InverseBeta_gen).c_str(),"GenMomentum",1000,-0.5,1.5);
	Calculate_Resolution( BetaTOF_RResolution_P, checkfile, treeMC, finalHistos,PRB.BetaBinsCent(),true);

	Resolution * BetaNaF_RResolution_P = new Resolution("BetaNaFvsR Resolution (P)",PRB,(IsPreselected + "&&" + IsProtonMC + "&&" + IsFromNaF   ).c_str(),("1/BetaRICH -" + InverseBeta_gen).c_str(),"GenMomentum",250,-0.05,0.15);
	Calculate_Resolution( BetaNaF_RResolution_P, checkfile, treeMC, finalHistos,PRB.BetaBinsCent(),true);

	Resolution * BetaAgl_RResolution_P = new Resolution("BetaAglvsR Resolution (P)",PRB,(IsPreselected + "&&" + IsProtonMC + "&&" + IsFromAgl   ).c_str(),("1/BetaRICH -" + InverseBeta_gen).c_str(),"GenMomentum",250,-0.075,0.075);
	Calculate_Resolution( BetaAgl_RResolution_P, checkfile, treeMC, finalHistos,PRB.BetaBinsCent(),true);

	Resolution * BetaTOF_RResolution_D = new Resolution("BetaTOFvsR Resolution (D)",PRB,(IsPreselected+ "&&" + IsDeutonMC + "&&" + IsOnlyFromToF).c_str(),("1/BetaHR -" + InverseBeta_gen).c_str(),"GenMomentum",1000,-0.5,1.5);
	Calculate_Resolution( BetaTOF_RResolution_D, checkfile, treeMC, finalHistos,PRB.BetaBinsCent(),true);

	Resolution * BetaNaF_RResolution_D = new Resolution("BetaNaFvsR Resolution (D)",PRB,(IsPreselected + "&&" + IsDeutonMC + "&&" + IsFromNaF   ).c_str(),("1/BetaRICH -" + InverseBeta_gen).c_str(),"GenMomentum",250,-0.05,0.15);
	Calculate_Resolution( BetaNaF_RResolution_D, checkfile, treeMC, finalHistos,PRB.BetaBinsCent(),true);

	Resolution * BetaAgl_RResolution_D = new Resolution("BetaAglvsR Resolution (D)",PRB,(IsPreselected + "&&" + IsDeutonMC + "&&" + IsFromAgl   ).c_str(),("1/BetaRICH -" + InverseBeta_gen).c_str(),"GenMomentum",250,-0.075,0.075);
	Calculate_Resolution( BetaAgl_RResolution_D, checkfile, treeMC, finalHistos,PRB.BetaBinsCent(),true);



	//beta resolution vs beta (ToF, NaF, Agl)
	Resolution * BetaTOFResolution_P = new Resolution("BetaTOFvsBeta Resolution (P)",ToFPB,(IsPreselected+ "&&" + IsProtonMC).c_str(),("1/BetaHR -" + InverseBeta_gen).c_str(),Beta_gen.c_str(),1000,-0.5,1.5);
	Calculate_Resolution( BetaTOFResolution_P, checkfile, treeMC, finalHistos);

	Resolution * BetaNaFResolution_P = new Resolution("BetaNaFvsBeta Resolution (P)",NaFPB,(IsPreselected + "&&" + IsProtonMC + "&&" + IsFromNaF ).c_str(),("1/BetaRICH -" + InverseBeta_gen).c_str(),Beta_gen.c_str(),250,-0.05,0.15);
	Calculate_Resolution( BetaNaFResolution_P, checkfile, treeMC, finalHistos);

	Resolution * BetaAglResolution_P = new Resolution("BetaAglvsBeta Resolution (P)",AglPB,(IsPreselected + "&&" + IsProtonMC + "&&" + IsFromAgl ).c_str(),("1/BetaRICH -" + InverseBeta_gen).c_str(),Beta_gen.c_str(),250,-0.075,0.075);
	Calculate_Resolution( BetaAglResolution_P, checkfile, treeMC, finalHistos);

	Resolution * BetaTOFResolution_D = new Resolution("BetaTOFvsBeta Resolution (D)",ToFDB,(IsPreselected+ "&&" + IsDeutonMC).c_str(),("1/BetaHR -" + InverseBeta_gen).c_str(),Beta_gen.c_str(),1000,-0.5,1.5);
	Calculate_Resolution( BetaTOFResolution_D, checkfile, treeMC, finalHistos);

	Resolution * BetaNaFResolution_D = new Resolution("BetaNaFvsBeta Resolution (D)",NaFDB,(IsPreselected + "&&" + IsDeutonMC + "&&" + IsFromNaF ).c_str(),("1/BetaRICH -" + InverseBeta_gen).c_str(),Beta_gen.c_str(),250,-0.05,0.15);
	Calculate_Resolution( BetaNaFResolution_D, checkfile, treeMC, finalHistos);

	Resolution * BetaAglResolution_D = new Resolution("BetaAglvsBeta Resolution (D)",AglDB,(IsPreselected + "&&" + IsDeutonMC + "&&" + IsFromAgl ).c_str(),("1/BetaRICH -" + InverseBeta_gen).c_str(),Beta_gen.c_str(),250,-0.075,0.075);
	Calculate_Resolution( BetaAglResolution_D, checkfile, treeMC, finalHistos);


	//E. dep. U. TOF vs beta (ToF)
	
	Resolution * EdepUTOFResolution_P = new Resolution("EdepUTOFvsBeta Resolution (P)",ToFPB,(IsPreselected+ "&&" + IsProtonMC).c_str(),"2/(TOFEndep[0]+TOFEndep[1])",Beta_gen.c_str(),1000,-1.5,1.5);
	Calculate_Resolution( EdepUTOFResolution_P, checkfile, treeMC, finalHistos);

	Resolution * EdepUTOFResolution_D = new Resolution("EdepUTOFvsBeta Resolution (D)",ToFPB,(IsPreselected+ "&&" + IsDeutonMC).c_str(),"2/(TOFEndep[0]+TOFEndep[1])",Beta_gen.c_str(),1000,-1.5,1.5);
	Calculate_Resolution( EdepUTOFResolution_D, checkfile, treeMC, finalHistos);


	//E. dep. L. TOF vs beta (ToF)
	
	Resolution * EdepLTOFResolution_P = new Resolution("EdepLTOFvsBeta Resolution (P)",ToFPB,(IsPreselected+ "&&" + IsProtonMC).c_str(),"2/(TOFEndep[2]+TOFEndep[3])",Beta_gen.c_str(),1000,-1.5,1.5);
	Calculate_Resolution( EdepLTOFResolution_P, checkfile, treeMC, finalHistos);

	Resolution * EdepLTOFResolution_D = new Resolution("EdepLTOFvsBeta Resolution (D)",ToFPB,(IsPreselected+ "&&" + IsDeutonMC).c_str(),"2/(TOFEndep[2]+TOFEndep[3])",Beta_gen.c_str(),1000,-1.5,1.5);
	Calculate_Resolution( EdepLTOFResolution_D, checkfile, treeMC, finalHistos);

	//E. dep. Inner Tracker vs beta (ToF)
	
	Resolution * EdepTrack_P = new Resolution("EdepTrackvsBeta Resolution (P)",ToFPB,(IsPreselected+ "&&" + IsProtonMC).c_str(),("7/"+EdepTrack).c_str(),Beta_gen.c_str(),500,-1.5,15);
	Calculate_Resolution( EdepTrack_P, checkfile, treeMC, finalHistos);

	Resolution * EdepTrack_D = new Resolution("EdepTrackvsBeta Resolution (D)",ToFPB,(IsPreselected+ "&&" + IsDeutonMC).c_str(),("7/"+EdepTrack).c_str(),Beta_gen.c_str(),500,-1.5,15);
	Calculate_Resolution( EdepTrack_D, checkfile, treeMC, finalHistos);





	return 0;
}


