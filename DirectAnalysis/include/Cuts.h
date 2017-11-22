#ifndef CUTS_H
#define CUTS_H

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <iterator>
#include "Commonglobals.cpp"
#include "Variables.hpp"

using namespace std;

std::vector<float> LatEdges={0.0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2};

TH1F* ProjectionXtoTH1F(TH2F* h2, string title, int binmin, int binmax);
TH1F* TH1DtoTH1F(TH1D* hd);

bool IsProtonMC    (Variables * vars){ return (vars->Massa_gen<1&&vars->Massa_gen>0);}
bool IsDeutonMC    (Variables * vars){ return (vars->Massa_gen<2&&vars->Massa_gen>1);}
bool IsHeliumMC    (Variables * vars){ return (vars->Massa_gen>2&&vars->Massa_gen>1);}

bool IsFragmentedPfromHeMC (Variables * vars) {return IsHeliumMC(vars)&&(GetPIDatL2(vars))==14&&(GetPIDatL3(vars))==14;}
bool IsFragmentedDfromHeMC (Variables * vars) {return IsHeliumMC(vars)&&(GetPIDatL2(vars))==45&&(GetPIDatL3(vars))==45;}
bool IsFragmentedTMC 	   (Variables * vars) {return IsHeliumMC(vars)&&(GetPIDatL2(vars))==46&&(GetPIDatL3(vars))==46;}
bool IsPureDMC 		   (Variables * vars) {return IsDeutonMC(vars)&&(GetPIDatL2(vars))==45&&(GetPIDatL3(vars))==45;}
bool IsFragmentedPfromDMC  (Variables * vars) {return IsDeutonMC(vars)&&(GetPIDatL2(vars))==14&&(GetPIDatL3(vars))==14;}

bool L1LooseCharge1(Variables * vars){ return (vars->qL1>0 && vars->qL1<1.6);} 
bool IsPrimary	   (Variables * vars){ return (vars->R>1.3*vars->Rcutoff); }
bool IsMC          (Variables * vars){ return (vars->Massa_gen>0);} 
bool IsData        (Variables * vars){ return (vars->Massa_gen==0);}
bool IsPreselectedHe (Variables * vars){ return (((int)vars->joinCutmask&187)==187&&(vars->EdepL1>0||vars->qL1>0)&&vars->R!=0)&&(vars->qL1>1.75&&vars->qL1<2.3);}
bool IsPreselectedInner (Variables * vars){ return (((int)vars->joinCutmask&187)==187&&(vars->EdepL1>0||vars->qL1>0)&&vars->R!=0);}
bool IsPreselected (Variables * vars){ return (((int)vars->joinCutmask&187)==187&&(vars->EdepL1>0||vars->qL1>0)&&vars->R!=0)&&L1LooseCharge1(vars);}
bool IsMinimumBias (Variables * vars){ return (((int)vars->joinCutmask&139)==139&&(vars->EdepL1>0||vars->qL1>0)&&vars->R!=0)&&L1LooseCharge1(vars);}

bool IsFromNaF     (Variables * vars){ return vars->IsFromNaF();}
bool IsFromAgl     (Variables * vars){ return vars->IsFromAgl();}
bool IsOnlyFromToF (Variables * vars){ return !((IsFromNaF(vars))||(IsFromAgl(vars)));}
bool InnerAndL1Charge2 (Variables * vars) { return (vars->qInner>1.8&&vars->qInner<2.5&&vars->qL1>1.8&&vars->qL1<2.5);}
bool ProtonsMassCut(Variables * vars){ return GetRecMassTOF(vars)>0.5&&GetRecMassTOF(vars)<1.5;}
bool DeutonsMassCut(Variables * vars){  if(IsFromNaF(vars)||IsFromAgl(vars))
						return GetRecMassRICH(vars)>1.6&&GetRecMassRICH(vars)<4.5;
				     	else
						return GetRecMassTOF(vars)>1.6&&GetRecMassTOF(vars)<4.5;
				     }
bool TemplatesMassCut(Variables * vars){  if(IsFromNaF(vars)||IsFromAgl(vars))
						return GetRecMassRICH(vars)>0.1&&GetRecMassRICH(vars)<4.5;
				     	else
						return GetRecMassTOF(vars)>0.1&&GetRecMassTOF(vars)<4.5;
				     }

bool ControlSampleMassCut(Variables * vars){  if(IsFromNaF(vars)||IsFromAgl(vars))
						return GetRecMassRICH(vars)>0.2;
				     	else
						return GetRecMassTOF(vars)>0.2;
				     }

bool IsTRDok(Variables * vars){	return (vars->EdepTRD>(1/EdepTRDbeta->Eval(vars->Beta))-100&&vars->EdepTRD<(1/EdepTRDbeta->Eval(vars->Beta))+100);}
bool IsNegativeCharged (Variables * vars) {return (vars->Beta>0&&vars->R<0);}


bool Qualitycut(Variables * vars, float cutvariable, float cutTOF, float cutNaF, float cutAgl){

        bool IsQual=false;
        if(IsOnlyFromToF(vars) && cutvariable<cutTOF)  IsQual=true;
        if(IsFromNaF(vars)     && cutvariable<cutNaF)  IsQual=true;
        if(IsFromAgl(vars)     && cutvariable<cutAgl)  IsQual=true;

        return IsQual;
}

bool QualChargeCut (Variables * vars){ return (vars->qInner>0.8&&vars->qInner<1.3&&vars->qUtof>0.8&&vars->qUtof<1.3&&vars->qLtof>0.8&&vars->qLtof<1.3);}

bool DistanceCut   (Variables * vars){ return QualChargeCut(vars);}//(Qualitycut(vars,vars->DistP,3,4,4)||Qualitycut(vars,vars->DistD,3,4,4));}
bool LikelihoodCut (Variables * vars){ return true;/*Qualitycut(vars,-vars->Likelihood,-0.83,-0.85,-0.85);*/}

bool IsGoodHe      (Variables * vars){ return (LikelihoodCut(vars) && vars->qInner>1.6 && vars->qInner<2.6 &&  vars->qUtof>1.6 && vars->qUtof<2.6 &&  vars->qLtof>1.6 && vars->qLtof<2.6);}

bool IsInLatZone   (Variables * vars, int lat) { return (vars->Latitude>=LatEdges[lat]&&vars->Latitude<LatEdges[lat+1]);}

bool ControlSample (Variables * vars) { return (IsPreselected(vars)&&vars->qInner>0.2&&vars->qInner<1.75&&ControlSampleMassCut(vars));}
bool PresControlSample (Variables * vars) { return (IsMinimumBias(vars)&&vars->qInner>0.2&&vars->qInner<1.75&&ControlSampleMassCut(vars));}

bool TofBetaSafetyCut (Variables * vars) {return vars->Beta<0.77 ;}
bool NafBetaSafetyCut (Variables * vars) {return vars->Beta<0.96 ;}
bool AglBetaSafetyCut (Variables * vars) {return vars->Beta<0.985;}




template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}


bool ApplyCuts(std::string cut, Variables * Vars){

	std::vector<std::string> spl = split(cut,'&');
	
	bool IsPassed = true;
	for(int i=0;i<spl.size();i++){
		if(spl[i]=="") 		     IsPassed=IsPassed;
		if(spl[i]=="IsProtonMC"	   ) IsPassed=IsPassed && IsProtonMC    (Vars);
		if(spl[i]=="IsDeutonMC"	   ) IsPassed=IsPassed && IsDeutonMC    (Vars);
		if(spl[i]=="IsHeliumMC"	   ) IsPassed=IsPassed && IsHeliumMC    (Vars);

		if(spl[i]=="IsFragmentedPfromHeMC") IsPassed=IsPassed && IsFragmentedPfromHeMC    (Vars);
		if(spl[i]=="IsFragmentedDfromHeMC") IsPassed=IsPassed && IsFragmentedDfromHeMC    (Vars);
		if(spl[i]=="IsFragmentedTMC") 	    IsPassed=IsPassed && IsFragmentedTMC    (Vars);
		if(spl[i]=="IsPureDMC") 	    IsPassed=IsPassed && IsPureDMC    (Vars);
		if(spl[i]=="IsFragmentedPfromDMC")  IsPassed=IsPassed && IsFragmentedPfromDMC    (Vars);

		if(spl[i]=="IsPrimary"	   ) IsPassed=IsPassed && IsPrimary     (Vars);
		if(spl[i]=="IsMC"          ) IsPassed=IsPassed && IsMC          (Vars);
		if(spl[i]=="IsData"	   ) IsPassed=IsPassed && IsData        (Vars);
		if(spl[i]=="IsPreselected" ) IsPassed=IsPassed && IsPreselected (Vars); 
		if(spl[i]=="IsPreselectedInner" ) IsPassed=IsPassed && IsPreselectedInner (Vars);
		if(spl[i]=="IsPreselectedHe" ) IsPassed=IsPassed && IsPreselectedHe (Vars);
		if(spl[i]=="IsMinimumBias" ) IsPassed=IsPassed && IsMinimumBias (Vars); 
		if(spl[i]=="IsOnlyFromToF" ) IsPassed=IsPassed && IsOnlyFromToF (Vars);
		if(spl[i]=="IsFromNaF"	   ) IsPassed=IsPassed && IsFromNaF     (Vars);
		if(spl[i]=="IsFromAgl"	   ) IsPassed=IsPassed && IsFromAgl     (Vars);
		if(spl[i]=="L1LooseCharge1") IsPassed=IsPassed && L1LooseCharge1(Vars);
		if(spl[i]=="InnerAndL1Charge2")   IsPassed=IsPassed && InnerAndL1Charge2(Vars);
		if(spl[i]=="DistanceCut")    IsPassed=IsPassed && DistanceCut(Vars);
		if(spl[i]=="LikelihoodCut")  IsPassed=IsPassed && LikelihoodCut(Vars);
		if(spl[i]=="QualChargeCut")  IsPassed=IsPassed && QualChargeCut(Vars);
		if(spl[i]=="ProtonsMassCut") IsPassed=IsPassed && ProtonsMassCut(Vars);
		if(spl[i]=="DeutonsMassCut") IsPassed=IsPassed && DeutonsMassCut(Vars);
		if(spl[i]=="TemplatesMassCut")IsPassed=IsPassed && TemplatesMassCut(Vars);
		if(spl[i]=="IsGoodHe")       IsPassed=IsPassed && IsGoodHe(Vars);
		if(spl[i]=="ControlSample")  IsPassed=IsPassed && ControlSample(Vars);
		if(spl[i]=="PresControlSample")  IsPassed=IsPassed && PresControlSample(Vars);
		if(spl[i]=="TofBetaSafetyCut")  IsPassed=IsPassed && TofBetaSafetyCut(Vars);
		if(spl[i]=="NafBetaSafetyCut")  IsPassed=IsPassed && NafBetaSafetyCut(Vars);
		if(spl[i]=="AglBetaSafetyCut")  IsPassed=IsPassed && AglBetaSafetyCut(Vars);
		if(spl[i]=="IsTRDok")  		IsPassed=IsPassed && IsTRDok(Vars);
		if(spl[i]=="IsNegativeCharged") IsPassed=IsPassed && IsNegativeCharged(Vars);
		
		for(int lat=0;lat<10;lat++)
			if(spl[i]==("IsInLatZone"+to_string(lat)).c_str())      IsPassed=IsPassed && IsInLatZone(Vars,lat);

	}

	return IsPassed;
}

int GetLatitude(Variables * vars){
	int latzone=-9;
	for (int lat=0;lat<10;lat++)
		if(fabs(vars->Latitude)>=LatEdges[lat]&&fabs(vars->Latitude)<LatEdges[lat+1])  latzone=lat;
		return latzone;		
}

TH1F* ProjectionXtoTH1F(TH2F* h2, string title, int binmin, int binmax) {
   TH1D* hd=h2->ProjectionX(title.data(),binmin, binmax);
   TH1F* hf=TH1DtoTH1F(hd);
   return hf;
}

TH1F* TH1DtoTH1F(TH1D* hd) {
   TH1F*  hf=(TH1F*)hd->Clone();
   hf->SetName(hd->GetTitle());
    return hf;
}

#endif
