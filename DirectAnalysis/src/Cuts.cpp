#include "Cuts.h"

using namespace std;

std::vector<float> LatEdges={0.0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2};

bool IsProtonMC    (Variables * vars){ return (vars->Massa_gen<1&&vars->Massa_gen>0);}
bool IsDeutonMC    (Variables * vars){ return (vars->Massa_gen<2&&vars->Massa_gen>1);}
bool IsHeliumMC    (Variables * vars){ return (vars->Massa_gen>3&&vars->Massa_gen>1);}
bool IsTritiumMC   (Variables * vars){ return (vars->Massa_gen>2&&vars->Massa_gen<3);}

bool IsFragmentedPfromHeMC (Variables * vars) {return IsHeliumMC(vars)&&(GetPIDatL2(vars))==14&&(GetPIDatL3(vars))==14;}
bool IsFragmentedDfromHeMC (Variables * vars) {return IsHeliumMC(vars)&&(GetPIDatL2(vars))==45&&(GetPIDatL3(vars))==45;}
bool IsFragmentedTMC 	   (Variables * vars) {return IsHeliumMC(vars)&&(GetPIDatL2(vars))==46&&(GetPIDatL3(vars))==46;}
bool IsPureDMC 		   (Variables * vars) {return IsDeutonMC(vars)&&(GetPIDatL2(vars))==45&&(GetPIDatL3(vars))==45;}
bool IsPurePMC 		   (Variables * vars) {return IsProtonMC(vars)&&(GetPIDatL2(vars))==14&&(GetPIDatL3(vars))==14;}
bool IsPureTMC 		   (Variables * vars) {return IsTritiumMC(vars)&&(GetPIDatL2(vars))==46&&(GetPIDatL3(vars))==46;}

bool IsFragmentedDMC 	   (Variables * vars) {return !(IsPureDMC(vars));}
bool IsFragmentedPMC 	   (Variables * vars) {return !(IsPurePMC(vars));}

bool IsFragmentedPfromDMC  (Variables * vars) {return IsDeutonMC(vars)&&(GetPIDatL2(vars))==14&&(GetPIDatL3(vars))==14;}
bool IsPrimary	   (Variables * vars){ return (vars->R>1.2*vars->Rcutoff_RTI && IsData(vars)); }
bool IsMC          (Variables * vars){ return (vars->Massa_gen>0);} 
bool IsData        (Variables * vars){ return (vars->Massa_gen==0);}
bool IsGoodL1Status (Variables * vars) {return (vars->qL1Status==0);}
bool IsGoodL2Status (Variables * vars) {return (vars->hitbits&0x2!=0) && (vars->qL2Status==0);}

bool HasL1 (Variables * vars) {return vars->qL1>0 && IsGoodL1Status(vars);}
bool HasL2 (Variables * vars) {return vars->qL2>0 && IsGoodL2Status(vars);} //(vars->hitbits&0x2!=0);}


bool IsL1Fiducial   (Variables * vars) {return ((vars->FiducialVolume)==255);}
bool IsGoodChiSquareX(Variables * vars) { return(vars->Chisquare<vars->Chi2Xcut->Eval(abs(vars->R)));}
bool IsGoodChiSquareY(Variables * vars) { return(vars->Chisquare_y<vars->Chi2Ycut->Eval(abs(vars->R)));}

bool IsStandardSel (Variables * vars) {return (vars->P_standard_sel&0x7FFF)==0;}
//fulltracker + bit6 + bit5 + bit4 + bit3 (vars->P_standard_sel&0x3E7B)==0;
//fulltracker + bit6 + bit5 + bit4 (vars->P_standard_sel&0x3E73)==0;
//fulltracker + bit6 + bit5 (vars->P_standard_sel&0x3E63)==0;
//fulltracker + bit6 (vars->P_standard_sel&0x3E43)==0;
//ultra-basic 	(vars->P_standard_sel&0x3403)==0;	
//fulltracker   (vars->P_standard_sel&0x3E03)==0;
//basic	        (vars->P_standard_sel&0x3C7F)==0;
//complete	(vars->P_standard_sel&0x7FFF)==0; }

// one bit at a time
bool IsStandardSel_v0 (Variables * vars) {return (vars->P_standard_sel&0x3403)==0&&vars->R>0;}
bool IsStandardSel_v1 (Variables * vars) {return (vars->P_standard_sel&0x3E03)==0&&vars->R>0;}
bool IsStandardSel_v2 (Variables * vars) {return (vars->P_standard_sel&0x3E07)==0&&vars->R>0;}
bool IsStandardSel_v3 (Variables * vars) {return (vars->P_standard_sel&0x3E0B)==0&&vars->R>0;}
bool IsStandardSel_v4 (Variables * vars) {return (vars->P_standard_sel&0x3E13)==0&&vars->R>0;}
bool IsStandardSel_v5 (Variables * vars) {return (vars->P_standard_sel&0x3E23)==0&&vars->R>0;}
bool IsStandardSel_v6 (Variables * vars) {return (vars->P_standard_sel&0x3E43)==0&&vars->R>0;}
bool IsStandardSel_v7 (Variables * vars) {return (vars->P_standard_sel&0x3E83)==0&&vars->R>0;}
bool IsStandardSel_v8 (Variables * vars) {return (vars->P_standard_sel&0x3F03)==0&&vars->R>0;}
bool IsStandardSel_v9 (Variables * vars) {return (vars->P_standard_sel&0x7E03)==0&&vars->R>0;}


//cascade
/*bool IsStandardSel_v0 (Variables * vars) {return (vars->P_standard_sel&0x3E03)==0;}
bool IsStandardSel_v1 (Variables * vars) {return (vars->P_standard_sel&0x3E07)==0;}
bool IsStandardSel_v2 (Variables * vars) {return (vars->P_standard_sel&0x3E0F)==0;}
bool IsStandardSel_v3 (Variables * vars) {return (vars->P_standard_sel&0x3E1F)==0;}
bool IsStandardSel_v4 (Variables * vars) {return (vars->P_standard_sel&0x3E3F)==0;}
bool IsStandardSel_v5 (Variables * vars) {return (vars->P_standard_sel&0x3E7F)==0;}
bool IsStandardSel_v6 (Variables * vars) {return (vars->P_standard_sel&0x3F83)==0;}
bool IsStandardSel_v7 (Variables * vars) {return (vars->P_standard_sel&0x3FFF)==0;}
bool IsStandardSel_v8 (Variables * vars) {return (vars->P_standard_sel&0x7FFF)==0;}
*/

bool Is2Tracks(Variables * vars) {return (vars->NTracks==2);}

bool Qualitycut(Variables * vars, float cutvariable, float cutTOF, float cutNaF, float cutAgl){

        bool IsQual=false;
        if(IsOnlyFromToF(vars) && cutvariable<cutTOF)  IsQual=true;
        if(IsFromNaF(vars)     && cutvariable<cutNaF)  IsQual=true;
        if(IsFromAgl(vars)     && cutvariable<cutAgl)  IsQual=true;

        return IsQual;
}

bool QualChargeCut (Variables * vars){ return (vars->qInner>0.8&&vars->qInner<1.3&&vars->qUtof>0.8&&vars->qUtof<1.3&&vars->qLtof>0.8&&vars->qLtof<1.3);}
bool QualChargeCut_notrack (Variables * vars){ return (vars->qUtof>0.8&&vars->qUtof<1.3&&vars->qLtof>0.8&&vars->qLtof<1.3);}


bool DistanceCut   (Variables * vars){ return QualChargeCut(vars);}//(Qualitycut(vars,vars->DistP,3,4,4)||Qualitycut(vars,vars->DistD,3,4,4));}

bool LikelihoodCut (Variables * vars){ return ((vars->BetaRICH_new<=0)||(vars->BetaRICH_new>0&&vars->BDTDiscr>0)); }


//baseline eff. corr
bool IsDownGoing    (Variables * vars)    	{ return (vars->Beta > 0.3); } 
bool IsGoodChi2	    (Variables * vars) 		{ return ( ((int)vars->joinCutmask&16)==16);}
bool IsLUT2         (Variables * vars)          { return (((int)vars->joinCutmask&4)==4);} 
bool IsPhysTrig     (Variables * vars)		{ return ((int)vars->joinCutmask&1)==1;}
bool IsGoodTrack    (Variables * vars) 		{return vars->R!=0 && HasL2(vars);}
bool IsCharge1Track (Variables * vars) 		{return (vars->qInner>0.7&&vars->qInner<1.5);}
bool IsCharge1TrackLoose (Variables * vars) 	{return (vars->qInner>0.3&&vars->qInner<1.8);}


//efficiency corrections
bool Is1TrTrack (Variables * vars)    { return  ( ((int)vars->NTracks)==1);}
bool IsMinTOF   (Variables * vars)    { return  ( ((int)vars->joinCutmask&2)==2);} 
bool IsCharge1UTOF (Variables * vars) {return (vars->qUtof>0.8&&vars->qUtof<1.3);}
bool IsCharge1LTOF (Variables * vars) {return (vars->qLtof>0.8&&vars->qLtof<1.3);}

//analysis selections
bool IsBaseline (Variables * vars){ return IsDownGoing(vars) && IsGoodTrack(vars) && IsGoodChi2(vars) && IsCharge1Track(vars);}
bool L1LooseCharge1(Variables * vars){ return (vars->qL1>0 && vars->qL1<2&&HasL1(vars));} 
bool IsCleaning	(Variables * vars) { return Is1TrTrack(vars)&&IsMinTOF(vars)&&IsCharge1UTOF(vars)&&IsCharge1LTOF(vars);  }
bool IsGoodTime (Variables * vars) { return ( ((int)vars->joinCutmask&32)==32);}
bool IsFromNaF_nosel     (Variables * vars){ return vars->IsFromNaF_nosel();}
bool IsFromAgl_nosel     (Variables * vars){ return vars->IsFromAgl_nosel();}
bool IsFromNaF     (Variables * vars){ return vars->IsFromNaF();}
bool IsFromAgl     (Variables * vars){ return vars->IsFromAgl();}
bool RICHBDTCut (Variables * vars){ return Qualitycut(vars,-vars->BDTDiscr,999999,-0.26,-0.25);  }
//////////////////////


//He fragm
bool IsPreselectedInner (Variables * vars){ return (((int)vars->joinCutmask&187)==187&&(vars->qL1>0)&&vars->R!=0&&HasL1(vars)&&HasL2(vars));}
bool IsPreselected (Variables * vars){ return (((int)vars->joinCutmask&187)==187&&(vars->qL1>0)&&L1LooseCharge1(vars)&&vars->R!=0&&HasL1(vars)&&HasL2(vars));}
bool IsPreselectedHe (Variables * vars){ return (((int)vars->joinCutmask&187)==187&&(vars->qL1>0)&&vars->R!=0&&(vars->qL1>1.75&&vars->qL1<2.3)&&HasL1(vars))&&HasL2(vars);}
bool IsPreselectedHeStep (Variables * vars,int step){
			float chargecut=1.65;
			return (((int)vars->joinCutmask&187)==187&&(vars->qL1>0)&&vars->R!=0&&(vars->qL1>(chargecut+0.05*step)&&vars->qL1<2.3)&&HasL1(vars)&&HasL2(vars));
}

bool IsOnlyFromToF (Variables * vars){ return !((IsFromNaF(vars))||(IsFromAgl(vars)));}
bool InnerAndL1Charge2 (Variables * vars) { return (vars->qInner>1.8&&vars->qInner<2.5&&vars->qL1>1.8&&vars->qL1<2.5);}
bool ProtonsMassCut(Variables * vars){ return GetRecMassTOF(vars)>0.5&&GetRecMassTOF(vars)<1.5;}
bool DeutonsMassCut(Variables * vars){  if(IsFromNaF(vars)||IsFromAgl(vars))
						return GetRecMassRICH(vars)>1.6&&GetRecMassRICH(vars)<4.5;
				     	else
						return GetRecMassTOF(vars)>1.75&&GetRecMassTOF(vars)<4.5;
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

bool IsNegativeCharged (Variables * vars) {return (vars->Beta>0&&vars->R<0);}
bool IsPositiveCharged (Variables * vars) {return (vars->Beta>0&&vars->R>0);}



bool IsGoodHe      (Variables * vars){ return (LikelihoodCut(vars) && vars->qInner>1.6 && vars->qInner<2.6 &&  vars->qUtof>1.6 && vars->qUtof<2.6 &&  vars->qLtof>1.6 && vars->qLtof<2.6);}

bool IsInLatZone   (Variables * vars, int lat) { return (vars->Latitude>=LatEdges[lat]&&vars->Latitude<LatEdges[lat+1]);}
bool IsHighEn	(Variables * vars) {return vars->R>30;}


bool TofBetaSafetyCut (Variables * vars) {return vars->Beta<0.8 ;}
bool NafBetaSafetyCut (Variables * vars) {return vars->BetaRICH_new<0.967 ;}
bool AglBetaSafetyCut (Variables * vars) {return vars->BetaRICH_new<0.99;}

bool IsClearQ1ExceptL2 (Variables * vars) {return (vars->qL2>0&&vars->qUtof>0.8&&vars->qUtof<1.3&&vars->qLtof>0.8&&vars->qLtof<1.3&&vars->qL1InnerNoL2>0.8&&vars->qL1InnerNoL2<1.3&&HasL1(vars)&&HasL2(vars));} 
bool IsClearQ2ExceptL2 (Variables * vars) {return (vars->qL2>0&&vars->qUtof>1.8&&vars->qUtof<2.3&&vars->qLtof>1.8&&vars->qLtof<2.3&&vars->qL1InnerNoL2>1.8&&vars->qL1InnerNoL2<2.3&&HasL1(vars)&&HasL2(vars));} 


//// Tracking Efficiency
//layerZ = [160.,52.985,29.185,25.215,1.685,-2.285,-25.215,-29.185,-136.0]
bool IsExtrapolInsideL9 (Variables * vars) {
	if(vars->entrypointcoo[0]==0 &&vars->entrypointcoo[1]==0) return true;
	float X_extrapol = vars->entrypointcoo[0] + (( -136.0-vars->entrypointcoo[2])/cos(vars->theta_track))*cos(vars->phi_track)*sin(vars->theta_track);
	float Y_extrapol = vars->entrypointcoo[1] + (( -136.0-vars->entrypointcoo[2])/cos(vars->theta_track))*sin(vars->phi_track)*sin(vars->theta_track);
	return (Y_extrapol>-45&&Y_extrapol<45);
}


bool IsExtrapolInsideL8 (Variables * vars) {
	if(vars->entrypointcoo[0]==0 &&vars->entrypointcoo[1]==0) return true;
	float X_extrapol = vars->entrypointcoo[0] + (( -29.185-vars->entrypointcoo[2])/cos(vars->theta_track))*cos(vars->phi_track)*sin(vars->theta_track);
	float Y_extrapol = vars->entrypointcoo[1] + (( -29.185-vars->entrypointcoo[2])/cos(vars->theta_track))*sin(vars->phi_track)*sin(vars->theta_track);
	return (Y_extrapol>-45&&Y_extrapol<45);
}

bool IsExtrapolInsideL1 (Variables * vars) {
	if(vars->entrypointcoo[0]==0 &&vars->entrypointcoo[1]==0) return true;
	float X_extrapol = vars->entrypointcoo[0] + (( 160.-vars->entrypointcoo[2])/cos(vars->theta_track))*cos(vars->phi_track)*sin(vars->theta_track);
	float Y_extrapol = vars->entrypointcoo[1] + (( 160.-vars->entrypointcoo[2])/cos(vars->theta_track))*sin(vars->phi_track)*sin(vars->theta_track);
	return (Y_extrapol>-45&&Y_extrapol<45);
}

bool IsGoodTOFStandaloneQ1(Variables * vars) {	
	return (vars->beta_SA>0 && vars->qUtof_SA > 0.8 && vars->qUtof_SA < 1.3 && vars->qLtof_SA > 0.8 && vars->qLtof_SA < 1.3 && vars-> qTrd_SA > 0.6 && vars-> qTrd_SA < 1.7);

} 

//// L1 pick-up efficicny

bool IsL1HitNearExtrapol (Variables * vars) {
	return (vars->hitdistfromint<3);
}


bool IsCleanL1Hit (Variables * vars) { 

	return ( vars->exthit_closest_status==0 && fabs(vars->exthit_closest_q -1)<0.4 );
}


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
		if(spl[i]=="IsPhysTrig"    ) IsPassed=IsPassed && IsPhysTrig    (Vars);
		if(spl[i]=="IsProtonMC"	   ) IsPassed=IsPassed && IsProtonMC    (Vars);
		if(spl[i]=="IsDeutonMC"	   ) IsPassed=IsPassed && IsDeutonMC    (Vars);
		if(spl[i]=="IsHeliumMC"	   ) IsPassed=IsPassed && IsHeliumMC    (Vars);
		if(spl[i]=="IsTritiumMC"   ) IsPassed=IsPassed && IsTritiumMC    (Vars);


		if(spl[i]=="IsFragmentedPfromHeMC") IsPassed=IsPassed && IsFragmentedPfromHeMC    (Vars);
		if(spl[i]=="IsFragmentedDfromHeMC") IsPassed=IsPassed && IsFragmentedDfromHeMC    (Vars);
		if(spl[i]=="IsFragmentedTMC") 	    IsPassed=IsPassed && IsFragmentedTMC    (Vars);
		if(spl[i]=="IsPureDMC") 	    IsPassed=IsPassed && IsPureDMC    (Vars);
		if(spl[i]=="IsPurePMC") 	    IsPassed=IsPassed && IsPurePMC    (Vars);
		if(spl[i]=="IsPureTMC") 	    IsPassed=IsPassed && IsPureTMC    (Vars);
		
		if(spl[i]=="IsFragmentedDMC") 	    IsPassed=IsPassed && IsFragmentedDMC    (Vars);
		if(spl[i]=="IsFragmentedPMC") 	    IsPassed=IsPassed && IsFragmentedPMC    (Vars);
		
		if(spl[i]=="IsFragmentedPfromDMC")  IsPassed=IsPassed && IsFragmentedPfromDMC    (Vars);
		if(spl[i]=="IsGoodL1Status")            IsPassed=IsPassed && IsGoodL1Status     (Vars);
		if(spl[i]=="IsGoodL2Status")            IsPassed=IsPassed && IsGoodL2Status     (Vars);	
		if(spl[i]=="HasL1")            IsPassed=IsPassed && HasL1     (Vars);
		if(spl[i]=="HasL2")            IsPassed=IsPassed && HasL2     (Vars);	
	
		if(spl[i]=="IsL1Fiducial")            IsPassed=IsPassed && IsL1Fiducial     (Vars);	
		if(spl[i]=="IsGoodChiSquareX")            IsPassed=IsPassed && IsGoodChiSquareX     (Vars);	
		if(spl[i]=="IsGoodChiSquareY")            IsPassed=IsPassed && IsGoodChiSquareY     (Vars);	

		if(spl[i]=="IsDownGoing" )            IsPassed=IsPassed && IsDownGoing    (Vars);
		if(spl[i]=="IsGoodTrack" )            IsPassed=IsPassed && IsGoodTrack    (Vars);	
		

		if(spl[i]=="IsClearQ1ExceptL2")   IsPassed=IsPassed && IsClearQ1ExceptL2     (Vars);
		if(spl[i]=="IsClearQ2ExceptL2")   IsPassed=IsPassed && IsClearQ2ExceptL2     (Vars);

		if(spl[i]=="IsPrimary"	   ) IsPassed=IsPassed && IsPrimary     (Vars);
		if(spl[i]=="IsMC"          ) IsPassed=IsPassed && IsMC          (Vars);
		if(spl[i]=="IsData"	   ) IsPassed=IsPassed && IsData        (Vars);
	
		if(spl[i]=="IsPreselectedInner" ) IsPassed=IsPassed && IsPreselectedInner (Vars);
		if(spl[i]=="IsPreselected" ) IsPassed=IsPassed && IsPreselected (Vars);
		if(spl[i]=="IsPreselectedHe" ) IsPassed=IsPassed && IsPreselectedHe (Vars);
		for(int j=0;j<10;j++) 
			if(spl[i]==("IsPreselectedHe"+to_string(j)) )  IsPassed=IsPassed && IsPreselectedHeStep (Vars,j);
		
		if(spl[i]=="IsBaseline" ) IsPassed=IsPassed && IsBaseline (Vars); 
		if(spl[i]=="IsMinTOF" )      IsPassed=IsPassed && IsMinTOF (Vars); 
		if(spl[i]=="IsStandardSel" )      IsPassed=IsPassed && IsStandardSel (Vars); 
	
		if(spl[i]=="IsStandardSel_v0" )      IsPassed=IsPassed && IsStandardSel_v0 (Vars); 
		if(spl[i]=="IsStandardSel_v1" )      IsPassed=IsPassed && IsStandardSel_v1 (Vars); 
		if(spl[i]=="IsStandardSel_v2" )      IsPassed=IsPassed && IsStandardSel_v2 (Vars); 
		if(spl[i]=="IsStandardSel_v3" )      IsPassed=IsPassed && IsStandardSel_v3 (Vars); 
		if(spl[i]=="IsStandardSel_v4" )      IsPassed=IsPassed && IsStandardSel_v4 (Vars); 
		if(spl[i]=="IsStandardSel_v5" )      IsPassed=IsPassed && IsStandardSel_v5 (Vars); 
		if(spl[i]=="IsStandardSel_v6" )      IsPassed=IsPassed && IsStandardSel_v6 (Vars); 
		if(spl[i]=="IsStandardSel_v7" )      IsPassed=IsPassed && IsStandardSel_v7 (Vars); 
		if(spl[i]=="IsStandardSel_v8" )      IsPassed=IsPassed && IsStandardSel_v8 (Vars); 
		if(spl[i]=="IsStandardSel_v9" )      IsPassed=IsPassed && IsStandardSel_v9 (Vars); 
	

		if(spl[i]=="IsCleaning" ) IsPassed=IsPassed && IsCleaning (Vars);

		if(spl[i]=="IsOnlyFromToF" ) IsPassed=IsPassed && IsOnlyFromToF (Vars);
		if(spl[i]=="IsFromNaF"	   ) IsPassed=IsPassed && IsFromNaF     (Vars);
		if(spl[i]=="IsFromAgl"	   ) IsPassed=IsPassed && IsFromAgl     (Vars);
		if(spl[i]=="IsFromNaF_nosel"	   ) IsPassed=IsPassed && IsFromNaF_nosel     (Vars);
		if(spl[i]=="IsFromAgl_nosel"	   ) IsPassed=IsPassed && IsFromAgl_nosel     (Vars);
		if(spl[i]=="RICHBDTCut"	   ) IsPassed=IsPassed && RICHBDTCut(Vars);

		if(spl[i]=="L1LooseCharge1") IsPassed=IsPassed && L1LooseCharge1(Vars);
		if(spl[i]=="InnerAndL1Charge2")   IsPassed=IsPassed && InnerAndL1Charge2(Vars);
		if(spl[i]=="DistanceCut")    IsPassed=IsPassed && DistanceCut(Vars);
		if(spl[i]=="LikelihoodCut")  IsPassed=IsPassed && LikelihoodCut(Vars);
		if(spl[i]=="QualChargeCut")  IsPassed=IsPassed && QualChargeCut(Vars);
		if(spl[i]=="QualChargeCut_notrack")  IsPassed=IsPassed && QualChargeCut_notrack(Vars);
		if(spl[i]=="ProtonsMassCut") IsPassed=IsPassed && ProtonsMassCut(Vars);
		if(spl[i]=="DeutonsMassCut") IsPassed=IsPassed && DeutonsMassCut(Vars);
		if(spl[i]=="TemplatesMassCut")IsPassed=IsPassed && TemplatesMassCut(Vars);
		if(spl[i]=="IsGoodHe")       IsPassed=IsPassed && IsGoodHe(Vars);
		if(spl[i]=="TofBetaSafetyCut")  IsPassed=IsPassed && TofBetaSafetyCut(Vars);
		if(spl[i]=="NafBetaSafetyCut")  IsPassed=IsPassed && NafBetaSafetyCut(Vars);
		if(spl[i]=="AglBetaSafetyCut")  IsPassed=IsPassed && AglBetaSafetyCut(Vars);
		if(spl[i]=="IsNegative") IsPassed=IsPassed && IsNegativeCharged(Vars);
		if(spl[i]=="IsPositive") IsPassed=IsPassed && IsPositiveCharged(Vars);
		
		if(spl[i]=="IsHighEn") IsPassed=IsPassed && IsHighEn(Vars);
		if(spl[i]=="IsLUT2") IsPassed=IsPassed && IsLUT2(Vars);
	
		if(spl[i]=="Is2Tracks") IsPassed=IsPassed && Is2Tracks(Vars);
	
		//eff. corr.
		if(spl[i]=="IsGoodChi2") IsPassed=IsPassed && IsGoodChi2(Vars);
		if(spl[i]=="IsGoodTime") IsPassed=IsPassed && IsGoodTime(Vars);
		if(spl[i]=="Is1TrTrack") IsPassed=IsPassed && Is1TrTrack(Vars);
		if(spl[i]=="IsCharge1Track") IsPassed=IsPassed && IsCharge1Track(Vars);
		if(spl[i]=="IsCharge1TrackLoose") IsPassed=IsPassed && IsCharge1TrackLoose(Vars);
	
		if(spl[i]=="IsCharge1UTOF") IsPassed=IsPassed && IsCharge1UTOF(Vars);
		if(spl[i]=="IsCharge1LTOF") IsPassed=IsPassed && IsCharge1LTOF(Vars);

		//tracking eff.	
  		if(spl[i]=="IsExtrapolInsideL8") IsPassed=IsPassed && IsExtrapolInsideL8(Vars);
		if(spl[i]=="IsExtrapolInsideL1") IsPassed=IsPassed && IsExtrapolInsideL1(Vars);
		if(spl[i]=="IsExtrapolInsideL9") IsPassed=IsPassed && IsExtrapolInsideL9(Vars);
		if(spl[i]=="IsGoodTOFStandaloneQ1") IsPassed=IsPassed && IsGoodTOFStandaloneQ1(Vars);

		//L1 pick-up eff.	
  		if(spl[i]=="IsL1HitNearExtrapol") IsPassed=IsPassed && IsL1HitNearExtrapol(Vars);
		if(spl[i]=="IsCleanL1Hit") IsPassed=IsPassed && IsCleanL1Hit(Vars);



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

