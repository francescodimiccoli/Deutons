#include "BetaSmearing.h"

float SmearBetaRICH(Variables * vars){
       
	float angle;
	float betasmear;
	float Beta = vars->BetaRICH_new;

        if(Beta<0.87) angle= acos(1/(1.25*Beta))*10e4;
        else          angle= acos(1/(1.15*Beta))*10e4;
        if(Beta<0.87) angle = angle + (-191.5) + Rand->Gaus(0,529.8);
        else          angle = angle + (-42.5) + Rand->Gaus(0,65.16);
	if(Beta<0.87) betasmear = 1/(1.25*cos(angle/10e4));
        else          betasmear = 1/(1.15*cos(angle/10e4));

	return betasmear;
}



float SmearBeta(Variables * vars){
	float Beta = vars->Beta;
        float time = 1.2/(Beta*3e-4);
        float tailcontrolfactor=110./90.;
	if(vars->R<2.7) tailcontrolfactor=1;//migration tail fixing

        time = time + (-21.2) + Rand->Gaus(0,tailcontrolfactor*98.8);
        return 1.2/(time*3e-4);

}

float GetSmearedBetaTOF  (Variables * vars) {   
	if(vars->Massa_gen>0) return SmearBeta(vars);
	else return GetBetaTOF(vars);
}

float GetSmearedBetaRICH  (Variables * vars) {   
	if(vars->Massa_gen>0) return SmearBetaRICH(vars);
	else return GetBetaRICH(vars);
}

float GetSmearedRecMassTOF(Variables * vars) { 
	if(vars->Massa_gen>0){		
		float smearbeta=GetSmearedBetaTOF(vars);   
		return (vars->R/smearbeta)*pow((1-pow(smearbeta,2)),0.5);
	}
	else return GetRecMassTOF(vars);
}

