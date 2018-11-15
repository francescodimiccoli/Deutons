#include "particle.h"


void  Particle::SetSlowDownModel(TF1 * SlowDownModel, float alpha_slowdown, float gamma_slowdown) { 
		slowdownmodel = SlowDownModel; 
		slowdownmodel->SetParameter(0,alpha_slowdown);
		slowdownmodel->SetParameter(1,gamma_slowdown);
		return;
		}    


float Particle::BetaTOIFromBeta (float beta){
	TF1 model("slowdown", "(x-[2]*(x-x*(1-[0]/x^[1]))) - [3]",beta*0.5,1) ;
	model.SetParameter(0,slowdownmodel->GetParameter(0));
	model.SetParameter(1,slowdownmodel->GetParameter(1));
	model.SetParameter(2,0.9/mass);
	model.SetParameter(3,beta);

	WrappedTF1 wf1(model);
	ROOT::Math::BrentRootFinder brf;

	brf.SetFunction( wf1, beta*0.5, 1 );
   	brf.Solve();
	float betaTOI = brf.Root();

	return betaTOI;
}	
 


void Particle::FillFromEk (float ek)
{
   ekin=ek;
   ekpermass=EkPerMassFromEk(ekin);
   etot=  EtotfromEk(ekin);
   beta = BetaFromEk (ekin);
   mom=MomFromEk (ekin);
   rig=RigFromMom (mom);

   beta_TOI = BetaTOIFromBeta(beta);
   ekin_TOI = EkFromBeta(beta_TOI);
   ekpermass_TOI=EkPerMassFromEk(ekin_TOI);
   etot_TOI=  EtotfromEk(ekin_TOI);
   mom_TOI=MomFromEk (ekin_TOI);
   rig_TOI=RigFromMom (mom_TOI);

}

void Particle::FillFromRig ( float r)
{
	rig=r;
	mom=MomFromRig (rig);
	ekin=EkFromMom (mom);
	ekpermass=EkPerMassFromEk(ekin);
	etot=EtotfromEk(ekin);  
	beta=BetaFromEk (ekin);

	beta_TOI = BetaTOIFromBeta(beta);
	if((beta_TOI-beta)>0.005){
		ekin_TOI = EkFromBeta(beta_TOI);
		ekpermass_TOI=EkPerMassFromEk(ekin_TOI);
		etot_TOI=  EtotfromEk(ekin_TOI);
		mom_TOI=MomFromEk (ekin_TOI);
		rig_TOI=RigFromMom (mom_TOI);
	}	
	else{
		rig_TOI=r;
		mom_TOI=MomFromRig (rig_TOI);
		ekin_TOI=EkFromMom (mom_TOI);
		ekpermass_TOI=EkPerMassFromEk(ekin_TOI);
		etot_TOI=EtotfromEk(ekin_TOI);  
		beta_TOI=BetaFromEk (ekin_TOI);
	}	
}

void Particle::FillFromBeta ( float Beta)
{
   beta = Beta;
   ekin = EkFromBeta(beta);
   ekpermass=EkPerMassFromEk(ekin);
   etot=  EtotfromEk(ekin);
   mom=MomFromEk (ekin);
   rig=RigFromMom (mom);

   beta_TOI = BetaTOIFromBeta(beta);
   ekin_TOI = EkFromBeta(beta_TOI);
   ekpermass_TOI=EkPerMassFromEk(ekin_TOI);
   etot_TOI=  EtotfromEk(ekin_TOI);
   mom_TOI=MomFromEk (ekin_TOI);
   rig_TOI=RigFromMom (mom_TOI);
}




void Particle::FillFromEkPerMass(float ekpermass)
{
   ekin=ekpermass * mass;
   FillFromEk(ekin);
}

