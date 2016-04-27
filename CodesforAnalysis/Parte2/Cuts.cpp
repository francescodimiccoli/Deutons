void Cuts(){
	//likelihood cut
	Likcut=false;
	if((((((int)Cutmask)>>11)==512||(((int)Cutmask)>>11)==0)&&-log(1-LDiscriminant)>2.6)||((!((((int)Cutmask)>>11)==512||(((int)Cutmask)>>11)==0))&&-log(1-LDiscriminant)>0.55)) Likcut=true;
	//for some reason, Distance cut doesn't reject pions, adding a low mass cut
	//mass cut
	float mass=0;
	float massRICH=0;
	bool masscut=false;
	if(BetaRICH>0) massRICH=(R/BetaRICH)*pow(1-pow(BetaRICH,2),0.5); else massRICH=999;
	mass=(R/Beta)*pow(1-pow(Beta,2),0.5); 
	if(mass>0.5&&massRICH>0.5) masscut=true;
	//distance cut
	Distcut=false;
	//if(fabs((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D))>0.2&&masscut) Distcut=true;
	if((Dist5D<4||Dist5D_P<4)) Distcut=true;
	//Helium rej. cut
	float EdepTOFud=(EdepTOFU+EdepTOFD)/2;
	Herejcut=false;
	if(fabs(EdepTrackbeta->Eval(Beta)-EdepTrack)/(pow(EdepTrackbeta->Eval(Beta),2)*etrack->Eval(Beta))<4||fabs(EdepTOFbeta->Eval(Beta)-EdepTOFud)/(pow(EdepTOFbeta->Eval(Beta),2)*etofu->Eval(Beta))<10)
	Herejcut=true;	
}


