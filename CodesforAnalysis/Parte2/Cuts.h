void Cuts(){
	//likelihood cut
	Likcut=false;
	if((((((int)Cutmask)>>11)==512||(((int)Cutmask)>>11)==0)&&-log(1-LDiscriminant)>2.6)||((!((((int)Cutmask)>>11)==512||(((int)Cutmask)>>11)==0))&&LDiscriminant>0.7)) Likcut=true;
	//distance cut
	Distcut=false;
	if(fabs((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D))>ddiscrcut) Distcut=true;
	//Helium rej. cut
	float EdepTOFud=(EdepTOFU+EdepTOFD)/2;
	Herejcut=false;
	if(fabs(EdepTrackbeta->Eval(Beta)-EdepTrack)/(pow(EdepTrackbeta->Eval(Beta),2)*etrack->Eval(Beta))<4.5||fabs(EdepTOFbeta->Eval(Beta)-EdepTOFud)/(pow(EdepTOFbeta->Eval(Beta),2)*etofu->Eval(Beta))<16)
	Herejcut=true;	
}


