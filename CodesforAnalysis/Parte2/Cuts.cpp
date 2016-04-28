void Cuts(){
	/////////////////////Analysis cuts//////////////////////////
	
	//likelihood cut
	Likcut=false;
	if((((((int)Cutmask)>>11)==512||(((int)Cutmask)>>11)==0)&&-log(1-LDiscriminant)>2.6)||((!((((int)Cutmask)>>11)==512||(((int)Cutmask)>>11)==0))&&-log(1-LDiscriminant)>0.55)) Likcut=true;
	
	//Distance cut
	Distcut=false;
	if((Dist5D<4||Dist5D_P<4)) Distcut=true;
	
	////////////////////////////////////////////////////////////
	
	
	/////////////////////Control sample cuts//////////////////////
	
	//Helium rej. cut
	float EdepTOFud=(EdepTOFU+EdepTOFD)/2;
	Herejcut=false;
	if(fabs(EdepTrackbeta->Eval(Beta)-EdepTrack)/(pow(EdepTrackbeta->Eval(Beta),2)*etrack->Eval(Beta))<4||fabs(EdepTOFbeta->Eval(Beta)-EdepTOFud)/(pow(EdepTOFbeta->Eval(Beta),2)*etofu->Eval(Beta))<10)
	Herejcut=true;	

	//Strong beta cut
	Betastrongcut = false;
	if(	(((((int)Cutmask)>>11)==512)&&BetaRICH<0.97 )
	      ||(((((int)Cutmask)>>11)== 0 )&&BetaRICH<0.985)	
	      ||((!((((int)Cutmask)>>11)==512||(((int)Cutmask)>>11)==0))&& Beta < 0.8)
	  ) 	
	Betastrongcut = true;
}


