void Cuts_Pre(){
	/////////////////////Control sample cuts//////////////////////
	
	//Helium rej. cut
	float EdepTOFud=(Tup.EdepTOFU+Tup.EdepTOFD)/2;
        Herejcut=false;
        if(fabs(EdepTrackbeta->Eval(Tup.Beta_pre)-Tup.EdepTrack)/(pow(EdepTrackbeta->Eval(Tup.Beta_pre),2)*etrack->Eval(Tup.Beta_pre))<4||fabs(EdepTOFbeta->Eval(Tup.Beta_pre)-EdepTOFud)/(pow(EdepTOFbeta->Eval(Tup.Beta_pre),2)*etofu->Eval(Tup.Beta_pre))<10)
        Herejcut=true;

	//Strong beta cut
	Betastrongcut = false;
        if(     (((((int)Tup.Cutmask)>>11)==512)&&Tup.BetaRICH<0.97 )
              ||(((((int)Tup.Cutmask)>>11)== 0 )&&Tup.BetaRICH<0.985)
              ||((!((((int)Tup.Cutmask)>>11)==512||(((int)Tup.Cutmask)>>11)==0))&& Tup.Beta_pre < 0.8)
          )
        Betastrongcut = true;

	return;
}





void Cuts(){
	/////////////////////Analysis cuts//////////////////////////
	
	//likelihood cut
	Likcut=false;
	if((((((int)Tup.Cutmask)>>11)==512||(((int)Tup.Cutmask)>>11)==0)&&-log(1-Tup.LDiscriminant)>2.6)||((!((((int)Tup.Cutmask)>>11)==512||(((int)Tup.Cutmask)>>11)==0))&&-log(1-Tup.LDiscriminant)>0.55)) Likcut=true;
	
	//Distance cut
	Distcut=false;
	if((Tup.Dist5D<4||Tup.Dist5D_P<4)) Distcut=true;
	
	////////////////////////////////////////////////////////////
	
	
	/////////////////////Control sample cuts//////////////////////
	
	//Helium rej. cut
	float EdepTOFud=(Tup.EdepTOFU+Tup.EdepTOFD)/2;
	Herejcut=false;
	if(fabs(EdepTrackbeta->Eval(Tup.Beta)-Tup.EdepTrack)/(pow(EdepTrackbeta->Eval(Tup.Beta),2)*etrack->Eval(Tup.Beta))<4||fabs(EdepTOFbeta->Eval(Tup.Beta)-EdepTOFud)/(pow(EdepTOFbeta->Eval(Tup.Beta),2)*etofu->Eval(Tup.Beta))<10)
	Herejcut=true;	

	//Strong beta cut
	Betastrongcut = false;
	if(	(((((int)Tup.Cutmask)>>11)==512)&&Tup.BetaRICH<0.97 )
	      ||(((((int)Tup.Cutmask)>>11)== 0 )&&Tup.BetaRICH<0.985)	
	      ||((!((((int)Tup.Cutmask)>>11)==512||(((int)Tup.Cutmask)>>11)==0))&& Tup.Beta < 0.8)
	  ) 	
	Betastrongcut = true;

	return;
}


