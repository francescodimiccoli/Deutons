void Cuts_Pre()
{
   /////////////////////Control sample cuts//////////////////////

   //Helium rej. cut
   float EdepTOFud=(Tup.EdepTOFU+Tup.EdepTOFD)/2;
   Herejcut=false;
   if(fabs(EdepTrackbeta->Eval(Tup.Beta_pre)-Tup.EdepTrack)/(pow(EdepTrackbeta->Eval(Tup.Beta_pre),2)*etrack->Eval(Tup.Beta_pre))<3||fabs(EdepTOFbeta->Eval(Tup.Beta_pre)-EdepTOFud)/(pow(EdepTOFbeta->Eval(Tup.Beta_pre),2)*etofu->Eval(Tup.Beta_pre))<10)
      Herejcut=true;	

   // L1 He cut
   IsHeL1=false;
   if(	  (cmask.isOnlyFromToF()&&((Tup.EdepL1-EdepL1beta->Eval(Tup.Beta_pre))/(pow(EdepL1beta->Eval(Tup.Beta_pre),2)*etrack->Eval(Tup.Beta_pre)) > 14.5))
	||(cmask.isFromNaF()&&((Tup.EdepL1-EdepL1beta->Eval(Tup.BetaRICH))/(pow(EdepL1beta->Eval(Tup.BetaRICH),2)*etrack->Eval(Tup.BetaRICH)) > 14.5))
	||(cmask.isFromAgl()&&((Tup.EdepL1-EdepL1beta->Eval(Tup.BetaRICH))/(pow(EdepL1beta->Eval(Tup.BetaRICH),2)*etrack->Eval(Tup.BetaRICH)) > 14.5)) 
	)	IsHeL1=true; 	


   //Strong beta cut
   Betastrongcut = false;
   if(     ((cmask.isFromNaF())&&Tup.BetaRICH<0.97 )
           ||(cmask.isFromAgl()&&Tup.BetaRICH<0.985 )
           ||((!(cmask.isFromNaF()||cmask.isFromAgl()))&& Tup.Beta_pre < 0.77)
     )
      Betastrongcut = true;

   ProtonsMassWindow = false;
   	if(Tup.Beta_pre<protons->Eval(Tup.R_pre)+0.1 && Tup.Beta_pre>protons->Eval(Tup.R_pre)-0.1) ProtonsMassWindow = true;	

  ProtonsMassThres = false;
        if(Tup.Beta_pre>protons->Eval(Tup.R_pre)) ProtonsMassThres = true;
   return;
}


bool Qualitycut(float cutvariable, float cutTOF, float cutNaF, float cutAgl){

	bool IsQual=false;	
	if(cmask.isOnlyFromToF() && cutvariable<cutTOF)  IsQual=true;
   	if(cmask.isFromNaF()	 && cutvariable<cutNaF)  IsQual=true;
   	if(cmask.isFromAgl()     && cutvariable<cutAgl)  IsQual=true;

	return IsQual;
}


void Cuts()
{
   /////////////////////Analysis cuts//////////////////////////

   //likelihood cut
   Likcut=false;
   Likcut=Qualitycut(log(1-Tup.LDiscriminant),-1,-2,-2.4);

   //Distance cut protons

   Distcut=false;
   if(Qualitycut(Tup.Dist5D_P,3,2,3)||Qualitycut(Tup.Dist5D,3,2,3)) Distcut=true;

   /////////////////////Control sample cuts//////////////////////

   //Helium rej. cut
   float EdepTOFud=(Tup.EdepTOFU+Tup.EdepTOFD)/2;
   Herejcut=false;
   if(fabs(EdepTrackbeta->Eval(Tup.Beta)-Tup.EdepTrack)/(pow(EdepTrackbeta->Eval(Tup.Beta),2)*etrack->Eval(Tup.Beta))<3||fabs(EdepTOFbeta->Eval(Tup.Beta)-EdepTOFud)/(pow(EdepTOFbeta->Eval(Tup.Beta),2)*etofu->Eval(Tup.Beta))<10)
      Herejcut=true;

	
	// L1 He cut
	IsHeL1=false;
   if(	(cmask.isOnlyFromToF()&&((Tup.EdepL1-EdepL1beta->Eval(Tup.Beta))/(pow(EdepL1beta->Eval(Tup.Beta),2)*etrack->Eval(Tup.Beta)) > 14.5))
        ||(cmask.isFromNaF()&&((Tup.EdepL1-EdepL1beta->Eval(Tup.BetaRICH))/(pow(EdepL1beta->Eval(Tup.BetaRICH),2)*etrack->Eval(Tup.BetaRICH)) > 14.5))
        ||(cmask.isFromAgl()&&((Tup.EdepL1-EdepL1beta->Eval(Tup.BetaRICH))/(pow(EdepL1beta->Eval(Tup.BetaRICH),2)*etrack->Eval(Tup.BetaRICH)) > 14.5))
        )       IsHeL1=true;


	//L1 He->P cut
	
	IsPfromHeL1=false;
   if(((Tup.EdepL1-EdepL1beta->Eval(Tup.Beta))/(pow(EdepL1beta->Eval(Tup.Beta),2)*etrack->Eval(Tup.Beta)) < 6)
        ||(cmask.isFromNaF()&&((Tup.EdepL1-EdepL1beta->Eval(Tup.BetaRICH))/(pow(EdepL1beta->Eval(Tup.BetaRICH),2)*etrack->Eval(Tup.BetaRICH)) < 6))
        ||(cmask.isFromAgl()&&((Tup.EdepL1-EdepL1beta->Eval(Tup.BetaRICH))/(pow(EdepL1beta->Eval(Tup.BetaRICH),2)*etrack->Eval(Tup.BetaRICH)) < 6))
        )       IsPfromHeL1=true;

	



   //Strong beta cut
   Betastrongcut = false;
   if(	((cmask.isFromNaF())&&Tup.BetaRICH<0.97 )
         ||(cmask.isFromAgl()&&Tup.BetaRICH<0.985)
         ||((!(cmask.isFromNaF()||cmask.isFromAgl()))&& Tup.Beta < 0.77)
     )
      Betastrongcut = true;

   ProtonsMassWindow = false;
        if(Tup.Beta<protons->Eval(Tup.R)+0.1 && Tup.Beta>protons->Eval(Tup.R)-0.1) ProtonsMassWindow = true;

   ProtonsMassThres = false;
        if(Tup.Beta>protons->Eval(Tup.R)) ProtonsMassThres = true;	
   return;
}
