#include "PlottingFunctions/DVSMCPreSeleff_Plot.h"


using namespace std;

DatavsMC * PreSel_DvsMC_P[3];


void Initialize_DVSMCpre(){ 

	for(int i=0;i<3;i++)
		PreSel_DvsMC_P[i]= new DatavsMC(("PreSel_DvsMC_P"+to_string(i)+"_").c_str(),11,1,3);
	return;
}

void DVSMCPreSeleff_D_Fill(int zona){

	//cuts
	if(!trgpatt.IsPhysical()||Tup.R_pre<=0||Tup.R_pre<SF*Tup.Rcutoff||!ProtonsMassWindow) return;
	if(!((Tup.R_pre>Rcut[zona]&&zona<10)||(zona==10)))  return;
	if(!Herejcut) return;
	//
	int Kbin;
	for(int S=0;S<3;S++){
		for(int cut=0;cut<3;cut++){
			if((Tup.EdepTOFU<EdepTOFbeta->Eval(Tup.Beta_pre)+0.5*(cut+1)&&Tup.EdepTOFU>EdepTOFbeta->Eval(Tup.Beta_pre)-0.5*(cut+1))) 
			{	
				//R bins
				Kbin = PRB.GetRBin(RUsed);
				if(cmask.notPassed(S)) (PreSel_DvsMC_P[S] -> DataEff[cut] -> beforeR) -> Fill(Kbin,zona);	
				if(cmask.passed(S))	      (PreSel_DvsMC_P[S] -> DataEff[cut] -> afterR ) -> Fill(Kbin,zona);     

				//Beta bins
				//ToF
				Kbin=ToFDB.GetBin(RUsed);	
				if(cmask.notPassed(S)) (PreSel_DvsMC_P[S] -> DataEff[cut] -> beforeTOF) -> Fill(Kbin,zona);
				if(cmask.passed(S))       (PreSel_DvsMC_P[S] -> DataEff[cut] -> afterTOF ) -> Fill(Kbin,zona);
				//NaF
				if(cmask.isFromNaF()) {	
					Kbin=NaFDB.GetBin(RUsed);
					if(cmask.notPassed(S)) (PreSel_DvsMC_P[S] -> DataEff[cut] -> beforeNaF) -> Fill(Kbin,zona);
					if(cmask.passed(S))       (PreSel_DvsMC_P[S] -> DataEff[cut] -> afterNaF ) -> Fill(Kbin,zona);
				}
				//Agl
				if(cmask.isFromAgl()) {
					Kbin=AglDB.GetBin(RUsed);
					if(cmask.notPassed(S)) (PreSel_DvsMC_P[S] -> DataEff[cut] -> beforeAgl) -> Fill(Kbin,zona);
					if(cmask.passed(S))       (PreSel_DvsMC_P[S] -> DataEff[cut] -> afterAgl ) -> Fill(Kbin,zona);
				}
			}
		}
	}
	return;

}
void DVSMCPreSeleff_Fill(){

	if(!trgpatt.IsPhysical()||Tup.R_pre<=0||!ProtonsMassWindow) return;
	if(!Herejcut) return;
	//cuts
	//
	int Kbin;
	for(int S=0;S<3;S++){
		if(Massa_gen<1) {
			for(int cut=0;cut<3;cut++)
				if((Tup.EdepTOFU<EdepTOFbeta->Eval(Tup.Beta_pre)+0.5*(cut+1)&&Tup.EdepTOFU>EdepTOFbeta->Eval(Tup.Beta_pre)-0.5*(cut+1)))  
				{

					//R bins
					Kbin = PRB.GetRBin(RUsed);	
					if(cmask.notPassed(S))		  (PreSel_DvsMC_P[S] -> MCEff[cut] -> beforeR) -> Fill(Kbin,Tup.mcweight);
					if(cmask.passed(S))                  (PreSel_DvsMC_P[S] -> MCEff[cut] -> afterR ) -> Fill(Kbin,Tup.mcweight);
					//Beta bins

					//ToF
					Kbin=ToFDB.GetBin(RUsed);	
					if(cmask.notPassed(S)) (PreSel_DvsMC_P[S] -> MCEff[cut] -> beforeTOF) -> Fill(Kbin,Tup.mcweight);
					if(cmask.passed(S))   (PreSel_DvsMC_P[S] -> MCEff[cut] -> afterTOF ) -> Fill(Kbin,Tup.mcweight);

					//NaF
					if(cmask.isFromNaF()) {	
						Kbin=NaFDB.GetBin(RUsed);	
						if(cmask.notPassed(S)) (PreSel_DvsMC_P[S] -> MCEff[cut] -> beforeNaF) -> Fill(Kbin,Tup.mcweight);
						if(cmask.passed(S))       (PreSel_DvsMC_P[S] -> MCEff[cut] -> afterNaF ) -> Fill(Kbin,Tup.mcweight);
					}
					//Agl
					if(cmask.isFromAgl()) {	
						Kbin=AglDB.GetBin(RUsed);
						if(cmask.notPassed(S)) (PreSel_DvsMC_P[S] -> MCEff[cut] -> beforeAgl) -> Fill(Kbin,Tup.mcweight);
						if(cmask.passed(S))       (PreSel_DvsMC_P[S] -> MCEff[cut] -> afterAgl ) -> Fill(Kbin,Tup.mcweight);
					}
				}
		}

	}
	return;
}


void DVSMCPreSeleff_Write(){

	for(int i=0;i<3;i++)
		PreSel_DvsMC_P[i] -> Write();
	return;
}


void DVSMCPreSeleff(string filename){
	
	   cout<<"******* Data vs MC:  PRESELECTIONS ********"<<endl;

	   cout<<"*** Reading  P1 file ****"<<endl;
	   TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");	

	   DatavsMC * PreSel_DvsMC_P[3]; 
	  
	   PreSel_DvsMC_P[0]= new DatavsMC(inputHistoFile,"PreSel_DvsMC_P0_",1,3);
	   PreSel_DvsMC_P[1]= new DatavsMC(inputHistoFile,"PreSel_DvsMC_P1_",1,3);
	   PreSel_DvsMC_P[2]= new DatavsMC(inputHistoFile,"PreSel_DvsMC_P2_",1,3);
	    

	   LATcorr * LATpreSelDATA = new LATcorr(inputHistoFile,"LATpreSelDATA"      ,"Results");


	   cout<<"******* Data vs MC:  PRESELECTIONS ********"<<endl;

	   TH1F* PreSel_Correction_R[3]  ; 
           TH1F* PreSel_Correction_TOF[3];
           TH1F* PreSel_Correction_NaF[3];
	   TH1F* PreSel_Correction_Agl[3];

           TH1F* PreSel_CorrectionFit_R[3]  ; 
           TH1F* PreSel_CorrectionFit_TOF[3];
           TH1F* PreSel_CorrectionFit_NaF[3];
	   TH1F* PreSel_CorrectionFit_Agl[3];


	   TH1F* EffData_R[3];  		
	   TH1F* EffData_TOF[3];  
	   TH1F* EffData_NaF[3];  
	   TH1F* EffData_Agl[3];  

	   TH1F* EffMC_R[3];  		
	   TH1F* EffMC_TOF[3];  
	   TH1F* EffMC_NaF[3];  
	   TH1F* EffMC_Agl[3];  


	   for(int sel=0;sel<3;sel++){
		   PreSel_DvsMC_P[sel] -> Assign_LatCorr( LATpreSelDATA   ->  LATcorrR_fit , 
				   LATpreSelDATA   ->  LATcorrR_fit ,
				   LATpreSelDATA   ->  LATcorrR_fit ,
				   LATpreSelDATA   ->  LATcorrR_fit );

		   PreSel_DvsMC_P[sel] ->Eval_DandMC_Eff(sel);  
		   PreSel_DvsMC_P[sel] ->Eval_Corrections();

		  cout<<"******* Data vs MC ********"<<endl;		  
   	
		   PreSel_Correction_R[sel]   =(TH1F*) PreSel_DvsMC_P[sel] -> GetCorrection_R()  ;
		   PreSel_Correction_TOF[sel] =(TH1F*) PreSel_DvsMC_P[sel] -> GetCorrection_TOF();
		   PreSel_Correction_NaF[sel] =(TH1F*) PreSel_DvsMC_P[sel] -> GetCorrection_NaF();
		   PreSel_Correction_Agl[sel] =(TH1F*) PreSel_DvsMC_P[sel] -> GetCorrection_Agl();

		   EffData_R[sel]     =(TH1F*) PreSel_DvsMC_P[sel] -> GetDataEff_R(0);
		   EffData_TOF[sel]   =(TH1F*) PreSel_DvsMC_P[sel] -> GetDataEff_TOF(0);
		   EffData_NaF[sel]   =(TH1F*) PreSel_DvsMC_P[sel] -> GetDataEff_NaF(0);
		   EffData_Agl[sel]   =(TH1F*) PreSel_DvsMC_P[sel] -> GetDataEff_Agl(0);


		   EffMC_R[sel]     =(TH1F*) PreSel_DvsMC_P[sel] -> GetMCEff_R(0);
		   EffMC_TOF[sel]   =(TH1F*) PreSel_DvsMC_P[sel] -> GetMCEff_TOF(0);
		   EffMC_NaF[sel]   =(TH1F*) PreSel_DvsMC_P[sel] -> GetMCEff_NaF(0);
		   EffMC_Agl[sel]   =(TH1F*) PreSel_DvsMC_P[sel] -> GetMCEff_Agl(0);


		   cout<<"******* SYST ERROR ********"<<endl;

		   PreSel_DvsMC_P[sel] ->Initialize_SystError();					   
		   PreSel_DvsMC_P[sel] ->Eval_SystError(); 

		   cout<<"******* FIT ********"<<endl;	
		   
		   PreSel_DvsMC_P[sel]->Eval_FittedCorrections();

		   PreSel_CorrectionFit_R[sel]   =(TH1F*) PreSel_DvsMC_P[sel] -> GetCorrection_R()  ;
		   PreSel_CorrectionFit_TOF[sel] =(TH1F*) PreSel_DvsMC_P[sel] -> GetCorrection_TOF();
		   PreSel_CorrectionFit_NaF[sel] =(TH1F*) PreSel_DvsMC_P[sel] -> GetCorrection_NaF();
		   PreSel_CorrectionFit_Agl[sel] =(TH1F*) PreSel_DvsMC_P[sel] -> GetCorrection_Agl();

	
	
	   PreSel_CorrectionFit_R[sel]    -> SetName(("PreSel_DvsMC_P"+to_string(sel)+ "_CorrectionR").c_str()  );
	   PreSel_CorrectionFit_TOF[sel]  -> SetName(("PreSel_DvsMC_P" +to_string(sel)+ "_CorrectionTOF").c_str() );
	   PreSel_CorrectionFit_NaF[sel]  -> SetName(("PreSel_DvsMC_P" +to_string(sel)+ "_CorrectionNaF").c_str() );
	   PreSel_CorrectionFit_Agl[sel]  -> SetName(("PreSel_DvsMC_P" +to_string(sel)+ "_CorrectionAgl").c_str() );


	   finalHistos.Add( PreSel_Correction_R[sel]   );
	   finalHistos.Add(PreSel_Correction_TOF[sel] );
	   finalHistos.Add( PreSel_Correction_NaF[sel] );
	   finalHistos.Add( PreSel_Correction_Agl[sel] );
	   }	
	   
	   finalHistos.writeObjsInFolder("Results");
	
	   cout<<"*** Plotting ...  ****"<<endl;


	   DVSMCPreSeleff_Plot(PreSel_Correction_R  ,
			   PreSel_Correction_TOF,
			   PreSel_Correction_NaF,
			   PreSel_Correction_Agl,

			   PreSel_CorrectionFit_R  , 	
			   PreSel_CorrectionFit_TOF, 	
			   PreSel_CorrectionFit_NaF, 
			   PreSel_CorrectionFit_Agl, 

			   EffData_R,  
			   EffData_TOF,
			   EffData_NaF,
			   EffData_Agl,

			   EffMC_R    ,
			   EffMC_TOF  ,
			   EffMC_NaF  ,
			   EffMC_Agl  



				   );

	   return;
}


