#include "EffCorr.h"


void EffCorr::Fill(TTree * treeMC,TTree * treeDT, Variables * vars, float (*discr_var) (Variables * vars),bool refill){

	EffMC      -> Fill(treeMC,vars,discr_var,refill);
	EffMC2     -> Fill(treeMC,vars,discr_var,refill);
	EffMCnopid -> Fill(treeMC,vars,discr_var,refill);
	EffData    -> Fill(treeDT,vars,discr_var,refill);
}

void EffCorr::Save(){
	EffMC  -> Save();
	EffMC2  -> Save();
	EffMCnopid  -> Save();
	EffData-> Save();
	EffData_glob-> Save();

}

void EffCorr::Eval_Efficiencies(){
	EffMC  -> Eval_Efficiency();
	EffMC2  -> Eval_Efficiency();
	EffMCnopid  -> Eval_Efficiency();
	EffData-> Eval_Efficiency();
}

void EffCorr::SaveResults(FileSaver finalhistos){
	EffMC  -> SaveResults(finalhistos);
	EffMC2  -> SaveResults(finalhistos);
	EffMCnopid  -> SaveResults(finalhistos);
	EffData-> SaveResults(finalhistos);

	for(int lat=0;lat<10;lat++) 
		if(LatCorrections[lat]){
			finalhistos.Add(LatCorrections[lat]); 	
			finalhistos.Add(LatEfficiencies[lat]);
			finalhistos.writeObjsInFolder((directory+"/"+basename).c_str());
		}
	finalhistos.Add(GlobalEfficiency);
	finalhistos.Add(GlobalCorrection);
 	finalhistos.Add(GlobalCorrection2);
 	finalhistos.Add(GlobalCorrectionnopid);

       finalhistos.writeObjsInFolder((directory+"/"+basename).c_str());	
}


void EffCorr::Eval_Corrections(){

	for(int lat=0;lat<10;lat++) {
		LatCorrections[lat] = ProjectionXtoTH1F((TH2F*)EffData->GetEfficiency(),(basename + "_Corr_lat"+to_string(lat)).c_str(),lat+1,lat+1);
		LatCorrections[lat]->SetName((basename + "_Corr_lat"+to_string(lat)).c_str());
		LatCorrections[lat]->Sumw2();
		LatEfficiencies[lat] = ProjectionXtoTH1F((TH2F*)EffData->GetEfficiency(),(basename + "_Corr_lat"+to_string(lat)).c_str(),lat+1,lat+1);
		LatEfficiencies[lat]->SetName((basename + "_Eff_lat"+to_string(lat)).c_str());
		LatEfficiencies[lat]->Sumw2();
		LatCorrections[lat]->Divide(EffMC->GetEfficiency());
	}

/*	TH1F * Global_Before=ProjectionXtoTH1F((TH2F*)EffData->GetBefore(),(basename + "_Corr_glob").c_str(),0,10);	
	TH1F * Global_After =ProjectionXtoTH1F((TH2F*)EffData->GetAfter() ,(basename + "_Corr_glob").c_str(),0,10);	

	GlobalEfficiency = (TH1F *) Global_After->Clone();
	GlobalEfficiency -> SetName((basename + "_Eff_glob").c_str());
	GlobalEfficiency->Sumw2();
	GlobalEfficiency->Divide(Global_Before);	
*/
	GlobalEfficiency = (TH1F *) EffData_glob->GetEfficiency()->Clone();
	GlobalCorrection = (TH1F *) GlobalEfficiency->Clone();
	GlobalCorrection -> SetName((basename + "_Corr_glob").c_str());
	GlobalCorrection->Sumw2();
	
	GlobalCorrection2=(TH1F *) EffData_glob->GetEfficiency()->Clone();
	GlobalCorrectionnopid=(TH1F *) EffData_glob->GetEfficiency()->Clone();
	GlobalCorrection2->SetName((basename + "_Corr_glob2").c_str());
	GlobalCorrection2->SetTitle((basename + "_Corr_glob2").c_str());
	GlobalCorrectionnopid->SetName((basename + "_Corr_globnopid").c_str());
	GlobalCorrectionnopid->SetTitle((basename + "_Corr_globnopid").c_str());

	GlobalCorrection->Divide(EffMC->GetEfficiency());	
	GlobalCorrection2->Divide(EffMC2->GetEfficiency());	
	GlobalCorrectionnopid->Divide(EffMCnopid->GetEfficiency());	
	
	return;
}



void EffCorr::SetToConstantValue(float value){
	if(!GlobalCorrection) return;
	for(int i=0; i<GlobalCorrection->GetNbinsX(); i++) {
		GlobalCorrection->SetBinContent(i+1,value);
		GlobalCorrection->SetBinError(i+1,0.01*value);
		GlobalCorrection2->SetBinContent(i+1,value);
		GlobalCorrection2->SetBinError(i+1,0.01*value);
		GlobalCorrectionnopid->SetBinContent(i+1,value);
		GlobalCorrectionnopid->SetBinError(i+1,0.01*value);
			
	}

}

void EffCorr::SetDefaultOutFile(FileSaver FinalHistos){
		EffMC      ->SetDefaultOutFile(FinalHistos);
		EffMC2     ->SetDefaultOutFile(FinalHistos);
		EffMCnopid ->SetDefaultOutFile(FinalHistos);
		EffData	   ->SetDefaultOutFile(FinalHistos);
		EffData_glob   ->SetDefaultOutFile(FinalHistos);
		 
}

