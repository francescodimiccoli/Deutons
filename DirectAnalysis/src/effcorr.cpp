#include "EffCorr.h"


void EffCorr::Fill(TTree * treeMC,TTree * treeDT, Variables * vars, float (*discr_var) (Variables * vars),bool refill){

	EffMC      -> Fill(treeMC,vars,discr_var,refill);
	EffMCpid2     -> Fill(treeMC,vars,discr_var,refill);
	EffMCpid -> Fill(treeMC,vars,discr_var,refill);
	EffData    -> Fill(treeDT,vars,discr_var,refill);
}

void EffCorr::Save(){
	EffMC  -> Save();
	EffMCpid2  -> Save();
	EffMCpid  -> Save();
	EffData-> Save();
	EffData_glob-> Save();

}

void EffCorr::Eval_Efficiencies(){
	EffMC  -> Eval_Efficiency();
	EffMCpid2  -> Eval_Efficiency();
	EffMCpid  -> Eval_Efficiency();
	if(!IsTrigEffCorr) {
		EffData-> Eval_Efficiency();
		EffData_glob-> Eval_Efficiency();
	}
	else{
		EffData-> Eval_TrigEfficiency();
		EffData_glob-> Eval_TrigEfficiency();
	}

}

void EffCorr::SaveResults(FileSaver finalhistos){
	EffMC  -> SaveResults(finalhistos);
	EffMCpid2  -> SaveResults(finalhistos);
	EffMCpid  -> SaveResults(finalhistos);
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
 	finalhistos.Add(GlobalCorrectionpid);
 	finalhistos.Add(Stat_Err);
 	finalhistos.Add(Syst_Err);



       finalhistos.writeObjsInFolder((directory+"/"+basename).c_str());	
}

void AddCorrSystError(TH1F * Correction) {
	for(int i=0; i<Correction->GetNbinsX(); i++) {
		float stat = Correction->GetBinError(i+1);
		float syst = fabs(Correction->GetBinContent(i+1) - 1)/pow(3,0.5);
		Correction->SetBinError(i+1, pow(pow(stat,2)+pow(syst,2),0.5));
	}
}

void EffCorr::Eval_Errors(){

	Stat_Err = (TH1F *) GlobalCorrection->Clone();
	Syst_Err = (TH1F *) GlobalCorrection->Clone();
	Stat_Err->SetName((basename + "Stat_Err").c_str());
	Stat_Err->SetTitle((basename +"Stat_Err").c_str());
	Syst_Err->SetName((basename + "Syst_Err").c_str());
	Syst_Err->SetTitle((basename +"Syst_Err").c_str());


	for(int i=0; i<GlobalCorrection->GetNbinsX(); i++) {
		float stat = GlobalCorrectionpid->GetBinError(i+1);
		float syst = fabs(GlobalCorrectionpid->GetBinContent(i+1) - 1)/pow(3,0.5);
		if(GlobalCorrection->GetBinContent(i+1)>0){
			Stat_Err->SetBinContent(i+1, stat / GlobalCorrectionpid->GetBinContent(i+1));
			Syst_Err->SetBinContent(i+1, syst / GlobalCorrectionpid->GetBinContent(i+1));
		}
		Stat_Err->SetBinError(i+1, 0);	
		Syst_Err->SetBinError(i+1, 0);	
	}

	AddCorrSystError(GlobalCorrection);
	AddCorrSystError(GlobalCorrection2);
	AddCorrSystError(GlobalCorrectionpid);

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

	cout<<EffData_glob->GetEfficiency()<<endl;
	GlobalEfficiency = (TH1F *) EffData_glob->GetEfficiency()->Clone();
	GlobalCorrection = (TH1F *) GlobalEfficiency->Clone();
	GlobalCorrection -> SetName((basename + "_Corr_glob").c_str());
	GlobalCorrection->Sumw2();
	
	GlobalCorrection2=(TH1F *) EffData_glob->GetEfficiency()->Clone();
	GlobalCorrectionpid=(TH1F *) EffData_glob->GetEfficiency()->Clone();
	GlobalCorrection2->SetName((basename + "_Corr_glob2").c_str());
	GlobalCorrection2->SetTitle((basename + "_Corr_glob2").c_str());
	GlobalCorrectionpid->SetName((basename + "_Corr_globpid").c_str());
	GlobalCorrectionpid->SetTitle((basename + "_Corr_globpid").c_str());


	GlobalCorrection->Divide(EffMC->GetEfficiency());	
	GlobalCorrection2->Divide(EffMCpid2->GetEfficiency());	
	GlobalCorrectionpid->Divide(EffMCpid->GetEfficiency());	

	Eval_Errors();
	return;
}



void EffCorr::SetToConstantValue(float value){
	if(!GlobalCorrection) return;
	for(int i=0; i<GlobalCorrection->GetNbinsX(); i++) {
		GlobalCorrection->SetBinContent(i+1,value);
		GlobalCorrection->SetBinError(i+1,0.01*value);
		GlobalCorrection2->SetBinContent(i+1,value);
		GlobalCorrection2->SetBinError(i+1,0.01*value);
		GlobalCorrectionpid->SetBinContent(i+1,value);
		GlobalCorrectionpid->SetBinError(i+1,0.01*value);
			
	}

}

void EffCorr::SetDefaultOutFile(FileSaver FinalHistos){
		EffMC      ->SetDefaultOutFile(FinalHistos);
		EffMCpid2     ->SetDefaultOutFile(FinalHistos);
		EffMCpid ->SetDefaultOutFile(FinalHistos);
		EffData	   ->SetDefaultOutFile(FinalHistos);
		EffData_glob   ->SetDefaultOutFile(FinalHistos);
		 
}

