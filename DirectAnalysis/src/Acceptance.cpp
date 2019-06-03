#include "Acceptance.h"
#include "histUtils.h"

void MCPar::Eval_trigrate(){
	std::vector<float> events;
	std::vector<float> triggers;
	std::ifstream infile;
	rundb* rdb= new rundb();
	int ret=rdb->readdb(("/cvmfs/ams.cern.ch/Offline/AMSDataDir/DataManagement/DataSetsDesc/"+filename).c_str());
	//rdb->Print();
    	rdb->Summary();	
	cout<<"********MC infos:*******"<<endl;
	cout<<"Trig. Rate: "<< rdb->GetTrigRate()<<endl; //rdb->GetTotEvents()/rdb->GetTotTrigg()<<endl;//rdb->GetTrigRate()<<endl;		
        cout<<"Total Gen: "<<rdb->GetTotTrigg()<<endl; 
	cout<<"Total events: "<<rdb->GetTotEvents()<<endl;
	cout<<"************************"<<endl;
        Trigrate=rdb->GetTrigRate();
	tot_ev = rdb->GetTotEvents();	
	tot_trig = rdb->GetTotTrigg();
}

void Acceptance::SetDefaultOutFile(FileSaver FinalHistos){
	FullSetEff      ->SetDefaultOutFile(FinalHistos);
	For_Acceptance     ->SetDefaultOutFile(FinalHistos);
}


void Acceptance::Set_MCPar(float rmin, float rmax, float Gen_factor, std::string Filename, float Art_ratio){
	cout<<"Setting MC parameters"<<endl;
	param.Rmin=rmin;
	param.Rmax=rmax;
	param.filename = Filename;
	param.gen_factor = Gen_factor;
	param.Eval_trigrate();
	param.art_ratio = Art_ratio;
}


void Acceptance::ApplyEfficCorr(EffCorr * Correction){
	if(Correction->GetGlobCorrection()) { 
		cout<<"Correction: "<<basename<<" "<<Correction->GetGlobCorrection()->GetEntries()<<endl;
		EfficiencyCorrections.push_back(Correction);
	}
	return;
}

void Acceptance::ApplyEfficFromData(EffCorr * Correction){
	if(Correction->GetGlobCorrection()) { 
		cout<<"Correction: "<<basename<<" "<<Correction->GetGlobCorrection()->GetEntries()<<endl;
		EfficiencyFromData.push_back(Correction);
	}
	return;
}

void Acceptance::Save(){
	FullSetEff  -> Save();
	For_Acceptance  -> Save();

}

void Acceptance::SaveResults(FileSaver finalhistos){
	FullSetEff  -> SaveResults(finalhistos);
	For_Acceptance  -> SaveResults(finalhistos);

	finalhistos.Add(Acc_StatErr);
	finalhistos.Add(Acc_SystErr);
 	finalhistos.Add(EffAcceptance ); 
 	finalhistos.Add(EffAcceptanceMC);
        finalhistos.writeObjsInFolder((directory+"/"+basename).c_str());	
}




