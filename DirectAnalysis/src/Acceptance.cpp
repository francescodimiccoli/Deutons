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
	if(Acc_StatErr)    finalhistos.Add(Acc_StatErr);
	if(Acc_SystErr)    finalhistos.Add(Acc_SystErr);
 	if(EffAcceptance)  finalhistos.Add(EffAcceptance ); 
 	if(EffAcceptanceMC)finalhistos.Add(EffAcceptanceMC);
        finalhistos.writeObjsInFolder((directory+"/"+basename).c_str());	
}



void Acceptance:: EvalEffAcc(){
	cout<<"********** MC flux reweighting **********"<<endl;
	Variables * vars = new Variables();
	cout<<"********** For_Acceptance **************"<<endl;
	float normalization = For_Acceptance->GetBefore()->GetEntries()*(pow(param.Trigrate,-1)) ;
	cout<<For_Acceptance->GetBefore()->GetEntries()<<" "<<normalization<<" "<<param.tot_trig<<endl;
	//float normalization = param.tot_trig/FRAC;
	
	//no_reweighting
	/*for(int i=0;i<bins.size();i++){
		if(EffAcc -> GetBinContent(i+1)>0){
			cout<<"********** MC flux reweighting **********"<<endl;
			float range;
			range = log(param.Rmax)-log(param.Rmin);
			
			float gen_Bins= normalization*(log(bins.RigBin(i+1))-log(bins.RigBin(i)))/range;  
			EffAcc -> SetBinContent(i+1,EffAcc -> GetBinContent(i+1)/gen_Bins); 
			EffAcc -> SetBinError(i+1,0);//pow(GenSpectrum -> GetBinContent(i+1),0.5)/gen_Bins); 			
			
		
		}
	}*/
	/////////////// 

	//with reweighting
	Histogram Spectrum = vars->reweighter.getTo();
	Histogram LogNorm  = vars->reweighter.getFrom();

	cout<<"********** Gen. Spectrum **********"<<endl;
	FullSetEff->GetBefore()->Reset();
	//total triggers in range
	for(int i=0;i<bins.size();i++){
		float bincontent = normalization*(log(bins.RigTOIBins()[i+1])-log(bins.RigTOIBins()[i]))/log(param.Rmax)-log(param.Rmin);
		float meanweight = 1;//Spectrum.integrate(bins.RigTOIBins()[i],bins.RigTOIBins()[i+1]) / LogNorm.integrate(bins.RigTOIBins()[i],bins.RigTOIBins()[i+1]);
		FullSetEff->GetBefore()->SetBinContent(i+1,param.art_ratio*bincontent*meanweight); //vars->reweighter.getWeight(bins.RigTOIBinsCent()[i]));
		FullSetEff->GetBefore()->SetBinError(i+1,0);
	}
	//Denominator->Scale(param.art_ratio*rangefactor*normalization/(Denominator->Integral()));

	cout<<"********** THE DIVISION **********"<<endl;
	FullSetEff->Eval_Efficiency();
	EffAcceptance= (TH1F*)FullSetEff->GetEfficiency()->Clone();
	EffAcceptanceMC= (TH1F*)FullSetEff->GetEfficiency()->Clone();	
	
	EffAcceptance->SetName((basename +"_Eff_Acceptance").c_str());
	EffAcceptance->SetTitle((basename +"_Eff_Acceptance").c_str());
	EffAcceptance -> Sumw2();
	EffAcceptance -> Scale(47.78/param.gen_factor);

	EffAcceptanceMC->SetName((basename +"_Eff_AcceptanceMC").c_str());
	EffAcceptanceMC->SetTitle((basename +"_Eff_AcceptanceMC").c_str());
	EffAcceptanceMC -> Sumw2();
	EffAcceptanceMC -> Scale(47.78/param.gen_factor);

	cout<<"*********** Efficiency corrections ************"<<endl;
	for(int i=0;i<EfficiencyCorrections.size();i++) {
		for(int j=0; j<EffAcceptance->GetNbinsX(); j++){
			if(EfficiencyCorrections[i]->IsEkin()) EffAcceptance -> SetBinContent( j+1, EffAcceptance -> GetBinContent(j+1)*EfficiencyCorrections[i]->GetCorrectionModel()->Eval(bins.ekpermassbincent_TOI[j]));
			else EffAcceptance -> SetBinContent( j+1, EffAcceptance -> GetBinContent(j+1)*EfficiencyCorrections[i]->GetCorrectionModel()->Eval(bins.rigbincent_TOI[j]) );
		}
	}	

	cout<<"*********** Eff. Acceptance parameters ************"<<endl;
	cout<<"TOT. Ev. gen: "<<normalization<<endl;
	cout<<"Trig. rate: "<<param.Trigrate<<endl;
}






