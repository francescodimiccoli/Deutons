#include "Acceptance.h"
#include "histUtils.h"

void MCPar::Eval_trigrate(){
	std::vector<int> runs;
	std::vector<float> events;
	std::vector<float> triggers;

	rundb* rdb= new rundb();
	int ret=rdb->readdb(("/cvmfs/ams.cern.ch/Offline/AMSDataDir/DataManagement/DataSetsDesc/"+filename).c_str());
	//rdb->Print();
    	rdb->Summary();	
	
	std::ifstream infile(runlist.c_str());
	std::string line;
	while (std::getline(infile, line))
		{
		std::istringstream iss(line);
		int start = line.find("//");
		int end =  line.find(".root");
		line.erase(0,start+2);
		line.erase(line.find(".root"),line.find(".root")+5);
		stringstream geek(line);
		int run = 0;
		geek >> run;
		runs.push_back(run);
		}
	for(int i=0;i<runs.size();i++) {
		rundb_el * elem = rdb->find(runs[i]);
		if(elem!=0){ 
			 //cout<<runs[i]<<endl;
			 events.push_back(elem->GetEv());
 			 triggers.push_back(elem->GetTrig());
			}
		else cout<<"NOT FOUND: "<<runs[i]<<endl;
	}
	tot_ev=0;
	tot_trig=0;
	Trigrate=0;
	for(int i=0; i<events.size(); i++)   tot_ev+=events[i];	
	for(int i=0; i<triggers.size(); i++) tot_trig+=triggers[i];	
	Trigrate = tot_ev/(double)tot_trig;

	cout<<"********MC infos:*******"<<endl;
	cout<<"*** MC DISTR: ****"<<endl;
	cout<<"MC name: "<<filename.c_str()<<endl;
	cout<<"Trig. Rate: "<< rdb->GetTrigRate()<<endl; 		
        cout<<"Total Gen: "<<rdb->GetTotTrigg()<<endl; 
	cout<<"Total events: "<<rdb->GetTotEvents()<<endl;
	cout<<"*** LIST: ****"<<endl;
	cout<<"Runs found in List: "<<events.size()<<" "<<triggers.size()<<endl;
	cout<<"Gen. in List: "<<tot_trig<<endl;
	cout<<"Events in List: "<<tot_ev<<endl;
	cout<<"Trig. Rate: "<<Trigrate<<endl;
	cout<<"************************"<<endl;
}

void Acceptance::SetDefaultOutFile(FileSaver FinalHistos){
	FullSetEff      ->SetDefaultOutFile(FinalHistos);
	FullSetEff_gen     ->SetDefaultOutFile(FinalHistos);
	For_Acceptance     ->SetDefaultOutFile(FinalHistos);
}


void Acceptance::Set_MCPar(float rmin, float rmax, float Gen_factor, std::string Filename, std::string Runlist, float Art_ratio){
	cout<<"Setting MC parameters"<<endl;
	param.Rmin=rmin;
	param.Rmax=rmax;
	param.filename = Filename;
	param.runlist = Runlist;
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
	FullSetEff  	-> Save();
	FullSetEff_gen  -> Save();
	For_Acceptance  -> Save();
      
	finalhistos.Add(Migr_rig);
        finalhistos.Add(Migr_beta);
       	finalhistos.Add(Migr_R);
        finalhistos.Add(Migr_B);
        finalhistos.writeObjsInFolder((directory+"/"+basename).c_str());	
}


void Acceptance::SaveResults(FileSaver finalhistos){
	FullSetEff  -> SaveResults(finalhistos);
	FullSetEff_gen  -> SaveResults(finalhistos);
	For_Acceptance  -> SaveResults(finalhistos);
 	if(EffAcceptance)  finalhistos.Add(EffAcceptance ); 
 	if(EffAcceptanceMC)finalhistos.Add(EffAcceptanceMC);
	
	if(Migr_rig) finalhistos.Add(Migr_rig);
	if(Migr_beta) finalhistos.Add(Migr_beta);
	if(Migr_R) finalhistos.Add(Migr_R);
	if(Migr_B) finalhistos.Add(Migr_B);
	
        finalhistos.writeObjsInFolder((directory+"/"+basename).c_str());	
}



void Acceptance:: EvalEffAcc(int timeindex,float SF){
	cout<<"********** MC flux reweighting **********"<<endl;
	Variables * vars = new Variables(timeindex);
	cout<<"********** For_Acceptance **************"<<endl;
	long int normalization = (For_Acceptance->GetBefore()->GetEntries()*20)/param.Trigrate;
	cout<<"Normalization Check:  "<<basename<<" "<<normalization<<" "<<param.tot_trig<<" "<<param.Trigrate<<endl;
	
	Histogram Spectrum = vars->reweighter.getTo();
	Histogram LogNorm  = vars->reweighter.getFrom();

	cout<<"********** Gen. Spectrum **********"<<endl;
	FullSetEff_gen->GetBefore()->Reset();
	//total triggers in range
	for(int i=0;i<bins.size();i++){
		float bincontent = normalization*(log(bins.RigTOIBins()[i+1])-log(bins.RigTOIBins()[i]))/log(param.Rmax)-log(param.Rmin);
		float meanweight = Spectrum.integrate(bins.RigTOIBins()[i],bins.RigTOIBins()[i+1])/LogNorm.integrate(bins.RigTOIBins()[i],bins.RigTOIBins()[i+1]);
		float MCweight = 1;
		//MCweight *= vars->GetCutoffCleaningWeight(bins.RigBins()[i],bins.RigTOIBins()[i],SF);
		MCweight *= vars->GetTimeDepWeight(bins.RigTOIBins()[i]);				      

	
		if(LogNorm.integrate(bins.RigTOIBins()[i],bins.RigTOIBins()[i+1])){
				FullSetEff_gen->GetBefore()->SetBinContent(i+1,bincontent);
				FullSetEff_gen->GetBefore()->SetBinError(i+1,0);
		}
	}
	
	cout<<"********** THE DIVISION **********"<<endl;
	FullSetEff_gen->Eval_Efficiency();
	EffAcceptance= (TH1F*)FullSetEff_gen->GetEfficiency()->Clone();
	EffAcceptanceMC= (TH1F*)FullSetEff_gen->GetEfficiency()->Clone();	
	EffAcceptance -> Sumw2();
	EffAcceptance -> Scale(47.78/param.gen_factor);
	EffAcceptanceMC -> Sumw2();
	EffAcceptanceMC -> Scale(47.78/param.gen_factor);

	EffAcceptance->SetName((basename +"_Eff_Acceptance").c_str());
	EffAcceptance->SetTitle((basename +"_Eff_Acceptance").c_str());
	EffAcceptanceMC->SetName((basename +"_Eff_AcceptanceMC").c_str());
	EffAcceptanceMC->SetTitle((basename +"_Eff_AcceptanceMC").c_str());
	cout<<"*********** Efficiency corrections ************"<<endl;
	for(int i=0;i<EfficiencyCorrections.size();i++) {
		for(int j=0; j<EffAcceptance->GetNbinsX(); j++){
			if(EfficiencyCorrections[i]->IsEkin()) 
				EffAcceptance -> SetBinContent( j+1, EffAcceptance -> GetBinContent(j+1)*EfficiencyCorrections[i]->GetCorrectionModel()->Eval(bins.ekpermassbincent[j]));
			else EffAcceptance -> SetBinContent( j+1, EffAcceptance -> GetBinContent(j+1)*EfficiencyCorrections[i]->GetCorrectionModel()->Eval(bins.rigbincent[j]) );
		}
	}	
	for(int j=0; j<EffAcceptance->GetNbinsX(); j++){
		float err =0;
		for(int i=0;i<EfficiencyCorrections.size();i++) {
			int bin=-1;
			bin = EfficiencyCorrections[i]->GetBins().GetRTOIBin(bins.RigBinCent(j));
			err+= pow(EfficiencyCorrections[i]->GetStat_Err()->GetBinContent(bin+1),2);
			err+= pow(EfficiencyCorrections[i]->GetSyst_Err()->GetBinContent(bin+1),2);
			err+= pow(EfficiencyCorrections[i]->GetSyst_Stat()->GetBinContent(bin+1),2);
		}	
		EffAcceptance -> SetBinError(j+1,pow(err,0.5)*EffAcceptance -> GetBinContent(j+1));
	}
	
}






