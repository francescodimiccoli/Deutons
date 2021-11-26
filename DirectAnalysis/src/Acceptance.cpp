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

		std::string namedir = filename;
		namedir.erase(namedir.find(".info"),namedir.find(".info")+5);

		std::istringstream iss(line);
		int start = line.find(namedir.c_str());
		int end =  line.find(".root");
		line.erase(0,start+namedir.size()+1);
		line.erase(line.find(".root"),line.find(".root")+5);
		stringstream geek(line);
		//cout<<"line: "<<line<<" "<<filename<<endl;
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


void Acceptance::Set_MCPar(float rmin, float rmax, float Gen_factor, std::string Filename, std::string Runlist, float compact_event_ratio, float Art_ratio){
	cout<<"Setting MC parameters"<<endl;
	param.Rmin=rmin;
	param.Rmax=rmax;
	param.filename = Filename;
	param.runlist = Runlist;
	param.gen_factor = Gen_factor;
	param.Eval_trigrate();
	param.art_ratio = Art_ratio;
	param.compact_event_ratio = compact_event_ratio;
}


void Acceptance::ApplyEfficCorr(EffCorr * Correction,float shift){
	if(Correction->GetGlobCorrection()) { 
		cout<<"Correction: "<<basename<<" "<<Correction->GetGlobCorrection()->GetEntries()<<endl;
		EfficiencyCorrections.push_back(Correction);
		EffCorrshift.push_back(shift);
	}
	return;
}

void Acceptance::ApplyEfficFromData(EffCorr * Correction,float shift){
	if(Correction->GetGlobCorrection()) { 
		cout<<"Correction: "<<basename<<" "<<Correction->GetGlobCorrection()->GetEntries()<<endl;
		EfficiencyFromData.push_back(Correction);
		EffDatashift.push_back(shift);
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
 	if(EffAcceptance_gen)  finalhistos.Add(EffAcceptance_gen ); 
 	if(EffAcceptanceMC_gen)finalhistos.Add(EffAcceptanceMC_gen);
	
	if(Migr_rig) finalhistos.Add(Migr_rig);
	if(Migr_beta) finalhistos.Add(Migr_beta);
	if(Migr_R) finalhistos.Add(Migr_R);
	if(Migr_B) finalhistos.Add(Migr_B);
	
        finalhistos.writeObjsInFolder((directory+"/"+basename).c_str());	
}


float GetAverage(TF1* Function, float xmin, float xmax) {

	float incr = (xmax-xmin)/10;
	float integral=0;

	for(int i=0;i<10;i++){
		integral+=Function->Eval(xmin+i*incr)*incr;
	}
	return integral/(xmax-xmin);
}

float GetAverage(TGraphErrors* Function, float xmin, float xmax) {

	float incr = (xmax-xmin)/10;
	float integral=0;

	for(int i=0;i<10;i++){
		integral+=Function->Eval(xmin+i*incr)*incr;
	}
	return integral/(xmax-xmin);
}


void Acceptance:: EvalEffAcc(int timeindex,float SF,bool IsHe){
	cout<<"********** MC flux reweighting **********"<<endl;
	Variables * vars = new Variables(timeindex);
	cout<<"********** For_Acceptance **************"<<endl;
	long int normalization = FRAC*(For_Acceptance->GetBefore()->GetEntries()/param.compact_event_ratio)/param.Trigrate;
	cout<<"Normalization Check:  "<<basename<<" "<<normalization<<" "<<param.tot_trig<<" "<<param.Trigrate<<endl;
	
	Histogram Spectrum;
	Histogram LogNorm;
	if(IsHe){  
		Spectrum = vars->reweighterHe.getTo();
		LogNorm  = vars->reweighter.getFrom();
	}
	else{  
		Spectrum = vars->reweighter.getTo();
		LogNorm  = vars->reweighter.getFrom();
	}
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
	FullSetEff->GetBefore()->Reset();
	//total triggers in range
	for(int i=0;i<bins.size();i++){
		float bincontent = normalization*(log(bins.RigTOIBins()[i+1])-log(bins.RigTOIBins()[i]))/log(param.Rmax)-log(param.Rmin);
		float meanweight = Spectrum.integrate(bins.RigTOIBins()[i],bins.RigTOIBins()[i+1])/LogNorm.integrate(bins.RigTOIBins()[i],bins.RigTOIBins()[i+1]);
		float MCweight = 1;
		//MCweight *= vars->GetCutoffCleaningWeight(bins.RigBins()[i],bins.RigTOIBins()[i],SF);
		MCweight *= vars->GetTimeDepWeight(bins.RigTOIBins()[i]);				      
	
		if(LogNorm.integrate(bins.RigTOIBins()[i],bins.RigTOIBins()[i+1])){
				FullSetEff->GetBefore()->SetBinContent(i+1,bincontent*meanweight*MCweight);
				FullSetEff->GetBefore()->SetBinError(i+1,0);
		}
	}
		
	cout<<"********** THE DIVISION **********"<<endl;
	FullSetEff_gen->Eval_Efficiency();
	FullSetEff->Eval_Efficiency();

	EffAcceptance= (TH1F*)FullSetEff->GetEfficiency()->Clone();
	EffAcceptanceMC= (TH1F*)FullSetEff->GetEfficiency()->Clone();	
	EffAcceptance -> Sumw2();
	EffAcceptance -> Scale(47.78/param.gen_factor);
	EffAcceptanceMC -> Sumw2();
	EffAcceptanceMC -> Scale(47.78/param.gen_factor);

	EffAcceptance->SetName((basename +"_Eff_Acceptance").c_str());
	EffAcceptance->SetTitle((basename +"_Eff_Acceptance").c_str());
	EffAcceptanceMC->SetName((basename +"_Eff_AcceptanceMC").c_str());
	EffAcceptanceMC->SetTitle((basename +"_Eff_AcceptanceMC").c_str());

	EffAcceptance_gen= (TH1F*)FullSetEff_gen->GetEfficiency()->Clone();
	EffAcceptanceMC_gen= (TH1F*)FullSetEff_gen->GetEfficiency()->Clone();	
	EffAcceptance_gen -> Sumw2();
	EffAcceptance_gen -> Scale(47.78/param.gen_factor);
	EffAcceptanceMC_gen -> Sumw2();
	EffAcceptanceMC_gen -> Scale(47.78/param.gen_factor);

	EffAcceptance_gen->SetName((basename +"_Eff_Acceptance_gen").c_str());
	EffAcceptance_gen->SetTitle((basename +"_Eff_Acceptance_gen").c_str());
	EffAcceptanceMC_gen->SetName((basename +"_Eff_AcceptanceMC_gen").c_str());
	EffAcceptanceMC_gen->SetTitle((basename +"_Eff_AcceptanceMC_gen").c_str());
	cout<<"*********** Efficiency corrections ************"<<endl;
	/*for(int i=0;i<EfficiencyCorrections.size();i++) {
		for(int j=0; j<EffAcceptance->GetNbinsX(); j++){
			if(EfficiencyCorrections[i]->IsEkin()) 
				EffAcceptance -> SetBinContent( j+1, EffAcceptance -> GetBinContent(j+1)*EfficiencyCorrections[i]->GetCorrectionModel()->Eval(bins.ekpermassbin_TOI[j]));
			else EffAcceptance -> SetBinContent( j+1, EffAcceptance -> GetBinContent(j+1)*EfficiencyCorrections[i]->GetCorrectionModel()->Eval(bins.rigbincent_TOI[j]) );
		}
	}*/	
	for(int i=0;i<EfficiencyCorrections.size();i++) {
		for(int j=0; j<EffAcceptance->GetNbinsX(); j++){
			if(EfficiencyCorrections[i]->IsEkin()){ 
				float correction;
				if(EfficiencyCorrections[i]->Get_TimeAvgCorrection())
					correction  = GetAverage(EfficiencyCorrections[i]->Get_TimeAvgCorrection(),bins.ekpermassbin_TOI[j],bins.ekpermassbin_TOI[j+1]);//EfficiencyCorrections[i]->Get_TimeAvgCorrection()->Eval(bins.ekpermassbincent_TOI[j] + EffCorrshift[i]); 
				else
					correction  = GetAverage(EfficiencyCorrections[i]->GetCorrectionModel(),bins.ekpermassbin_TOI[j],bins.ekpermassbin_TOI[j+1]);//EfficiencyCorrections[i]->GetCorrectionModel()->Eval(bins.ekpermassbincent_TOI[j] + EffCorrshift[i]);
				EffAcceptance -> SetBinContent( j+1, EffAcceptance -> GetBinContent(j+1)*correction);
				EffAcceptance_gen -> SetBinContent( j+1, EffAcceptance_gen -> GetBinContent(j+1)*correction);
				}
			else { 
				float correction;
				if(EfficiencyCorrections[i]->Get_TimeAvgCorrection())
					correction  = GetAverage(EfficiencyCorrections[i]->Get_TimeAvgCorrection(),bins.rigbin_TOI[j],bins.rigbin_TOI[j+1]);//EfficiencyCorrections[i]->Get_TimeAvgCorrection()->Eval(bins.rigbincent_TOI[j]);
				else	
					correction  = GetAverage(EfficiencyCorrections[i]->GetCorrectionModel(),bins.rigbin_TOI[j],bins.rigbin_TOI[j+1]);//EfficiencyCorrections[i]->GetCorrectionModel()->Eval(bins.rigbincent_TOI[j]);
				EffAcceptance -> SetBinContent( j+1, EffAcceptance -> GetBinContent(j+1)*correction);
				EffAcceptance_gen -> SetBinContent( j+1, EffAcceptance_gen -> GetBinContent(j+1)*correction);
			     }
		}
	}	

	for(int i=0;i<EfficiencyFromData.size();i++) {
		for(int j=0; j<EffAcceptance->GetNbinsX(); j++){
		cout<<"EFF FROM DATA: "<<EfficiencyFromData[i]->GetName()<<" "<<EfficiencyFromData[i]->IsEkin()<<" "<<bins.ekpermassbincent_TOI[j]<<" "<<EffDatashift[i] <<endl;
			if(EfficiencyFromData[i]->IsEkin()){ 
				float correction;
					correction= GetAverage(EfficiencyFromData[i]->GetDataEffModel(),bins.ekpermassbin_TOI[j],bins.ekpermassbin_TOI[j+1]); //EfficiencyFromData[i]->GetDataEffModel()->Eval(bins.ekpermassbincent_TOI[j]+EffDatashift[i]);
				if(EfficiencyCorrections[i]->Get_TimeAvgCorrection()) 
				//	correction*=EfficiencyCorrections[i]->Get_TimeAvgCorrection()->Eval(bins.ekpermassbincent_TOI[j] + EffCorrshift[i])/
				//	EfficiencyCorrections[i]->GetCorrectionModel()->Eval(bins.ekpermassbincent_TOI[j] + EffCorrshift[i]);	
					correction*=GetAverage(EfficiencyCorrections[i]->Get_TimeAvgCorrection(),bins.ekpermassbin_TOI[j],bins.ekpermassbin_TOI[j+1])/
						    GetAverage(EfficiencyCorrections[i]->GetCorrectionModel(),bins.ekpermassbin_TOI[j],bins.ekpermassbin_TOI[j+1]);
				EffAcceptance -> SetBinContent( j+1, EffAcceptance -> GetBinContent(j+1)*correction);
				EffAcceptance_gen -> SetBinContent( j+1, EffAcceptance_gen -> GetBinContent(j+1)*correction);
				}

			else {
				float correction;
					correction=GetAverage(EfficiencyFromData[i]->GetDataEffModel(),bins.rigbin_TOI[j],bins.rigbin_TOI[j+1]);//EfficiencyFromData[i]->GetDataEffModel()->Eval(bins.rigbincent_TOI[j]);
				if(EfficiencyCorrections[i]->Get_TimeAvgCorrection()) 
				//	correction*=EfficiencyCorrections[i]->Get_TimeAvgCorrection()->Eval(bins.rigbincent_TOI[j])/
				//	            EfficiencyCorrections[i]->GetCorrectionModel()->Eval(bins.rigbincent_TOI[j]);	
					correction*=GetAverage(EfficiencyCorrections[i]->Get_TimeAvgCorrection(),bins.rigbin_TOI[j],bins.rigbin_TOI[j+1])/
						    GetAverage(EfficiencyCorrections[i]->GetCorrectionModel(),bins.rigbin_TOI[j],bins.rigbin_TOI[j+1]);	
					EffAcceptance -> SetBinContent( j+1, EffAcceptance -> GetBinContent(j+1)*correction );
					EffAcceptance_gen -> SetBinContent( j+1, EffAcceptance_gen -> GetBinContent(j+1)*correction );
			     }		
		}
	}	



	for(int j=0; j<EffAcceptance->GetNbinsX(); j++){
		float err =0;
		for(int i=0;i<EfficiencyCorrections.size();i++) {
			int bin=-1;
			if(EfficiencyCorrections[i]->IsEkin()) bin = EfficiencyCorrections[i]->GetGlobCorrection()->FindBin(bins.ekpermassbin_TOI[j]);
			else bin = EfficiencyCorrections[i]->GetGlobCorrection()->FindBin(bins.rigbin_TOI[j]);
	
			if(EfficiencyCorrections[i]->GetStat_Err()->GetBinContent(bin+1)<0.1) err+= pow(EfficiencyCorrections[i]->GetStat_Err()->GetBinContent(bin+1),2);
			if(EfficiencyCorrections[i]->GetSyst_Err()->GetBinContent(bin+1)<0.1) err+= pow(EfficiencyCorrections[i]->GetSyst_Err()->GetBinContent(bin+1),2);
			if(EfficiencyCorrections[i]->GetSyst_Stat()->GetBinContent(bin+1)<0.1) err+= pow(EfficiencyCorrections[i]->GetSyst_Stat()->GetBinContent(bin+1),2);
		}	
		for(int i=0;i<EfficiencyFromData.size();i++) {
			int bin=-1;
			if(EfficiencyFromData[i]->IsEkin()) bin = EfficiencyFromData[i]->GetGlobCorrection()->FindBin(bins.ekpermassbin_TOI[j]);
			else bin = EfficiencyFromData[i]->GetGlobCorrection()->FindBin(bins.rigbin_TOI[j]);
				
			if(EfficiencyFromData[i]->GetStat_Err()->GetBinContent(bin+1)<0.1) err+= pow(EfficiencyFromData[i]->GetStat_Err()->GetBinContent(bin+1),2);
			if(EfficiencyFromData[i]->GetSyst_Err()->GetBinContent(bin+1)<0.1) err+= pow(EfficiencyFromData[i]->GetSyst_Err()->GetBinContent(bin+1),2);
			if(EfficiencyFromData[i]->GetSyst_Stat()->GetBinContent(bin+1)<0.1) err+= pow(EfficiencyFromData[i]->GetSyst_Stat()->GetBinContent(bin+1),2);
		}	
		err+=pow(0.02,2); //fragmentation syst error
		EffAcceptance -> SetBinError(j+1,pow(err,0.5)*EffAcceptance -> GetBinContent(j+1));
		EffAcceptance_gen -> SetBinError(j+1,pow(err,0.5)*EffAcceptance_gen -> GetBinContent(j+1));
	}
	
}






