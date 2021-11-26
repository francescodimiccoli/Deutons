#include "EffCorr.h"


void EffCorr::Fill(TTree * treeMC,TTree * treeDT, Variables * vars, float (*discr_var) (Variables * vars),bool refill){

	EffMC      -> Fill(treeMC,vars,discr_var,refill);
	EffMCpid2     -> Fill(treeMC,vars,discr_var,refill);
	EffMCpid -> Fill(treeMC,vars,discr_var,refill);
	EffData    -> Fill(treeDT,vars,discr_var,refill);
}

void EffCorr::Save(){
	EffMC  -> Save();
	EffMCpid  -> Save();
	EffData-> Save();
	EffData_glob-> Save();

}

void EffCorr::Eval_Efficiencies(){

	
	if(!IsTrigEffCorr) {
		EffData-> Eval_Efficiency();
		EffData_glob-> Eval_Efficiency();
		EffMC     -> Eval_Efficiency();
		EffMCpid  -> Eval_Efficiency();
	}
	else{
		EffData	    -> Eval_TrigEfficiency();
		EffData_glob-> Eval_TrigEfficiency();

		EffMC     ->Eval_TrigEfficiency();
                EffMCpid  ->Eval_TrigEfficiency();

	}

}

void EffCorr::SaveResults(FileSaver finalhistos){
	EffMC  -> SaveResults(finalhistos);
	EffMCpid  -> SaveResults(finalhistos);
	EffData-> SaveResults(finalhistos);
	EffData_glob-> SaveResults(finalhistos);


	finalhistos.Add(GlobalEfficiency);
	finalhistos.Add(GlobalEfficiency_plot);
	finalhistos.Add(GlobalCorrection);
 	finalhistos.Add(GlobalCorrectionpid);
	if(GlobalCorrection_timeavg) finalhistos.Add(GlobalCorrection_timeavg);

	 finalhistos.Add(Stat_Err);
 	finalhistos.Add(Syst_Err);
	finalhistos.Add(Syst_Stat);
	finalhistos.Add(CorrectionModel);
	finalhistos.Add(DataEffModel);




       finalhistos.writeObjsInFolder((directory+"/"+basename).c_str());	
}

void AddCorrSystError(TH1F * Correction,TH1F * Syst_Err, TH1F * Syst_Stat) {
	for(int i=0; i<Correction->GetNbinsX(); i++) {
		float stat = Correction->GetBinError(i+1);
		float syst = Syst_Err->GetBinContent(i+1)*Correction->GetBinContent(i+1);
		float syst_stat = Syst_Stat->GetBinContent(i+1)*Correction->GetBinContent(i+1);
		Correction->SetBinError(i+1, pow(pow(stat,2)+pow(syst,2)+pow(syst_stat,2),0.5));
	}
}

void EffCorr::Eval_Errors(){

	Stat_Err = (TH1F *) GlobalCorrection->Clone();
	Syst_Err = (TH1F *) GlobalCorrection->Clone();
	Syst_Stat = (TH1F *) GlobalCorrection->Clone();
	Stat_Err->SetName((basename + "Stat_Err").c_str());
	Stat_Err->SetTitle((basename +"Stat_Err").c_str());
	Syst_Err->SetName((basename + "Syst_Err").c_str());
	Syst_Err->SetTitle((basename +"Syst_Err").c_str());
	Syst_Stat->SetName((basename + "Syst_Stat").c_str());
	Syst_Stat->SetTitle((basename +"Syst_Stat").c_str());


	TH1F * mceff = (TH1F*)EffMC->GetEfficiency();
	for(int i=0; i<GlobalCorrection->GetNbinsX(); i++) {
		cout<<"ECCO "<<EffMC->GetEfficiency()->GetNbinsX()<<" "<<Stat_Err<<" "<<Syst_Err<<" "<<Syst_Stat<<endl;	
		cout<<i<<endl;
		float stat = mceff->GetBinError(i+1); 
		float syst = 0.2*fabs(EffData_glob->GetEfficiency()->GetBinContent(i+1)-EffData_glob->GetEfficiency()->GetBinContent(i));   //fabs(GlobalCorrection->GetBinContent(i+1) - 1);
		if(GlobalCorrection->GetBinContent(i+1)>0){
			Stat_Err->SetBinContent(i+1, stat / GlobalCorrection->GetBinContent(i+1));
			Syst_Err->SetBinContent(i+1, syst / GlobalCorrection->GetBinContent(i+1));
			Syst_Stat->SetBinContent(i+1, 0);

			if(syst_stat>0){ 
				if(ForEffCorr.RigBinCent(i+1)<14)		
					Syst_Stat->SetBinContent(i+1,syst_stat->GetBinContent( syst_stat->FindBin(ForEffCorr.RigBinCent(i+1) )));
				else Syst_Stat->SetBinContent(i+1,syst_stat->GetBinContent(syst_stat->FindBin(14))); 	
			}

		}
		Stat_Err->SetBinError(i+1, 0);	
		Syst_Err->SetBinError(i+1, 0);	
		Syst_Stat->SetBinError(i+1, 0);	
		
	}
	cout<<"ECCO "<<EffMC<<" "<<Stat_Err<<" "<<Syst_Err<<" "<<Syst_Stat<<endl;	
	
	AddCorrSystError(GlobalCorrection,Syst_Err,Syst_Stat);
	AddCorrSystError(GlobalCorrectionpid,Syst_Err,Syst_Stat);
	
}

void EffCorr::Eval_Corrections(float shift){

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
	
	GlobalCorrectionpid=(TH1F *) EffData_glob->GetEfficiency()->Clone();
	GlobalCorrectionpid->SetName((basename + "_Corr_globpid").c_str());
	GlobalCorrectionpid->SetTitle((basename + "_Corr_globpid").c_str());

	for(int i=0; i<GlobalCorrection->GetNbinsX();i++) GlobalCorrection->SetBinError(i+1,0);
	for(int i=0; i<GlobalCorrection->GetNbinsX();i++) GlobalCorrectionpid->SetBinError(i+1,0);

	GlobalCorrection->Divide(EffMC->GetEfficiency());	
	GlobalCorrectionpid->Divide(EffMCpid->GetEfficiency());	

	if(avg_time){
		TH1F * timeavg;
		timeavg = ProjectionXtoTH1F(avg_time,(basename + "_Corr_avg").c_str(),avg_time->GetYaxis()->FindBin(time),avg_time->GetYaxis()->FindBin(time));
		if(IsEkinCorrection) timeavg = ConvertBinnedHisto(timeavg,"Global Correction;Ekin [GeV/n];Correction",Bins,true);
		else timeavg = ConvertBinnedHisto(timeavg,"Global Correction;Ekin [GeV/n];Correction",Bins,false);

		GlobalCorrection_timeavg = new TGraphErrors(timeavg);
		GlobalCorrection_timeavg->SetName((basename + "_Corr_avgtime").c_str());
	}

	for(int i=0; i<GlobalCorrection->GetNbinsX();i++) {
		if(GlobalCorrection->GetBinError(i+1)>=0.1 || GlobalEfficiency->GetBinContent(i+1)<=0.015){
			GlobalCorrection->SetBinContent(i+1,0);
			GlobalCorrectionpid->SetBinContent(i+1,0);
			GlobalCorrection->SetBinError(i+1,2);
			GlobalCorrectionpid->SetBinError(i+1,2);
	
		}

	}
	if(IsEkinCorrection){
		GlobalCorrection = ConvertBinnedHisto(GlobalCorrection,"Global Correction;Ekin [GeV/n];Correction",Bins,true);
		GlobalEfficiency_plot = ConvertBinnedHisto(GlobalEfficiency,"Data Efficiency;Ekin [GeV/n];Correction",Bins,true);
		}
	else{
		GlobalCorrection = ConvertBinnedHisto(GlobalCorrection,"Global Correction;R [GV];Efficiency",Bins,false);
		GlobalEfficiency_plot = ConvertBinnedHisto(GlobalEfficiency,"Data Efficiency;R [GV];Efficiency",Bins,false);
	}
	ModelWithSpline(shift);
	
	Eval_Errors();

	return;
}



void EffCorr::SetToConstantValue(float value){
	if(!GlobalCorrection) return;
	for(int i=0; i<GlobalCorrection->GetNbinsX(); i++) {
		GlobalCorrection->SetBinContent(i+1,value);
		GlobalCorrection->SetBinError(i+1,0.01*value);
		GlobalCorrectionpid->SetBinContent(i+1,value);
		GlobalCorrectionpid->SetBinError(i+1,0.01*value);
			
	}

}

void EffCorr::SetDefaultOutFile(FileSaver FinalHistos){
		EffMC      ->SetDefaultOutFile(FinalHistos);
		EffMCpid ->SetDefaultOutFile(FinalHistos);
		EffData	   ->SetDefaultOutFile(FinalHistos);
		EffData_glob   ->SetDefaultOutFile(FinalHistos);
		 
}

double FitFunc(double *x, double *p){
	float value;
	for(int i=0;i<24;i++){
		if(x[0]>=4*i&&x[0]<4*(i+1)) value = p[i];
	}
	return value;
}

double Model::Function(double *x, double *p){
	double X[nodes];
	double Y[nodes];
	for(int i=0;i<nodes;i++){
		X[i]=xs[i];
		Y[i]=p[i];
	}
	TSpline3 * spline = new TSpline3("spline",X,Y,nodes);
	return spline->Eval(x[0]); 
}

void EffCorr::ModelWithSpline(float shift){

	for(int i=0;i<GlobalCorrection->GetNbinsX();i++) if(((fabs(GlobalCorrection->GetBinContent(i+1)-GlobalCorrection->GetBinContent(i))>0.05&&fabs(GlobalCorrection->GetBinContent(i+1)-GlobalCorrection->GetBinContent(i+2))>0.05) || GlobalCorrection->GetBinError(i+1)>0.4)&&i>10) { GlobalCorrection->SetBinContent(i+1,(GlobalCorrection->GetBinContent(i) + GlobalCorrection->GetBinContent(i-1)) /2); GlobalCorrection->SetBinError(i+1,0.01); i++;}; 	

	int nodes = 11;//Bins.size()/4+1;
	double p[nodes];
	double eff[nodes];
	float fibonacci[11]={0,0.5, 1., 2., 3., 5., 8., 15., 23., 30.,37};
	int j=0;
	for(int i=0;i<nodes;i++){
		p[j]=GlobalCorrection->GetBinContent(GlobalCorrection->FindBin(shift+fibonacci[i]));
		eff[j]=EffData_glob->GetEfficiency()->GetBinContent(GlobalCorrection->FindBin(shift+fibonacci[i]));
		cout<<"NODINODI: "<<shift+fibonacci[i]<<" "<<p[j]<<endl;
		j++;
	}

	cout<<"******** "<<basename<<" ***********"<<endl;
	std::vector<float> spline_x;
	for(int i=0;i<nodes;i++) {
		spline_x.push_back(shift + fibonacci[i]);
		cout<<spline_x[i]<<endl;
	}
	Model *M = new Model(spline_x);

	CorrectionModel = new TF1((basename+"_Corr").c_str(),M,&Model::Function,0,50,nodes,"Correction_Model","Function");
 	DataEffModel = new TF1((basename+"_Eff").c_str(),M,&Model::Function,0,50,nodes,"DataEff_Model","Function");
 
	for(int i=0;i<nodes;i++) {
		CorrectionModel->SetParameter(i,p[i]);
		CorrectionModel->SetParLimits(i,0.8*p[i],1.2*p[i]);
		DataEffModel->SetParameter(i,eff[i]);
		DataEffModel->SetParLimits(i,0.8*eff[i],1.2*eff[i]);
	}
	float endrange=9999;
	for(int i=GlobalCorrection->GetNbinsX();i>0;i--) if(!(GlobalCorrection->GetBinContent(i+1)>0)) endrange = GlobalCorrection->GetBinLowEdge(i+1); else break;

	GlobalCorrection->Fit((basename+"_Corr").c_str(),"","M",shift,endrange-0.5); 
	GlobalEfficiency_plot->Fit((basename+"_Eff").c_str(),"","M",shift,endrange-0.5); 

}
