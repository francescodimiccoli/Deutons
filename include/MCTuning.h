int frac=2;

class Tuning{

	private:
	std::vector<std::vector<TH1F *>> DistribMC;
	std::vector<TH1F *> OriginalMC;
	std::vector<TH1F *> DistribDataPrimSec;
	std::vector<TH1F *> DistribData;
	std::vector<TH1F *> CutoffFilters;
	std::vector<std::vector<TH1F *>> Residuals;
	std::vector<TH2F *> ChiSquare;
	int steps;
	float sigma;
	float shift;
	Binning bins;
	TRandom3 * Rand;
	std::vector<int> bestChi2;
	TH1F * shiftbest;
	TH1F * sigmabest;
	TH1F * Exposure_Time;
	bool CutoffFilterMode=false;

	public:
	Tuning(Binning Bins, int Steps, float Sigma, float Shift,TH1F * ExposureTime=0x0){

		for(int bin=0;bin<Bins.size();bin++){
			DistribMC.push_back(std::vector<TH1F *>());
			Residuals.push_back(std::vector<TH1F *>());
			bestChi2.push_back(0);
			for(int i=0;i<Steps*Steps;i++){
				TH1F * Distrib  = new TH1F(("Distrib_"+to_string(bin)+"_"+to_string(i)).c_str(),("Distrib_"+to_string(bin)+"_"+to_string(i)).c_str(),300,0.4,1.5);
				TH1F * Residual = new TH1F(("Residual_"+to_string(bin)+"_"+to_string(i)).c_str(),("Residual_"+to_string(bin)+"_"+to_string(i)).c_str(),300,0.4,1.5);
				DistribMC[bin].push_back(Distrib);
				Residuals[bin].push_back(Residual);
			}	
			DistribData.push_back(new TH1F(("DistribData_"+to_string(bin)).c_str(),("Distrib_"+to_string(bin)).c_str(),300,0.4,1.5));
			DistribDataPrimSec.push_back(new TH1F(("DistribDataPrimSec_"+to_string(bin)).c_str(),("Distrib_"+to_string(bin)).c_str(),300,0.4,1.5));
			OriginalMC.push_back( new TH1F(("OriginalMC_"+to_string(bin)).c_str(),("OriginalMC_"+to_string(bin)).c_str(),300,0.4,1.5));
			ChiSquare.push_back(new TH2F(("ChiSquare_"+to_string(bin)).c_str(),("ChiSquare_"+to_string(bin)).c_str(),Steps,0,2*Sigma,Steps,-Shift,Shift));
		}
	
	        Exposure_Time=(TH1F*)ExposureTime;
		shiftbest = new TH1F("Shift Best","Shift Best",Bins.size(),0,Bins.size());
		sigmabest = new TH1F("Sigma Best","Sigma Best",Bins.size(),0,Bins.size());
	

		Rand= new TRandom3(time(0));
		steps=Steps;
		sigma=Sigma;
		shift=Shift;
		bins=Bins;
	}	

	Tuning(FileSaver File,Binning Bins, int Steps, float Sigma, float Shift,TH1F * ExposureTime=0x0){
		
		TFile * file = File.GetFile();
		for(int bin=0;bin<Bins.size();bin++){
			DistribMC.push_back(std::vector<TH1F *>());
			Residuals.push_back(std::vector<TH1F *>());
			bestChi2.push_back(0);
			for(int i=0;i<Steps*Steps;i++){
				TH1F * Distrib  = (TH1F*)  file->Get(("Bin " +to_string(bin) + "/MC Distributions/Distrib_"+to_string(bin)+"_"+to_string(i)).c_str());
				TH1F * Residual  = (TH1F*)  file->Get(("Bin " +to_string(bin) + "/Residuals/Residual_"+to_string(bin)+"_"+to_string(i)).c_str());
				DistribMC[bin].push_back(Distrib); 		
				Residuals[bin].push_back(Residual);
			}
			DistribData.push_back((TH1F*)  file->Get(("Bin " +to_string(bin) + "/Data Distribution/DistribData_"+to_string(bin)).c_str()));
			DistribDataPrimSec.push_back((TH1F*)  file->Get(("Bin " +to_string(bin) + "/Data Distribution/DistribDataPrimSec_"+to_string(bin)).c_str()));
			CutoffFilters.push_back((TH1F*)  file->Get(("Bin " +to_string(bin) + "/Data Distribution/CutoffFilter_"+to_string(bin)).c_str()));
			OriginalMC.push_back((TH1F*)  file->Get(("Bin " +to_string(bin) + "/Results/OriginalMC_"+to_string(bin)).c_str()));
			ChiSquare.push_back(new TH2F(("ChiSquare_"+to_string(bin)).c_str(),("ChiSquare_"+to_string(bin)).c_str(),Steps,0,2*Sigma,Steps,-Shift,Shift));
		}
		Exposure_Time=(TH1F*)  file->Get("Exposure Time");
		shiftbest = new TH1F("Shift Best","Shift Best",Bins.size(),0,Bins.size());
		sigmabest = new TH1F("Sigma Best","Sigma Best",Bins.size(),0,Bins.size());
	
		Rand= new TRandom3(time(0));
		steps=Steps;
		sigma=Sigma;
		shift=Shift;
		bins=Bins;
	}	
	
	void UseCutoffFilterMode() {CutoffFilterMode=true; return;};
	float SmearBeta(float Beta, float stepsigma, float stepshift);
	float SmearR(float R,float stepshift);
	void Fill(TNtuple * treeMC,TNtuple * treeDT, Variables * vars );
	void FillEventByEvent(std::vector<std::vector<TH1F *>> Histos, float var, float discr_var, bool CUT, float weight);
	void FillEventByEvent(std::vector<TH1F *> Histo, float var, float discr_var, bool CUT, float weight);	
	void Normalize();
	void EvalCutoffFilters();
	void Save(FileSaver finalhisto,bool recreate=false);
	void SaveResults(FileSaver finalhisto);
	void EvalResiduals();
	Binning  GetBinning() {return bins;};


	TH1F * GetOriginalMC(int bin){return OriginalMC[bin];}
	TH1F * GetDataPrim(int bin){return DistribData[bin];}
	TH1F * GetDataPrimSec(int bin){return DistribDataPrimSec[bin];}

	TH1F * GetBestChiSquare(int bin){return DistribMC[bin][bestChi2[bin]];}
	TH1F * GetBestChiSquareResiduals(int bin){return Residuals[bin][bestChi2[bin]];}
	TH1F * GetOriginalMCResiduals(int bin){return Residuals[bin][steps/2];}	
};

TSpline3 * ExtractCutoffWeight(TH1F * ExposureTime){

	ExposureTime->Scale(1/ExposureTime->GetBinContent(ExposureTime->GetMaximumBin()));
	
	double x[ExposureTime->GetNbinsX()];
	double y[ExposureTime->GetNbinsX()];

	for(int i=0;i<ExposureTime->GetNbinsX();i++){
		x[i]=ExposureTime->GetBinCenter(i+1);
		y[i]=ExposureTime->GetBinContent(i+1);
	}	
	
	TSpline3 * CutoffWeight = new TSpline3("CutoffWeight",x,y,ExposureTime->GetNbinsX());
	CutoffWeight->SetName("CutoffWeight");
	return CutoffWeight;

}


void Tuning::Normalize(){

	for(int bin=0;bin<bins.size();bin++){

		float integral=DistribData[bin]->Integral();
		DistribData[bin]->Sumw2();
		DistribData[bin]->Scale(1/integral);
	
		integral=DistribDataPrimSec[bin]->Integral();
		DistribDataPrimSec[bin]->Sumw2();
		DistribDataPrimSec[bin]->Scale(1/integral);

		integral=OriginalMC[bin]->Integral();
		OriginalMC[bin]->Sumw2();
		OriginalMC[bin]->Scale(1/integral);

		for(int i=0;i<DistribMC[bin].size();i++){
			float integral=DistribMC[bin][i]->Integral();
			DistribMC[bin][i]->Sumw2();
			DistribMC[bin][i]->Scale(1/integral);
		}
	}
	return;

}


void Tuning::Save(FileSaver finalhisto,bool recreate){
	
	for(int bin=0;bin<bins.size();bin++){
		finalhisto.Add(DistribData[bin]);
		finalhisto.Add(DistribDataPrimSec[bin]);
		finalhisto.Add(CutoffFilters[bin]);
		finalhisto.writeObjsInFolder(("Bin " + to_string(bin) + "/Data Distribution").c_str(),recreate);
		

		for(int i=0;i<DistribMC[bin].size();i++)
			finalhisto.Add(DistribMC[bin][i]);
		finalhisto.writeObjsInFolder(("Bin " + to_string(bin) + "/MC Distributions").c_str(),recreate);
	}
	return;
}

void Tuning::SaveResults(FileSaver finalhisto){

	for(int bin=0;bin<bins.size();bin++){

		for(int i=0;i<Residuals[bin].size();i++)
			finalhisto.Add(Residuals[bin][i]);
		finalhisto.writeObjsInFolder(("Bin " + to_string(bin) + "/Residuals").c_str());
		finalhisto.Add(ChiSquare[bin]);
		finalhisto.writeObjsInFolder(("Bin " + to_string(bin) + "/ChiSquare").c_str());

		finalhisto.Add(GetOriginalMC(bin));
		finalhisto.Add(GetDataPrim(bin));
		finalhisto.Add(GetDataPrimSec(bin));
		

		TH1F * Best = GetBestChiSquare(bin);
		Best->SetName("Best Fit MC");
		finalhisto.Add(Best);

		if(CutoffFilterMode){
			TH1F * Best_Filtered = (TH1F *)GetBestChiSquare(bin) ->Clone();
			Best_Filtered->Sumw2();
			Best_Filtered->Multiply(CutoffFilters[bin]);
			Best_Filtered->SetName("Best Fit MC (Cutoff filtered)");
			Best_Filtered->Scale(1/Best_Filtered->Integral());
			finalhisto.Add(Best_Filtered);
		}
		TH1F * ResidualOriginal = (TH1F *)GetOriginalMCResiduals(bin)->Clone();
		ResidualOriginal->SetName("ResidualOriginal");
		TH1F * ResidualBest = (TH1F *)GetBestChiSquareResiduals(bin)->Clone();
                ResidualBest->SetName("ResidualBest");

		finalhisto.Add(ResidualOriginal);
		finalhisto.Add(ResidualBest);
		finalhisto.writeObjsInFolder(("Bin " + to_string(bin) + "/Results").c_str());
		
		finalhisto.Add(shiftbest);
		finalhisto.Add(sigmabest);
		finalhisto.writeObjsInFolder("Best Fit Parameters");
	}
	return;	
}

void Tuning::EvalCutoffFilters(){
        for(int bin=0;bin<bins.size();bin++){
		TH1F * CutoffFilter=(TH1F *)DistribData[bin]->Clone();
		CutoffFilter->Sumw2();
		CutoffFilter->SetName(("CutoffFilter_"+to_string(bin)).c_str());
		CutoffFilter->Divide((TH1F*)DistribDataPrimSec[bin]->Clone());
		CutoffFilters.push_back(CutoffFilter);	
	}
}

void Tuning::Fill(TNtuple * treeMC,TNtuple * treeDT, Variables * vars ){

        vars->ReadAnalysisBranches(treeDT);

        for(int i=0;i<treeDT->GetEntries()/frac;i++){
                vars->AnalysisVariablseReset();
                UpdateProgressBar(i, treeDT->GetEntries()/frac);
                treeDT->GetEvent(i);
		FillEventByEvent( DistribData, vars->R, vars->Beta,ApplyCuts("ControlSample&LikelihoodCut&QualChargeCut&IsPrimary",vars),1);
       		FillEventByEvent( DistribDataPrimSec, vars->R, vars->Beta,ApplyCuts("ControlSample&LikelihoodCut&QualChargeCut",vars),1);
        
	 }

	EvalCutoffFilters();	

        vars->ReadAnalysisBranches(treeMC);

	TSpline3* W;
	float totalweight=vars->mcweight;
        if(Exposure_Time)   W = ExtractCutoffWeight(Exposure_Time);


        for(int i=0;i<treeMC->GetEntries()/frac;i++){
                vars->AnalysisVariablseReset();
                UpdateProgressBar(i, treeMC->GetEntries()/frac);
                treeMC->GetEvent(i);
                float totalweight=vars->mcweight;
		if(W&&(!CutoffFilterMode)) totalweight*=W->Eval(vars->R);
		if(IsProtonMC(vars)) {
			FillEventByEvent( DistribMC, vars->R, vars->Beta,ApplyCuts("ControlSample&LikelihoodCut&QualChargeCut",vars),totalweight);
        		FillEventByEvent( OriginalMC, vars->R, vars->Beta,ApplyCuts("ControlSample&LikelihoodCut&QualChargeCut",vars),vars->mcweight);
			}
	}

        return;
}

float Tuning::SmearR(float R,float stepshift){
        float B=0.8;
        float L=1.2;
        float sagitta=(0.0375*B/R)*pow(L,2)*1e3;
        sagitta = sagitta + (1/sagitta)*Rand->Gaus(0,(float)((2*shift/steps)*stepshift));
        return (0.0375*B/sagitta)*pow(L,2)*1e3;
}



float Tuning::SmearBeta(float Beta, float stepsigma, float stepshift){

	float time = 1.2/(Beta*3e-4);
	float shiftstart=-shift;
	time = time + (shiftstart+(2*shift/(float)steps)*stepshift) + Rand->Gaus(0,(float)((2*sigma/steps)*stepsigma));
	return 1.2/(time*3e-4);

}

void Tuning::FillEventByEvent(std::vector<std::vector<TH1F *>> Histos, float var, float beta, bool CUT, float weight){

	int kbin;
	for(int i=0;i<steps;i++)
		for(int j=0;j<steps;j++){	
			float betasmear= SmearBeta(beta,(float)i,(float)j);
			float Rsmear = var;//SmearR(var,(float)j);
			
			if(bins.IsUsingBetaEdges()) kbin =  bins.GetBin(betasmear);
			else kbin =  bins.GetBin(Rsmear);

			float mass = Rsmear/betasmear * pow((1-pow(betasmear,2)),0.5);
			if(/*mass>0.4&&mass<4.3&&*/CUT&&kbin>0)
				Histos[kbin][i*steps+j]->Fill(betasmear,weight);
			}
	return;
}

void Tuning::FillEventByEvent(std::vector<TH1F *> Histo, float var, float beta, bool CUT, float weight){
	
	int kbin;

	if(bins.IsUsingBetaEdges()) kbin =  bins.GetBin(beta);
	else kbin =  bins.GetBin(var);


	float mass=var/beta * pow((1-pow(beta,2)),0.5);
	if(CUT&&kbin>0/*&&mass>0.5&&mass<1.3*/){
		Histo[kbin]->Fill(beta,weight);
	}
	return;

}

void Tuning::EvalResiduals(){

	for(int bin=0;bin<bins.size();bin++){
		
		TH1F * Data;
		if(CutoffFilterMode) Data=(TH1F*)DistribDataPrimSec[bin]->Clone();
		else 		     Data=(TH1F*)DistribData[bin]->Clone(); 
	
		for(int i=0;i<Residuals[bin].size();i++){
			for(int j=0;j<DistribMC[bin][i]->GetNbinsX();j++){
				float sigma = pow((pow(DistribMC[bin][i]->GetBinError(j+1),2)+pow(Data->GetBinError(j+1),2)),0.5);
				if(sigma>0){
					Residuals[bin][i]->SetBinContent(j+1,(DistribMC[bin][i]->GetBinContent(j+1)-Data->GetBinContent(j+1))/sigma);
					Residuals[bin][i]->SetBinError(j+1,1);
				}
			}		}

		for(int i=0;i<steps;i++){
			for(int j=0;j<steps;j++){
				float chi=0;
				for(int b=0;b<Residuals[bin][i*steps+j]->GetNbinsX();b++) 
					{chi+=pow(Residuals[bin][i*steps+j]->GetBinContent(b+1),2);}
				ChiSquare[bin]->SetBinContent(i+1,j+1,pow(chi,0.5));
			}
		}
	
		float DoF=0;	
		for(int b=0;b<Residuals[bin][50]->GetNbinsX();b++)
                        if(Residuals[bin][50]->GetBinContent(b+1)>0) DoF=DoF+1;


		//chisquare minimization
                ChiSquare[bin]->Scale(1/DoF);
	
		float chimin=999999;
		for(int i=0;i<steps;i++){
			for(int j=0;j<steps;j++){
				float chi=ChiSquare[bin]->GetBinContent(i+1,j+1);
				if(pow(chi,0.5)<chimin)
					{ bestChi2[bin]=i*steps+j; shiftbest->SetBinContent(bin+1,(-shift+(2*shift/steps)*(j+1))); sigmabest->SetBinContent(bin+1,((2*sigma/steps)*(i+1))); chimin=pow(chi,0.5);}
			}
		}

		float BinSTD=0;
		for(int i=0;i<steps;i++){
                        for(int j=0;j<steps;j++){	
			BinSTD+=pow(ChiSquare[bin]->GetBinContent(i+1,j+1)-chimin,2);
			}
		}
		BinSTD=pow(BinSTD,0.5)/(steps*steps);
                shiftbest->SetBinError(bin+1,1.5/BinSTD);
		sigmabest->SetBinError(bin+1,1.5/BinSTD);
		//	
	
	}

	return;
}



	


