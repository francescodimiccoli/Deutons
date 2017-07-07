using namespace std;

// not used
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




struct TFit {
   TH1F * Templ_P ;
   TH1F * Templ_D ;

   TH1F * Data;
   TH1F * DataPrim;	


   float wheightP,wheightD,wheightHe;
   TFractionFitter *Tfit;
   int Tfit_outcome;
   
   float DCounts=0;
   float PCounts=0;
   float ChiSquare=0;			
   float StatErr=0;	

   TH1F * Templ_DPrim;
   float DCountsPrim=0;		

   TFit(){}	
   TFit(TH1F * templ_P, TH1F * templ_D, TH1F * data, TH1F * dataPrim) { Templ_P= templ_P; Templ_D=templ_D; Data=data; DataPrim=dataPrim; }

};


struct Systpar{
	int steps;
	float sigma;
	float shift;
};

struct BestChi {
	int i=0;
	int j=0;
	float chimin=0;
	void FindMinimum(TH2F * Histo) {
		float Best = 9999999;
		for(int x=0;x<Histo->GetNbinsX();x++)
			for(int y=0;y<Histo->GetNbinsY();y++)
				if(Histo->GetBinContent(x+1,y+1)<Best){
					Best=Histo->GetBinContent(x+1,y+1);
				        i=x;
					j=y;
					chimin=Best;	
				}
	}
};

class TemplateFIT {

	private:
	std::vector<std::vector<std::vector<TFit *>>> fits;
	std::vector<BestChi *> BestChiSquare;	
	std::vector<TH1F *> TransferFunction;


	std::vector<TH2F *> DCountsSpread;
	std::vector<TH1F *> WeightedDCounts;
	std::vector<TH2F *> TFitChisquare;

	TH1F * StatError;
	TH1F * SystError;
	TH1F * HeContError;

	TH1F * ProtonCounts;
	TH1F * DeuteronCounts;

	TH1F * ProtonCountsPrim;
	TH1F * DeuteronCountsPrim;


	TH1F * BestChiSquares;
	TH1F * OriginalChiSquares;

	TH1F * BestFitSigma;
	TH1F * BestFitShift;

	Binning bins;
        std::string var;
        std::string cut;
	std::string cutprimary;
	std::string discr_var;

	std::string basename;
	std::string hecontname;

	TH1F * Exposure_Time;

	Systpar systpar;
	float fitrangemin;
	float fitrangemax;
	bool fitDisabled=false;
	bool isrich=false;

	public:	
	//standard constructor
	TemplateFIT(std::string Basename,std::string HeContname,Binning Bins, std::string Cut, int Nbins, float Xmin, float Xmax, bool IsRich=false ,int steps=11,float sigma=50,float shift=40,TH1F * ExposureTime=0x0){
		
		for(int bin=0;bin<Bins.size();bin++){
			fits.push_back(std::vector<std::vector<TFit *>>());
			for(int i=0;i<steps;i++){
				fits[bin].push_back(std::vector<TFit *>());
				for(int j=0;j<steps;j++){

					TFit * fit = new TFit;
					string named    =Basename + "_Data_" +to_string(bin)+" "+to_string(i)+" "+to_string(j);
					string namedprim=Basename + "_DataPrim_" +to_string(bin)+" "+to_string(i)+" "+to_string(j);
					string nameP    =Basename + "_MCP_"      +to_string(bin)+" "+to_string(i)+" "+to_string(j);
					string nameD    =Basename + "_MCD_"      +to_string(bin)+" "+to_string(i)+" "+to_string(j);

					fit->Templ_P =  new TH1F(nameP.c_str(),nameP.c_str(),Nbins,Xmin,Xmax);
					fit->Templ_D =  new TH1F(nameD.c_str(),nameD.c_str(),Nbins,Xmin,Xmax);
					fit->Data    =  new TH1F(named.c_str(),named.c_str(),Nbins,Xmin,Xmax);
					fit->DataPrim=  new TH1F(namedprim.c_str(),namedprim.c_str(),Nbins,Xmin,Xmax);
					fits[bin][i].push_back(fit);
				}
			}
		}	
		basename=Basename;
		hecontname=HeContname;
		cut = Cut;
		cutprimary=Cut+"&IsPrimary";
		bins=Bins;
		
		StatError  = new TH1F("StatError","StatError",bins.size(),0,bins.size()) ;
        	SystError  = new TH1F("SystError","SystError",bins.size(),0,bins.size()) ;

		ProtonCounts    = new TH1F("Proton Counts","Proton Counts",bins.size(),0,bins.size()) ;
        	DeuteronCounts  = new TH1F("Deuteron Counts","Deuteron Counts",bins.size(),0,bins.size()) ;

		ProtonCountsPrim    = new TH1F("Primary Proton Counts","Primary Proton Counts",bins.size(),0,bins.size()) ;
        	DeuteronCountsPrim  = new TH1F("Primary Deuteron Counts","Primary Deuteron Counts",bins.size(),0,bins.size()) ;
	
	
		BestFitSigma  = new TH1F("Best FIt Sigma","Best FIt Sigma",bins.size(),0,bins.size()) ;
        	BestFitShift  = new TH1F("Best Fit Shift","Best FIt Shift",bins.size(),0,bins.size()) ;

		BestChiSquares     = new TH1F("Best ChiSquare","Best ChiSquare",bins.size(),0,bins.size()) ;
        	OriginalChiSquares = new TH1F("Original ChiSquare","Original CHiSquare",bins.size(),0,bins.size()) ;

		Exposure_Time=(TH1F*)ExposureTime;

		isrich = IsRich;

		systpar.sigma=sigma;
		systpar.shift=shift;
		systpar.steps=steps;

		fitrangemin=0.6;
		fitrangemax=3;	
	}

	//reading constructor

	TemplateFIT(FileSaver  File, std::string Basename,std::string HeContname,Binning Bins, bool IsRich=false, int steps=11,float sigma=50,float shift=40,TH1F * ExposureTime=0x0){

		TFile * file = File.GetFile();

		for(int bin=0;bin<Bins.size();bin++){
			fits.push_back(std::vector<std::vector<TFit *>>());
			for(int i=0;i<steps;i++){
				fits[bin].push_back(std::vector<TFit *>());
				for(int j=0;j<steps;j++){

					TFit * fit = new TFit;
					string named    =Basename + "/Bin "+ to_string(bin)+"/Data/" + Basename + "_Data_" +to_string(bin)+" "+to_string(0)+" "+to_string(5);
					string namedprim=Basename + "/Bin "+ to_string(bin)+"/Data/" + Basename + "_DataPrim_" +to_string(bin)+" "+to_string(0)+" "+to_string(5);
					string nameP    =Basename + "/Bin "+ to_string(bin)+"/TemplateP/" + Basename + "_MCP_"      +to_string(bin)+" "+to_string(i)+" "+to_string(j);
					string nameD    =Basename + "/Bin "+ to_string(bin)+"/TemplateD/" + Basename + "_MCD_"      +to_string(bin)+" "+to_string(i)+" "+to_string(j);

					fit->Templ_P =  (TH1F *)file->Get(nameP.c_str());
					fit->Templ_D =  (TH1F *)file->Get(nameD.c_str());
					fit->Data    =  (TH1F *)file->Get(named.c_str());
					fit->DataPrim=  (TH1F *)file->Get(namedprim.c_str());
					fits[bin][i].push_back(fit);
				}
			}
		}	



		basename=Basename;
		hecontname=HeContname;

		bins=Bins; 

		StatError  = new TH1F("StatError","StatError",bins.size(),0,bins.size()) ;
		SystError  = new TH1F("SystError","SystError",bins.size(),0,bins.size()) ;

		ProtonCounts    = new TH1F("Proton Counts","Proton Counts",bins.size(),0,bins.size()) ;
		DeuteronCounts  = new TH1F("Deuteron Counts","Deuteron Counts",bins.size(),0,bins.size()) ;
	
		ProtonCountsPrim    = new TH1F("Primary Proton Counts","Primary Proton Counts",bins.size(),0,bins.size()) ;
        	DeuteronCountsPrim  = new TH1F("Primary Deuteron Counts","Primary Deuteron Counts",bins.size(),0,bins.size()) ;
	

		BestFitSigma  = new TH1F("Best Fit Sigma","Best Fit Sigma",bins.size(),0,bins.size()) ;
		BestFitShift  = new TH1F("Best Fit Shift","Best Fit Shift",bins.size(),0,bins.size()) ;

		BestChiSquares     = new TH1F("Best ChiSquare","Best ChiSquare",bins.size(),0,bins.size()) ;
        	OriginalChiSquares = new TH1F("Original ChiSquare","Original CHiSquare",bins.size(),0,bins.size()) ;

		isrich=IsRich;

		Exposure_Time=(TH1F*)ExposureTime;

		systpar.sigma=sigma;
		systpar.shift=shift;
		systpar.steps=steps;

		fitrangemin=0.6;
		fitrangemax=3;	
	


	}

	void Eval_TransferFunction();
	void Fill(TNtuple * treeMC,TNtuple * treeDT, Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars) );

	float SmearBeta(float Beta, float stepsigma );
	float SmearBetaRICH(float Beta, float stepsigma);
	float SmearR(float R,float stepshift);

	void FillEventByEventData(float var, float discr_var, bool CUT, bool CUTPRIM, float weight);
	void FillEventByEventMC(float var, float discr_var, bool CUTP, bool CUTD, float weight);

	void Do_TemplateFIT(TFit * Fit);
	void ExtractCounts(FileSaver finalhisto,FileSaver finalResults);
	void EvalFinalParameters();
	void EvalFinalErrors(FileSaver finalResults);
	void CalculateFinalPDCounts();
	void Save(FileSaver finalhisto,bool recreate=false);	
	void SaveFitResults(FileSaver finalhisto);

	void SetSystematicParameters(int steps,float sigma,float shift){ systpar.steps=steps; systpar.shift=shift; systpar.sigma=sigma; return;};
	void SetFitRange(float min, float max){fitrangemin=min; fitrangemax=max; return;}

	void DisableFit(){fitDisabled=true;}

	std::string GetName(){return basename;}
	TH1F * GetStatError(){ return StatError;}
	TH1F * GetSystError() { return SystError;}	
	Binning  GetBinning() {return bins;}


	TH2F * GetDCountsSpread(int bin)	{ return DCountsSpread[bin];}
	TH2F * GetChiSquareSpread(int bin)      { return TFitChisquare[bin];}	
	TH1F * GetWeightedDCounts(int bin)     { return WeightedDCounts[bin];}

};

void TemplateFIT::Eval_TransferFunction(){
	for(int bin=0;bin<fits.size();bin++){
		TH1F * transferfunction = (TH1F *) fits[bin][0][0]->DataPrim->Clone();
		transferfunction->Sumw2();
		transferfunction->Divide(fits[bin][0][0]->Data);
		transferfunction->Smooth();
		TransferFunction.push_back(transferfunction);
	}
	return;
}

void TemplateFIT::Fill(TNtuple * treeMC,TNtuple * treeDT, Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars) ){

	cout<<basename.c_str()<<" Filling ... (Data)"<< endl;
	vars->ReadAnalysisBranches(treeDT);

	for(int i=0;i<treeDT->GetEntries();i++){
		vars->AnalysisVariablseReset();		
		UpdateProgressBar(i, treeDT->GetEntries());
		treeDT->GetEvent(i);
		FillEventByEventData(var(vars), discr_var(vars),ApplyCuts(cut,vars),ApplyCuts(cutprimary,vars),vars->mcweight);
	}
	
	cout<<basename.c_str()<<" Filling ... (MC Protons)"<< endl;
	vars->ReadAnalysisBranches(treeMC);

	for(int i=0;i<treeMC->GetEntries();i++){
		vars->AnalysisVariablseReset();		
		UpdateProgressBar(i, treeMC->GetEntries());
		treeMC->GetEvent(i);
		std::string cutP=cut+"&IsProtonMC";
		std::string cutD=cut+"&IsDeutonMC";
		FillEventByEventMC(vars->R, discr_var(vars),ApplyCuts(cutP,vars),ApplyCuts(cutD,vars),vars->mcweight);
	}

	return;
}


void TemplateFIT::FillEventByEventData(float var, float discr_var, bool CUT, bool CUTPRIM, float weight){

	int kbin;
	kbin = 	bins.GetBin(discr_var);
	if(CUT&&kbin>0){
		for(int i=0;i<systpar.steps;i++)
			for(int j=0;j<systpar.steps;j++){
				fits[kbin][i][j]->Data->Fill(var,weight);		
				if(CUTPRIM) fits[kbin][i][j]->DataPrim->Fill(var,weight);
				}
	}
	return;	

}


float TemplateFIT::SmearBetaRICH(float Beta, float stepsigma){
	float angle= acos(1/(1.15*Beta))*10e4;
	
	angle = angle + Rand->Gaus(0,(float)((2*systpar.sigma/systpar.steps)*stepsigma));
	return 1/(1.15*cos(angle/10e4));
}


float TemplateFIT::SmearBeta(float Beta, float stepsigma){

	float time = 1.2/(Beta*3e-4);
	time = time + Rand->Gaus(0,(float)((2*systpar.sigma/systpar.steps)*stepsigma));
	return 1.2/(time*3e-4);

}

float TemplateFIT::SmearR(float R,float stepshift){
	float B=0.8;
	float L=1.2;
	float sagitta=(0.0375*B/R)*pow(L,2)*1e3;
	sagitta = sagitta + (1/sagitta)*Rand->Gaus(0,(float)((2*systpar.shift/systpar.steps)*stepshift));
	return (0.0375*B/sagitta)*pow(L,2)*1e3;

}

void TemplateFIT::FillEventByEventMC(float var, float discr_var, bool CUTP, bool CUTD, float weight){

	if((CUTP||CUTD)){
		for(int i=0;i<systpar.steps;i++)
			for(int j=0;j<systpar.steps;j++){
			
				float betasmear;
				if(!isrich) betasmear = SmearBeta(discr_var,(float)i);
				else 	   betasmear = SmearBetaRICH(discr_var,(float)i);
				
				float rsmear= SmearR(var,(float)j);	
				int kbin =  bins.GetBin(betasmear);
				float mass = rsmear/betasmear * pow((1-pow(betasmear,2)),0.5);
				if(CUTP&&kbin>0) fits[kbin][i][j]->Templ_P->Fill(mass,weight);		
				if(CUTD&&kbin>0) fits[kbin][i][j]->Templ_D->Fill(mass,weight);
			
			}
	}
	return;	

}


void TemplateFIT::Save(FileSaver finalhisto,bool recreate){

	for(int bin=0;bin<bins.size();bin++){ 
		finalhisto.Add(fits[bin][0][5]->Data);
		finalhisto.Add(fits[bin][0][5]->DataPrim);

		finalhisto.writeObjsInFolder((basename + "/Bin "+ to_string(bin)+"/Data").c_str(),recreate);
	}

	for(int bin=0;bin<bins.size();bin++){ 
		for(int i=0;i<systpar.steps;i++)
                        for(int j=0;j<systpar.steps;j++){

				finalhisto.Add(fits[bin][i][j]->Templ_P);
		}
		finalhisto.writeObjsInFolder((basename + "/Bin "+ to_string(bin)+"/TemplateP").c_str(),recreate);
	}

	for(int bin=0;bin<bins.size();bin++){ 
		for(int i=0;i<systpar.steps;i++)
                        for(int j=0;j<systpar.steps;j++){

				finalhisto.Add(fits[bin][i][j]->Templ_D);
		}
		finalhisto.writeObjsInFolder((basename + "/Bin "+ to_string(bin)+"/TemplateD").c_str(),recreate);
	}



	return;
}

void TemplateFIT::SaveFitResults(FileSaver finalhisto){


	for(int bin=0;bin<bins.size();bin++){
	TH1F * OriginalP=(TH1F*)fits[bin][0][0]->Templ_P->Clone();
	TH1F * BestP=(TH1F*)fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->Templ_P->Clone();
	OriginalP->SetName("Original Proton MC ");
	BestP->SetName("Best #chi^{2} Mod. Proton MC ");
	finalhisto.Add(OriginalP);
	finalhisto.Add(BestP);
	finalhisto.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesP/Bin"+to_string(bin)).c_str());
	
	TH1F * OriginalD=(TH1F*)fits[bin][0][0]->Templ_D->Clone();
	TH1F * BestD=(TH1F*)fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->Templ_D->Clone();
	OriginalD->SetName("Original Deuton MC ");
	BestD->SetName("Best #chi^{2} Mod. Deuton MC ");
	finalhisto.Add(OriginalD);
	finalhisto.Add(BestD);
	finalhisto.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesD/Bin"+to_string(bin)).c_str());
	
	}

	for(int bin=0;bin<bins.size();bin++){ 
		finalhisto.Add(fits[bin][0][0]->Data);
		finalhisto.Add(fits[bin][0][0]->DataPrim);
	finalhisto.writeObjsInFolder((basename+"/Fit Results/Data/Bin"+to_string(bin)).c_str());	
	}
	
	for(int bin=0;bin<bins.size();bin++){
                for(int i=0;i<systpar.steps;i++)
                        for(int j=0;j<systpar.steps;j++){
				if(!(i==0&&i==5))finalhisto.Add(fits[bin][i][j]->Templ_P);
			}	
	finalhisto.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesP/Bin"+to_string(bin)).c_str());
	}

	for(int bin=0;bin<bins.size();bin++){
                for(int i=0;i<systpar.steps;i++)
                        for(int j=0;j<systpar.steps;j++){
				if(!(i==0&&i==5))finalhisto.Add(fits[bin][i][j]->Templ_D);
				}
	finalhisto.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesD/Bin"+to_string(bin)).c_str());
	}
	
	for(int bin=0;bin<bins.size();bin++){
		if(BestChiSquare[bin]&&fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->Tfit){
			TH1F * FIT=(TH1F*)fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->Tfit-> GetPlot();
			if(FIT){
			FIT->SetName(("Fraction Fit bin" + to_string(bin)).c_str());
			finalhisto.Add(FIT);
			finalhisto.writeObjsInFolder((basename+"/Fit Results/FractionFits/Bin"+to_string(bin)).c_str());
			}
		}
	}


	for(int i=0;i<bins.size();i++) 
                finalhisto.Add(TransferFunction[i]);
        finalhisto.writeObjsInFolder((basename + "/Fit Results/TrasnferFunctions/").c_str());

	for(int i=0;i<bins.size();i++){
		finalhisto.Add(DCountsSpread[i]);
		finalhisto.writeObjsInFolder((basename+"/Fit Results/Spreads/DCounts").c_str());
	}
	for(int i=0;i<bins.size();i++){
		finalhisto.Add(TFitChisquare[i]);
		finalhisto.writeObjsInFolder((basename+"/Fit Results/Spreads/ChiSquare").c_str());	
	}
	for(int i=0;i<bins.size();i++){
		finalhisto.Add(WeightedDCounts[i]);
		finalhisto.writeObjsInFolder((basename+"/Fit Results/Spreads/Weighted D counts").c_str());	
	}


	finalhisto.Add(BestFitSigma);
	finalhisto.Add(BestFitShift);
	finalhisto.Add(StatError);
	finalhisto.Add(SystError);
	finalhisto.Add(HeContError);
	finalhisto.Add(ProtonCounts);
	finalhisto.Add(DeuteronCounts);
	finalhisto.Add(ProtonCountsPrim);
	finalhisto.Add(DeuteronCountsPrim);

	finalhisto.Add(BestChiSquares   ); 
	finalhisto.Add(OriginalChiSquares);

	finalhisto.writeObjsInFolder((basename+"/Fit Results/").c_str());


}


void TemplateFIT::Do_TemplateFIT(TFit * Fit){

	TObjArray *Tpl;
	Tpl = new TObjArray(2);
	Tpl -> Add( Fit ->  Templ_P );
	Tpl -> Add( Fit ->  Templ_D );
	Fit -> Tfit = new TFractionFitter(Fit -> Data, Tpl ,"q");
	
	float min=fitrangemin;
	float max=fitrangemax;

	Fit -> Tfit -> SetRangeX(Fit -> Data -> FindBin(min), Fit -> Data -> FindBin(max));
	
	bool fitcondition = (Fit -> Data->Integral()>5000)&&(Fit -> Templ_P->Integral()>1000) &&(Fit -> Templ_D->Integral()>500);

	if(fitcondition) { 
	
		Fit -> Tfit -> Constrain(0, 0.0001 ,1 );
                Fit -> Tfit -> Constrain(1, 0.0001 ,1 );
			 
		if(Fit -> Tfit ) Fit -> Tfit_outcome = Fit -> Tfit -> Fit();

		for(int fit_attempt=0; fit_attempt<20; fit_attempt++) {
			cout<<fit_attempt<<endl;
			if(Fit -> Tfit_outcome == 0) break;
			else {
				cout<<fit_attempt<<endl;
				Fit -> Tfit_outcome = Fit -> Tfit -> Fit();
			}
		}

		if(Fit -> Tfit_outcome==0){
			TH1F * Result = (TH1F *) Fit-> Tfit -> GetPlot();
			float itot= Result->Integral();
			double w1,e1 = 0;
			double w2,e2 = 0;
			Fit -> Tfit ->GetResult(0,w1,e1);
			Fit -> Tfit ->GetResult(1,w2,e2);
			float i1 = Fit-> Templ_P  ->Integral(Fit->Templ_P -> FindBin(min), Fit->Templ_P -> FindBin(max));
			float i2 = Fit-> Templ_D  ->Integral(Fit->Templ_D -> FindBin(min), Fit->Templ_D -> FindBin(max));
			Fit ->wheightP= w1*itot/i1;
			Fit ->wheightD= w2*itot/i2;
			cout<<w1<<" "<<w2<<endl;

				Fit ->  Templ_P  -> Scale(Fit ->wheightP);
			Fit ->  Templ_D  -> Scale(Fit ->wheightD);

			float Cov01=0;

			Cov01= Fit -> Tfit->GetFitter()->GetCovarianceMatrixElement(0,1);
			Fit -> StatErr = pow(pow(w2*e2,2)+pow(w1*e1,2)-2*Cov01*w1*w2,0.5);

			Fit -> ChiSquare = Fit -> Tfit -> GetChisquare()/(float) (Fit ->  Templ_P ->GetNbinsX()*1.3);
			Fit -> DCounts = Fit ->  Templ_D -> Integral();
			Fit -> PCounts = Fit ->  Templ_P -> Integral();
		}
		else{
			Fit ->wheightP= 0;
			Fit ->wheightD= 0;
		}
	}
	return;
}

float EvalFitProbability(float chi){
	TF1 * CumulativeChiSquare=new TF1("Chi","exp(-x)*x^-0.5",0.001,10);
	return CumulativeChiSquare->Integral(chi,100);
}


void TemplateFIT::ExtractCounts(FileSaver finalhisto, FileSaver finalResults){
	
	Eval_TransferFunction();

	if(!fitDisabled){

	for(int bin=0;bin<bins.size();bin++){

		// Template Fits
		for(int sigma=0;sigma<systpar.steps;sigma++){
                        for(int shift=0;shift<systpar.steps;shift++){
                                cout<<fits[bin][sigma][shift]<<endl;
                                Do_TemplateFIT(fits[bin][sigma][shift]);

				fits[bin][sigma][shift]->Templ_DPrim=(TH1F*) fits[bin][sigma][shift]->Templ_D->Clone();
				fits[bin][sigma][shift]->Templ_DPrim->Multiply(TransferFunction[bin]);
				fits[bin][sigma][shift]->DCountsPrim = fits[bin][sigma][shift]->Templ_DPrim->Integral();

			}
		}

		// Histograms for systematic error evaluation
		TH2F * dcountsspread = new TH2F(("DCountsSpread Bin " +to_string(bin)).c_str(),("DCountsSpread Bin " +to_string(bin)).c_str(),systpar.steps,0,2*systpar.sigma,systpar.steps,-systpar.shift,systpar.shift);
        	TH2F * tfitchisquare = new TH2F(("ChiSquare Bin " +to_string(bin)).c_str(),("ChiSquare Bin " +to_string(bin)).c_str(),systpar.steps,0,2*systpar.sigma,systpar.steps,-systpar.shift,systpar.shift);
		TH1F * weighteddcounts  = new TH1F(("Weighted Counts Bin " +to_string(bin)).c_str(),("Weighted Counts Bin " +to_string(bin)).c_str(),35,0.4*fits[bin][0][0]->DCounts,1.7*fits[bin][0][0]->DCounts);
		BestChi * MinimumChi = new BestChi();
		

		for(int sigma=0;sigma<systpar.steps;sigma++){
			for(int shift=0;shift<systpar.steps;shift++){
				dcountsspread->SetBinContent(sigma+1,shift+1,fits[bin][sigma][shift]->DCounts);
				if(fits[bin][sigma][shift]->ChiSquare>0)
						tfitchisquare->SetBinContent(sigma+1,shift+1,fits[bin][sigma][shift]->ChiSquare);
				else  tfitchisquare->SetBinContent(sigma+1,shift+1,500);

				if(	fits[bin][sigma][shift]->DCounts>0.1*fits[bin][0][0]->DCounts &&
					fits[bin][sigma][shift]->DCounts<2*fits[bin][0][0]->DCounts  )

				weighteddcounts -> Fill(fits[bin][sigma][shift]->DCounts,EvalFitProbability(fits[bin][sigma][shift]->ChiSquare));
			}
		}

		MinimumChi->FindMinimum(tfitchisquare);	
		

		BestChiSquare.push_back(MinimumChi);			
		DCountsSpread.push_back(dcountsspread);
		TFitChisquare.push_back(tfitchisquare);
		WeightedDCounts.push_back(weighteddcounts);


	}

	EvalFinalParameters();
	EvalFinalErrors(finalResults);
	CalculateFinalPDCounts();	
	}
	
	return;
}

void TemplateFIT::EvalFinalParameters(){
	for(int bin=0;bin<bins.size();bin++){
		float BinSTD=0;
		for(int i=0;i<systpar.steps;i++){
			for(int j=0;j<systpar.steps;j++){
				BinSTD+=pow(TFitChisquare[bin]->GetBinContent(i+1,j+1)-BestChiSquare[bin]->chimin,2);
			}
		}
		BinSTD=pow(BinSTD,0.5)/(systpar.steps*systpar.steps);

		BestFitSigma->SetBinContent(bin+1,(float)((2*systpar.sigma/systpar.steps)*BestChiSquare[bin]->i));
		BestFitShift->SetBinContent(bin+1,-systpar.shift+(2*systpar.shift/(float)systpar.steps)*BestChiSquare[bin]->j);

		
	
		if(BinSTD>0){
		BestFitSigma->SetBinError(bin+1,pow(pow(4/BinSTD,2)+pow(1.5*(2*systpar.sigma/(float)systpar.steps),2),0.5));
		BestFitShift->SetBinError(bin+1,pow(pow(4/BinSTD,2)+pow(1.5*(2*systpar.shift/(float)systpar.steps),2),0.5));
		}

		if(BestChiSquare[bin]->chimin<100) 	       BestChiSquares     ->SetBinContent(bin+1,BestChiSquare[bin]->chimin);
                if(TFitChisquare[bin]->GetBinContent(1,6)<100) OriginalChiSquares ->SetBinContent(bin+1,TFitChisquare[bin]->GetBinContent(1,6));
		BestChiSquares     ->SetBinError(bin+1,0.25);	
                OriginalChiSquares ->SetBinError(bin+1,0.25);
	}
	return;
}

void TemplateFIT::EvalFinalErrors(FileSaver finalResults){

	if(finalResults.GetFile())
		HeContError = (TH1F *) finalResults.Get(("HeliumFragmentation/"+hecontname+"/"+hecontname+"_Eff").c_str());

	for(int bin=0;bin<bins.size();bin++){

		if(fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->DCounts>0){
			StatError -> SetBinContent(bin+1,fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->StatErr);
			SystError -> SetBinContent(bin+1,WeightedDCounts[bin]->GetStdDev()/fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->DCounts);
			StatError -> SetBinError(bin+1,0);
			SystError -> SetBinError(bin+1,0);
		}

	}
	return;

}


void TemplateFIT::CalculateFinalPDCounts(){
	
	if(!HeContError) HeContError = new TH1F("","",bins.size(),0,bins.size());
	
		HeContError-> Smooth(1);
		StatError -> Smooth(1);
		SystError -> Smooth(1);
		
	for(int bin=0;bin<bins.size();bin++){
		
		float toterr= pow(pow(StatError -> GetBinContent(bin+1),2) + pow(SystError -> GetBinContent(bin+1),2) + pow(2*HeContError -> GetBinContent(bin+1),2) ,0.5);
		
		ProtonCounts->SetBinContent(bin+1,fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->PCounts);
		ProtonCounts->SetBinError(bin+1,toterr*fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->PCounts);

		DeuteronCounts->SetBinContent(bin+1,fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->DCounts);
                DeuteronCounts->SetBinError(bin+1,toterr*fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->DCounts);		
		
		ProtonCountsPrim->SetBinContent(bin+1,fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->PCounts);
		ProtonCountsPrim->SetBinError(bin+1,toterr*fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->PCounts);

		DeuteronCountsPrim->SetBinContent(bin+1,fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->DCountsPrim);
                DeuteronCountsPrim->SetBinError(bin+1,toterr*fits[bin][BestChiSquare[bin]->i][BestChiSquare[bin]->j]->DCountsPrim);		
		
		


	}
		

}

