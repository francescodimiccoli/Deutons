using namespace std;



TSpline3 * Model_Histo(TH1F * Histo);
TF1 * WeightGausswithSpline(TSpline3 * Spline, float mean,float mu, float sigma);
TH1F * InfinitesimalTermForConvolution(TH1F * Histo, TSpline3 * Spline, float mean, float mu, float sigma,int steps=500);
TH1F * ConvolveWithGaus(TH1F * Histo, float mu, float sigma,int steps=500);
TH1F * MultiplyWithGaus(TH1F * Original, float mean,float sigma);
float GetSigmatoConvolve(float distort_factor, TH1F * Original);
float GetSigmatoMultiply(float distort_factor, TH1F * Original);
TH1F * Distort_Histo(TH1F * Original,float distort_factor, float shift_factor=0);
std::vector< std::vector<TH1F *> > Multiple_Distortions(TH1F * Original,float distort_factor, float shift_factor, int steps=5);

			 



struct TFit {
   TH1F * Templ_P ;
   TH1F * Templ_D ;
   TH1F * Templ_He;
   TH1F * Data;
   float wheightP,wheightD,wheightHe;
   TFractionFitter *Tfit;
   int Tfit_outcome;
   
   float DCounts=0;
   float PCounts=0;
   float ChiSquare=0;			
   float StatErr=0;	

};



class TemplateFIT {

	private:
	std::vector<TFit *> fits;
	std::vector<TH1F *> TemplateP;
	std::vector<TH1F *> TemplateD;
	std::vector<TH1F *> Data;

	std::vector<TH2F *> DCountsSpread;
	std::vector<TH1F *> WeightedDCounts;
	std::vector<TH2F *> TFitChisquare;

	TH1F * StatError;
	TH1F * SystError;

	TH1F * ProtonCounts;
	TH1F * DeuteronCounts;

	Binning bins;
        std::string var;
        std::string cut;
	std::string discr_var;

	std::string basename;

	public:	
	//standard constructor
	TemplateFIT(std::string Basename,Binning Bins, std::string Cut, int Nbins, float Xmin, float Xmax){
		
		for(int i=0;i<Bins.size();i++){
				string name=Basename + "_Data_" + to_string(i);
                                TH1F * Histo = new TH1F((name).c_str(),(name).c_str(),Nbins,Xmin,Xmax);
                                Data.push_back(Histo);
                        }
		for(int i=0;i<Bins.size();i++){
				string name=Basename + "_MCP_" + to_string(i);
                                TH1F * Histo = new TH1F((name).c_str(),(name).c_str(),Nbins,Xmin,Xmax);
                                TemplateP.push_back(Histo);
                        }
		for(int i=0;i<Bins.size();i++){
				string name=Basename + "_MCD_" + to_string(i);
                                TH1F * Histo = new TH1F((name).c_str(),(name).c_str(),Nbins,Xmin,Xmax);
                                TemplateD.push_back(Histo);
                        }
	
		basename=Basename;
		cut = Cut;
		bins=Bins;
		
		StatError  = new TH1F("StatError","StatError",bins.size(),0,bins.size()) ;
        	SystError  = new TH1F("SystError","SystError",bins.size(),0,bins.size()) ;

		ProtonCounts    = new TH1F("Proton Counts","Proton Counts",bins.size(),0,bins.size()) ;
        	DeuteronCounts  = new TH1F("Deuteron Counts","Deuteron Counts",bins.size(),0,bins.size()) ;


	}

	//reading constructor
	
	TemplateFIT(FileSaver  File, std::string Basename,Binning Bins){
		std::string path = Basename + "/";
                        for(int i=0;i<Bins.size();i++){
                                TH1F * Histo = (TH1F*)  File.Get((path+ "Data/" + Basename  +"_Data_" + to_string(i)).c_str());
                                Data.push_back(Histo);
				cout<<Data[i]->GetEntries()<<endl;
                                }
			 for(int i=0;i<Bins.size();i++){
                                TH1F * Histo = (TH1F*)  File.Get((path+ "TemplateP/" + Basename  + "_MCP_" + to_string(i)).c_str());
                                TemplateP.push_back(Histo);
                                }
			 for(int i=0;i<Bins.size();i++){
                                TH1F * Histo = (TH1F*)  File.Get((path+ "TemplateD/" + Basename  + "_MCD_" + to_string(i)).c_str());
                                TemplateD.push_back(Histo);
			 }

			 basename=Basename;
			 bins=Bins; 

			 StatError  = new TH1F("StatError","StatError",bins.size(),0,bins.size()) ;
			 SystError  = new TH1F("SystError","SystError",bins.size(),0,bins.size()) ;

			ProtonCounts    = new TH1F("Proton Counts","Proton Counts",bins.size(),0,bins.size()) ;
        		DeuteronCounts  = new TH1F("Deuteron Counts","Deuteron Counts",bins.size(),0,bins.size()) ;


	}


	void Fill(TNtuple * treeMC,TNtuple * treeDT, Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars) );
	void FillEventByEvent(std::vector<TH1F *> Histos, float var, float discr_var, bool CUT,float weight=1);

	void Do_TemplateFIT(TFit * Fit);
	void ExtractCounts(FileSaver finalhisto);
	void CalculateFinalPDCounts();
	void Save(FileSaver finalhisto,bool recreate=false);	
	void SaveFitResults(FileSaver finalhisto);

	std::string GetName(){return basename;}
	TH1F * GetStatError(){ return StatError;}
	TH1F * GetSystError() { return SystError;}	
	Binning  GetBinning() {return bins;}


	TH2F * GetDCountsSpread(int bin)	{ return DCountsSpread[bin];}
	TH2F * GetChiSquareSpread(int bin)      { return TFitChisquare[bin];}	
	TH1F * GetWeightedDCounts(int bin)     { return WeightedDCounts[bin];}

};



void TemplateFIT::Fill(TNtuple * treeMC,TNtuple * treeDT, Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars) ){

	cout<<basename.c_str()<<" Filling ... (Data)"<< endl;
	vars->ReadAnalysisBranches(treeDT);

	for(int i=0;i<treeDT->GetEntries();i++){
		vars->AnalysisVariablseReset();		
		UpdateProgressBar(i, treeDT->GetEntries());
		treeDT->GetEvent(i);
		FillEventByEvent( Data, var(vars), discr_var(vars),ApplyCuts(cut,vars),vars->mcweight);
	}

	cout<<basename.c_str()<<" Filling ... (MC Protons)"<< endl;
	vars->ReadAnalysisBranches(treeMC);

	for(int i=0;i<treeMC->GetEntries();i++){
		vars->AnalysisVariablseReset();		
		UpdateProgressBar(i, treeMC->GetEntries());
		treeMC->GetEvent(i);
		if(IsProtonMC(vars)) FillEventByEvent( TemplateP, var(vars), discr_var(vars),ApplyCuts(cut,vars),vars->mcweight);
		if(IsDeutonMC(vars)) FillEventByEvent( TemplateD, var(vars), discr_var(vars),ApplyCuts(cut,vars),vars->mcweight);
	}

	return;
}


void TemplateFIT::FillEventByEvent(std::vector<TH1F *> Histos, float var, float discr_var, bool CUT, float weight){

	int kbin;
	kbin = 	bins.GetBin(discr_var);
	if(CUT&&kbin>0){ 
		Histos[kbin]->Fill(var,weight);		
	}
	return;	

}

void TemplateFIT::Save(FileSaver finalhisto,bool recreate){

	for(int i=0;i<bins.size();i++) 
		finalhisto.Add(Data[i]);
	finalhisto.writeObjsInFolder((basename + "/Data").c_str(),recreate);

	for(int i=0;i<bins.size();i++) 
		finalhisto.Add(TemplateP[i]);
	finalhisto.writeObjsInFolder((basename + "/TemplateP").c_str(),recreate);

	for(int i=0;i<bins.size();i++) 
		finalhisto.Add(TemplateD[i]);
	finalhisto.writeObjsInFolder((basename + "/TemplateD").c_str(),recreate);

	return;
}

void TemplateFIT::SaveFitResults(FileSaver finalhisto){
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

	finalhisto.Add(StatError);
	finalhisto.Add(SystError);
	finalhisto.Add(ProtonCounts);
	finalhisto.Add(DeuteronCounts);

	finalhisto.writeObjsInFolder((basename+"/Fit Results/").c_str());
}


void TemplateFIT::Do_TemplateFIT(TFit * Fit){

	TObjArray *Tpl;
	Tpl = new TObjArray(2);
	Tpl -> Add( Fit ->  Templ_P );
	Tpl -> Add( Fit ->  Templ_D );
	Fit -> Tfit = new TFractionFitter(Fit -> Data, Tpl ,"q");
	float highP=1;
	float lowP=0.5;
	float highD=0.25;
	float lowD=0.00001;
	float min=0.7;
	float max=2.5;
	Fit -> Tfit -> SetRangeX(Fit -> Data -> FindBin(min), Fit -> Data -> FindBin(max));
	
	bool fitcondition = (Fit -> Data->Integral()>5000)&&(Fit -> Templ_P->Integral()>1000) &&(Fit -> Templ_D->Integral()>500);

	if(fitcondition) {  
		Fit -> Tfit_outcome = Fit -> Tfit -> Fit();

		for(int fit_attempt=0; fit_attempt<20; fit_attempt++) {
			cout<<fit_attempt<<endl;
			if(Fit -> Tfit_outcome == 0) break;
			else {
				Fit -> Tfit -> Constrain(0, lowP+(float)fit_attempt/1000 ,highP );
				Fit -> Tfit -> Constrain(1, lowD-(float)fit_attempt/10000 ,highD-fit_attempt/10000 );
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

			Fit -> ChiSquare = Fit -> Tfit -> GetChisquare()/(float) Fit ->  Templ_P ->GetNbinsX();
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


TSpline3 * Model_Histo(TH1F * Histo){
	int nbins =Histo->GetNbinsX();
	double X[nbins];
	double Y[nbins];

	for(int i=0;i<Histo->GetNbinsX();i++){
		X[i]=Histo->GetBinCenter(i+1);
		Y[i]=Histo->GetBinContent(i+1);
	}

	TSpline3 * Model = new TSpline3("Model",X,Y,nbins);
	return Model;
}


TF1 * WeightGausswithSpline(TSpline3 * Spline, float mean,float mu, float sigma){
	TF1 * f1 = new TF1("f1","gaus",0,10);
	f1->SetParameters(Spline->Eval(mean),mean-mu,sigma);
	return f1;
}


TH1F * InfinitesimalTermForConvolution(TH1F * Histo, TSpline3 * Spline, float mean, float mu, float sigma,int steps){
        float min = Histo->GetXaxis()->GetBinLowEdge(1);
        float max = Histo->GetXaxis()->GetBinUpEdge(Histo->GetNbinsX());
        TH1F * Term = new TH1F("","",steps, min, max);
        TF1 * f1 = WeightGausswithSpline(Spline,mean,mu,sigma);
        for (int i=0;i<steps;i++){
                Term->SetBinContent(i+1,f1->Eval(Term->GetXaxis()->GetBinCenter(i)));
        }
        return Term;
}

TH1F * ConvolveWithGaus(TH1F * Histo, float mu, float sigma,int steps){

        float min = Histo->GetXaxis()->GetBinLowEdge(1);
        float max = Histo->GetXaxis()->GetBinUpEdge(Histo->GetNbinsX());
        int nbins =Histo->GetNbinsX();
        float original_area=Histo->Integral();

        TSpline3 * Model =  Model_Histo(Histo);
        TH1F * Term = InfinitesimalTermForConvolution(Histo,Model,Histo->GetBinCenter(1),mu,sigma,steps);

        for(int i=1;i<steps;i++)
                Term->Add(InfinitesimalTermForConvolution(Histo,Model,(Histo->GetBinCenter(i)),mu,sigma,steps));

        Term -> Rebin(steps/float(nbins));
        Term -> Scale(original_area/Term->Integral());
        return Term;
}

TH1F * MultiplyWithGaus(TH1F * Original, float mean,float sigma){

        TH1F * Histo = (TH1F*) Original->Clone();
        int nbins =Histo->GetNbinsX();
        float original_area=Histo->Integral();

        TF1 * f1 = new TF1("f1","gaus",0,10);
        f1->SetParameters(1,mean,sigma);

        for(int i=0;i<nbins;i++)
                Histo->SetBinContent(i+1,Histo->GetBinContent(i+1)*f1->Eval(Histo->GetBinCenter(i+1)));

        Histo -> Scale(original_area/Histo->Integral());
        return Histo;
}


float GetSigmatoConvolve(float distort_factor, TH1F * Original){

        float sigma_or = Original->GetStdDev();
        float sigma_f  = sigma_or+sigma_or*distort_factor;

        float sigma_c = pow(pow(sigma_f,2) - pow(sigma_or,2) , 0.5);

        return sigma_c;
}

float GetSigmatoMultiply(float distort_factor, TH1F * Original){

        float sigma_or = Original->GetStdDev();
        float sigma_f  = sigma_or-sigma_or*distort_factor;

        float sigma_c = 1/pow(1/pow(sigma_f,2) - pow(1/sigma_or,2) , 0.5);

        return sigma_c;
}




TH1F * Distort_Histo(TH1F * Original,float distort_factor, float shift_factor){

        TH1F * Distorted;

        if(distort_factor>0)        Distorted = ConvolveWithGaus(Original,0.0,GetSigmatoConvolve(distort_factor,Original));
        else  if(distort_factor<0)  Distorted = MultiplyWithGaus(Original,Original->GetMean(),GetSigmatoMultiply(fabs(distort_factor),Original));
        else                        Distorted = (TH1F*) Original ->Clone();

        TH1F * Distorted_Shifted = (TH1F*) Distorted ->Clone();
        if(fabs(shift_factor-0)>0.005) Distorted_Shifted  = ConvolveWithGaus(Distorted,shift_factor*Original->GetMean(),0.01);
        cout<<"shift: "<<shift_factor<<" sigma: "<<distort_factor<<endl;

        return Distorted_Shifted;

}




std::vector< std::vector<TH1F *> > Multiple_Distortions(TH1F * Original,float distort_factor, float shift_factor, int steps){

        std::vector<std::vector<TH1F *>> Collection;
        float distort_step = 2*distort_factor/float(steps-1);
        float shift_step = 2*shift_factor/float(steps-1);
        int c=0;
        int color_step = (50)/steps/steps;
        for(int i=0;i<steps;i++){
                Collection.push_back(std::vector<TH1F *>());
                for(int j=0;j<steps;j++){
                        TH1F * Distorted = Distort_Histo(Original,-distort_factor+i*distort_step,-shift_factor+j*shift_step);

                        string str1 = to_string ((-shift_factor+j*shift_step)*100);
                        string str2 = to_string ((-distort_factor+i*distort_step)*100);
                        str1.erase ( str1.find_last_not_of('0') + 1, std::string::npos );
                        str2.erase ( str2.find_last_not_of('0') + 1, std::string::npos );
                        Distorted->SetTitle(("Shift: " +str1 + "0% , Sigma: " + str2+"0%").c_str());
                        Distorted->SetLineWidth(2);
                        Collection[i].push_back(Distorted);
                        c+=color_step;
                }
        }
        return Collection;
}


void TemplateFIT::ExtractCounts(FileSaver finalhisto){



	for(int bin=0;bin<bins.size();bin++){
		std::vector< std::vector<TH1F *> > Collection = Multiple_Distortions((TH1F*)TemplateP[bin],0.07,0.02,5);
		cout<<"Bin: "<<bin<<endl;


		// Main FIT
		TFit * fit = new TFit;
		fit->Templ_P= (TH1F*)TemplateP[bin]->Clone();
		fit->Templ_D= (TH1F*)TemplateD[bin]->Clone();
		fit->Data   = (TH1F*)Data[bin]->Clone();
		Do_TemplateFIT(fit);

		fit->Templ_P->SetName( "Original Proton MC");
		fit->Templ_D->SetName( "Original Deuteron MC");
		TH1F * Result = (TH1F*)fit->Tfit->GetPlot();
		if(Result) Result->SetName("Fit Result");
		
		fits.push_back(fit);

		float temp = fit->DCounts;


		finalhisto.Add(fit->Templ_P);
		finalhisto.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesP/Bin"+to_string(bin)).c_str());
		finalhisto.Add(fit->Templ_D);
		finalhisto.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesD/Bin"+to_string(bin)).c_str());	
		if(Result){
		finalhisto.Add(Result);
		finalhisto.writeObjsInFolder((basename+"/Fit Results/FractionFit/Bin"+to_string(bin)).c_str());	
		}		

		// FITS for systematic error
		TH2F * dcountsspread = new TH2F(("DCountsSpread Bin " +to_string(bin)).c_str(),("DCountsSpread Bin " +to_string(bin)).c_str(),5,-7,7,5,-2,2);
        	TH2F * tfitchisquare = new TH2F(("ChiSquare Bin " +to_string(bin)).c_str(),("ChiSquare Bin " +to_string(bin)).c_str(),5,-7,7,5,-2,2);
		TH1F * weighteddcounts  = new TH1F(("Weighted Counts Bin " +to_string(bin)).c_str(),("Weighted Counts Bin " +to_string(bin)).c_str(),50,fits[bin]->DCounts - 0.5*fits[bin]->DCounts,fits[bin]->DCounts + 0.5*fits[bin]->DCounts);

		for(int sigma=0;sigma<Collection.size();sigma++){
			for(int shift=0;shift<Collection[sigma].size();shift++){

				TFit * fit = new TFit;	
				fit->Templ_P= (TH1F*)Collection[sigma][shift]->Clone();
				fit->Templ_D= (TH1F*)TemplateD[bin]->Clone();
				fit->Data   = (TH1F*)Data[bin]->Clone();
				Do_TemplateFIT(fit);

				fit->Templ_P->SetName( Collection[sigma][shift]->GetTitle());
				fit->Templ_D->SetName( Collection[sigma][shift]->GetTitle()); 	

				finalhisto.Add(fit->Templ_P);
				finalhisto.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesP/Bin"+to_string(bin)).c_str());
				finalhisto.Add(fit->Templ_D);
				finalhisto.writeObjsInFolder((basename+"/Fit Results/ScaledTemplatesD/Bin"+to_string(bin)).c_str());			

				
				dcountsspread->SetBinContent(sigma+1,shift+1,fit->DCounts);
				if(fit->ChiSquare>0)
						tfitchisquare->SetBinContent(sigma+1,shift+1,fit->ChiSquare);
				else  tfitchisquare->SetBinContent(sigma+1,shift+1,500);
				weighteddcounts -> Fill(dcountsspread->GetBinContent(sigma+1,shift+1),1/tfitchisquare->GetBinContent(sigma+1,shift+1));

			}

		}


		DCountsSpread.push_back(dcountsspread);
                TFitChisquare.push_back(tfitchisquare);
		WeightedDCounts.push_back(weighteddcounts);

		if(temp>0){
			StatError -> SetBinContent(bin+1,fits[bin]->StatErr);
			SystError -> SetBinContent(bin+1,WeightedDCounts[bin]->GetStdDev()/temp);
			StatError -> SetBinError(bin+1,0);
			SystError -> SetBinError(bin+1,0);
		}

	}

	CalculateFinalPDCounts();	
	return;
}

void TemplateFIT::CalculateFinalPDCounts(){
	for(int bin=0;bin<bins.size();bin++){
		
		float toterr= pow(pow(StatError -> GetBinContent(bin+1),2) + pow(SystError -> GetBinContent(bin+1),2) ,0.5);
		
		ProtonCounts->SetBinContent(bin+1,fits[bin]->PCounts);
		ProtonCounts->SetBinError(bin+1,toterr*fits[bin]->PCounts);

		DeuteronCounts->SetBinContent(bin+1,fits[bin]->DCounts);
                DeuteronCounts->SetBinError(bin+1,toterr*fits[bin]->DCounts);		
		
	}
		

}

