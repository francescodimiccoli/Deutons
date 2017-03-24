class Resolution{

	private:
		std::vector<TH1F*> Histos;
		Binning bins;
		std::string var;
		std::string cut;
		std::string discr_var;
		std::string basename;
		int nbins;
		float xmin,xmax;
		TH1F * Sigmas_Histo;
		TH1F * Resolutions_Histo;
		TH1F * Means_Histo;

	public:
		//standard constructor
		Resolution(std::string Basename,Binning Bins, std::string Cut, std::string Var, std::string Discr_var, int Nbins, float Xmin, float Xmax){
			for(int i=0;i<Bins.size();i++){
				TH1F * Histo = new TH1F((Basename + to_string(i)).c_str(),(Basename + to_string(i)).c_str(),Nbins,Xmin,Xmax);
				Histos.push_back(Histo);	
			}
			basename=Basename;	
			var=Var;
			bins=Bins;
			xmin=Xmin;
			xmax=Xmax;
			nbins=Nbins;
			if(Cut.size()>0) cut=Cut; else cut="1>0";
			discr_var=Discr_var;	

			Sigmas_Histo      = new TH1F((Basename + "_sigmas").c_str(),(Basename + "_sigmas").c_str(),Bins.size(),0,Bins.size()); 		
			Resolutions_Histo = new TH1F((Basename + "_reso"  ).c_str(),(Basename + "_reso"  ).c_str(),Bins.size(),0,Bins.size());
			Means_Histo 	  = new TH1F((Basename + "_means" ).c_str(),(Basename + "_means" ).c_str(),Bins.size(),0,Bins.size());
		}

		//reading constructors
		Resolution(FileSaver File,std::string Basename, Binning Bins){
			std::string path = Basename + "/";
			for(int i=0;i<Bins.size();i++){
				TH1F * Histo = (TH1F*)  File.Get((path+Basename + to_string(i)).c_str());
				Histos.push_back(Histo);
				basename=Basename;
				bins=Bins;
			}
			nbins=Bins.size();

			Sigmas_Histo      = (TH1F*)  File.Get((path+Basename + "/Fit Results/" + Basename+ "_sigmas" ).c_str());
			Resolutions_Histo = (TH1F*)  File.Get((path+Basename + "/Fit Results/" + Basename+ "_sigmas" ).c_str());
			Means_Histo       = (TH1F*)  File.Get((path+Basename + "/Fit Results/" + Basename+ "_sigmas" ).c_str());	

			if(!Sigmas_Histo||!Resolutions_Histo||!Means_Histo){					
				Sigmas_Histo      = new TH1F((Basename + "_sigmas").c_str(),(Basename + "_sigmas").c_str(),Bins.size(),0,Bins.size());	
				Resolutions_Histo = new TH1F((Basename + "_reso"  ).c_str(),(Basename + "_reso"  ).c_str(),Bins.size(),0,Bins.size());
				Means_Histo 	  = new TH1F((Basename + "_means" ).c_str(),(Basename + "_means" ).c_str(),Bins.size(),0,Bins.size());
			}
		}

		void Fill(TTree * tree);	
		void Save(FileSaver finalhisto,bool recreate=false);
		void Normalize();
		void Eval_Resolution();
	bool CheckHistos();
	std::string GetName() {return basename;}
	Binning GetBinning(){return bins;}


	TH1F * Get_Resolutions() {return Resolutions_Histo;}
	TH1F * Get_Means() {return Means_Histo;}
	TH1F * Get_Sigmas() {return Sigmas_Histo;}

	TSpline3* ModelSigmaWithSpline();	
	TF1* ModelSigmaWithPoly();
};


void Resolution::Fill(TTree * tree){

	cout<<basename.c_str()<<" Filling ..."<< endl;
	
	std::string histo = "htemp("+to_string(nbins)+","+to_string(xmin)+","+to_string(xmax)+")";
	cout<<histo.c_str()<<endl;

	for(int i=0;i<bins.size();i++){
		cout<<(" Filling Bin " + to_string(i)).c_str()<< endl;
		std::string discr_cut = discr_var + ">" + to_string(bins.GetBinLowEdge(i)) + "&&" + discr_var + "<=" + to_string(bins.GetBinLowEdge(i+1));   
		cout<<(var + ">>" + histo).c_str()<<" " <<(cut + "&&" + discr_cut).c_str()<<endl;
		tree->Draw((var + ">>" + histo).c_str(), (cut + "&&" + discr_cut).c_str() );
	
		Histos[i] = (TH1F*) gDirectory->Get("htemp");
		Histos[i]->SetName((basename + to_string(i)).c_str());
	}
	
	return;
}


void Resolution::Save(FileSaver finalhisto,bool recreate){

	for(int i=0;i<bins.size();i++) 
		finalhisto.Add(Histos[i]);
	finalhisto.writeObjsInFolder(basename.c_str(),recreate);
	return;
}



void Resolution::Normalize(){

	for(int i=0;i<Histos.size();i++){
		float integral=Histos[i]->Integral();
		Histos[i]->Sumw2();
		Histos[i]->Scale(1/integral);
	}
	return;

}


void Resolution::Eval_Resolution(){

	for(int i=0;i<Histos.size();i++){
		double FitRangeEdges[2];
		double quantiles[2] = {0.025,0.8};
		Histos[i]->GetQuantiles(2,FitRangeEdges,quantiles);

		TF1 * fitfunc = new TF1("fitfunc","gaus",-1,1);
		fitfunc->SetParameter(0,Histos[i]->GetBinContent(Histos[i]->GetMaximumBin()));
		fitfunc->SetParameter(1,Histos[i]->GetMean());
		fitfunc->SetParameter(2,Histos[i]->GetRMS());

		Histos[i]->Fit("fitfunc","","",FitRangeEdges[0],FitRangeEdges[1]);
	
		float mean  	  = bins.GetBinCenter(i);
		float meanshift   = fitfunc->GetParameter(1)/(1/mean);
		float meanshifterr= fitfunc->GetParError(1)/(1/mean);
		
		float sigma 	  = pow(mean,2)*fitfunc->GetParameter(2);
		float sigmaerr    = pow(mean,2)*fitfunc->GetParError(2);	

		Sigmas_Histo->SetBinContent(i+1,sigma);
		Sigmas_Histo->SetBinError(i+1,sigmaerr);		

		Means_Histo->SetBinContent(i+1,meanshift);
		Means_Histo->SetBinError(i+1,meanshifterr);

		Resolutions_Histo->SetBinContent(i+1,sigma/mean);
		Resolutions_Histo->SetBinError(i+1,sigmaerr/mean);				
	}	
	return;
}

TSpline3* Resolution::ModelSigmaWithSpline(){
		
	TSpline3* Model;
	if(Sigmas_Histo->Integral()==0)	cout<<"******** ERROR: Sigma values seems not to be yet calculated: returning NULL **********"<<endl;
	else{
		FitFunction * Fit = new FitFunction(Sigmas_Histo,1);
		Fit->FitValues();
		TH1F * FittedValues = (TH1F *) Fit->ReturnFittedValues();	
		double x[nbins]={0};
		double y[nbins]={0};
		for(int j =0;j<nbins;j++){
			x[j]=bins.GetBinCenter(j);
			y[j]=FittedValues->GetBinContent(j+1);	
		}

		Model = new TSpline3((basename+"_Model").c_str(),x,y,nbins);
		Model -> SetName((basename+"_Model").c_str());
	}
	return Model;
}


TF1* Resolution::ModelSigmaWithPoly(){
		
	TF1* Model=new TF1((basename+"_Model").c_str(),"pol3");
	TGraphErrors *Graph=new TGraphErrors();
	if(Sigmas_Histo->Integral()==0)	cout<<"******** ERROR: Sigma values seems not to be yet calculated: returning NULL **********"<<endl;
	else{
		double x[nbins]={0};
		double y[nbins]={0};
		double x_err[nbins]={0};
                double y_err[nbins]={0};

		for(int j =0;j<nbins;j++){
			x[j]=bins.GetBinCenter(j);
			y[j]=Sigmas_Histo->GetBinContent(j+1);	
			y_err[j]=Sigmas_Histo->GetBinError(j+1);	
			Graph->SetPoint(j,x[j],y[j]);
			Graph->SetPointError(j,x_err[j],y_err[j]);
		}

		Graph -> Fit((basename+"_Model").c_str(),"","",(float)x[0],(float)x[nbins]);	
		Model -> SetName((basename+"_Model").c_str());
	}
	return Model;
}


bool Resolution::CheckHistos(){
	bool check=true;
	for(int i=0;i<bins.size();i++){
		if(!Histos[i])
			cout<<"**** BIN "<<i<<" : Histogram not found on File"<<endl;
		check=(check && Histos[i]);
	}		
	return check;
}
