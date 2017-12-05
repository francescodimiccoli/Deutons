#ifndef FITFUNCTION_H
#define FITFUNCTION_H


using namespace std;

class FitFunction {

	private: 
	TH1 * RoughValues;
	TH1 * FittedValues;
	int n_smooth=0;
	int mc_types=0;

	public:
	FitFunction( TH1 * values, int nsmooth=1){
		RoughValues = (TH1 *) values -> Clone();
		FittedValues= (TH1 *) values -> Clone(); 	
	
		n_smooth = nsmooth;
		mc_types = RoughValues -> GetNbinsY();	
	}	
		

	void FitValues();
	TH1 * ReturnRoughValues() {return RoughValues;}
	TH1 * ReturnFittedValues() {return FittedValues;}


};


float FitError(TH1F * fit,TH1F * values,float bins,int param){

        float fitX2=0;
        for(int i=1;i<bins;i++) {
                if((fit->GetBinContent(i+1)-values->GetBinContent(i+1))/values->GetBinError(i+1)<10)
			fitX2+=pow((fit->GetBinContent(i+1)-values->GetBinContent(i+1))/values->GetBinError(i+1),2);
                }
        fitX2=(fitX2)/(bins-param);

        float sigmamean=0;
        for(int i=1;i<bins;i++)
                sigmamean+=pow(values->GetBinError(i+1),2);
        sigmamean=pow(sigmamean,0.5)/(bins-param);

        float errorefit=fitX2*sigmamean;
        return errorefit;
}




void FitFunction::FitValues(){

	if(mc_types==1){

		FittedValues -> Smooth(n_smooth);

		float fiterror = FitError((TH1F*)FittedValues,(TH1F*)RoughValues,FittedValues ->GetNbinsX(),0);	

		for(int nbins = 0; nbins < FittedValues ->GetNbinsX(); nbins ++){
			FittedValues -> SetBinError(nbins +1, fiterror);
		}
	}
	
	return;

}






TSpline3* ModelWithSpline(TH1F * Histo,std::string basename,Binning Bins){

        TSpline3* Model=new TSpline3 ();
        if(Histo->Integral()==0)        cout<<"******** ERROR: Exp. values seems not to be yet calculated: returning NULL **********"<<endl;
        else{
                int bins = Histo->GetNbinsX();
                FitFunction * Fit = new FitFunction(Histo,2);
                Fit->FitValues();
                TH1F * FittedValues = (TH1F *) Fit->ReturnFittedValues();
                int nbinsnotzero=0;
                for(int i=0;i<bins;i++) if(Histo->GetBinContent(i+1)>0) nbinsnotzero++;
                double x[nbinsnotzero]={0};
                double y[nbinsnotzero]={0};
                int i=0;
                for(int j =0;j<bins;j++){
                      if(Histo->GetBinContent(j+1)>0){
                        x[i]=Bins.GetBinCenter(j);
                        y[i]=FittedValues->GetBinContent(j+1);
                        if(i==0) y[i] = Histo->GetBinContent(j+1); // otherwise the first bin is smoothed with 0
			i++;
                        }
                }

                Model = new TSpline3((basename).c_str(),x,y,nbinsnotzero);
                Model -> SetName((basename).c_str());
        }
        return Model;
}


TF1* ModelWithPoly(TH1F * Histo,std::string basename,Binning Bins){

        TF1* Model=new TF1((basename).c_str(),"pol4");
        TGraphErrors *Graph=new TGraphErrors();
        int bins = Histo->GetNbinsX();
        if(Histo->Integral()==0)        cout<<"******** ERROR: Exp. values seems not to be yet calculated: returning NULL **********"<<endl;
        else{
                double x[bins]={0};
                double y[bins]={0};
                double x_err[bins]={0};
                double y_err[bins]={0};

                for(int j =0;j<bins;j++){
                        x[j]=Bins.GetBinCenter(j);
                        y[j]=Histo->GetBinContent(j+1);
                        y_err[j]=Histo->GetBinError(j+1);
                        Graph->SetPoint(j,x[j],y[j]);
                        Graph->SetPointError(j,x_err[j],y_err[j]);
                }
		int fit_attempt=0;
		while(Model->GetParameter(0)==0){
                	Graph -> Fit((basename).c_str(),"","",(float)x[fit_attempt],(float)x[bins-fit_attempt]);
			fit_attempt++;
			cout<<0.01*(float)fit_attempt<<endl;
		}
                Model -> SetName((basename).c_str());
        }
        return Model;
}
#endif



