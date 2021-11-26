#include "Unfolding.h"
#include "TVirtualFitter.h"

int slicenormalizex(TH2F * h){
        for(int i=0;i<h->GetNbinsX();i++){
		TH1D * proj = h->ProjectionY("proj",i+1,i+1);
		for(int j=0;j<h->GetNbinsY();j++){
		//	cout<<"proj: "<<proj->Integral()<<" "<<h->GetEntries()<<endl;
			if(proj->Integral()>0){ 
				h->SetBinContent(i+1,j+1, h->GetBinContent(i+1,j+1)/proj->Integral());
				h->SetBinError(i+1,j+1, h->GetBinError(i+1,j+1)/proj->Integral());
			}
			else{
				h->SetBinContent(i+1,j+1,0);
				h->SetBinError(i+1,j+1,0);
			}
		}

	}
	return 0;
}


double splineIntegral(TSpline3* spline,float xmin,float xmax) {
                double integral=0;
                float dx = (xmax-xmin)/10;
                for(int i=0;i<10;i++){
                        integral+=spline->Eval(xmin+i*dx)*dx;
                }
                return integral;
};

double Folder::foldwithmatrix(TSpline3* flux, float x){
	float result=0;
	int bin = Migr_matr->GetXaxis()->FindBin(x);
	float delta = Migr_matr->GetXaxis()->GetBinLowEdge(bin+1)-Migr_matr->GetXaxis()->GetBinLowEdge(bin);
	float expo_x = Expo->Eval(Migr_matr->GetXaxis()->GetBinCenter(bin));
	float acce_x = Acce->Eval(Migr_matr->GetXaxis()->GetBinCenter(bin));
	for(int i=0; i< Migr_matr->GetNbinsX();i++) {

			float expo_i = Expo->Eval(Migr_matr->GetXaxis()->GetBinCenter(i+1));
			float delta_i =  Migr_matr->GetXaxis()->GetBinLowEdge(i+2)-Migr_matr->GetXaxis()->GetBinLowEdge(i+1);
			result += (splineIntegral(flux,Migr_matr->GetXaxis()->GetBinLowEdge(i+1),Migr_matr->GetXaxis()->GetBinLowEdge(i+2)))*expo_x*acce_x*Migr_matr->GetBinContent(i+1,bin);
	}
	return result/(delta*expo_x*acce_x);
}


double Folder::Fold(double *x, double *p){
	double X[knots];
	double Y[knots];
	for(int i=0;i<knots;i++) {X[i]=0; X[i] = GetXknot(i);  }
	for(int i=0;i<knots;i++) {Y[i]=0; Y[i] = p[i]; }
	TSpline3 * foldedCounts = new TSpline3("foldedCounts",X,Y,knots);
	return foldwithmatrix(foldedCounts,x[0]);
}

float Folder::GetXknot(int indx){
	float logmin = log(0.1);
	float logmax = log(10);
	float binstep = (logmax-logmin)  / knots;
	float edges[knots];
	for(int i=0;i<knots;i++) edges[i]=(0.1*exp(i*binstep))/8;
	float min = Migr_matr->GetXaxis()->GetBinLowEdge(1)+offset;
	float max = Migr_matr->GetXaxis()->GetBinLowEdge(Migr_matr->GetNbinsX());
	//return min +(max-min)*edges[indx];
	float center;
	if(indx<(knots-1)) center = fit_min + offset + ((fit_max-fit_min)/(knots-1))*indx;
	else center = fit_max+ (max-fit_max)/6;

	return Migr_matr->GetXaxis()->GetBinCenter(Migr_matr->GetXaxis()->FindBin(center));
}

void Folder::SetKnots(double *y, std::string name) {
	double x[knots];
	for(int i=0;i<knots;i++) {x[i]=0; x[i] = GetXknot(i);  }
	UnfoldedCounts = new TSpline3(name.c_str(),x,y,knots);
}


UnfoldRes Unfold(TH2* migr_matr_norm, TH1 * measured_R, TGraph * expo, TGraph* acce, float fit_min,float fit_max,int knots,float offset){
         	UnfoldRes Res;
	        Folder * F = new Folder(migr_matr_norm,measured_R, expo, acce, fit_min,fit_max,knots,offset);
                int iter=100;
                TF1 * f[iter];

                float params[iter][F->GetNknots()];
                float chi2fit[iter];
                float value = 9999999999999;
                int min=0;
                for(int i=0;i<iter;i++) {chi2fit[i]=value;}

                int counter = 0;
                for(int att=0;att<iter;att++){
                        f[att]=new TF1(("f"+to_string(att)).c_str(),F,&Folder::Fold,0,20,F->GetNknots(),"Folder","Fold");
                        for(int i=0; i<F->GetNknots();i++){
                                if(F->GetXknot(i)<=measured_R->GetXaxis()->GetBinLowEdge(measured_R->GetNbinsX())){
		                      f[att]->SetParLimits(i,0.7*measured_R->GetBinContent(measured_R->GetXaxis()->FindBin(F->GetXknot(i))),
                                                       1.2*measured_R->GetBinContent(measured_R->GetXaxis()->FindBin(F->GetXknot(i))));}
                                else f[att]->SetParLimits(i,0.4*measured_R->GetBinContent(measured_R->GetNbinsX()),
                                                       0.9*measured_R->GetBinContent(measured_R->GetNbinsX()));
				}


                        for(int i=0; i<F->GetNknots();i++){
                                if(att==0) f[att]->SetParameter(i,measured_R->GetBinContent(measured_R->GetXaxis()->FindBin(F->GetXknot(i))));
                                if(att>0){
                                        f[att]->SetParameter(i,f[att-1]->GetParameter(i));
                                }
                        }

                        if(att>2){
                                if(fabs(chi2fit[att-2]-chi2fit[att-1])<0.1) counter++;
				else counter=0;
				//if(chi2fit[att-1]>1000) counter = 11;
                                if(counter>=10){
                                        for(int i=0; i<F->GetNknots();i++){
           						double limmin,limmax;
							f[att]->GetParLimits(i,limmin,limmax);		
		                                     f[att]->SetParameter(i,F->GetRandomGen()->Uniform(limmin,limmax));
					}
                                        counter=0;
                                }
                        }
			TVirtualFitter::Fitter(measured_R)->SetMaxIterations(100000);
                        measured_R->Fit(("f"+to_string(att)).c_str(),"NQWLM0","",fit_min,fit_max);
                        chi2fit[att] =f[att]->GetChisquare();
                        cout<<att<<" "<<chi2fit[att]<<" "<<counter<<" "<<knots<<endl;
        			                
        
                        if(chi2fit[att]<value){
                                value=chi2fit[att];
                                min = att;
                        }
                }
                cout<<"*********** FINAL FIT ********"<<endl;
                cout<<"minimum chi2: "<<min<<" attempt: "<<chi2fit[min]/((measured_R->GetNbinsX()-F->GetNknots()+1))<<endl; 
		measured_R->Fit(("f"+to_string(min)).c_str(),"WLM0","",fit_min,fit_max);
             	cout<<f[min]->GetChisquare()/(measured_R->GetNbinsX()-F->GetNknots())<<endl;         
                double y[F->GetNknots()];
        
                for(int i=0; i<F->GetNknots();i++){
                        y[i]=f[min]->GetParameter(i);
                }
                F->SetKnots(y,("Spline "+to_string(min)).c_str());
       		 
		Res.spline = F->GetSpline();
		Res.funct = f[min];
          
	        return Res;

}

