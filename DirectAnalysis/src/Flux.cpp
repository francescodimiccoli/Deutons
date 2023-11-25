#include "Flux.h"
#include "FluxRebin.h"
#include "histUtils.h"

#include "TRandom.h"
#include "TH1D.h"
#include "TCanvas.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "TVirtualFitter.h"

TH1F * Twentisevenify_(TH1F* input){

	std::string name = input->GetName();
	TH1F * output = (TH1F*)input->Clone((name+"_27").c_str());
	for(int i=0;i<output->GetNbinsX();i++){
		output->SetBinContent(i+1,input->GetBinContent(i+1)*pow(input->GetBinCenter(i+1),2.7));
		output->SetBinError(i+1,input->GetBinError(i+1)*pow(input->GetBinCenter(i+1),2.7));
	}

	return output;
}

void Regularize_Counts(TH1 * flux){
	TF1* trend = new TF1("trend","[0]*x^[1]",0.1,25);
	flux->Fit("trend","","",0,10);
	for(int i=0;i<flux->GetNbinsX();i++){
		flux->SetBinContent(i+1,trend->Eval(flux->GetBinCenter(i+1)));
		flux->SetBinError(i+1,flux->GetBinError(i+1)/5.);
	}
	return;
}

void Flux::AverageCountsWithOther(FileSaver FileRes, std::string namecounts, std::string nameexposure,float toll){
	TFile * fileres = FileRes.GetFile();

	Counts_ForAverage = (TH1F *) fileres->Get((namecounts).c_str())->Clone("Counts_ForAverage")  ;
	Exposure_ForAverage  = (TH1F *) fileres->Get((nameexposure).c_str())->Clone("Exposure_ForAverage");

	Counts_ForAverage->Divide(Exposure_ForAverage);


	for(int i=0;i<Counts_ForAverage->GetNbinsX();i++){
		Counts_ForAverage->SetBinContent(i+1,Counts_ForAverage->GetBinContent(i+1)/(Counts_ForAverage->GetBinLowEdge(i+2)-Counts_ForAverage->GetBinLowEdge(i+1)));
		Counts_ForAverage->SetBinError(i+1,Counts_ForAverage->GetBinError(i+1)/(Counts_ForAverage->GetBinLowEdge(i+2)-Counts_ForAverage->GetBinLowEdge(i+1)));
	}

	
	TH1F * This 	= ConvertBinnedHisto(Counts,"This",bins,true);
	TH1F * Exposure_this     =  ConvertBinnedHisto(ExposureTime,"Exposure_this",bins,true);

	This->Sumw2();
        This -> Divide(Exposure_this);
	
	for(int i=0;i<This->GetNbinsX();i++){	
		This->SetBinError(i+1,This->GetBinError(i+1)/(bins.EkPerMasTOIBins()[i+1]-bins.EkPerMasTOIBins()[i]));
		This->SetBinContent(i+1,This->GetBinContent(i+1)/(bins.EkPerMasTOIBins()[i+1]-bins.EkPerMasTOIBins()[i]));
	}

	TGraphErrors * RatioOther_g= new TGraphErrors(Counts_ForAverage);

	for(int i=0;i<This->GetNbinsX();i++){
		if(This->GetBinCenter(i+1)>Counts_ForAverage->GetBinLowEdge(1)){
			float ave=(This->GetBinContent(i+1) + RatioOther_g->Eval(This->GetBinCenter(i+1)))/2;
			ave/=This->GetBinContent(i+1) ;
			if(ave>1-toll&&ave<1+toll) Counts->SetBinContent(i+1,ave*Counts->GetBinContent(i+1));	
		}
	}
	return;

}

void Flux::Unfold_Counts(float fit_min, float fit_max,int knots,float offset,bool regularize){
	if(fit_min==0 && fit_max==0) {
		Unfolding_factor = (TH1F*) Counts_density_unfolded->Clone("Unfolding_factor");
        	Unfolding_factor->Reset();
		return;
	}

	std::vector<TH1D *> random_meas;
        std::vector<TSpline3 *> spline_random;
        std::vector<TF1 *> funct_random;
        int iter =3;
        for(int i=0; i<iter;i++) {
 	
		Counts_density->Multiply(ExposureTime_plot_rig);
               	Counts_density->Multiply(Eff_Acceptance_plot_rig);
	
		TH1D * tmp =(TH1D*) Counts_density->Clone("");
		if(i!=0){
			tmp->Reset();
			tmp->FillRandom(Counts_density,Counts_density->Integral());
		}	

		for(int j=0;j<Counts_density->GetNbinsX();j++) tmp->SetBinError(j+1,Counts_density->GetBinError(j+1));

		Counts_density->Divide(ExposureTime_plot_rig);
		Counts_density->Divide(Eff_Acceptance_plot_rig);
		tmp->Divide(ExposureTime_plot_rig);
		tmp->Divide(Eff_Acceptance_plot_rig);
	/*	if(regularize) Regularize_Counts(tmp);
		else{
			tmp->Smooth();
			tmp->Smooth();
			tmp->Smooth();
		}*/
		random_meas.push_back(tmp);
	}

	TGraph * expo = new TGraph(ExposureTime_plot_rig);
	TGraph * acce = new TGraph(Eff_Acceptance_plot_rig);
	
	for(int i=0;i<iter;i++) {
		UnfoldRes Res;// = Unfold(migr_matr,random_meas[i], expo,acce, fit_min, fit_max,knots,offset);	
		spline_random.push_back(Res.spline);
		funct_random.push_back(Res.funct);
	}

	Unfolding_factor = (TH1F*) random_meas[0]->Clone("Unfolding_factor");
	Unfolding_factor->Sumw2();
	Unfolding_factor_raw = (TH1F*) random_meas[0]->Clone("Unfolding_factor_raw");
	Unfolding_factor_raw->Sumw2();
	
        /*for(int i=0;i< Unfolding_factor->GetNbinsX();i++) {
                float mean=0;
                for(int j=0;j<spline_random.size();j++) {
                        mean+=splineIntegral(spline_random[j],Unfolding_factor->GetBinLowEdge(i+1),Unfolding_factor->GetBinLowEdge(i+2))/(bins.RigTOIBins()[i+1]-bins.RigTOIBins()[i]);
                }
                float devst=0;
                for(int j=0;j<spline_random.size();j++) {
                        float value = splineIntegral(spline_random[j],Unfolding_factor->GetBinLowEdge(i+1),Unfolding_factor->GetBinLowEdge(i+2))/(bins.RigTOIBins()[i+1]-bins.RigTOIBins()[i]) ;
                        devst+=pow(value - mean/spline_random.size(),2);
                }
              //  devst=pow(devst/spline_random.size(),0.5);
		if(Counts_density->GetBinContent(i+1)>0){
			Unfolding_factor_raw->SetBinContent(i+1, (mean/spline_random.size()) /Unfolding_factor->GetBinContent(i+1));
	                Unfolding_factor_raw->SetBinError(i+1,(sqrt(devst)/spline_random.size())/Unfolding_factor->GetBinContent(i+1)   );

		}
		else{
                        Unfolding_factor_raw->SetBinContent(i+1,0);
                        Unfolding_factor_raw->SetBinError(i+1,0);
                }
        }


	TF1 * modelunf = new TF1("modelunf","pol2",0,30);
	Unfolding_factor_raw->Fit("modelunf","w0");
	
	for(int i=0;i< Unfolding_factor->GetNbinsX();i++) {
		Unfolding_factor->SetBinContent(i+1,modelunf->Eval(Unfolding_factor->GetBinCenter(i+1)));
		Unfolding_factor->SetBinError(i+1,Unfolding_factor_raw->GetBinError(i+1));
	
	}	

	for(int i=0;i<Counts_density_unfolded->GetNbinsX();i++)
		{ 
			if(!forcefolded) Counts_density_unfolded->SetBinContent(i+1,Counts_density_unfolded->GetBinContent(i+1)*Unfolding_factor->GetBinContent(i+1));
			float error = pow(Counts_density->GetBinError(i+1)/Counts_density->GetBinContent(i+1),2);
			error+=pow(Unfolding_factor->GetBinError(i+1)/Unfolding_factor->GetBinContent(i+1),2);
			Counts_density_unfolded->SetBinError(i+1,pow(error,0.5)*Counts_density_unfolded->GetBinContent(i+1));
		}


	test = new TCanvas(("Test Unfolding"+basename).c_str());
	test->cd();
	TH1F * Counts_density_draw= (TH1F*)Counts_density->Clone("Counts density");
	Counts_density_draw->SetLineColor(2);
	Counts_density_draw->SetLineWidth(2);
	//TH2F * frame=new TH2F("frame","frame",1000, 0.5*Counts_density_draw->GetBinLowEdge(1),1.5*Counts_density_draw->GetBinLowEdge(Counts_density_draw->GetNbinsX()),1000,0.01*Counts_density_draw->GetBinContent(Counts_density_draw->GetMaximumBin()),1.5*Counts_density_draw->GetBinContent(Counts_density_draw->GetMaximumBin()) );
	Counts_density_draw->Draw();
	for(int i=0; i<iter;i++) random_meas[i]->Draw("same");
	for(int i=0;i<iter;i++) {spline_random[i]->SetMarkerStyle(8); spline_random[i]->Draw("PCsame");}
	for(int i=0;i<iter;i++) funct_random[i]->Draw("same");

*/


	return;

}

Double_t SimpleUnfModel(double *x, double *p){
        if(x[0]<=p[3])
        return p[0]*x[0]+p[1];
        else return p[2]*x[0]+(p[0]-p[2])*p[3]+p[1];
}


TH1F* RegularizeRooUnfold2(TH1F * unf_factor){

	float fitrange = (unf_factor->GetBinLowEdge(unf_factor->GetNbinsX()) - unf_factor->GetBinLowEdge(1))/4;
	float fitrangemin = unf_factor->GetBinLowEdge(1)+fitrange;
	float fitrangemax = unf_factor->GetBinLowEdge(unf_factor->GetNbinsX()-1) - fitrange;


	TH1F* clone2 =(TH1F*)unf_factor->Clone("RooUnfolding_factor");
	
	if(unf_factor->GetNbinsX()<6) {
		for(int i=0;i<unf_factor->GetNbinsX();i++) {
			clone2->SetBinContent(i+1,1);   
			clone2->SetBinError(i+1,0);
		}
	return clone2;}

	//mitigate end range
	TF1 * model = new TF1("model","[0]*x+[1]",0,20);
	model->SetParLimits(0,-0.01/(4*fitrange),0.01/(4*fitrange));
	unf_factor->Fit("model","WN","",fitrangemin,fitrangemax);

	for(int i=0;i<unf_factor->GetNbinsX();i++) {
		if(unf_factor->GetBinCenter(i+1)>fitrangemax){
			clone2->SetBinContent(i+1,model->Eval(unf_factor->GetBinCenter(i+1)));
		}
	}	

	//smoothing con media mobile
	for(int i=0;i<unf_factor->GetNbinsX();i++) {
		if(i==0){
		clone2->SetBinContent(i+1,0.33*(unf_factor->GetBinContent(i+1)+unf_factor->GetBinContent(i+2)+unf_factor->GetBinContent(i+3)));
		clone2->SetBinError(i+1,0.13*unf_factor->GetBinError(i+1));
		}
		else if(i>=unf_factor->GetNbinsX()-2){
		clone2->SetBinContent(i+1,0.33*(unf_factor->GetBinContent(i)+unf_factor->GetBinContent(i-1)+unf_factor->GetBinContent(i-2)));
		clone2->SetBinError(i+1,0.13*unf_factor->GetBinError(i+1));
		}
		else{
		clone2->SetBinContent(i+1,0.33*(unf_factor->GetBinContent(i)+unf_factor->GetBinContent(i+1)+unf_factor->GetBinContent(i+2)));
		clone2->SetBinError(i+1,0.13*unf_factor->GetBinError(i+1));
		}
	}


	clone2->Smooth();


	//shift correction

	TF1* normalizeto1 = new TF1("normalizeto1","pol0",0,20);
	clone2->Fit("normalizeto1","WN","",fitrangemin,20);
	for(int i=0;i<unf_factor->GetNbinsX();i++) {
                if(clone2->GetBinCenter(i+1)>0){
                        clone2->SetBinContent(i+1,clone2->GetBinContent(i+1)/normalizeto1->Eval(unf_factor->GetBinCenter(i+1)));
	                unf_factor->SetBinContent(i+1,unf_factor->GetBinContent(i+1)/normalizeto1->Eval(unf_factor->GetBinCenter(i+1)));
	        }
        }


	return clone2;
}


TH1F* RegularizeRooUnfold(TH1F * unf_factor){

	float fitrange = (unf_factor->GetBinLowEdge(unf_factor->GetNbinsX()) - unf_factor->GetBinLowEdge(1))/4;
	float fitrangemin = unf_factor->GetBinLowEdge(1)+fitrange;
	float fitrangemax = unf_factor->GetBinLowEdge(unf_factor->GetNbinsX()-1) - fitrange;

	TH1F* clone =(TH1F*)unf_factor->Clone("RooUnfolding_factor");
	TH1F* clone2 =(TH1F*)unf_factor->Clone("RooUnfolding_factor");

	if(unf_factor->GetNbinsX()<=6) {
		for(int i=0;i<unf_factor->GetNbinsX();i++) {
			clone2->SetBinContent(i+1,1);   
			clone2->SetBinError(i+1,0.01);
		}
	return clone2;}


	TF1* SimpleModel = new TF1("SimpleModel",SimpleUnfModel,0,20,4);	

	SimpleModel->SetParLimits(0,0,0.5);
	SimpleModel->SetParLimits(1,-0.4,2);
	SimpleModel->SetParLimits(2,0,0.005);
	SimpleModel->SetParLimits(3, clone->GetBinLowEdge(2) + fitrange,clone->GetBinLowEdge(2) + (3/2.)*fitrange);
	
	clone->Fit("SimpleModel","","N",clone->GetBinLowEdge(2),clone->GetBinLowEdge(clone->GetNbinsX()-2));
	
	for(int i=0;i<unf_factor->GetNbinsX();i++) {
		double x[1] = {clone->GetBinCenter(i+1)};
		double c1[1];
		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(1,1,x,c1,0.68);
		clone2->SetBinContent(i+1,SimpleModel->Eval(unf_factor->GetBinCenter(i+1)));
	}

	//shift correction

	TF1* normalizeto1 = new TF1("normalizeto1","pol0",0,20);
	clone2->Fit("normalizeto1","WN","",fitrangemin,20);
	for(int i=0;i<unf_factor->GetNbinsX();i++) {
                if(clone2->GetBinCenter(i+1)>0){
                        clone2->SetBinContent(i+1,clone2->GetBinContent(i+1)/normalizeto1->Eval(unf_factor->GetBinCenter(i+1)));
			clone2->SetBinError(i+1,0.1*fabs(clone2->GetBinContent(i+1)-1) );
	                unf_factor->SetBinContent(i+1,unf_factor->GetBinContent(i+1)/normalizeto1->Eval(unf_factor->GetBinCenter(i+1)));
	        }
        }

	clone2->Smooth(2);
	return clone2;

}


void Flux::Roounfold(int iterations, bool usepdf){

	Counts_density_roounf = (TH1F*) Counts_density->Clone();
	Counts_density_roounf -> Reset();
	Counts_density_roounf ->SetName("Counts_density_roounf"); 

	RooUnfolding_factor = (TH1F*) Counts_density_Ekin->Clone();
	RooUnfolding_factor -> Reset();
	RooUnfolding_factor->SetName("RooUnfolding_factor"); 

	RooUnfolding_Err = (TH1F*) Counts_density->Clone();
	RooUnfolding_Err -> Reset();
	RooUnfolding_Err->SetName("RooUnfolding_Err"); 

	//default values for NO UNFOLDING
	for(int i=0;i< RooUnfolding_factor->GetNbinsX();i++){
		RooUnfolding_factor->SetBinContent(i+1,1);
		RooUnfolding_factor->SetBinError(i+1,0.01);
		RooUnfolding_Err->SetBinContent(i+1,0.01);
		RooUnfolding_Err->SetBinError(i+1,0.0);
	}

	if(!rooUnfolding) return;         

	TH1F * tmp = (TH1F*) roounfold_meas->Clone("tmp");
	tmp->Reset();

	ExposureTime_plot     =  ConvertBinnedHisto(ExposureTime,"ExposureTime Ekin",bins,true);
	Eff_Acceptance_plot = ConvertBinnedHisto(Eff_Acceptance,"Eff. Acceptance Ekin",bins,true);

	//	Counts_density_Ekin->Divide(ExposureTime_plot);
	Counts_density_Ekin=(TH1F*)Eff_Acceptance_plot->Clone();

	if(usepdf) {
		for(int i=-1;i<tmp->GetNbinsX()+1;i++){
			roounfold_meas->SetBinContent(i+1,roounfold_meas->GetBinContent(i+1)/(tmp->GetBinLowEdge(i+2)-tmp->GetBinLowEdge(i+1) ));
			roounfold_meas->SetBinError(i+1,roounfold_meas->GetBinError(i+1)/(tmp->GetBinLowEdge(i+2)-tmp->GetBinLowEdge(i+1) ));

			roounfold_true->SetBinContent(i+1,roounfold_true->GetBinContent(i+1)/(tmp->GetBinLowEdge(i+2)-tmp->GetBinLowEdge(i+1) ));
			roounfold_true->SetBinError(i+1,roounfold_true->GetBinError(i+1)/(tmp->GetBinLowEdge(i+2)-tmp->GetBinLowEdge(i+1) ));
		}
	}	

	//"extend" range of measured histo using MC	
	if(Counts_density_Ekin->GetNbinsX()>5)
		for(int i=0;i<roounfold_meas->GetNbinsX();i++) { 

			TF1 * p = new TF1("","pol1",0,20);

			Counts_density_Ekin->Fit(p,"","",Counts_density_Ekin->GetBinLowEdge(1),Counts_density_Ekin->GetBinLowEdge(3));

			if(tmp->GetBinLowEdge(i+1)<Counts_density_Ekin->GetBinLowEdge(1)){
				tmp->SetBinContent(i+1,roounfold_meas->GetBinContent(i+1)/roounfold_meas->GetBinContent(roounfold_meas->FindBin(Counts_density_Ekin->GetBinCenter(1)))*Counts_density_Ekin->GetBinContent(1));
				tmp->SetBinError(i+1,Counts_density_Ekin->GetBinError(1));
			}
			else if(tmp->GetBinLowEdge(i+1)<Counts_density_Ekin->GetBinLowEdge(Counts_density_Ekin->GetNbinsX())){

				tmp->SetBinContent(i+1, Counts_density_Ekin->GetBinContent(Counts_density_Ekin->FindBin(tmp->GetBinCenter(i+1))) ); 
				tmp->SetBinError(i+1, Counts_density_Ekin->GetBinError(Counts_density_Ekin->FindBin(tmp->GetBinCenter(i+1))));
			}
			else {

				/*			TF1 * ff = new TF1("ff","[0]*x^[1]",0,50);
							Counts_density_Ekin->Fit(ff,"","w",Counts_density_Ekin->GetBinCenter(Counts_density_Ekin->GetNbinsX()-4),Counts_density_Ekin->GetBinCenter(Counts_density_Ekin->GetNbinsX()));
							tmp->SetBinContent(i+1,ff->Eval(tmp->GetBinCenter(i+1)));
							tmp->SetBinError(i+1,Counts_density_Ekin->GetBinError(Counts_density_Ekin->GetNbinsX()));
							*/
				tmp->SetBinContent(i+1,roounfold_meas->GetBinContent(i+1)/roounfold_meas->GetBinContent(roounfold_meas->FindBin(Counts_density_Ekin->GetBinCenter(Counts_density_Ekin->GetNbinsX())))*Counts_density_Ekin->GetBinContent(Counts_density_Ekin->GetNbinsX()));
				tmp->SetBinError(i+1,Counts_density_Ekin->GetBinError(Counts_density_Ekin->GetNbinsX()));


			}
		}

	if(usepdf){
		for(int i=-1;i<tmp->GetNbinsX()+1;i++){
			tmp->SetBinContent(i+1,tmp->GetBinContent(i+1)/(tmp->GetBinLowEdge(i+2)-tmp->GetBinLowEdge(i+1) ));
			tmp->SetBinError(i+1,tmp->GetBinError(i+1)/(tmp->GetBinLowEdge(i+2)-tmp->GetBinLowEdge(i+1) ));
		}
	}

	//	Counts_density_Ekin->Multiply(ExposureTime_plot);
	//     	Counts_density_Ekin->Multiply(Eff_Acceptance_plot);

	tmp->Smooth();

	if(usepdf){	for(int i=-2;i<migr_matr_rootunfold->GetNbinsX()+2;i++)
		for(int j=-1;j<migr_matr_rootunfold->GetNbinsX()+1;j++){
			float x1 = migr_matr_rootunfold->GetXaxis()->GetBinLowEdge(i+1);
			float x2 = migr_matr_rootunfold->GetXaxis()->GetBinLowEdge(i+2);
			migr_matr_rootunfold->SetBinContent(i+1,j+1,migr_matr_rootunfold->GetBinContent(i+1,j+1)/(x2-x1));
			migr_matr_rootunfold->SetBinError(i+1,j+1,roounfold_meas->GetBinError(i+1,j+1)/(x2-x1));
		}

		//	roounfold_true->Scale(roounfold_meas->Integral()/roounfold_true->Integral());
	}
	RooUnfoldResponse response (roounfold_meas,roounfold_true,migr_matr_rootunfold);

	cout << "==================================== UNFOLD ===================================" << endl;
	if(!usepdf)	response.UseOverflow();        




	RooUnfoldBayes   unfold (&response, tmp, iterations); 


	TH1D* hReco= (TH1D*) unfold.Hreco();	
	cout<<" Unfolding successful"<< hReco->GetNbinsX()<<" "<<Counts_density->GetNbinsX()<<" "<<Counts_density_Ekin->GetNbinsX()<<endl;		


	Counts_density_hRecoroounf = (TH1F*) hReco->Clone("RooUnfolding_Posterior");
	Prior = (TH1F*)tmp->Clone("RooUnfolding_Prior");

	//	RooUnfolding_factor_raw= (TH1F*) RooUnfolding_factor->Clone("RooUnfolding_factor_raw");	

	RooUnfolding_factor_raw= (TH1F*) Counts_density_hRecoroounf ->Clone("RooUnfolding_factor_raw");	
	RooUnfolding_factor_raw->Divide(tmp);


	RooUnfolding_factor = (TH1F*) RooUnfolding_factor_raw->Clone("RooUnfolding_factor");
	RooUnfolding_factor->SetBinContent(1,0.9*RooUnfolding_factor_raw->GetBinContent(2));
	RooUnfolding_factor->Smooth();

	RooUnfolding_factor_g= new TGraphErrors(RooUnfolding_factor);
	RooUnfolding_factor_g->SetName("RooUnfolding_factor_g");


	for(int i=0;i<RooUnfolding_factor->GetNbinsX();i++) {
		RooUnfolding_factor->SetBinContent(i+1,RooUnfolding_factor_g->Eval(RooUnfolding_factor->GetBinCenter(i+1)));
		RooUnfolding_factor->SetBinError(i+1,0.01);


         	float unferr =  fabs(RooUnfolding_factor->GetBinContent(i+1) -1)/sqrt(12);
		RooUnfolding_Err->SetBinContent(i+1, unferr);
	}	

/*	for(int i=0;i<RooUnfolding_factor->GetNbinsX();i++) {
		RooUnfolding_factor->SetBinContent(i+1,RooUnfolding_factor_raw->GetBinContent(RooUnfolding_factor_raw->FindBin(Counts_density_Ekin->GetBinCenter(i+1))));
		//		 float unf_err =  fabs(RooUnfolding_factor->GetBinContent(i+1) -1)*0.1;
		if(i==0) RooUnfolding_factor->SetBinContent(i+1,0.9*RooUnfolding_factor_raw->GetBinContent(i+2));
	
		RooUnfolding_factor->SetBinError(i+1,RooUnfolding_factor_raw->GetBinError(RooUnfolding_factor_raw->FindBin(Counts_density_Ekin->GetBinCenter(i+1))));		
	}
*/

	// EXTERNAL UNFOLDING
/*	if(Unfolding_factor_timeavg && Unfolding_factor_timeavg_err)
	{

		for(int i=0;i<RooUnfolding_factor->GetNbinsX();i++) {
			RooUnfolding_factor->SetBinContent(i+1,Unfolding_factor_timeavg->Eval(RooUnfolding_factor->GetBinCenter(i+1)));
			RooUnfolding_factor->SetBinError(i+1,Unfolding_factor_timeavg_err->Eval(RooUnfolding_factor->GetBinCenter(i+1)));


			float unferr =  fabs(RooUnfolding_factor->GetBinContent(i+1) -1)/sqrt(12);
			RooUnfolding_Err->SetBinContent(i+1, unferr);
		}
	}

*/
//

	
	Counts_density_Ekin_unf = (TH1F*) Counts_density_Ekin->Clone("Counts density Ekin unfolded");
	Counts_density_Ekin_unf -> Multiply(RooUnfolding_factor);

	for(int i=0;i<Counts_density_roounf->GetNbinsX();i++){
		int ekinbin= RooUnfolding_factor->FindBin(Counts_density_Ekin_unf -> GetBinCenter(i+1))-1 ;
	//	  Counts_density_roounf->SetBinContent(i+1,Counts_density->GetBinContent(i+1)*RooUnfolding_factor->GetBinContent(ekinbin));		
		  Counts_density_roounf->SetBinContent(i+1,Counts_density->GetBinContent(i+1)*RooUnfolding_factor_g->Eval(Counts_density_Ekin_unf -> GetBinCenter(i+1)));		
		  Counts_density_roounf->SetBinError(i+1,Counts_density->GetBinError(i+1)*RooUnfolding_factor->GetBinContent(ekinbin));		
	}


	return; 

}



void Flux::Eval_Flux(float corr_acc, float fit_min, float fit_max,int knots, float offset,bool regularize){


	// COUNTS
	cout<<endl;
	cout<<"********************************************"<<endl;

	if(!(Counts>0)) { Counts=new TH1F("Counts_dummy","Counts_dummy",bins.size(),0,bins.size());
		Counts->FillRandom("gaus",2000); 	}

	cout<<"Counts "<<Counts->GetName()<<" "<<Counts->GetEntries()<<endl;	

	Counts_density 		= ConvertBinnedHisto(Counts,"Counts density",bins,false);
	Counts_density_unfolded = ConvertBinnedHisto(Counts,"Counts density Unfolded",bins,false);
	Counts_density_Ekin 	= ConvertBinnedHisto(Counts,"Counts density Ekin",bins,true);
	Counts_density_Ekin_unf	= ConvertBinnedHisto(Counts,"Counts density Ekin unfolded",bins,true);

	Counts_density 		-> Sumw2();
	Counts_density_unfolded -> Sumw2();
	Counts_density_Ekin 	-> Sumw2();
	Counts_density_Ekin_unf	-> Sumw2();
	Counts_density27= Twentisevenify_(Counts_density);
	DeltaE=(TH1F*)Counts_density->Clone("DeltaE"); 	
	DeltaE->Reset();

	
	
	// ROOUNFOLD
	if(Counts>0) Roounfold(4,Usepdf);	

	// DeltaE
	for(int i=0;i<bins.size();i++){
		Counts_density->SetBinError(i+1,Counts->GetBinError(i+1)/(bins.RigTOIBins()[i+1]-bins.RigTOIBins()[i]));
                Counts_density->SetBinContent(i+1,Counts->GetBinContent(i+1)/(bins.RigTOIBins()[i+1]-bins.RigTOIBins()[i]));	
		Counts_density_unfolded->SetBinError(i+1,Counts->GetBinError(i+1)/(bins.RigTOIBins()[i+1]-bins.RigTOIBins()[i]));
                Counts_density_unfolded->SetBinContent(i+1,Counts->GetBinContent(i+1)/(bins.RigTOIBins()[i+1]-bins.RigTOIBins()[i]));	
	if(Counts_density_roounf){
			Counts_density_roounf->SetBinError(i+1,Counts_density_roounf->GetBinError(i+1)/(bins.RigTOIBins()[i+1]-bins.RigTOIBins()[i]));
			Counts_density_roounf->SetBinContent(i+1,Counts_density_roounf->GetBinContent(i+1)/(bins.RigTOIBins()[i+1]-bins.RigTOIBins()[i]));	
		}
	Counts_density_Ekin->SetBinError(i+1,Counts->GetBinError(i+1)/(bins.EkPerMasTOIBins()[i+1]-bins.EkPerMasTOIBins()[i]));
                Counts_density_Ekin->SetBinContent(i+1,Counts->GetBinContent(i+1)/(bins.EkPerMasTOIBins()[i+1]-bins.EkPerMasTOIBins()[i]));
		Counts_density_Ekin_unf->SetBinError(i+1,Counts->GetBinError(i+1)/(bins.EkPerMasTOIBins()[i+1]-bins.EkPerMasTOIBins()[i]));
                Counts_density_Ekin_unf->SetBinContent(i+1,Counts->GetBinContent(i+1)/(bins.EkPerMasTOIBins()[i+1]-bins.EkPerMasTOIBins()[i]));
	      

		Counts_density27->SetBinContent(i+1,Counts_density27->GetBinContent(i+1)/(bins.RigTOIBins()[i+1]-bins.RigTOIBins()[i]));
		Counts_density27->SetBinError(i+1,Counts_density27->GetBinError(i+1)/(bins.RigTOIBins()[i+1]-bins.RigTOIBins()[i]));

		DeltaE->SetBinContent(i+1,(bins.RigTOIBins()[i+1]-bins.RigTOIBins()[i]));
	}
	      	
	// EXPOSURE TIME	
	if(ExposureTime) {
		ExposureTime_plot_rig =  ConvertBinnedHisto(ExposureTime,"ExposureTime R",bins,false);
		ExposureTime_plot     =  ConvertBinnedHisto(ExposureTime,"ExposureTime Ekin",bins,true);
	

		//////////////////////// ????????????????????????????????	
		//ExposureTime_plot_rig->Scale(1.2369646);
                //ExposureTime_plot->Scale(1.2369646);
    		//////////////////////// ????????????????????????????????
	
		Counts_density_Ekin -> Divide(ExposureTime_plot);
        	Counts_density_Ekin_unf -> Divide(ExposureTime_plot);
                Counts_density -> Divide(ExposureTime_plot_rig);
	        Counts_density_unfolded -> Divide(ExposureTime_plot_rig);
	        if(Counts_density_roounf)  Counts_density_roounf -> Divide(ExposureTime_plot_rig);
		}


	
	if(Eff_Acceptance) {

		Eff_Acceptance_plot = ConvertBinnedHisto(Eff_Acceptance,"Eff. Acceptance Ekin",bins,true);
		Eff_Acceptance_plot_rig = ConvertBinnedHisto(Eff_Acceptance,"Eff. Acceptance R",bins,false);

		Eff_Acceptance_plot->Scale(corr_acc);
		Eff_Acceptance_plot_rig->Scale(corr_acc);
		
		Counts_density -> Divide(Eff_Acceptance_plot_rig);
	        Counts_density_unfolded -> Divide(Eff_Acceptance_plot_rig);
	
		Unfold_Counts(fit_min,fit_max,knots,offset,regularize);	

		Counts_density -> Multiply(Eff_Acceptance_plot_rig);
	        Counts_density_unfolded -> Multiply(Eff_Acceptance_plot_rig);

	}
	

	if(RooUnfolding_factor)  Counts_density_Ekin_unf->Multiply(RooUnfolding_factor);


	FluxEstim = (TH1F *) Counts_density_Ekin->Clone();
	FluxEstim -> SetName((basename+"_Flux").c_str());
	FluxEstim -> SetTitle((basename+"_Flux").c_str());
	FluxEstim -> Sumw2();

	FluxEstim_ekin_unf = (TH1F *) Counts_density_Ekin_unf->Clone();
	FluxEstim_ekin_unf -> SetName((basename+"_Flux_ekin_unf").c_str());
	FluxEstim_ekin_unf -> SetTitle((basename+"_Flux_ekin_unf").c_str());
	FluxEstim_ekin_unf -> Sumw2();

	FluxEstim_rig = (TH1F*) Counts_density ->Clone();
	FluxEstim_rig -> SetName((basename+"_Flux_rig").c_str());
	FluxEstim_rig -> SetTitle((basename+"_Flux_rig").c_str());
	FluxEstim_rig -> Sumw2();
	
	if(Counts_density_roounf )	FluxEstim_unf = (TH1F*) Counts_density_roounf ->Clone();  //default: roounfold bayes
	else FluxEstim_unf = (TH1F*) Counts_density_unfolded ->Clone();  

	FluxEstim_unf -> SetName((basename+"_Flux_unf").c_str());
	FluxEstim_unf -> SetTitle((basename+"_Flux_unf").c_str());
	FluxEstim_unf -> Sumw2();
	

	//ACCEPTANCE
	if(Eff_Acceptance) {

		Eff_Acceptance_plot = ConvertBinnedHisto(Eff_Acceptance,"Eff. Acceptance Ekin",bins,true);
		Eff_Acceptance_plot_rig = ConvertBinnedHisto(Eff_Acceptance,"Eff. Acceptance R",bins,false);

		Eff_Acceptance_plot->Scale(corr_acc);
		Eff_Acceptance_plot_rig->Scale(corr_acc);
		
	
	//	ModelAcceptanceWithSpline(0);
		
		FluxEstim -> Divide(Eff_Acceptance_plot);
		FluxEstim_ekin_unf -> Divide(Eff_Acceptance_plot);
		FluxEstim_rig -> Divide(Eff_Acceptance_plot_rig);
		FluxEstim_unf -> Divide(Eff_Acceptance_plot_rig);
	
	}		

/*	for(int i=0;i<FluxEstim_rig->GetNbinsX();i++){
			FluxEstim ->SetBinContent(i+1,FluxEstim->GetBinContent(i+1)/AcceptanceModel->Eval(bins.ekpermassbincent_TOI[i]));
			FluxEstim ->SetBinError(i+1,FluxEstim->GetBinError(i+1)/AcceptanceModel->Eval(bins.ekpermassbincent_TOI[i]));
			FluxEstim_rig ->SetBinContent(i+1,FluxEstim_rig->GetBinContent(i+1)/AcceptanceModel->Eval(bins.RigTOIBinsCent()[i]));
			FluxEstim_rig ->SetBinError(i+1,FluxEstim_rig->GetBinError(i+1)/AcceptanceModel->Eval(bins.RigTOIBinsCent()[i]));
			FluxEstim_unf ->SetBinContent(i+1,FluxEstim_unf->GetBinContent(i+1)/AcceptanceModel->Eval(bins.RigTOIBinsCent()[i]));
			FluxEstim_unf ->SetBinError(i+1,FluxEstim_unf->GetBinError(i+1)/AcceptanceModel->Eval(bins.RigTOIBinsCent()[i]));
		}
	}
*/
	Eval_Errors();

	return;		
}

void Flux::SaveResults(FileSaver finalhistos,std::string bn){
	std::cout<<"Results: "<<FluxEstim<<" "<<ExposureTime<<" "<<Eff_Acceptance<<std::endl;
	if(bn!="") basename=bn; 

	if(FluxEstim) finalhistos.Add(FluxEstim); 	
	if(FluxEstim_ekin_unf) finalhistos.Add(FluxEstim_ekin_unf); 	
	if(FluxEstim_rig) finalhistos.Add(FluxEstim_rig);
	if(FluxEstim_unf) finalhistos.Add(FluxEstim_unf);
	if(FluxEstim_rig) finalhistos.Add(new TGraphErrors(Twentisevenify_(FluxEstim_rig)));
	if(FluxEstim_unf) finalhistos.Add(new TGraphErrors(Twentisevenify_(FluxEstim_unf)));
	 	
	if(ExposureTime) { 
				ExposureTime = ConvertBinnedHisto(ExposureTime,"ExposureTime R",bins,false);
				finalhistos.Add(ExposureTime);
				finalhistos.Add(ExposureTime_plot);
				finalhistos.Add(ExposureTime_plot_rig);
	}

	if(Eff_Acceptance) { finalhistos.Add(Eff_Acceptance);
		finalhistos.Add(Eff_Acceptance_plot);	
		finalhistos.Add(Eff_Acceptance_plot_rig);	
	}

	if(MC_Acceptance) { 
		cout<<" MC ACCEPTANCE "<<MC_Acceptance<<endl;
		MC_Acceptance=ConvertBinnedHisto(MC_Acceptance,"MC Acceptance",bins,false);
		finalhistos.Add(MC_Acceptance);
	}


	if(Counts){ finalhistos.Add(Counts);
		    finalhistos.Add(Counts_density);
		    finalhistos.Add(Counts_density27);
		    finalhistos.Add(Counts_density_unfolded);
		    if(Counts_density_roounf) finalhistos.Add(Counts_density_roounf);
		    if(Counts_density_hRecoroounf) finalhistos.Add(Counts_density_hRecoroounf);
		    if(Prior) finalhistos.Add(Prior);
		    finalhistos.Add(Counts_density_Ekin);
	}
        if(DeltaE) finalhistos.Add(DeltaE);
	if(Unfolding_factor_raw) finalhistos.Add(Unfolding_factor_raw);
	if(Unfolding_factor) finalhistos.Add(Unfolding_factor);
	if(RooUnfolding_factor) finalhistos.Add(RooUnfolding_factor);
	if(RooUnfolding_factor_raw) finalhistos.Add(RooUnfolding_factor_raw);
	if(RooUnfolding_factor_g) finalhistos.Add(RooUnfolding_factor_g);


	if(Counts_ForAverage) finalhistos.Add(Counts_ForAverage);
	if(Exposure_ForAverage) finalhistos.Add(Exposure_ForAverage);

	if(test) finalhistos.Add(test);
	if(migr_matr) finalhistos.Add(migr_matr);

	if(Counts_statErr) {
		Counts_statErr = ConvertBinnedHisto(Counts_statErr,"Counts_Stat_Error",bins,false);
		finalhistos.Add(Counts_statErr);
	}
	if(Counts_systErr) {
		Counts_systErr = ConvertBinnedHisto(Counts_systErr,"Counts_Syst_Error",bins,false);
		finalhistos.Add(Counts_systErr);

	}
	if(Counts_systErr_unb) {
		Counts_systErr_unb = ConvertBinnedHisto(Counts_systErr_unb,"Counts_Syst_Error_unb",bins,false);
		finalhistos.Add(Counts_systErr_unb);
	}
	if(Acc_Err) {
		Acc_Err = ConvertBinnedHisto(Acc_Err,"Acceptance_TOT_Error",bins,false);
		finalhistos.Add(Acc_Err);

	}
	if(Acc_ErrStat) {
		Acc_ErrStat = ConvertBinnedHisto(Acc_ErrStat,"Acceptance_Stat_Error",bins,false);
		finalhistos.Add(Acc_ErrStat);

	}

	if(Unfolding_Err) {
		Unfolding_Err = ConvertBinnedHisto(Unfolding_Err,"Unfolding_Syst_Error",bins,false);
		finalhistos.Add(Unfolding_Err);

	}

	if(RooUnfolding_Err) {
		RooUnfolding_Err = ConvertBinnedHisto(RooUnfolding_Err,"RooUnfolding_Syst_Error",bins,false);
		finalhistos.Add(RooUnfolding_Err);
	}
	if(TOT_statErr) {
		TOT_statErr = ConvertBinnedHisto(TOT_statErr,"TOT_Stat_Error",bins,false);
		finalhistos.Add(TOT_statErr);
	}
	
	if(FluxEstim_unf){
		FluxEstim_unf_stat = (TH1F*) FluxEstim_unf->Clone();
	 	FluxEstim_unf_stat -> SetName((basename+"_Flux_unf_stat").c_str());
        	FluxEstim_unf_stat -> SetTitle((basename+"_Flux_unf_stat").c_str());
        	FluxEstim_unf_stat -> Sumw2();

		for(int i=0;i<FluxEstim_unf_stat->GetNbinsX();i++){
			float error=0;
			if(Counts_statErr) error+=pow(Counts_statErr->GetBinContent(i+1),2);
			//if(Counts_systErr) error+=pow(Counts_systErr->GetBinContent(i+1),2);
			if(Acc_ErrStat) 	   error+=pow(Acc_ErrStat->GetBinContent(i+1),2);
			//if(Unfolding_Err)  error+=pow(Unfolding_Err->GetBinContent(i+1),2);
			FluxEstim_unf_stat->SetBinError(i+1, pow(error,0.5)*FluxEstim_unf_stat->GetBinContent(i+1));
		}
		finalhistos.Add(FluxEstim_unf_stat);

	}
	
	if(FluxEstim_rig){
		FluxEstim_rig_stat = (TH1F*) FluxEstim_rig->Clone();
	 	FluxEstim_rig_stat -> SetName((basename+"_Flux_rig_stat").c_str());
        	FluxEstim_rig_stat -> SetTitle((basename+"_Flux_rig_stat").c_str());
        	FluxEstim_rig_stat -> Sumw2();

		for(int i=0;i<FluxEstim_rig_stat->GetNbinsX();i++){
			float error=0;
			if(Counts_statErr) error+=pow(Counts_statErr->GetBinContent(i+1),2);
			//if(Counts_systErr) error+=pow(Counts_systErr->GetBinContent(i+1),2);
			if(Acc_ErrStat) 	   error+=pow(Acc_ErrStat->GetBinContent(i+1),2);
		//	if(Unfolding_Err)  error+=pow(Unfolding_Err->GetBinContent(i+1),2);
			FluxEstim_rig_stat->SetBinError(i+1, pow(error,0.5)*FluxEstim_rig_stat->GetBinContent(i+1));
		}
		finalhistos.Add(FluxEstim_rig_stat);

	}
	


	finalhistos.writeObjsInFolder(("Fluxes/"+basename).c_str());
}


void Flux::Eval_ExposureTime(Variables * vars, TTree * treeDT,FileSaver finalhistos,bool refill){
		
	if(!ExposureTime||refill){
	
		ExposureTime = new TH1F(exposurename.c_str(),exposurename.c_str(),bins.size(),0,bins.size());

		int ActualTime=0;
		vars->ReadBranches(treeDT);
		for(int i=0;i<treeDT->GetEntries()/FRAC;i++){
			UpdateProgressBar(i, treeDT->GetEntries()/FRAC);
			treeDT->GetEvent(i);
			vars->Update();
			if((int)vars->U_time!=ActualTime) {
				UpdateZoneLivetime(vars->Livetime,vars->Rcutoff_IGRFRTI,ExposureTime,bins);
				ActualTime=vars->U_time;
			}
		}
		finalhistos.Add(ExposureTime); 	
		finalhistos.writeObjsInFolder(("Fluxes/"+basename).c_str());

	}
	return;
}


TH1F * Flux::Eval_FluxRatio(Flux * Denominator,std::string name){

        TH1F * Numerator = (TH1F*) FluxEstim -> Clone();
        Numerator->SetName(name.c_str());
        Numerator->SetTitle(name.c_str());
        Numerator->Sumw2();

        Numerator->Divide(Denominator->GetFlux());
        for(int i=0;i<Numerator->GetNbinsX();i++){
                float A = FluxEstim -> GetBinContent(i+1);
                float B = Denominator->GetFlux() -> GetBinContent(i+1);
                float sigA = FluxEstim -> GetBinError(i+1);
                float sigB = Denominator->GetFlux() -> GetBinError(i+1);
                Numerator->SetBinError(i+1,pow(pow(sigA/A,2)+pow(sigB/B,2)+2*sigA/A*sigB/B,0.5)*A/B);
        }

        return Numerator;
}

void Flux::Eval_Errors(){
	if(Counts>0){
		if(!Counts_statErr){
			Counts_statErr = (TH1F*) Counts->Clone("Counts_statErr");
			Counts_systErr = (TH1F*) Counts->Clone("Counts_systErr");
			Counts_systErr_unb = (TH1F*) Counts->Clone("Counts_systErr_unb");
			Counts_statErr->Reset();	
			Counts_systErr->Reset();	
			Counts_systErr_unb->Reset();	


			for(int i=0;i<Counts->GetNbinsX();i++){
				if(Counts->GetBinContent(i+1)){	
					Counts_statErr->SetBinContent(i+1,Counts->GetBinError(i+1)/Counts->GetBinContent(i+1));
					Counts_systErr->SetBinContent(i+1,0);
					Counts_systErr_unb->SetBinContent(i+1,0);
				}
			}
		}
		/*
		else{
			Counts_systErr = (TH1F*) Counts->Clone("Counts_systErr");
			Counts_systErr->Reset();	

			for(int i=0;i<Counts->GetNbinsX();i++){
			if(Counts->GetBinContent(i+1)>0)
				Counts_systErr->SetBinContent(i+1,pow(fabs( (pow(Counts->GetBinError(i+1)/Counts->GetBinContent(i+1),2) -
							          pow(Counts_statErr->GetBinContent(i+1),2))),0.5));
								
			}
		}*/

		//protection against super stat error in counts
		for(int i=0;i<Counts->GetNbinsX();i++){
			Counts_statErr->SetBinContent(i+1,min(Counts_statErr->GetBinContent(i+1),0.025));
		}
		TOT_statErr = (TH1F*) Counts_statErr->Clone("TOT_statErr");
	}

	if(Counts_systErr&&Counts_systErr_unb) for(int i=0;i<Counts_systErr_unb->GetNbinsX();i++) 
		if(Counts_systErr_unb->GetBinContent(i+1)-Counts_systErr->GetBinContent(i+1)>0) 
		   Counts_systErr_unb->SetBinContent(i+1,Counts_systErr_unb->GetBinContent(i+1)-Counts_systErr->GetBinContent(i+1));
		else
		   Counts_systErr_unb->SetBinContent(i+1,Counts_systErr_unb->GetBinContent(i+1));

	if(Eff_Acceptance>0){
		Acc_Err =(TH1F*) EffAcceptance->GetSyst_Err();
		Acc_Err->SetName("Acc_Err");
		Acc_ErrStat=(TH1F*)EffAcceptance->GetStat_Err();
		//correcting the fragmentation error 
		for(int i=0;i<Acc_Err->GetNbinsX();i++){
			float err= Acc_Err->GetBinContent(i+1);
			err=pow(err,2);
			err-=0.02*0.02;
			err+=pow(systfragm,2);
			err=sqrt(err);
			Acc_Err->SetBinContent(i+1,err);
			Acc_Err->SetBinError(i+1,0);
		}

		for(int j=0;j<TOT_statErr->GetNbinsX();j++) 
			TOT_statErr->SetBinContent(j+1,pow(pow(Acc_ErrStat->GetBinContent(j+1),2)+
						       pow(RooUnfolding_factor->GetBinError(j+1)/RooUnfolding_factor->GetBinContent(j+1),2)+ 
						       pow(TOT_statErr->GetBinContent(j+1),2),0.5));

	}
	if(Unfolding_factor>0){
		Unfolding_Err = (TH1F*) Unfolding_factor->Clone("Unfolding_Err");	
		Unfolding_Err->Reset();
		for(int i=0;i<Unfolding_factor->GetNbinsX();i++){
                                Unfolding_Err->SetBinContent(i+1,Unfolding_factor->GetBinError(i+1)/Unfolding_factor->GetBinContent(i+1));
                        }

	}

	for(int i=0;i<FluxEstim_rig->GetNbinsX();i++){
		float error=0;
		if(Counts_statErr) error+=pow(Counts_statErr->GetBinContent(i+1),2);
		if(Counts_systErr) error+=pow(Counts_systErr->GetBinContent(i+1),2);
		if(Counts_systErr_unb) error+=pow(Counts_systErr_unb->GetBinContent(i+1),2);
		if(Acc_ErrStat)            error+=pow(Acc_ErrStat->GetBinContent(i+1),2);
		if(Acc_Err)               error+=pow(Acc_Err->GetBinContent(i+1),2);
		FluxEstim_rig->SetBinError(i+1, pow(error,0.5)*FluxEstim_rig->GetBinContent(i+1));
	}

	for(int i=0;i<FluxEstim_unf->GetNbinsX();i++){
		float error=0;
		if(Counts_statErr) error+=pow(Counts_statErr->GetBinContent(i+1),2);
		if(Counts_systErr) error+=pow(Counts_systErr->GetBinContent(i+1),2);
		if(Counts_systErr_unb) error+=pow(Counts_systErr_unb->GetBinContent(i+1),2);
		if(Acc_ErrStat)            error+=pow(Acc_ErrStat->GetBinContent(i+1),2);
		if(Acc_Err)            error+=pow(Acc_Err->GetBinContent(i+1),2);
		if(RooUnfolding_Err)  error+=pow(RooUnfolding_Err->GetBinContent(i+1),2);
		FluxEstim_unf->SetBinError(i+1, pow(error,0.5)*FluxEstim_unf->GetBinContent(i+1));
	}


}

void Flux::ModelAcceptanceWithSpline(float shift){


	int nodes = bins.size()/4+2;
	double p[nodes];
	double x[nodes];

	for(int i=0;i<nodes;i++){
		x[i]=bins.RigTOIBins()[0]+((bins.RigTOIBins()[bins.size()] - bins.RigTOIBins()[0])/(nodes-1))*i;
		p[i]=Eff_Acceptance_plot_rig->GetBinContent(Eff_Acceptance_plot_rig->FindBin(x[i]));
	}
	x[nodes-1]=bins.RigTOIBins()[bins.size()];
	p[nodes-1]=Eff_Acceptance_plot_rig->GetBinContent(Eff_Acceptance_plot_rig->FindBin(x[nodes-1])-1);
	
	cout<<"******** Model Acc. "<<basename<<" ***********"<<endl;
	std::vector<float> spline_x;
	for(int i=0;i<nodes;i++) {
		spline_x.push_back(shift + x[i]);
		cout<<x[i]<<" "<<p[i]<<endl;
	}
	Model *M = new Model(spline_x);

	AcceptanceModel = new TF1(basename.c_str(),M,&Model::Function,shift,50,nodes,"Model","Function");
 
	for(int i=0;i<nodes;i++) {
		AcceptanceModel->SetParameter(i,p[i]);
		AcceptanceModel->SetParLimits(i,0.8*p[i],1.2*p[i]);
		if(i==nodes-1)	AcceptanceModel->SetParLimits(i,0.95*p[i],1.02*p[i]);
	
	}	

	Eff_Acceptance_plot_rig->Fit(basename.c_str(),"","M",x[0],x[nodes-1]);

}



void Flux::Average_with_another(TH1F * ratio){

	TGraphErrors * RT = new TGraphErrors(ratio);
}
