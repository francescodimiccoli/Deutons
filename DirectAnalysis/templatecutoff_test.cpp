float bins[18] = {

  0.53,0.5432,0.5568,0.5708, 0.585,0.5996,0.6146,  0.63,0.6458,0.6619,0.6784,0.6954,0.7128,0.7306,0.7489,0.7676,0.7868,0.8065,
};


float GetCutoffWeight(TSpline3 * ExposureTime, float particle_m, float beta, float m){
        float weight =0;
        if(m<=particle_m) weight = ExposureTime->Eval(m*beta*sqrt(1/(1-beta*beta)));
        else weight = ExposureTime->Eval(particle_m*beta*sqrt(1/(1-beta*beta)));
        if (weight<0) return 0;
        else return weight;
}


void ReweightHisto(TH1F * Histo,TSpline3 * ExposureTime,float particle_m, float beta){
	for(int i=0;i<Histo->GetNbinsX();i++){
		Histo->SetBinContent(i+1,Histo->GetBinContent(i+1)*GetCutoffWeight(ExposureTime,particle_m,beta,Histo->GetBinLowEdge(i)));
		Histo->SetBinError(i+1,sqrt(Histo->GetBinContent(i+1)*GetCutoffWeight(ExposureTime,particle_m,beta,Histo->GetBinCenter(i+1))));
	}
	return;
}

void ScaleErrors(TH1F * Histo, float scale){

	for(int i=0;i<Histo->GetNbinsX();i++){
		Histo->SetBinError(i+1,Histo->GetBinError(i+1)*scale);
	}
	

} 

void Do_TemplateFit(TH1F * Data, TH1F* modelP, TH1F* modelD, TH1F* modelPsec);

float GetChiSquare(TH1 * Result, TH1 * Data, float min, float max){

	int binmin = Data->FindBin(min);
	int binmax = Data->FindBin(max);
	float chi = 0;
	float err = 0;
	for(int i=binmin; i<binmax; i++){
		if(Data->GetBinContent(i+1)>0&&Result->GetBinContent(i+1)>0){
			err = pow(pow(Data->GetBinError(i+1),2) + pow(Result->GetBinError(i+1),2),0.5);
			chi += pow((Data->GetBinContent(i+1) - Result->GetBinContent(i+1)),2)/pow(err,2);
		}
	}
	chi /= ((binmax-binmin)-3);
	cout<<"chi: "<<chi<<endl;
	return chi;

}


int templatecutoff_test(){

	TFile *file = TFile::Open("AnalysisFiles/1306886400-1383264000/Result.root_Results");
	cout<<file<<endl;	
	TH1F * DataPrim    = (TH1F *) file->Get("TOFfits/Fit Results/Data/Bin3/TOFfits_DataPrim_3 0 3");
	TH1F * Data        = (TH1F *) file->Get("TOFfits/Fit Results/Data/Bin3/TOFfits_Data_3 0 3");
	TH1F * Templ_P     = (TH1F *) file->Get("TOFfits/Fit Results/ScaledTemplatesP/Bin3/Best #chi^{2} Mod. Proton MC");
	TH1F * Templ_D     = (TH1F *) file->Get("TOFfits/Fit Results/ScaledTemplatesD/Bin3/Best #chi^{2} Mod. Deuton MC");
	TH1F * Templ_He    = (TH1F *) file->Get("TOFfits/Fit Results/ScaledTemplatesHe/Bin3/Best #chi^{2} Tritium MC");
	TSpline3 * ExposureTime = (TSpline3 *) file->Get("TOFfits/ExposureModel/Spline3");	

	TH1F * Templ_PPrim     = (TH1F *)Templ_P->Clone(); 
	TH1F * Templ_PSec      = (TH1F *)Templ_P->Clone(); 
	TH1F * Templ_DPrim     = (TH1F *)Templ_D->Clone(); 
	TH1F * Templ_DPrim     = (TH1F *)Templ_D->Clone(); 
	TH1F * Templ_HePrim    = (TH1F *)Templ_He->Clone(); 



	ReweightHisto(Templ_PPrim,ExposureTime,0.938,bins[3]);
	ReweightHisto(Templ_PSec,ExposureTime,3,bins[3]);
	ReweightHisto(Templ_DPrim,ExposureTime,1.875,bins[3]);
	ReweightHisto(Templ_HePrim,ExposureTime,3,bins[3]);

	Templ_PPrim->SetLineColor(2); 
	Templ_PSec->SetLineColor(6); 
	Templ_DPrim ->SetLineColor(4); 
	Templ_He->SetLineColor(3); 
	Templ_PPrim->SetMarkerColor(2); 
	Templ_PSec->SetMarkerColor(6); 
	Templ_DPrim ->SetMarkerColor(4); 
	Templ_He->SetMarkerColor(3); 

	float norm = DataPrim->GetBinContent(DataPrim->GetMaximumBin())/Templ_PPrim->GetBinContent(Templ_PPrim->GetMaximumBin());

	Templ_PPrim ->Scale(norm);
	Templ_PSec  ->Scale(norm);
	Templ_DPrim ->Scale(norm);
	Templ_HePrim->Scale(norm);


	Templ_PPrim_copy =(TH1F *)Templ_PPrim->Clone(); 
        Templ_PSec_copy  =(TH1F *)Templ_PSec->Clone(); 
	Templ_DPrim_copy =(TH1F *)Templ_DPrim->Clone(); 
	Templ_HePrim_copy=(TH1F *)Templ_HePrim->Clone(); 

	TCanvas * c1 = new TCanvas("Templates");
	c1->cd();
	Templ_PSec_copy->Draw("histo");
	Templ_PPrim_copy->Draw("same,histo");
	Templ_DPrim_copy ->Draw("same,histo");
	Templ_HePrim_copy->Draw("same,histo");
	gPad->SetLogy();
	gPad->SetLogx();
	
	TH1F *Result =  Do_TemplateFit(DataPrim,Templ_PPrim,Templ_DPrim,Templ_PSec);

	TH1F* Sum = (TH1F*) Templ_PPrim->Clone();
	Sum->Add(Templ_PSec);
	Sum->Add(Templ_DPrim);
	Sum->SetLineColor(1);
	
	cout<<"*** CHISQUARE *****"<<endl;
	cout<<GetChiSquare(Sum,DataPrim,0,4)<<endl;


	TCanvas * c2 = new TCanvas("Fit");
	c2->cd();
	DataPrim->Draw("Psame");
	Templ_PPrim->Draw("same,histo");
	Templ_PSec->Draw("same,histo");
	Templ_DPrim ->Draw("same,histo");
	Templ_He->Draw("same,histo");
	Sum->Draw("same,histo");
	gPad->SetLogy();
	gPad->SetLogx();
	


	return 0;
}


TH1F * Do_TemplateFit(TH1F * Data, TH1F* modelP, TH1F* modelD, TH1F* modelPsec){

	Tpl = new TObjArray(3);
	Tpl->Add(modelP);
	Tpl->Add(modelPsec);
	Tpl->Add(modelD);
	TFractionFitter * Tfit = new TFractionFitter(Data, Tpl ,"q");

	Tfit -> Constrain(1, 0.0 ,1);
	Tfit -> Constrain(2, 0.0,1);
	Tfit -> Constrain(3, 0.001 ,0.2);

	int fit_outcome = Tfit->Fit();
	
	if(fit_outcome==0){
			TH1F * Result = (TH1F *) Tfit -> GetPlot();
			float itot= Result->Integral();
			double w1,e1 = 0;
			double w2,e2 = 0;
			double w3,e3 = 0;
		
			Tfit ->GetResult(0,w1,e1);
			Tfit ->GetResult(1,w2,e2);
			Tfit ->GetResult(2,w3,e3);
			
			float i1=1;
			float i2=1;
			float i3=1;
			
			i1 = modelP  ->Integral();
		 	i2 = modelPsec  ->Integral();
			i3 = modelD->Integral();
			
			cout<<w1<<" "<<w2<<" "<<w3<<endl;
 
			modelP->Scale(w1*itot/i1);	
			modelPsec->Scale(w2*itot/i2);	
			modelD->Scale(w3*itot/i3);	

			return Result;
	}
	return (TH1F*) Data->Clone();
}
