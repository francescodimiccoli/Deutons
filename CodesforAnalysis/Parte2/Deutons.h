
TemplateFIT * FitTOF_Dbins	= new TemplateFIT("FitTOF_Dbins",18,0,3,11);
TemplateFIT * FitNaF_Dbins	= new TemplateFIT("FitNaF_Dbins",18,0,3,11);
TemplateFIT * FitAgl_Dbins	= new TemplateFIT("FitAgl_Dbins",18,0,3,11);

TemplateFIT * FitTOF_Pbins	= new TemplateFIT("FitTOF_Pbins",18,0,3,11);
TemplateFIT * FitNaF_Pbins	= new TemplateFIT("FitNaF_Pbins",18,0,3,11);
TemplateFIT * FitAgl_Pbins	= new TemplateFIT("FitAgl_Pbins",18,0,3,11);



void DeutonsMC_Fill(TNtuple *ntupla, int l){
	int k = ntupla->GetEvent(l);
	if(Beta<=0||R<=0) return;
	float mass = 0;
	if(!(Likcut&&Distcut)) return;	

	for(int m=0;m<nbinsToF;m++){ //TOF
		
		mass = ((R/Beta)*pow((1-pow(Beta,2)),0.5));
		if(Var>BetaD[m]&&Var<=BetaD[m+1]){ 
			if(Massa_gen<1&&Massa_gen>0.5)	FitTOF_Dbins -> TemplateP -> Fill(mass,m);
			if(Massa_gen<2&&Massa_gen>1.5)	FitTOF_Dbins -> TemplateD -> Fill(mass,m);	
			if(Massa_gen<4&&Massa_gen>2.5)	FitTOF_Dbins -> TemplateHe-> Fill(mass,m);
		}
		if(Var>BetaP[m]&&Var<=BetaP[m+1]) {
			if(Massa_gen<1&&Massa_gen>0.5)	FitTOF_Pbins -> TemplateP -> Fill(mass,m);
			if(Massa_gen<2&&Massa_gen>1.5)	FitTOF_Pbins -> TemplateD -> Fill(mass,m);	
			if(Massa_gen<4&&Massa_gen>2.5)	FitTOF_Pbins -> TemplateHe-> Fill(mass,m);

		}
	}
	for(int m=0;m<nbinsNaF;m++) { //NaF
		
		if((((int)Cutmask)>>11)==512){
			mass = ((R/BetaRICH)*pow((1-pow(BetaRICH,2)),0.5));
			if(Var2>BetaNaFD[m]&&Var2<=BetaNaFD[m+1]) {
				if(Massa_gen<1&&Massa_gen>0.5)	FitNaF_Dbins -> TemplateP -> Fill(mass,m);	
				if(Massa_gen<2&&Massa_gen>1.5)	FitNaF_Dbins -> TemplateD -> Fill(mass,m);
				if(Massa_gen<4&&Massa_gen>2.5)	FitNaF_Dbins -> TemplateHe-> Fill(mass,m);
			}
			if(Var2>BetaNaFP[m]&&Var2<=BetaNaFP[m+1]) {
				if(Massa_gen<1&&Massa_gen>0.5)	FitNaF_Pbins -> TemplateP -> Fill(mass,m);	
				if(Massa_gen<2&&Massa_gen>1.5)	FitNaF_Pbins -> TemplateD -> Fill(mass,m);
				if(Massa_gen<4&&Massa_gen>2.5)	FitNaF_Pbins -> TemplateHe-> Fill(mass,m);
			}

		}
	}
	for(int m=0;m<nbinsAgl;m++) {	//Agl
		if((((int)Cutmask)>>11)==0){
			mass = ((R/BetaRICH)*pow((1-pow(BetaRICH,2)),0.5));
			if(Var2>BetaAglD[m]&&Var2<=BetaAglD[m+1]) {
				if(Massa_gen<1&&Massa_gen>0.5)	FitAgl_Dbins -> TemplateP -> Fill(mass,m);	
				if(Massa_gen<2&&Massa_gen>1.5)	FitAgl_Dbins -> TemplateD -> Fill(mass,m);
				if(Massa_gen<4&&Massa_gen>2.5)	FitAgl_Dbins -> TemplateHe-> Fill(mass,m);
			}
			if(Var2>BetaAglP[m]&&Var2<=BetaAglP[m+1]) {
				if(Massa_gen<1&&Massa_gen>0.5)	FitAgl_Pbins -> TemplateP -> Fill(mass,m);	
				if(Massa_gen<2&&Massa_gen>1.5)	FitAgl_Pbins -> TemplateD -> Fill(mass,m);
				if(Massa_gen<4&&Massa_gen>2.5)	FitAgl_Pbins -> TemplateHe-> Fill(mass,m);
			}

		}

	}
}			


void DeutonsDATA_Fill(TNtuple *ntupla, int l,int zona){
	int k = ntupla->GetEvent(l);
	if(Beta<=0||R<=0) return;
	float mass = 0;
	if(!(Likcut&&Distcut)) return;	

	for(int m=0;m<nbinsToF;m++){ //TOF
		mass = ((R/Beta)*pow((1-pow(Beta,2)),0.5));
		if(Var>BetaD[m]&&Var<=BetaD[m+1]){ 
			if(R>1.2*Rcutoff)	FitTOF_Dbins -> Data_Prim -> Fill(mass,m);
			((TH3*)FitTOF_Dbins -> Data_Geomag) -> Fill(mass,m,zona);
		}
		if(Var>BetaP[m]&&Var<=BetaP[m+1]) {
			if(R>1.2*Rcutoff)	FitTOF_Pbins -> Data_Prim -> Fill(mass,m);

		}
	}
	for(int m=0;m<nbinsNaF;m++){//NaF
		
		if((((int)Cutmask)>>11)==512){
			mass = ((R/BetaRICH)*pow((1-pow(BetaRICH,2)),0.5));
			if(Var2>BetaNaFD[m]&&Var2<=BetaNaFD[m+1]) {
				if(R>1.2*Rcutoff)	FitNaF_Dbins -> Data_Prim -> Fill(mass,m);	
				((TH3*)FitNaF_Dbins -> Data_Geomag) -> Fill(mass,m,zona);
			}
			if(Var2>BetaNaFP[m]&&Var2<=BetaNaFP[m+1]) {
				if(R>1.2*Rcutoff)	FitNaF_Pbins -> Data_Prim -> Fill(mass,m);	
			}

		}
			}
	for(int m=0;m<nbinsAgl;m++){ //Agl
		if((((int)Cutmask)>>11)==0){
			mass = ((R/BetaRICH)*pow((1-pow(BetaRICH,2)),0.5));
			if(Var2>BetaAglD[m]&&Var2<=BetaAglD[m+1]) {
				if(R>1.2*Rcutoff)	FitAgl_Dbins -> Data_Prim -> Fill(mass,m);	
				((TH3*)FitAgl_Dbins -> Data_Geomag) -> Fill(mass,m,zona);
			}
			if(Var2>BetaAglP[m]&&Var2<=BetaAglP[m+1]) {
				if(R>1.2*Rcutoff)	FitAgl_Pbins -> Data_Prim -> Fill(mass,m);	
			}

		}

	}
	return;
}		


void DeutonsMC_Write(){

	FitTOF_Dbins -> Write();
	FitNaF_Dbins -> Write();
	FitAgl_Dbins -> Write();
	
	FitTOF_Pbins -> Write();
	FitNaF_Pbins -> Write();
	FitAgl_Pbins -> Write();
	return;
}


void DeutonsTemplFits(){
	string nomefile="../Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
	TFile * file1 = TFile::Open(nomefile.c_str(),"READ");

	TemplateFIT * FitTOF_Dbins	= new TemplateFIT(file1,"FitTOF_Dbins",0,3,11);
	TemplateFIT * FitNaF_Dbins	= new TemplateFIT(file1,"FitNaF_Dbins",0,3,11);
	TemplateFIT * FitAgl_Dbins	= new TemplateFIT(file1,"FitAgl_Dbins",0,3,11);

	TemplateFIT * FitTOF_Pbins	= new TemplateFIT(file1,"FitTOF_Pbins",0,3,11);
	TemplateFIT * FitNaF_Pbins	= new TemplateFIT(file1,"FitNaF_Pbins",0,3,11);
	TemplateFIT * FitAgl_Pbins	= new TemplateFIT(file1,"FitAgl_Pbins",0,3,11);

	cout<<"******************** DEUTONS TEMPlATE FITS ************************"<<endl;

	FitTOF_Dbins -> TemplateFits(); 
	FitNaF_Dbins -> TemplateFits();
	FitAgl_Dbins -> TemplateFits();

	FitTOF_Pbins -> TemplateFits();
	FitNaF_Pbins -> TemplateFits();
	FitAgl_Pbins -> TemplateFits();

	cout<<"*** Updating P1 file ****"<<endl;

        nomefile="../Histos/"+mese+"/"+mese+"_"+frac+"_P1.root";
        file1 = TFile::Open(nomefile.c_str(),"UPDATE");

	FitTOF_Dbins -> DCounts -> Write ("D_FluxCounts_TOF");
	FitNaF_Dbins -> DCounts -> Write ("D_FluxCounts_NaF");
        FitAgl_Dbins -> DCounts -> Write ("D_FluxCounts_Agl");
                       
        FitTOF_Pbins -> PCounts -> Write ("P_FluxCounts_TOF");
	FitNaF_Pbins -> PCounts -> Write ("P_FluxCounts_NaF");
        FitAgl_Pbins -> PCounts -> Write ("P_FluxCounts_Agl");

	file1 -> Write();
	file1 -> Close(); 	
	
	TCanvas *c30TOF[12][nbinsbeta];
	TCanvas *c30_NaF[12][nbinsbeta];
	TCanvas *c30_Agl[12][nbinsbeta];
	




	return;


}
/*


TCanvas *c30[12][nbinsbeta];
TCanvas *c30_bis[12][nbinsbeta];
TCanvas *c30_tris[12][nbinsbeta];
TH3F * DCountsgeoTOF = new TH3F("DCountsgeoTOF","DCountsgeoTOF",nbinsbeta,0,nbinsbeta,12,0,12,2,0,2);
TH3F * DCountsgeoNaF = new TH3F("DCountsgeoNaF","DCountsgeoNaF",nbinsbeta,0,nbinsbeta,12,0,12,2,0,2);
TH3F * DCountsgeoAgl = new TH3F("DCountsgeoAgl","DCountsgeoAgl",nbinsbeta,0,nbinsbeta,12,0,12,2,0,2);
TH2F * PCountsTOF = new TH2F("PCountsTOF","PCountsTOF",nbinsbeta,0,nbinsbeta,2,0,2);
TH2F * PCountsNaF = new TH2F("PCountsNaF","PCountsNaF",nbinsbeta,0,nbinsbeta,2,0,2);
TH2F * PCountsAgl = new TH2F("PCountsAgl","PCountsAgl",nbinsbeta,0,nbinsbeta,2,0,2);

void DeutonsTemplFits(TFile * file1){
	TH2F* DTemplatesTOF=(TH2F*)file1->Get("DTemplatesTOF");
	TH2F* PTemplatesTOF=(TH2F*)file1->Get("PTemplatesTOF");
	TH2F* HeTemplatesTOF=(TH2F*)file1->Get("HeTemplatesTOF");
	TH2F* DTemplatesNaF=(TH2F*)file1->Get("DTemplatesNaF");
	TH2F* PTemplatesNaF=(TH2F*)file1->Get("PTemplatesNaF");
	TH2F* HeTemplatesNaF=(TH2F*)file1->Get("HeTemplatesNaF");
	TH2F* DTemplatesAgl=(TH2F*)file1->Get("DTemplatesAgl");
	TH2F* PTemplatesAgl=(TH2F*)file1->Get("PTemplatesAgl");
	TH2F* HeTemplatesAgl=(TH2F*)file1->Get("HeTemplatesAgl");
	TH2F* DTemplatesTOF2=(TH2F*)file1->Get("DTemplatesTOF2");
	TH2F* PTemplatesTOF2=(TH2F*)file1->Get("PTemplatesTOF2");
	TH2F* HeTemplatesTOF2=(TH2F*)file1->Get("HeTemplatesTOF2");
	TH2F* DTemplatesNaF2=(TH2F*)file1->Get("DTemplatesNaF2");
	TH2F* PTemplatesNaF2=(TH2F*)file1->Get("PTemplatesNaF2");
	TH2F* HeTemplatesNaF2=(TH2F*)file1->Get("HeTemplatesNaF2");
	TH2F* DTemplatesAgl2=(TH2F*)file1->Get("DTemplatesAgl2");
	TH2F* PTemplatesAgl2=(TH2F*)file1->Get("PTemplatesAgl2");
	TH2F* HeTemplatesAgl2=(TH2F*)file1->Get("HeTemplatesAgl2");
	TH3F* DhistosgeoTOF=(TH3F*)file1->Get("DhistosgeoTOF");
	TH3F* DhistosgeoNaF=(TH3F*)file1->Get("DhistosgeoNaF");
	TH3F* DhistosgeoAgl=(TH3F*)file1->Get("DhistosgeoAgl");
	TH2F* DhistosTOF=(TH2F*)file1->Get("DhistosTOF");
	TH2F* DhistosNaF=(TH2F*)file1->Get("DhistosNaF");
	TH2F* DhistosAgl=(TH2F*)file1->Get("DhistosAgl");
	TH2F* DhistosTOF2=(TH2F*)file1->Get("DhistosTOF2");
	TH2F* DhistosNaF2=(TH2F*)file1->Get("DhistosNaF2");
	TH2F* DhistosAgl2=(TH2F*)file1->Get("DhistosAgl2");
	
	cout<<"******************** Mass TOF TEMPLATE FITS *******************"<<endl;
	TH1F *DataMassTOF[12][nbinsbeta];
	TH1F *DTemplTOF[nbinsbeta];
	TH1F *PTemplTOF[nbinsbeta];
	TH1F *HeTemplTOF[nbinsbeta];
	TH1F *DTemplTOF2[nbinsbeta];
        TH1F *PTemplTOF2[nbinsbeta];
        TH1F *HeTemplTOF2[nbinsbeta];	
	
	string nome;
	for(int m=0;m<nbinsbeta;m++){
		DTemplTOF[m]=new TH1F("","",50,0,3);
		PTemplTOF[m]=new TH1F("","",50,0,3);
		HeTemplTOF[m]=new TH1F("","",50,0,3);
		DTemplTOF2[m]=new TH1F("","",50,0,3);
                PTemplTOF2[m]=new TH1F("","",50,0,3);
                HeTemplTOF2[m]=new TH1F("","",50,0,3);
		for(int l=0;l<12;l++){
			DataMassTOF[l][m]=new TH1F("","",50,0,3);	
		}	
	}
	for(int m=0;m<nbinsbeta;m++){
		for(int i=0;i<50;i++){
		DTemplTOF[m]->SetBinContent(i+1,DTemplatesTOF->GetBinContent(i+1,m+1));
		PTemplTOF[m]->SetBinContent(i+1,PTemplatesTOF->GetBinContent(i+1,m+1));
		HeTemplTOF[m]->SetBinContent(i+1,HeTemplatesTOF->GetBinContent(i+1,m+1));
		DTemplTOF2[m]->SetBinContent(i+1,DTemplatesTOF2->GetBinContent(i+1,m+1));
                PTemplTOF2[m]->SetBinContent(i+1,PTemplatesTOF2->GetBinContent(i+1,m+1));
                HeTemplTOF2[m]->SetBinContent(i+1,HeTemplatesTOF2->GetBinContent(i+1,m+1));
		for(int l=1;l<11;l++)  DataMassTOF[l][m]->SetBinContent(i+1,DhistosgeoTOF->GetBinContent(l+1,i+1,m+1));
		DataMassTOF[0][m]->SetBinContent(i+1,DhistosTOF2->GetBinContent(i+1,m+1)); 
		DataMassTOF[11][m]->SetBinContent(i+1,DhistosTOF->GetBinContent(i+1,m+1));
		}	
	}
	
	
	TFractionFitter * fitT[nbinsbeta][12]={{NULL}};
	TObjArray *Tpl[nbinsbeta][12]={{NULL}};
	int cut=17;
	int s1[nbinsbeta][12]={{0}};
	float Err[nbinsbeta][12]={{0}};
	TH1F *PTemplTOFW[nbinsbeta][12];
	TH1F *DTemplTOFW[nbinsbeta][12];
	TH1F *HeTemplTOFW[nbinsbeta][12];
	for(int i=0;i<12;i++) for(int m=0;m<nbinsbeta;m++)c30[i][m]=new TCanvas(); 	
		bool He=true;
	for(int i=0;i<12;i++) for(int m=0;m<nbinsbeta;m++) {
		c30[i][m]->cd();
		PTemplTOFW[m][i]=new TH1F("","",50,0,3);
		DTemplTOFW[m][i]=new TH1F("","",50,0,3);
		HeTemplTOFW[m][i]=new TH1F("","",50,-0,3);
		PTemplTOFW[m][i]->SetFillStyle(3001);
		DTemplTOFW[m][i]->SetFillStyle(3001);
		HeTemplTOFW[m][i]->SetFillStyle(3001);
		PTemplTOFW[m][i]->SetFillColor(2);
		DTemplTOFW[m][i]->SetFillColor(4);
		HeTemplTOFW[m][i]->SetFillColor(3);
		gPad->SetLogy();

		TH1F *Result;
		THStack *Stack=new THStack("","");
		Tpl[m][i] = new TObjArray(2);
		if(i!=0){
			cut=17;
			Tpl[m][i]->Add(PTemplTOF[m]);
			Tpl[m][i]->Add(DTemplTOF[m]);
			if(He) Tpl[m][i]->Add(HeTemplTOF[m]);
		}
		else{
			cut=0;
			Tpl[m][i]->Add(PTemplTOF2[m]);
                        Tpl[m][i]->Add(DTemplTOF2[m]);
                        if(He) Tpl[m][i]->Add(HeTemplTOF2[m]);
		}
		fitT[m][i] = new TFractionFitter( DataMassTOF[i][m], Tpl[m][i],"q");
		if(He){
			fitT[m][i]->Constrain(0,0.0,1.0);
			fitT[m][i]->Constrain(1,0.0,1.0);
			if(i==0) fitT[m][i]->Constrain(1,0.0,0.01);
			fitT[m][i]->Constrain(2,0.0,0.01);
			fitT[m][i]->SetRangeX(cut,50);
		}
		else {
			fitT[m][i]->Constrain(0,0,1);
			fitT[m][i]->Constrain(1,0,1);
			if(i==0) fitT[m][i]->Constrain(1,0.0,0.1);
			fitT[m][i]->SetRangeX(cut,50);
		}	
		if(DataMassTOF[i][m]->Integral(0,50)>50) s1[m][i]=fitT[m][i]->Fit();
		else s1[m][i]=1;
		if(i!=0){
			cut=17;
			for(int z=0;z<5;z++){
				if(s1[m][i]==-1){
					cut+=1;
					fitT[m][i]->SetRangeX(cut,50);
					s1[m][i]=fitT[m][i]->Fit();
				}
				if(s1[m][i]==0) break;
			}	
		}
		double w1,w2,w3=0;
		double e1,e2,e3=0;
		if(s1[m][i]==0){
			fitT[m][i]->GetResult(0,w1,e1);
			fitT[m][i]->GetResult(1,w2,e2);
			if(He) fitT[m][i]->GetResult(2,w3,e3);
			Result = (TH1F*) fitT[m][i]->GetPlot();
			if(i!=0){
				float itot= Result->Integral();
				float i1 = PTemplTOF[m]->Integral(cut,50);
				float i2 = DTemplTOF[m]->Integral(cut,50);
				float i3=0;
				if(He) i3 = HeTemplTOF[m]->Integral();
				if(i1>0) for(int j=0; j<PTemplTOF[m]->GetNbinsX();j++) PTemplTOFW[m][i]->SetBinContent(j,w1*PTemplTOF[m]->GetBinContent(j)/i1*itot);
				if(i2>0) for(int j=0; j<DTemplTOF[m]->GetNbinsX();j++) DTemplTOFW[m][i]->SetBinContent(j,w2*DTemplTOF[m]->GetBinContent(j)/i2*itot);
				if(i3>0) for(int j=0; j<HeTemplTOF[m]->GetNbinsX();j++) HeTemplTOFW[m][i]->SetBinContent(j,w3*HeTemplTOF[m]->GetBinContent(j)/i3*itot);
			}
			else{
				float itot= Result->Integral();
                                float i1 = PTemplTOF2[m]->Integral(cut,50);
                                float i2 = DTemplTOF2[m]->Integral(cut,50);
                                float i3=0;
                                if(He) i3 = HeTemplTOF2[m]->Integral();
                                if(i1>0) for(int j=0; j<PTemplTOF2[m]->GetNbinsX();j++) PTemplTOFW[m][i]->SetBinContent(j,w1*PTemplTOF2[m]->GetBinContent(j)/i1*itot);
                                if(i2>0) for(int j=0; j<DTemplTOF2[m]->GetNbinsX();j++) DTemplTOFW[m][i]->SetBinContent(j,w2*DTemplTOF2[m]->GetBinContent(j)/i2*itot);
                                if(i3>0) for(int j=0; j<HeTemplTOF2[m]->GetNbinsX();j++) HeTemplTOFW[m][i]->SetBinContent(j,w3*HeTemplTOF2[m]->GetBinContent(j)/i3*itot);	
			}
			Stack->Add(DTemplTOFW[m][i]);
			Stack->Add(PTemplTOFW[m][i]);
			if(He) Stack->Add(HeTemplTOFW[m][i]);
			DataMassTOF[i][m]->SetMarkerStyle(8);
			Stack->Draw();
			DataMassTOF[i][m]->Draw("epsame");
			Result->SetLineColor(5);
			Result->SetLineWidth(2);
			Result->Draw("SAME");
			float Cov01=fitT[m][i]->GetFitter()->GetCovarianceMatrixElement(0,1);
			float Cov02=fitT[m][i]->GetFitter()->GetCovarianceMatrixElement(0,2);
			float Cov12=fitT[m][i]->GetFitter()->GetCovarianceMatrixElement(1,2);
			float Sigma=pow((pow(e2/w2,2)+pow(e1/w1,2))/2,0.5);
			Err[m][i]= Sigma*DTemplTOFW[m][i]->Integral();
			DCountsgeoTOF->SetBinContent(m+1,i+1,0,DTemplTOFW[m][i]->Integral(0,50));
			DCountsgeoTOF->SetBinContent(m+1,i+1,1,Err[m][i]);	
			if(i==0){
				PCountsTOF->SetBinContent(m+1,0,PTemplTOFW[m][0]->Integral(0,50));
                        	PCountsTOF->SetBinContent(m+1,1,Err[m][i]);
			}
		}
		else{
			DataMassTOF[i][m]->SetMarkerStyle(8);
			DataMassTOF[i][m]->Draw("ep");
			 if(i!=0){PTemplTOF[m]->Draw("same");
                                DTemplTOF[m]->Draw("same");
                                if(He) HeTemplTOF[m]->Draw("same");
                                }
                        if(i==0) {PTemplTOF2[m]->Draw("same");
                                 DTemplTOF2[m]->Draw("same");
                                if(He) HeTemplTOF2[m]->Draw("same");
                                }
			if(i==0){
				PCountsTOF->SetBinContent(m+1,0,DataMassTOF[0][m]->Integral(0,21));
				PCountsTOF->SetBinContent(m+1,1,DataMassTOF[0][m]->Integral(0,21)/2);
			}
		}

	}

	cout<<"******************** MASS NaF TEMPLATE FITS *******************"<<endl;
        TH1F *DataMassNaF[12][nbinsbeta];
        TH1F *DTemplNaF[nbinsbeta];
        TH1F *PTemplNaF[nbinsbeta];
        TH1F *HeTemplNaF[nbinsbeta];
	TH1F *DTemplNaF2[nbinsbeta];
        TH1F *PTemplNaF2[nbinsbeta];
        TH1F *HeTemplNaF2[nbinsbeta];
        for(int m=0;m<nbinsbeta;m++){
                DTemplNaF[m]=new TH1F("","",50,0,3);
                PTemplNaF[m]=new TH1F("","",50,0,3);
                HeTemplNaF[m]=new TH1F("","",50,0,3);
                DTemplNaF2[m]=new TH1F("","",50,0,3);
                PTemplNaF2[m]=new TH1F("","",50,0,3);
                HeTemplNaF2[m]=new TH1F("","",50,0,3);
		for(int l=0;l<12;l++){
                        DataMassNaF[l][m]=new TH1F("","",50,0,3);
                }
        }
        for(int m=0;m<nbinsbeta;m++){
                for(int i=0;i<50;i++){
                DTemplNaF[m]->SetBinContent(i+1,DTemplatesNaF->GetBinContent(i+1,m+1));
                PTemplNaF[m]->SetBinContent(i+1,PTemplatesNaF->GetBinContent(i+1,m+1));
                HeTemplNaF[m]->SetBinContent(i+1,HeTemplatesNaF->GetBinContent(i+1,m+1));
		DTemplNaF2[m]->SetBinContent(i+1,DTemplatesNaF2->GetBinContent(i+1,m+1));
                PTemplNaF2[m]->SetBinContent(i+1,PTemplatesNaF2->GetBinContent(i+1,m+1));
                HeTemplNaF2[m]->SetBinContent(i+1,HeTemplatesNaF2->GetBinContent(i+1,m+1));
                for(int l=1;l<11;l++)  DataMassNaF[l][m]->SetBinContent(i+1,DhistosgeoNaF->GetBinContent(l+1,i+1,m+1));
                DataMassNaF[11][m]->SetBinContent(i+1,DhistosNaF->GetBinContent(i+1,m+1));
		DataMassNaF[0][m]->SetBinContent(i+1,DhistosNaF2->GetBinContent(i+1,m+1));
		}
        }
	TFractionFitter * fitTNaF[nbinsbeta][12]={{NULL}};
	TObjArray *TplNaF[nbinsbeta][12]={{NULL}};
	int cutNaF=17;
	int s1NaF[nbinsbeta][12]={{0}};
	float ErrNaF[nbinsbeta][12]={{0}};
	TH1F *PTemplNaFW[nbinsbeta][12];
	TH1F *DTemplNaFW[nbinsbeta][12];
	TH1F *HeTemplNaFW[nbinsbeta][12];
	for(int i=0;i<12;i++) for(int m=0;m<nbinsbeta;m++) c30_bis[i][m]=new TCanvas(); 	
		He=false;
	for(int i=0;i<12;i++) for(int m=0;m<nbinsbeta;m++) {
		c30_bis[i][m]->cd();
		PTemplNaFW[m][i]=new TH1F("","",50,0,3);
		DTemplNaFW[m][i]=new TH1F("","",50,0,3);
		HeTemplNaFW[m][i]=new TH1F("","",50,0,3);
		PTemplNaFW[m][i]->SetFillStyle(3001);
		DTemplNaFW[m][i]->SetFillStyle(3001);
		HeTemplNaFW[m][i]->SetFillStyle(3001);
		PTemplNaFW[m][i]->SetFillColor(2);
		DTemplNaFW[m][i]->SetFillColor(4);
		HeTemplNaFW[m][i]->SetFillColor(3);
		gPad->SetLogy();

		TH1F *Result;
		THStack *Stack=new THStack("","");
		TplNaF[m][i] = new TObjArray(2);
		if(i!=0){
			cutNaF=17;
			TplNaF[m][i]->Add(PTemplNaF[m]);
			TplNaF[m][i]->Add(DTemplNaF[m]);
			if(He) TplNaF[m][i]->Add(HeTemplNaF[m]);
		}
		else{
			cutNaF=0;		
			TplNaF[m][i]->Add(PTemplNaF2[m]);
               		TplNaF[m][i]->Add(DTemplNaF2[m]);
                	if(He) TplNaF[m][i]->Add(HeTemplNaF2[m]);
		}
		fitTNaF[m][i] = new TFractionFitter( DataMassNaF[i][m], TplNaF[m][i],"q");
		if(He){
			fitTNaF[m][i]->Constrain(0,0.5,1.0);
			fitTNaF[m][i]->Constrain(1,0.0,1.0);
			if(i==0) fitTNaF[m][i]->Constrain(1,0.0,0.1);
			fitTNaF[m][i]->Constrain(2,0.0,1.0);
		}
		else {
			fitTNaF[m][i]->Constrain(0,0,1);
			fitTNaF[m][i]->Constrain(1,0,1);
			if(i==0) fitTNaF[m][i]->Constrain(1,0.0,0.1);
			fitTNaF[m][i]->SetRangeX(cutNaF,50);
		}
		if((DataMassNaF[i][m]->Integral(25,50)>50&&i!=0)||i==0) s1NaF[m][i]=fitTNaF[m][i]->Fit();
		else s1NaF[m][i]=1;
		double w1,w2,w3=0;
		double e1,e2,e3=0;
		if(s1NaF[m][i]==0){
			fitTNaF[m][i]->GetResult(0,w1,e1);
			fitTNaF[m][i]->GetResult(1,w2,e2);
			if(He) fitTNaF[m][i]->GetResult(2,w3,e3);
			Result = (TH1F*) fitTNaF[m][i]->GetPlot();
			if(i!=0){
				float itot= Result->Integral();
				float i1 = PTemplNaF[m]->Integral(cutNaF,50);
				float i2 = DTemplNaF[m]->Integral(cutNaF,50);
				float i3=0;
				if(He) i3 = HeTemplNaF[m]->Integral();
				if(i1>0) for(int j=0; j<PTemplNaF[m]->GetNbinsX();j++) PTemplNaFW[m][i]->SetBinContent(j,w1*PTemplNaF[m]->GetBinContent(j)/i1*itot);
				if(i2>0) for(int j=0; j<DTemplNaF[m]->GetNbinsX();j++) DTemplNaFW[m][i]->SetBinContent(j,w2*DTemplNaF[m]->GetBinContent(j)/i2*itot);
				if(i3>0) for(int j=0; j<HeTemplNaF[m]->GetNbinsX();j++) HeTemplNaFW[m][i]->SetBinContent(j,w3*HeTemplNaF[m]->GetBinContent(j)/i3*itot);
			}
			else{
				float itot= Result->Integral();
				float i1 = PTemplNaF2[m]->Integral(cutNaF,50);
                                float i2 = DTemplNaF2[m]->Integral(cutNaF,50);
                                float i3=0;
                                if(He) i3 = HeTemplNaF2[m]->Integral();
				if(i1>0) for(int j=0; j<PTemplNaF2[m]->GetNbinsX();j++) PTemplNaFW[m][i]->SetBinContent(j,w1*PTemplNaF2[m]->GetBinContent(j)/i1*itot);
                                if(i2>0) for(int j=0; j<DTemplNaF2[m]->GetNbinsX();j++) DTemplNaFW[m][i]->SetBinContent(j,w2*DTemplNaF2[m]->GetBinContent(j)/i2*itot);
                                if(i3>0) for(int j=0; j<HeTemplNaF2[m]->GetNbinsX();j++) HeTemplNaFW[m][i]->SetBinContent(j,w3*HeTemplNaF2[m]->GetBinContent(j)/i3*itot);
			}
			Stack->Add(PTemplNaFW[m][i]);
			Stack->Add(DTemplNaFW[m][i]);
			if(He) Stack->Add(HeTemplNaFW[m][i]);
			DataMassNaF[i][m]->SetMarkerStyle(8);
			Stack->Draw();
			DataMassNaF[i][m]->Draw("epsame");
			Result->SetLineColor(5);
			Result->SetLineWidth(2);
			Result->Draw("SAME");
			
			float Cov01=fitTNaF[m][i]->GetFitter()->GetCovarianceMatrixElement(0,1);
			float Cov02=fitTNaF[m][i]->GetFitter()->GetCovarianceMatrixElement(0,2);
			float Cov12=fitTNaF[m][i]->GetFitter()->GetCovarianceMatrixElement(1,2);
			float Sigma=pow((pow(e2/w2,2)+pow(e1/w1,2))/2,0.5);
			ErrNaF[m][i]= Sigma*DTemplNaFW[m][i]->Integral();
			DCountsgeoNaF->SetBinContent(m+1,i+1,0,DTemplNaFW[m][i]->Integral(0,50));
			DCountsgeoNaF->SetBinContent(m+1,i+1,1,ErrNaF[m][i]);	
			if(i==0){
                                PCountsNaF->SetBinContent(m+1,0,PTemplNaFW[m][0]->Integral(0,50));
                                PCountsNaF->SetBinContent(m+1,1,Err[m][i]);
                        }

		}
		else{
			DataMassNaF[i][m]->SetMarkerStyle(8);
			DataMassNaF[i][m]->Draw("ep");
			 if(i!=0){PTemplNaF[m]->Draw("same");
                                DTemplNaF[m]->Draw("same");
                                if(He) HeTemplNaF[m]->Draw("same");
                                }
                        if(i==0) {PTemplNaF2[m]->Draw("same");
                                 DTemplNaF2[m]->Draw("same");
                                if(He) HeTemplNaF2[m]->Draw("same");
                                PCountsNaF->SetBinContent(m+1,0,DataMassNaF[i][m]->Integral(5,25));
				}

		}

	}
	cout<<"******************** MASS Agl TEMPLATE FITS *******************"<<endl;
        TH1F *DataMassAgl[12][nbinsbeta];
        TH1F *DTemplAgl[nbinsbeta];
        TH1F *PTemplAgl[nbinsbeta];
        TH1F *HeTemplAgl[nbinsbeta];
	TH1F *DTemplAgl2[nbinsbeta];
        TH1F *PTemplAgl2[nbinsbeta];
        TH1F *HeTemplAgl2[nbinsbeta];
        for(int m=0;m<nbinsbeta;m++){
                DTemplAgl[m]=new TH1F("","",50,0,3);
                PTemplAgl[m]=new TH1F("","",50,0,3);
                HeTemplAgl[m]=new TH1F("","",50,0,3);
		DTemplAgl2[m]=new TH1F("","",50,0,3);
                PTemplAgl2[m]=new TH1F("","",50,0,3);
                HeTemplAgl2[m]=new TH1F("","",50,0,3);
                for(int l=0;l<12;l++){
                        DataMassAgl[l][m]=new TH1F("","",50,0,3);
                }
        }
        for(int m=0;m<nbinsbeta;m++){
                for(int i=0;i<50;i++){
                DTemplAgl[m]->SetBinContent(i+1,DTemplatesAgl->GetBinContent(i+1,m+1));
                PTemplAgl[m]->SetBinContent(i+1,PTemplatesAgl->GetBinContent(i+1,m+1));
                HeTemplAgl[m]->SetBinContent(i+1,HeTemplatesAgl->GetBinContent(i+1,m+1));
                DTemplAgl2[m]->SetBinContent(i+1,DTemplatesAgl2->GetBinContent(i+1,m+1));
                PTemplAgl2[m]->SetBinContent(i+1,PTemplatesAgl2->GetBinContent(i+1,m+1));
                HeTemplAgl2[m]->SetBinContent(i+1,HeTemplatesAgl2->GetBinContent(i+1,m+1));
		for(int l=1;l<11;l++)  DataMassAgl[l][m]->SetBinContent(i+1,DhistosgeoAgl->GetBinContent(l+1,i+1,m+1));
                DataMassAgl[11][m]->SetBinContent(i+1,DhistosAgl->GetBinContent(i+1,m+1));
		DataMassAgl[0][m]->SetBinContent(i+1,DhistosAgl2->GetBinContent(i+1,m+1));
		}
        }

        TFractionFitter * fitTAgl[nbinsbeta][12]={{NULL}};
        TObjArray *TplAgl[nbinsbeta][12]={{NULL}};
        int cutAgl=16;
        int s1Agl[nbinsbeta][12]={{0}};
        float ErrAgl[nbinsbeta][12]={{0}};
        TH1F *PTemplAglW[nbinsbeta][12];
        TH1F *DTemplAglW[nbinsbeta][12];
        TH1F *HeTemplAglW[nbinsbeta][12];
        for(int i=0;i<12;i++) for(int m=0;m<nbinsbeta;m++) c30_tris[i][m]=new TCanvas();
                He=true;
        for(int i=0;i<12;i++) for(int m=0;m<nbinsbeta;m++) {
                c30_tris[i][m]->cd();
                cutAgl=16;
		PTemplAglW[m][i]=new TH1F("","",50,0,3);
                DTemplAglW[m][i]=new TH1F("","",50,0,3);
                HeTemplAglW[m][i]=new TH1F("","",50,0,3);
                PTemplAglW[m][i]->SetFillStyle(3001);
                DTemplAglW[m][i]->SetFillStyle(3001);
                HeTemplAglW[m][i]->SetFillStyle(3001);
                PTemplAglW[m][i]->SetFillColor(2);
                DTemplAglW[m][i]->SetFillColor(4);
                HeTemplAglW[m][i]->SetFillColor(3);
                gPad->SetLogy();

                TH1F *Result;
                THStack *Stack=new THStack("","");
                TplAgl[m][i] = new TObjArray(2);
                if(i!=0){
			cutAgl=16;
			TplAgl[m][i]->Add(PTemplAgl[m]);
                	TplAgl[m][i]->Add(DTemplAgl[m]);
                	if(He) TplAgl[m][i]->Add(HeTemplAgl[m]);
                }
		else{
			cutAgl=0;
			TplAgl[m][i]->Add(PTemplAgl2[m]);
                        TplAgl[m][i]->Add(DTemplAgl2[m]);
                        if(He) TplAgl[m][i]->Add(HeTemplAgl2[m]);
		}
		fitTAgl[m][i] = new TFractionFitter( DataMassAgl[i][m], TplAgl[m][i],"q");
                if(He){
                        fitTAgl[m][i]->Constrain(0,0.5,1.0);
                        fitTAgl[m][i]->Constrain(1,0.0,1.0);
                        if(i==0) fitTAgl[m][i]->Constrain(1,0.0,0.1);
			fitTAgl[m][i]->Constrain(2,0.0,0.01);
                	fitTAgl[m][i]->SetRangeX(cutAgl,50);
		}
                else {
                        fitTAgl[m][i]->Constrain(0,0,1);
                        fitTAgl[m][i]->Constrain(1,0,1);
                        if(i==0) fitTAgl[m][i]->Constrain(1,0.0,0.1);
			fitTAgl[m][i]->SetRangeX(cutAgl,50);
                }
                if((DataMassAgl[i][m]->Integral(25,50)>50&&i!=0)||i==0) s1Agl[m][i]=1;//fitTAgl[m][i]->Fit();
                else s1Agl[m][i]=1;
		for(int z=0;z<5;z++){
			if(s1Agl[m][i]==-1){
				cutAgl+=1;
				fitTAgl[m][i]->SetRangeX(cutAgl,50);
				s1Agl[m][i]=fitTAgl[m][i]->Fit();
			}
			if(s1Agl[m][i]==0) break;
		}
		double w1,w2,w3=0;
                double e1,e2,e3=0;
                if(s1Agl[m][i]==0){
                        fitTAgl[m][i]->GetResult(0,w1,e1);
                        fitTAgl[m][i]->GetResult(1,w2,e2);
                        if(He) fitTAgl[m][i]->GetResult(2,w3,e3);
                        Result = (TH1F*) fitTAgl[m][i]->GetPlot();
                        if(i!=0){
				float itot= Result->Integral();
                        	float i1 = PTemplAgl[m]->Integral(cutAgl,50);
                        	float i2 = DTemplAgl[m]->Integral(cutAgl,50);
                        	float i3=0;
                        	if(He) i3 = HeTemplAgl[m]->Integral();
                        	if(i1>0) for(int j=0; j<PTemplAgl[m]->GetNbinsX();j++) PTemplAglW[m][i]->SetBinContent(j,w1*PTemplAgl[m]->GetBinContent(j)/i1*itot);
                        	if(i2>0) for(int j=0; j<DTemplAgl[m]->GetNbinsX();j++) DTemplAglW[m][i]->SetBinContent(j,w2*DTemplAgl[m]->GetBinContent(j)/i2*itot);
                        	if(i3>0) for(int j=0; j<HeTemplAgl[m]->GetNbinsX();j++) HeTemplAglW[m][i]->SetBinContent(j,w3*HeTemplAgl[m]->GetBinContent(j)/i3*itot);
                        }
			else{
				float itot= Result->Integral();
                                float i1 = PTemplAgl2[m]->Integral(cutAgl,50);
                                float i2 = DTemplAgl2[m]->Integral(cutAgl,50);
                                float i3=0;
                                if(He) i3 = HeTemplAgl2[m]->Integral();
                                if(i1>0) for(int j=0; j<PTemplAgl[m]->GetNbinsX();j++) PTemplAglW[m][i]->SetBinContent(j,w1*PTemplAgl2[m]->GetBinContent(j)/i1*itot);
                                if(i2>0) for(int j=0; j<DTemplAgl[m]->GetNbinsX();j++) DTemplAglW[m][i]->SetBinContent(j,w2*DTemplAgl2[m]->GetBinContent(j)/i2*itot);
                                if(i3>0) for(int j=0; j<HeTemplAgl[m]->GetNbinsX();j++) HeTemplAglW[m][i]->SetBinContent(j,w3*HeTemplAgl2[m]->GetBinContent(j)/i3*itot);
			}
			Stack->Add(PTemplAglW[m][i]);
                        Stack->Add(DTemplAglW[m][i]);
                        if(He) Stack->Add(HeTemplAglW[m][i]);
                        DataMassAgl[i][m]->SetMarkerStyle(8);
                        Stack->Draw();
                        DataMassAgl[i][m]->Draw("epsame");
                        Result->SetLineColor(5);
                        Result->SetLineWidth(2);
                        Result->Draw("SAME");
                        float Cov01=fitTAgl[m][i]->GetFitter()->GetCovarianceMatrixElement(0,1);
                        float Cov02=fitTAgl[m][i]->GetFitter()->GetCovarianceMatrixElement(0,2);
                        float Cov12=fitTAgl[m][i]->GetFitter()->GetCovarianceMatrixElement(1,2);
                        float Sigma=pow((pow(e2/w2,2)+pow(e1/w1,2))/2,0.5);
                        ErrAgl[m][i]= Sigma*DTemplAglW[m][i]->Integral();
                        DCountsgeoAgl->SetBinContent(m+1,i+1,0,DTemplAglW[m][i]->Integral(0,50));
			DCountsgeoAgl->SetBinContent(m+1,i+1,1,ErrAgl[m][i]);
                	if(i==0){
                                PCountsAgl->SetBinContent(m+1,0,PTemplAglW[m][0]->Integral(0,50));
                                PCountsAgl->SetBinContent(m+1,1,Err[m][i]);
                        }

		}
                else{
                        DataMassAgl[i][m]->SetMarkerStyle(8);
                        DataMassAgl[i][m]->Draw("ep");
                        if(i!=0){PTemplAgl[m]->Draw("same");
                        	DTemplAgl[m]->Draw("same");
                        	if(He) HeTemplAgl[m]->Draw("same");
				}
			if(i==0) {PTemplAgl2[m]->Draw("same");
                                 DTemplAgl2[m]->Draw("same");
                                if(He) HeTemplAgl2[m]->Draw("same");
				PCountsAgl->SetBinContent(m+1,0,DataMassAgl[i][m]->Integral(0,25));
				}
                }

        }

	cout<<"********** TEMPL: FIT RESULTS ********************"<<endl;
	cout<<"**** TOF ********"<<endl;
	for(int i=0;i<12;i++){
		for(int m=0;m<nbinsbeta;m++) cout<<s1[m][i]<<" ";
		cout<<endl;
	}
	cout<<"**** NaF ********"<<endl;
        for(int i=0;i<12;i++){
                for(int m=0;m<nbinsbeta;m++) cout<<s1NaF[m][i]<<" ";
                cout<<endl;
        }
	cout<<"**** Agl ********"<<endl;
        for(int i=0;i<12;i++){
                for(int m=0;m<nbinsbeta;m++) cout<<s1Agl[m][i]<<" ";
                cout<<endl;
        }

}*/
