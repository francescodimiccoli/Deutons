TCanvas *c30[12][18];
TCanvas *c30_bis[12][18];
TCanvas *c30_tris[12][18];

TH2F *DTemplatesTOF_Dist=new TH2F("DTemplatesTOF_Dist","DTemplatesTOF_Dist",50,-1,1,18,0,18);
TH2F *PTemplatesTOF_Dist=new TH2F("PTemplatesTOF_Dist","PTemplatesTOF_Dist",50,-1,1,18,0,18);
TH2F *HeTemplatesTOF_Dist=new TH2F("HeTemplatesTOF_Dist","HeTemplatesTOF_Dist",50,-1,1,18,0,18);

TH2F *DTemplatesNaF=new TH2F("DTemplatesNaF","DTemplatesNaF",50,0,3,18,0,18);
TH2F *PTemplatesNaF=new TH2F("PTemplatesNaF","PTemplatesNaF",50,0,3,18,0,18);
TH2F *HeTemplatesNaF=new TH2F("HeTemplatesNaF","HeTemplatesNaF",50,0,3,18,0,18);

TH2F *DTemplatesAgl=new TH2F("DTemplatesAgl","DTemplatesAgl",50,0,3,18,0,18);
TH2F *PTemplatesAgl=new TH2F("PTemplatesAgl","PTemplatesAgl",50,0,3,18,0,18);
TH2F *HeTemplatesAgl=new TH2F("HeTemplatesAgl","HeTemplatesAgl",50,0,3,18,0,18);

TH2F *DTemplatesTOF_Dist2=new TH2F("DTemplatesTOF_Dist2","DTemplatesTOF_Dist2",50,-1,1,18,0,18);
TH2F *PTemplatesTOF_Dist2=new TH2F("PTemplatesTOF_Dist2","PTemplatesTOF_Dist2",50,-1,1,18,0,18);
TH2F *HeTemplatesTOF_Dist2=new TH2F("HeTemplatesTOF_Dist2","HeTemplatesTOF_Dist2",50,-1,1,18,0,18);

TH2F *DTemplatesNaF2=new TH2F("DTemplatesNaF2","DTemplatesNaF2",50,0,3,18,0,18);
TH2F *PTemplatesNaF2=new TH2F("PTemplatesNaF2","PTemplatesNaF2",50,0,3,18,0,18);
TH2F *HeTemplatesNaF2=new TH2F("HeTemplatesNaF2","HeTemplatesNaF2",50,0,3,18,0,18);

TH2F *DTemplatesAgl2=new TH2F("DTemplatesAgl2","DTemplatesAgl2",50,0,3,18,0,18);
TH2F *PTemplatesAgl2=new TH2F("PTemplatesAgl2","PTemplatesAgl2",50,0,3,18,0,18);
TH2F *HeTemplatesAgl2=new TH2F("HeTemplatesAgl2","HeTemplatesAgl2",50,0,3,18,0,18);

TH3F *DhistosgeoTOF_Dist=new TH3F("DhistosgeoTOF_Dist","DhistosgeoTOF_Dist",11,0,11,50,-1,1,18,0,18);
TH3F *DhistosgeoNaF=new TH3F("DhistosgeoNaF","DhistosgeoNaF",11,0,11,50,0,3,18,0,18);
TH3F *DhistosgeoAgl=new TH3F("DhistosgeoAgl","DhistosgeoAgl",11,0,11,50,0,3,18,0,18);

TH2F *DhistosTOF_Dist=new TH2F("DhistosTOF_Dist","DhistosTOF_Dist",50,-1,1,18,0,18);
TH2F *DhistosNaF=new TH2F("DhistosNaF","DhistosNaF",50,0,3,18,0,18);
TH2F *DhistosAgl=new TH2F("DhistosAgl","DhistosAgl",50,0,3,18,0,18);

TH2F *DhistosTOF_Dist2=new TH2F("DhistosTOF_Dist2","DhistosTOF_Dist2",50,-1,1,18,0,18);
TH2F *DhistosNaF2=new TH2F("DhistosNaF2","DhistosNaF2",50,0,3,18,0,18);
TH2F *DhistosAgl2=new TH2F("DhistosAgl2","DhistosAgl2",50,0,3,18,0,18);

TH3F * DCountsgeoTOF_Dist = new TH3F("DCountsgeoTOF_Dist","DCountsgeoTOF_Dist",18,0,18,12,0,12,2,0,2);
TH3F * DCountsgeoNaF = new TH3F("DCountsgeoNaF","DCountsgeoNaF",18,0,18,12,0,12,2,0,2);
TH3F * DCountsgeoAgl = new TH3F("DCountsgeoAgl","DCountsgeoAgl",18,0,18,12,0,12,2,0,2);

TH2F * PCountsTOF_Dist = new TH2F("PCountsTOF_Dist","PCountsTOF_Dist",18,0,18,2,0,2);
TH2F * PCountsNaF = new TH2F("PCountsNaF","PCountsNaF",18,0,18,2,0,2);
TH2F * PCountsAgl = new TH2F("PCountsAgl","PCountsAgl",18,0,18,2,0,2);


void DeutonsMC_Fill(TNtuple *ntupla, int l){
        int k = ntupla->GetEvent(l);
	if(Beta<=0||R<=0) return;
        float mm=0;
	if(Likcut&&R<4/*&&Distcut*/){
		if(Massa_gen<1&&Massa_gen>0.5) 
			for(int m=0;m<18;m++)  {if(Var>BetaD[m]&&Var<=BetaD[m+1]) PTemplatesTOF_Dist->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),m); ;
						if(Beta>Betabins[m]&&Beta<=Betabins[m+1]) PTemplatesTOF_Dist2->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),m);} 
		if(Massa_gen<2&&Massa_gen>1.5)
                        for(int m=0;m<18;m++)  {if(Var>BetaD[m]&&Var<=BetaD[m+1]) DTemplatesTOF_Dist->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),m);	
						if(Beta>Betabins[m]&&Beta<=Betabins[m+1]) DTemplatesTOF_Dist2->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),m);}
		 if(Massa_gen<4&&Massa_gen>2.5)
                        for(int m=0;m<18;m++)  {if(Var>BetaD[m]&&Var<=BetaD[m+1]) HeTemplatesTOF_Dist->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),m);
						if(Beta>Betabins[m]&&Beta<=Betabins[m+1]) HeTemplatesTOF_Dist2->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),m);}
	}
	if(Likcut&&R<20&&Distcut){
		mm=1;
		//if(Dist5D>distcut) mm=2;
		if(Massa_gen<1&&Massa_gen>0.5)
			for(int m=0;m<18;m++) 
			{if((((int)Cutmask)>>11)==512&&Var2>BetaNaFD[m]&&Var2<=BetaNaFD[m+1]) PTemplatesNaF->Fill(mm*((R/BetaRICH)*pow((1-pow(BetaRICH,2)),0.5)),m);
				if((((int)Cutmask)>>11)==512&&BetaRICH>BetabinsNaF[m]&&BetaRICH<=BetabinsNaF[m+1]) PTemplatesNaF2->Fill(mm*((R/BetaRICH)*pow((1-pow(BetaRICH,2)),0.5)),m);}	
		if(Massa_gen<2&&Massa_gen>1.5)
			for(int m=0;m<18;m++) 
			{if((((int)Cutmask)>>11)==512&&Var2>BetaNaFD[m]&&Var2<=BetaNaFD[m+1]) DTemplatesNaF->Fill(mm*((R/BetaRICH)*pow((1-pow(BetaRICH,2)),0.5)),m);
				if((((int)Cutmask)>>11)==512&&BetaRICH>BetabinsNaF[m]&&BetaRICH<=BetabinsNaF[m+1]) DTemplatesNaF2->Fill(mm*((R/BetaRICH)*pow((1-pow(BetaRICH,2)),0.5)),m);}	
		if(Massa_gen<4&&Massa_gen>2.5)
			for(int m=0;m<18;m++) 
			{if((((int)Cutmask)>>11)==512&&Var2>BetaNaFD[m]&&Var2<=BetaNaFD[m+1]) HeTemplatesNaF->Fill(mm*((R/BetaRICH)*pow((1-pow(BetaRICH,2)),0.5)),m);
				if((((int)Cutmask)>>11)==512&&BetaRICH>BetabinsNaF[m]&&BetaRICH<=BetabinsNaF[m+1]) HeTemplatesNaF2->Fill(mm*((R/BetaRICH)*pow((1-pow(BetaRICH,2)),0.5)),m);}
	}
	if(Likcut&&R<40&&Distcut){
		mm=1;
                //if(Dist5D>distcut) mm=2;
		if(Massa_gen<1&&Massa_gen>0.5)
			for(int m=0;m<18;m++) 
			{if((((int)Cutmask)>>11)==0&&Var2>BetaAglD[m]&&Var2<=BetaAglD[m+1]) {PTemplatesAgl->Fill(mm*((R/BetaRICH)*pow((1-pow(BetaRICH,2)),0.5)),m);}
				if((((int)Cutmask)>>11)==0&&BetaRICH>BetabinsAgl[m]&&BetaRICH<=BetabinsAgl[m+1]) PTemplatesAgl2->Fill(mm*((R/BetaRICH)*pow((1-pow(BetaRICH,2)),0.5)),m);}
		if(Massa_gen<2&&Massa_gen>1.5)
			for(int m=0;m<18;m++) 
			{if((((int)Cutmask)>>11)==0&&Var2>BetaAglD[m]&&Var2<=BetaAglD[m+1]) DTemplatesAgl->Fill(mm*((R/BetaRICH)*pow((1-pow(BetaRICH,2)),0.5)),m);
				if((((int)Cutmask)>>11)==0&&BetaRICH>BetabinsAgl[m]&&BetaRICH<=BetabinsAgl[m+1]) DTemplatesAgl2->Fill(mm*((R/BetaRICH)*pow((1-pow(BetaRICH,2)),0.5)),m);}
		if(Massa_gen<4&&Massa_gen>2.5)
			for(int m=0;m<18;m++) 
			{if((((int)Cutmask)>>11)==0&&Var2>BetaAglD[m]&&Var2<=BetaAglD[m+1]) HeTemplatesAgl->Fill(mm*((R/BetaRICH)*pow((1-pow(BetaRICH,2)),0.5)),m);
				if((((int)Cutmask)>>11)==0&&BetaRICH>BetabinsAgl[m]&&BetaRICH<=BetabinsAgl[m+1]) HeTemplatesAgl2->Fill(mm*((R/BetaRICH)*pow((1-pow(BetaRICH,2)),0.5)),m);}
	}

	
        return;
}

void DeutonsDATA_Fill(TNtuple *ntupla, int l,int zona){
        int k = ntupla->GetEvent(l);
	if(Beta<=0||R<=0) return;
	float mm=0;
	if(Likcut&&R<4/*&&Distcut*/)
		 for(int m=0;m<18;m++)  {if(Var>BetaD[m]&&Var<=BetaD[m+1]) {
						DhistosgeoTOF_Dist->Fill(zona,(Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),m);
						if(R>1.2*Rcutoff) DhistosTOF_Dist->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),m);
						}
					 if(Beta>Betabins[m]&&Beta<=Betabins[m+1]) if(R>1.2*Rcutoff) DhistosTOF_Dist2->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),m);	 
					}
	if(Likcut&&R<20&&Distcut){
		mm=1;
                //if(Dist5D>distcut) mm=2; 
		for(int m=0;m<18;m++) 
				if((((int)Cutmask)>>11)==512)
					{if(Var2>BetaNaFD[m]&&Var2<=BetaNaFD[m+1]) {
						DhistosgeoNaF->Fill(zona,mm*((R/BetaRICH)*pow((1-pow(BetaRICH,2)),0.5)),m);
						if(R>1.2*Rcutoff) DhistosNaF->Fill(mm*((R/BetaRICH)*pow((1-pow(BetaRICH,2)),0.5)),m);
						}
					if(BetaRICH>BetabinsNaF[m]&&BetaRICH<=BetabinsNaF[m+1]) if(R>1.2*Rcutoff) DhistosNaF2->Fill(mm*((R/BetaRICH)*pow((1-pow(BetaRICH,2)),0.5)),m);
					}
	}
	if(Likcut&&R<40&&Distcut){
                mm=1;
                //if(Dist5D>distcut) mm=2; 
		for(int m=0;m<18;m++) 
				if((((int)Cutmask)>>11)==0)	
					{if(Var2>BetaAglD[m]&&Var2<=BetaAglD[m+1]) {
						DhistosgeoAgl->Fill(zona,mm*((R/BetaRICH)*pow((1-pow(BetaRICH,2)),0.5)),m);
						if(R>1.2*Rcutoff) DhistosAgl->Fill(mm*((R/BetaRICH)*pow((1-pow(BetaRICH,2)),0.5)),m);
						}	
					 if(BetaRICH>BetabinsAgl[m]&&BetaRICH<=BetabinsAgl[m+1]) if(R>1.2*Rcutoff) DhistosAgl2->Fill(mm*((R/BetaRICH)*pow((1-pow(BetaRICH,2)),0.5)),m);
					}
	}
}

void DeutonsMC_Copy(TFile * file){
	DTemplatesTOF_Dist=(TH2F*)file->Get("DTemplatesTOF_Dist");
	PTemplatesTOF_Dist=(TH2F*)file->Get("PTemplatesTOF_Dist");
	HeTemplatesTOF_Dist=(TH2F*)file->Get("HeTemplatesTOF_Dist");
	DTemplatesNaF=(TH2F*)file->Get("DTemplatesNaF");
        PTemplatesNaF=(TH2F*)file->Get("PTemplatesNaF");
        HeTemplatesNaF=(TH2F*)file->Get("HeTemplatesNaF");
	DTemplatesAgl=(TH2F*)file->Get("DTemplatesAgl");
        PTemplatesAgl=(TH2F*)file->Get("PTemplatesAgl");
        HeTemplatesAgl=(TH2F*)file->Get("HeTemplatesAgl");
	DTemplatesTOF_Dist2=(TH2F*)file->Get("DTemplatesTOF_Dist2");
        PTemplatesTOF_Dist2=(TH2F*)file->Get("PTemplatesTOF_Dist2");
        HeTemplatesTOF_Dist2=(TH2F*)file->Get("HeTemplatesTOF_Dist2");
        DTemplatesNaF2=(TH2F*)file->Get("DTemplatesNaF2");
        PTemplatesNaF2=(TH2F*)file->Get("PTemplatesNaF2");
        HeTemplatesNaF2=(TH2F*)file->Get("HeTemplatesNaF2");
        DTemplatesAgl2=(TH2F*)file->Get("DTemplatesAgl2");
        PTemplatesAgl2=(TH2F*)file->Get("PTemplatesAgl2");
        HeTemplatesAgl2=(TH2F*)file->Get("HeTemplatesAgl2");
	DhistosgeoTOF_Dist=(TH3F*)file->Get("DhistosgeoTOF_Dist");
	DhistosgeoNaF=(TH3F*)file->Get("DhistosgeoNaF");
	DhistosgeoAgl=(TH3F*)file->Get("DhistosgeoAgl");
	DhistosTOF_Dist=(TH2F*)file->Get("DhistosTOF_Dist");
        DhistosNaF=(TH2F*)file->Get("DhistosNaF");
        DhistosAgl=(TH2F*)file->Get("DhistosAgl");
	DhistosTOF_Dist2=(TH2F*)file->Get("DhistosTOF_Dist2");
        DhistosNaF2=(TH2F*)file->Get("DhistosNaF2");
        DhistosAgl2=(TH2F*)file->Get("DhistosAgl2");	
}

void DeutonsMC_Write(){
        DTemplatesTOF_Dist->Write();
        PTemplatesTOF_Dist->Write();
        HeTemplatesTOF_Dist->Write();
        DTemplatesNaF->Write();
        PTemplatesNaF->Write();
        HeTemplatesNaF->Write();
        DTemplatesAgl->Write();
        PTemplatesAgl->Write();
        HeTemplatesAgl->Write();
        DTemplatesTOF_Dist2->Write();
        PTemplatesTOF_Dist2->Write();
        HeTemplatesTOF_Dist2->Write();
        DTemplatesNaF2->Write();
        PTemplatesNaF2->Write();
        HeTemplatesNaF2->Write();
        DTemplatesAgl2->Write();
        PTemplatesAgl2->Write();
        HeTemplatesAgl2->Write();
        DhistosgeoTOF_Dist->Write();
        DhistosgeoNaF->Write();
        DhistosgeoAgl->Write();
        DhistosTOF_Dist->Write();
        DhistosNaF->Write();
        DhistosAgl->Write();
        DhistosTOF_Dist2->Write();
        DhistosNaF2->Write();
        DhistosAgl2->Write();
}

void DeutonsTemplFits(TFile * file1){
	
	cout<<"******************** DISTANCE TEMPLATE FITS *******************"<<endl;
	TH1F *DataDistTOF[12][18];
	TH1F *DTemplTOF[18];
	TH1F *PTemplTOF[18];
	TH1F *HeTemplTOF[18];
	TH1F *DTemplTOF2[18];
        TH1F *PTemplTOF2[18];
        TH1F *HeTemplTOF2[18];	
	
	string nome;
	for(int m=0;m<18;m++){
		DTemplTOF[m]=new TH1F("","",50,-1,1);
		PTemplTOF[m]=new TH1F("","",50,-1,1);
		HeTemplTOF[m]=new TH1F("","",50,-1,1);
		DTemplTOF2[m]=new TH1F("","",50,-1,1);
                PTemplTOF2[m]=new TH1F("","",50,-1,1);
                HeTemplTOF2[m]=new TH1F("","",50,-1,1);
		for(int l=0;l<12;l++){
			DataDistTOF[l][m]=new TH1F("","",50,-1,1);	
		}	
	}
	for(int m=0;m<18;m++){
		for(int i=0;i<50;i++){
		DTemplTOF[m]->SetBinContent(i+1,DTemplatesTOF_Dist->GetBinContent(i+1,m+1));
		PTemplTOF[m]->SetBinContent(i+1,PTemplatesTOF_Dist->GetBinContent(i+1,m+1));
		HeTemplTOF[m]->SetBinContent(i+1,HeTemplatesTOF_Dist->GetBinContent(i+1,m+1));
		DTemplTOF2[m]->SetBinContent(i+1,DTemplatesTOF_Dist2->GetBinContent(i+1,m+1));
                PTemplTOF2[m]->SetBinContent(i+1,PTemplatesTOF_Dist2->GetBinContent(i+1,m+1));
                HeTemplTOF2[m]->SetBinContent(i+1,HeTemplatesTOF_Dist2->GetBinContent(i+1,m+1));
		for(int l=1;l<11;l++)  DataDistTOF[l][m]->SetBinContent(i+1,DhistosgeoTOF_Dist->GetBinContent(l+1,i+1,m+1));
		DataDistTOF[0][m]->SetBinContent(i+1,DhistosTOF_Dist2->GetBinContent(i+1,m+1)); 
		DataDistTOF[11][m]->SetBinContent(i+1,DhistosTOF_Dist->GetBinContent(i+1,m+1));
		}	
	}
	
	/////////////////// TEMP
	TH1F *DataDistTOF_cut[12][18];
        TH1F *DTemplTOF_cut[18];
        TH1F *PTemplTOF_cut[18];
        TH1F *HeTemplTOF_cut[18];
	TH1F *DTemplTOF2_cut[18];
        TH1F *PTemplTOF2_cut[18];
        TH1F *HeTemplTOF2_cut[18];

		
	for(int m=0;m<18;m++){
                DTemplTOF_cut[m]=new TH1F("","",50,-1,1);
                PTemplTOF_cut[m]=new TH1F("","",50,-1,1);
                HeTemplTOF_cut[m]=new TH1F("","",50,-1,1);
                DTemplTOF2_cut[m]=new TH1F("","",50,-1,1);
                PTemplTOF2_cut[m]=new TH1F("","",50,-1,1);
                HeTemplTOF2_cut[m]=new TH1F("","",50,-1,1);
		for(int l=0;l<12;l++){
                        DataDistTOF_cut[l][m]=new TH1F("","",50,-1,1);
                }
        }

	for(int m=0;m<18;m++){
                for(int i=0;i<50;i++){
				if(i>29){	
				DTemplTOF_cut[m]->SetBinContent(i+1,DTemplatesTOF_Dist->GetBinContent(i+1,m+1));
				PTemplTOF_cut[m]->SetBinContent(i+1,PTemplatesTOF_Dist->GetBinContent(i+1,m+1));
				HeTemplTOF_cut[m]->SetBinContent(i+1,HeTemplatesTOF_Dist->GetBinContent(i+1,m+1));
				for(int l=1;l<11;l++)  DataDistTOF_cut[l][m]->SetBinContent(i+1,DhistosgeoTOF_Dist->GetBinContent(l+1,i+1,m+1));
				DataDistTOF_cut[11][m]->SetBinContent(i+1,DhistosTOF_Dist->GetBinContent(i+1,m+1));
			}
                }
        }
	for(int m=0;m<18;m++){
                for(int i=0;i<50;i++){
                                if(i<16||i>35){
					DTemplTOF2_cut[m]->SetBinContent(i+1,DTemplatesTOF_Dist2->GetBinContent(i+1,m+1));
                                        PTemplTOF2_cut[m]->SetBinContent(i+1,PTemplatesTOF_Dist2->GetBinContent(i+1,m+1));
                                        HeTemplTOF2_cut[m]->SetBinContent(i+1,HeTemplatesTOF_Dist2->GetBinContent(i+1,m+1));
					DataDistTOF_cut[0][m]->SetBinContent(i+1,DhistosTOF_Dist2->GetBinContent(i+1,m+1));
			}
		}
	}	
	///////////////////////////
	
	TFractionFitter * fitT[18][12]={{NULL}};
	TObjArray *Tpl[18][12]={{NULL}};
	int cut=25*(1+ddiscrcut);
	int s1[18][12]={{0}};
	float Err[18][12]={{0}};
	TH1F *PTemplTOFW[18][12];
	TH1F *DTemplTOFW[18][12];
	TH1F *HeTemplTOFW[18][12];
	for(int i=0;i<12;i++) for(int m=0;m<18;m++)c30[i][m]=new TCanvas(); 	
		bool He=false;
	for(int i=0;i<12;i++) for(int m=0;m<18;m++) {
		c30[i][m]->cd();
		PTemplTOFW[m][i]=new TH1F("","",50,-1,1);
		DTemplTOFW[m][i]=new TH1F("","",50,-1,1);
		HeTemplTOFW[m][i]=new TH1F("","",50,-1,1);
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
			cut=25*(1+ddiscrcut);;
			Tpl[m][i]->Add(PTemplTOF_cut[m]);
			Tpl[m][i]->Add(DTemplTOF_cut[m]);
			if(He) Tpl[m][i]->Add(HeTemplTOF_cut[m]);
		}
		else{
			cut=0;
			Tpl[m][i]->Add(PTemplTOF2_cut[m]);
                        Tpl[m][i]->Add(DTemplTOF2_cut[m]);
                        if(He) Tpl[m][i]->Add(HeTemplTOF2_cut[m]);
		}
		fitT[m][i] = new TFractionFitter( DataDistTOF_cut[i][m], Tpl[m][i],"q");
		if(He){
			fitT[m][i]->Constrain(0,0.0,1.0);
			fitT[m][i]->Constrain(1,0.0,1.0);
			fitT[m][i]->Constrain(2,0.0,1.0);
			fitT[m][i]->SetRangeX(cut,50);
		}
		else {
			fitT[m][i]->Constrain(0,0,1);
			fitT[m][i]->Constrain(1,0,1);
			fitT[m][i]->SetRangeX(cut,50);
		}	
		if(DataDistTOF_cut[i][m]->Integral(0,50)>50) s1[m][i]=fitT[m][i]->Fit();
		else s1[m][i]=1;
		if(i!=0){
			cut=25*(1+ddiscrcut);
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
				float i1 = PTemplTOF_cut[m]->Integral(cut,50);
				float i2 = DTemplTOF_cut[m]->Integral(cut,50);
				float i3=0;
				if(He) i3 = HeTemplTOF_cut[m]->Integral();
				if(i1>0) for(int j=0; j<PTemplTOF[m]->GetNbinsX();j++) PTemplTOFW[m][i]->SetBinContent(j,w1*PTemplTOF[m]->GetBinContent(j)/i1*itot);
				if(i2>0) for(int j=0; j<DTemplTOF[m]->GetNbinsX();j++) DTemplTOFW[m][i]->SetBinContent(j,w2*DTemplTOF[m]->GetBinContent(j)/i2*itot);
				if(i3>0) for(int j=0; j<HeTemplTOF[m]->GetNbinsX();j++) HeTemplTOFW[m][i]->SetBinContent(j,w3*HeTemplTOF[m]->GetBinContent(j)/i3*itot);
			}
			else{
				float itot= Result->Integral();
                                float i1 = PTemplTOF2_cut[m]->Integral(cut,50);
                                float i2 = DTemplTOF2_cut[m]->Integral(cut,50);
                                float i3=0;
                                if(He) i3 = HeTemplTOF2_cut[m]->Integral();
                                if(i1>0) for(int j=0; j<PTemplTOF2[m]->GetNbinsX();j++) PTemplTOFW[m][i]->SetBinContent(j,w1*PTemplTOF2[m]->GetBinContent(j)/i1*itot);
                                if(i2>0) for(int j=0; j<DTemplTOF2[m]->GetNbinsX();j++) DTemplTOFW[m][i]->SetBinContent(j,w2*DTemplTOF2[m]->GetBinContent(j)/i2*itot);
                                if(i3>0) for(int j=0; j<HeTemplTOF2[m]->GetNbinsX();j++) HeTemplTOFW[m][i]->SetBinContent(j,w3*HeTemplTOF2[m]->GetBinContent(j)/i3*itot);	
			}
			Stack->Add(DTemplTOFW[m][i]);
			Stack->Add(PTemplTOFW[m][i]);
			if(He) Stack->Add(HeTemplTOFW[m][i]);
			DataDistTOF[i][m]->SetMarkerStyle(8);
			Stack->Draw();
			DataDistTOF[i][m]->Draw("epsame");
			Result->SetLineColor(5);
			Result->SetLineWidth(2);
			Result->Draw("SAME");
			float Cov01=fitT[m][i]->GetFitter()->GetCovarianceMatrixElement(0,1);
			float Cov02=fitT[m][i]->GetFitter()->GetCovarianceMatrixElement(0,2);
			float Cov12=fitT[m][i]->GetFitter()->GetCovarianceMatrixElement(1,2);
			float Sigma=pow((pow(e2/w2,2)+pow(e1/w1,2))/2,0.5);
			Err[m][i]= Sigma*DTemplTOFW[m][i]->Integral();
			DCountsgeoTOF_Dist->SetBinContent(m+1,i+1,0,DTemplTOFW[m][i]->Integral(cut,50));
			DCountsgeoTOF_Dist->SetBinContent(m+1,i+1,1,Err[m][i]);	
			if(i==0){
				PCountsTOF_Dist->SetBinContent(m+1,0,PTemplTOFW[m][0]->Integral(0,(int)(26*(1-ddiscrcut))));
                        	PCountsTOF_Dist->SetBinContent(m+1,1,Err[m][i]);
			}
		}
		else{
			DataDistTOF[i][m]->SetMarkerStyle(8);
			DataDistTOF[i][m]->Draw("ep");
			 if(i!=0){PTemplTOF[m]->Draw("same");
                                DTemplTOF[m]->Draw("same");
                                if(He) HeTemplTOF[m]->Draw("same");
                                }
                        if(i==0) {PTemplTOF2[m]->Draw("same");
                                 DTemplTOF2[m]->Draw("same");
                                if(He) HeTemplTOF2[m]->Draw("same");
                                }

		}

	}

	cout<<"******************** MASS NaF TEMPLATE FITS *******************"<<endl;
        TH1F *DataMassNaF[12][18];
        TH1F *DTemplNaF[18];
        TH1F *PTemplNaF[18];
        TH1F *HeTemplNaF[18];
	TH1F *DTemplNaF2[18];
        TH1F *PTemplNaF2[18];
        TH1F *HeTemplNaF2[18];
        for(int m=0;m<18;m++){
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
        for(int m=0;m<18;m++){
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
	TFractionFitter * fitTNaF[18][12]={{NULL}};
	TObjArray *TplNaF[18][12]={{NULL}};
	int cutNaF=17;
	int s1NaF[18][12]={{0}};
	float ErrNaF[18][12]={{0}};
	TH1F *PTemplNaFW[18][12];
	TH1F *DTemplNaFW[18][12];
	TH1F *HeTemplNaFW[18][12];
	for(int i=0;i<12;i++) for(int m=0;m<18;m++) c30_bis[i][m]=new TCanvas(); 	
		He=false;
	for(int i=0;i<12;i++) for(int m=0;m<18;m++) {
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
			TplNaF[m][i]->Add(PTemplNaF[m]);
			TplNaF[m][i]->Add(DTemplNaF[m]);
			if(He) TplNaF[m][i]->Add(HeTemplNaF[m]);
		}
		else{
			TplNaF[m][i]->Add(PTemplNaF2[m]);
               		TplNaF[m][i]->Add(DTemplNaF2[m]);
                	if(He) TplNaF[m][i]->Add(HeTemplNaF2[m]);
		}
		fitTNaF[m][i] = new TFractionFitter( DataMassNaF[i][m], TplNaF[m][i],"q");
		if(He){
			fitTNaF[m][i]->Constrain(0,0.5,1.0);
			fitTNaF[m][i]->Constrain(1,0.0,1.0);
			fitTNaF[m][i]->Constrain(2,0.0,1.0);
		}
		else {
			fitTNaF[m][i]->Constrain(0,0,1);
			fitTNaF[m][i]->Constrain(1,0,1);
			fitTNaF[m][i]->SetRangeX(cutNaF,40);
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
                                PCountsNaF->SetBinContent(m+1,0,/*PTemplNaFW[m][0]*/DataMassNaF[i][m]->Integral(5,25));
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
        TH1F *DataMassAgl[12][18];
        TH1F *DTemplAgl[18];
        TH1F *PTemplAgl[18];
        TH1F *HeTemplAgl[18];
	TH1F *DTemplAgl2[18];
        TH1F *PTemplAgl2[18];
        TH1F *HeTemplAgl2[18];
        for(int m=0;m<18;m++){
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
        for(int m=0;m<18;m++){
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

        TFractionFitter * fitTAgl[18][12]={{NULL}};
        TObjArray *TplAgl[18][12]={{NULL}};
        int cutAgl=16;
        int s1Agl[18][12]={{0}};
        float ErrAgl[18][12]={{0}};
        TH1F *PTemplAglW[18][12];
        TH1F *DTemplAglW[18][12];
        TH1F *HeTemplAglW[18][12];
        for(int i=0;i<12;i++) for(int m=0;m<18;m++) c30_tris[i][m]=new TCanvas();
                He=false;
        for(int i=0;i<12;i++) for(int m=0;m<18;m++) {
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
			TplAgl[m][i]->Add(PTemplAgl[m]);
                	TplAgl[m][i]->Add(DTemplAgl[m]);
                	if(He) TplAgl[m][i]->Add(HeTemplAgl[m]);
                }
		else{
			TplAgl[m][i]->Add(PTemplAgl2[m]);
                        TplAgl[m][i]->Add(DTemplAgl2[m]);
                        if(He) TplAgl[m][i]->Add(HeTemplAgl2[m]);
		}
		fitTAgl[m][i] = new TFractionFitter( DataMassAgl[i][m], TplAgl[m][i],"q");
                if(He){
                        fitTAgl[m][i]->Constrain(0,0.5,1.0);
                        fitTAgl[m][i]->Constrain(1,0.0,1.0);
                        fitTAgl[m][i]->Constrain(2,0.0,1.0);
                }
                else {
                        fitTAgl[m][i]->Constrain(0,0,1);
                        fitTAgl[m][i]->Constrain(1,0,1);
                        fitTAgl[m][i]->SetRangeX(cutAgl,40);
                }
                if((DataMassAgl[i][m]->Integral(25,50)>50&&i!=0)||i==0) s1Agl[m][i]=fitTAgl[m][i]->Fit();
                else s1Agl[m][i]=1;
		for(int z=0;z<5;z++){
			if(s1Agl[m][i]==-1){
				cutAgl+=1;
				fitTAgl[m][i]->SetRangeX(cutAgl,40);
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
                        DCountsgeoAgl->SetBinContent(m+1,i+1,0,DTemplAglW[m][i]->Integral(5,50));
			DCountsgeoAgl->SetBinContent(m+1,i+1,1,ErrAgl[m][i]);
                	if(i==0){
                                PCountsAgl->SetBinContent(m+1,0,/*PTemplAglW[m][0]*/DataMassAgl[i][m]->Integral(5,25));
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
		for(int m=0;m<18;m++) cout<<s1[m][i]<<" ";
		cout<<endl;
	}
	cout<<"**** NaF ********"<<endl;
        for(int i=0;i<12;i++){
                for(int m=0;m<18;m++) cout<<s1NaF[m][i]<<" ";
                cout<<endl;
        }
	cout<<"**** Agl ********"<<endl;
        for(int i=0;i<12;i++){
                for(int m=0;m<18;m++) cout<<s1Agl[m][i]<<" ";
                cout<<endl;
        }

}
