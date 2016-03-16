
TH2F *DTemplatesTOF_Dist=new TH2F("DTemplatesTOF_Dist","DTemplatesTOF_Dist",50,-1,1,18,0,18);
TH2F *PTemplatesTOF_Dist=new TH2F("PTemplatesTOF_Dist","PTemplatesTOF_Dist",50,-1,1,18,0,18);
TH2F *HeTemplatesTOF_Dist=new TH2F("HeTemplatesTOF_Dist","HeTemplatesTOF_Dist",50,-1,1,18,0,18);

TH2F *DTemplatesNaF_Dist=new TH2F("DTemplatesNaF_Dist","DTemplatesNaF_Dist",50,-1,1,18,0,18);
TH2F *PTemplatesNaF_Dist=new TH2F("PTemplatesNaF_Dist","PTemplatesNaF_Dist",50,-1,1,18,0,18);
TH2F *HeTemplatesNaF_Dist=new TH2F("HeTemplatesNaF_Dist","HeTemplatesNaF_Dist",50,-1,1,18,0,18);

TH2F *DTemplatesAgl_Dist=new TH2F("DTemplatesAgl_Dist","DTemplatesAgl_Dist",50,-1,1,18,0,18);
TH2F *PTemplatesAgl_Dist=new TH2F("PTemplatesAgl_Dist","PTemplatesAgl_Dist",50,-1,1,18,0,18);
TH2F *HeTemplatesAgl_Dist=new TH2F("HeTemplatesAgl_Dist","HeTemplatesAgl_Dist",50,-1,1,18,0,18);

TH2F *DTemplatesTOF_Dist2=new TH2F("DTemplatesTOF_Dist2","DTemplatesTOF_Dist2",50,-1,1,18,0,18);
TH2F *PTemplatesTOF_Dist2=new TH2F("PTemplatesTOF_Dist2","PTemplatesTOF_Dist2",50,-1,1,18,0,18);
TH2F *HeTemplatesTOF_Dist2=new TH2F("HeTemplatesTOF_Dist2","HeTemplatesTOF_Dist2",50,-1,1,18,0,18);

TH2F *DTemplatesNaF_Dist2=new TH2F("DTemplatesNaF_Dist2","DTemplatesNaF_Dist2",50,-1,1,18,0,18);
TH2F *PTemplatesNaF_Dist2=new TH2F("PTemplatesNaF_Dist2","PTemplatesNaF_Dist2",50,-1,1,18,0,18);
TH2F *HeTemplatesNaF_Dist2=new TH2F("HeTemplatesNaF_Dist2","HeTemplatesNaF_Dist2",50,-1,1,18,0,18);

TH2F *DTemplatesAgl_Dist2=new TH2F("DTemplatesAgl_Dist2","DTemplatesAgl_Dist2",50,-1,1,18,0,18);
TH2F *PTemplatesAgl_Dist2=new TH2F("PTemplatesAgl_Dist2","PTemplatesAgl_Dist2",50,-1,1,18,0,18);
TH2F *HeTemplatesAgl_Dist2=new TH2F("HeTemplatesAgl_Dist2","HeTemplatesAgl_Dist2",50,-1,1,18,0,18);

TH3F *DhistosgeoTOF_Dist=new TH3F("DhistosgeoTOF_Dist","DhistosgeoTOF_Dist",11,0,11,50,-1,1,18,0,18);
TH3F *DhistosgeoNaF_Dist=new TH3F("DhistosgeoNaF_Dist","DhistosgeoNaF_Dist",11,0,11,50,-1,1,18,0,18);
TH3F *DhistosgeoAgl_Dist=new TH3F("DhistosgeoAgl_Dist","DhistosgeoAgl_Dist",11,0,11,50,-1,1,18,0,18);

TH2F *DhistosTOF_Dist=new TH2F("DhistosTOF_Dist","DhistosTOF_Dist",50,-1,1,18,0,18);
TH2F *DhistosNaF_Dist=new TH2F("DhistosNaF_Dist","DhistosNaF_Dist",50,-1,1,18,0,18);
TH2F *DhistosAgl_Dist=new TH2F("DhistosAgl_Dist","DhistosAgl_Dist",50,-1,1,18,0,18);

TH2F *DhistosTOF_Dist2=new TH2F("DhistosTOF_Dist2","DhistosTOF_Dist2",50,-1,1,18,0,18);
TH2F *DhistosNaF_Dist2=new TH2F("DhistosNaF_Dist2","DhistosNaF_Dist2",50,-1,1,18,0,18);
TH2F *DhistosAgl_Dist2=new TH2F("DhistosAgl_Dist2","DhistosAgl_Dist2",50,-1,1,18,0,18);


void DeutonsMC_Dist_Fill(TNtuple *ntupla, int l){
        int k = ntupla->GetEvent(l);
	if(Beta<=0||R<=0) return;
        float mm=0;
	if(Likcut&&R<4&&Distcut){
		if(Massa_gen<1&&Massa_gen>0.5) 
			for(int m=0;m<18;m++)  {if(Var>BetaD[m]&&Var<=BetaD[m+1]) PTemplatesTOF_Dist->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),m+0.1); ;
						if(Var>BetaP[m]&&Var<=BetaP[m+1]) PTemplatesTOF_Dist2->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),m+0.1);} 
		if(Massa_gen<2&&Massa_gen>1.5)
                        for(int m=0;m<18;m++)  {if(Var>BetaD[m]&&Var<=BetaD[m+1]) DTemplatesTOF_Dist->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),m+0.1);	
						if(Var>BetaP[m]&&Var<=BetaP[m+1]) DTemplatesTOF_Dist2->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),m+0.1);}
		 if(Massa_gen<4&&Massa_gen>2.5)
                        for(int m=0;m<18;m++)  {if(Var>BetaD[m]&&Var<=BetaD[m+1]) HeTemplatesTOF_Dist->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),m+0.1);
						if(Var>BetaP[m]&&Var<=BetaP[m+1]) HeTemplatesTOF_Dist2->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),m+0.1);}
	}
	if(Likcut&&R<20&&Distcut){
		mm=1;
		//if(Dist5D>distcut) mm=2;
		if(Massa_gen<1&&Massa_gen>0.5)
			for(int m=0;m<18;m++) 
			{if((((int)Cutmask)>>11)==512&&Var2>BetaNaFD[m]&&Var2<=BetaNaFD[m+1]) PTemplatesNaF_Dist->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),m+0.1);
			 if((((int)Cutmask)>>11)==512&&Var2>BetaNaFP[m]&&Var2<=BetaNaFP[m+1])	 PTemplatesNaF_Dist2->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),m+0.1);}	
		if(Massa_gen<2&&Massa_gen>1.5)
			for(int m=0;m<18;m++) 
			{if((((int)Cutmask)>>11)==512&&Var2>BetaNaFD[m]&&Var2<=BetaNaFD[m+1]) DTemplatesNaF_Dist->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),m+0.1);
			 if((((int)Cutmask)>>11)==512&&Var2>BetaNaFP[m]&&Var2<=BetaNaFP[m+1]) DTemplatesNaF_Dist2->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),m+0.1);}	
		if(Massa_gen<4&&Massa_gen>2.5)
			for(int m=0;m<18;m++) 
			{if((((int)Cutmask)>>11)==512&&Var2>BetaNaFD[m]&&Var2<=BetaNaFD[m+1]) HeTemplatesNaF_Dist->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),m+0.1);
			if((((int)Cutmask)>>11)==512&&Var2>BetaNaFP[m]&&Var2<=BetaNaFP[m+1]) HeTemplatesNaF_Dist2->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),m+0.1);}
	}
	if(Likcut&&R<40&&Distcut){
		mm=1;
                //if(Dist5D>distcut) mm=2;
		if(Massa_gen<1&&Massa_gen>0.5)
			for(int m=0;m<18;m++) 
			{if((((int)Cutmask)>>11)==0&&Var2>BetaAglD[m]&&Var2<=BetaAglD[m+1]) {PTemplatesAgl_Dist->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),m+0.1); }
			if((((int)Cutmask)>>11)==0&&Var2>BetaAglP[m]&&Var2<=BetaAglP[m+1]) PTemplatesAgl_Dist2->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),m+0.1);}
		if(Massa_gen<2&&Massa_gen>1.5)
			for(int m=0;m<18;m++) 
			{if((((int)Cutmask)>>11)==0&&Var2>BetaAglD[m]&&Var2<=BetaAglD[m+1]) DTemplatesAgl_Dist->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),m+0.1);
			if((((int)Cutmask)>>11)==0&&Var2>BetaAglP[m]&&Var2<=BetaAglP[m+1]) DTemplatesAgl_Dist2->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),m+0.1);}
		if(Massa_gen<4&&Massa_gen>2.5)
			for(int m=0;m<18;m++) 
			{if((((int)Cutmask)>>11)==0&&Var2>BetaAglD[m]&&Var2<=BetaAglD[m+1]) HeTemplatesAgl_Dist->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),m+0.1);
			if((((int)Cutmask)>>11)==0&&Var2>BetaAglP[m]&&Var2<=BetaAglP[m+1]) HeTemplatesAgl_Dist2->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),m+0.1);}
	}

	
        return;
}

void DeutonsDATA_Dist_Fill(TNtuple *ntupla, int l,int zona){
        int k = ntupla->GetEvent(l);
	if(Beta<=0||R<=0) return;
	float mm=0;
	if(Likcut&&R<4&&Distcut)
		 for(int m=0;m<18;m++)  {if(Var>BetaD[m]&&Var<=BetaD[m+1]) {
						DhistosgeoTOF_Dist->Fill(zona,(Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),m+0.1);
						if(R>1.2*Rcutoff) DhistosTOF_Dist->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),m+0.1);
						}
					 if(Var>BetaP[m]&&Var<=BetaP[m+1]) if(R>1.2*Rcutoff) DhistosTOF_Dist2->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),m+0.1);	 
					}
	if(Likcut&&R<20&&Distcut){
		mm=1;
                //if(Dist5D>distcut) mm=2; 
		for(int m=0;m<18;m++) 
				if((((int)Cutmask)>>11)==512)
					{if(Var2>BetaNaFD[m]&&Var2<=BetaNaFD[m+1]) {
						DhistosgeoNaF_Dist->Fill(zona,(Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),m+0.1);
						if(R>1.2*Rcutoff) DhistosNaF_Dist->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),m+0.1);
						}
					if(Var2>BetaNaFP[m]&&Var2<=BetaNaFP[m+1]) if(R>1.2*Rcutoff) DhistosNaF_Dist2->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),m+0.1);
					}
	}
	if((((int)Cutmask)>>11)==0&&Dist5D+Dist5D_P>20) return;
	if(Likcut&&R<40&&Distcut){
                mm=1;
		for(int m=0;m<18;m++) 
				if((((int)Cutmask)>>11)==0)	
					{
					if(Var2>BetaAglD[m]&&Var2<=BetaAglD[m+1]) {
						DhistosgeoAgl_Dist->Fill(zona,(Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),m+0.1);
						if(R>1.2*Rcutoff) DhistosAgl_Dist->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),m+0.1);
						}	
					 if(Var2>BetaAglP[m]&&Var2<=BetaAglP[m+1]) if(R>1.2*Rcutoff) DhistosAgl_Dist2->Fill((Dist5D_P-Dist5D)/(Dist5D_P+Dist5D),m+0.1);
					}
	}
}

void DeutonsMC_Dist_Copy(TFile * file){
	DTemplatesTOF_Dist=(TH2F*)file->Get("DTemplatesTOF_Dist");
	PTemplatesTOF_Dist=(TH2F*)file->Get("PTemplatesTOF_Dist");
	HeTemplatesTOF_Dist=(TH2F*)file->Get("HeTemplatesTOF_Dist");
	DTemplatesNaF_Dist=(TH2F*)file->Get("DTemplatesNaF_Dist");
        PTemplatesNaF_Dist=(TH2F*)file->Get("PTemplatesNaF_Dist");
        HeTemplatesNaF_Dist=(TH2F*)file->Get("HeTemplatesNaF_Dist");
	DTemplatesAgl_Dist=(TH2F*)file->Get("DTemplatesAgl_Dist");
        PTemplatesAgl_Dist=(TH2F*)file->Get("PTemplatesAgl_Dist");
        HeTemplatesAgl_Dist=(TH2F*)file->Get("HeTemplatesAgl_Dist");
	DTemplatesTOF_Dist2=(TH2F*)file->Get("DTemplatesTOF_Dist2");
        PTemplatesTOF_Dist2=(TH2F*)file->Get("PTemplatesTOF_Dist2");
        HeTemplatesTOF_Dist2=(TH2F*)file->Get("HeTemplatesTOF_Dist2");
        DTemplatesNaF_Dist2=(TH2F*)file->Get("DTemplatesNaF_Dist2");
        PTemplatesNaF_Dist2=(TH2F*)file->Get("PTemplatesNaF_Dist2");
        HeTemplatesNaF_Dist2=(TH2F*)file->Get("HeTemplatesNaF_Dist2");
        DTemplatesAgl_Dist2=(TH2F*)file->Get("DTemplatesAgl_Dist2");
        PTemplatesAgl_Dist2=(TH2F*)file->Get("PTemplatesAgl_Dist2");
        HeTemplatesAgl_Dist2=(TH2F*)file->Get("HeTemplatesAgl_Dist2");
	DhistosgeoTOF_Dist=(TH3F*)file->Get("DhistosgeoTOF_Dist");
	DhistosgeoNaF_Dist=(TH3F*)file->Get("DhistosgeoNaF_Dist");
	DhistosgeoAgl_Dist=(TH3F*)file->Get("DhistosgeoAgl_Dist");
	DhistosTOF_Dist=(TH2F*)file->Get("DhistosTOF_Dist");
        DhistosNaF_Dist=(TH2F*)file->Get("DhistosNaF_Dist");
        DhistosAgl_Dist=(TH2F*)file->Get("DhistosAgl_Dist");
	DhistosTOF_Dist2=(TH2F*)file->Get("DhistosTOF_Dist2");
        DhistosNaF_Dist2=(TH2F*)file->Get("DhistosNaF_Dist2");
        DhistosAgl_Dist2=(TH2F*)file->Get("DhistosAgl_Dist2");	
}

void DeutonsMC_Dist_Write(){
        DTemplatesTOF_Dist->Write();
        PTemplatesTOF_Dist->Write();
        HeTemplatesTOF_Dist->Write();
        DTemplatesNaF_Dist->Write();
        PTemplatesNaF_Dist->Write();
        HeTemplatesNaF_Dist->Write();
        DTemplatesAgl_Dist->Write();
        PTemplatesAgl_Dist->Write();
        HeTemplatesAgl_Dist->Write();
        DTemplatesTOF_Dist2->Write();
        PTemplatesTOF_Dist2->Write();
        HeTemplatesTOF_Dist2->Write();
        DTemplatesNaF_Dist2->Write();
        PTemplatesNaF_Dist2->Write();
        HeTemplatesNaF_Dist2->Write();
        DTemplatesAgl_Dist2->Write();
        PTemplatesAgl_Dist2->Write();
        HeTemplatesAgl_Dist2->Write();
        DhistosgeoTOF_Dist->Write();
        DhistosgeoNaF_Dist->Write();
        DhistosgeoAgl_Dist->Write();
        DhistosTOF_Dist->Write();
        DhistosNaF_Dist->Write();
        DhistosAgl_Dist->Write();
        DhistosTOF_Dist2->Write();
        DhistosNaF_Dist2->Write();
        DhistosAgl_Dist2->Write();
}



TCanvas *c40[12][18];
TCanvas *c40_bis[12][18];
TCanvas *c40_tris[12][18];

TH3F * DCountsgeoTOF_Dist = new TH3F("DCountsgeoTOF_Dist","DCountsgeoTOF_Dist",18,0,18,12,0,12,2,0,2);
TH3F * DCountsgeoNaF_Dist = new TH3F("DCountsgeoNaF_Dist","DCountsgeoNaF_Dist",18,0,18,12,0,12,2,0,2);
TH3F * DCountsgeoAgl_Dist = new TH3F("DCountsgeoAgl_Dist","DCountsgeoAgl_Dist",18,0,18,12,0,12,2,0,2);
TH2F * PCountsTOF_Dist = new TH2F("PCountsTOF_Dist","PCountsTOF_Dist",18,0,18,2,0,2);
TH2F * PCountsNaF_Dist = new TH2F("PCountsNaF_Dist","PCountsNaF_Dist",18,0,18,2,0,2);
TH2F * PCountsAgl_Dist = new TH2F("PCountsAgl_Dist","PCountsAgl_Dist",18,0,18,2,0,2);


void DeutonsTemplFits_Dist(TFile * file1){
	cout<<"******************** Dist TOF TEMPLATE FITS *******************"<<endl;
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
	
	TFractionFitter * fitT[18][12]={{NULL}};
	TObjArray *Tpl[18][12]={{NULL}};
	int cut1=21;
	int cut2=29;
	int s1[18][12]={{0}};
	float Err[18][12]={{0}};
	TH1F *PTemplTOFW[18][12];
	TH1F *DTemplTOFW[18][12];
	TH1F *HeTemplTOFW[18][12];
	for(int i=0;i<12;i++) for(int m=0;m<18;m++)c40[i][m]=new TCanvas(); 	
		bool He=true;
	for(int i=0;i<12;i++) for(int m=0;m<18;m++) {
		c40[i][m]->cd();
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
			Tpl[m][i]->Add(PTemplTOF[m]);
			Tpl[m][i]->Add(DTemplTOF[m]);
			if(He) Tpl[m][i]->Add(HeTemplTOF[m]);
		}
		else{	
			Tpl[m][i]->Add(PTemplTOF2[m]);
                        Tpl[m][i]->Add(DTemplTOF2[m]);
                        if(He) Tpl[m][i]->Add(HeTemplTOF2[m]);
		}	
		fitT[m][i] = new TFractionFitter( DataDistTOF[i][m], Tpl[m][i],"Q");
		if(He){
			fitT[m][i]->Constrain(0,0.0,1.0);
			fitT[m][i]->Constrain(1,0.0,0.1);
			fitT[m][i]->Constrain(2,0.0,0.01);
			fitT[m][i]->SetRangeX(0,50);
		}
		else {
			fitT[m][i]->Constrain(0,0,1);
			fitT[m][i]->Constrain(1,0,1);
			fitT[m][i]->SetRangeX(0,50);
		}
		cut1=21;
        	cut2=29;
		for(int z=cut1;z<=cut2;z++) fitT[m][i]->ExcludeBin(z);
		if(DataDistTOF[i][m]->Integral(0,50)>50) s1[m][i]=fitT[m][i]->Fit();
		else s1[m][i]=1;
		if(s1[m][i]==-1){
                        for(int z=0;z<5;z++){
                                if(s1[m][i]==-1){
                                        cut2+=1;
					cut1-=1;
                                        fitT[m][i]->ExcludeBin(cut1);
					fitT[m][i]->ExcludeBin(cut2);
					s1[m][i]=fitT[m][i]->Fit();
                                }
                                if(s1[m][i]==0) {s1[m][i]=2;break;}
                        }
                }
		
		double w1,w2,w3=0;
		double e1,e2,e3=0;
		if(s1[m][i]==0||s1[m][i]==2){
			fitT[m][i]->GetResult(0,w1,e1);
			fitT[m][i]->GetResult(1,w2,e2);
			if(He) fitT[m][i]->GetResult(2,w3,e3);
			Result = (TH1F*) fitT[m][i]->GetPlot();
			if(!Result) continue;	
			if(i!=0){
				
				float i1 = PTemplTOF[m]->Integral(0,cut1) + PTemplTOF[m]->Integral(cut2,50);
                                float i2 = DTemplTOF[m]->Integral(0,cut1) + DTemplTOF[m]->Integral(cut2,50);
                                float i3=0;
                                float itot= Result->Integral();
                                if(He) i3 = HeTemplTOF[m]->Integral(0,cut1) + HeTemplTOF[m]->Integral(cut2,50);
				if(i1>0) for(int j=0; j<PTemplTOF[m]->GetNbinsX();j++) PTemplTOFW[m][i]->SetBinContent(j,w1*PTemplTOF[m]->GetBinContent(j)/i1*itot);
				if(i2>0) for(int j=0; j<DTemplTOF[m]->GetNbinsX();j++) DTemplTOFW[m][i]->SetBinContent(j,w2*DTemplTOF[m]->GetBinContent(j)/i2*itot);
				if(i3>0) for(int j=0; j<HeTemplTOF[m]->GetNbinsX();j++) HeTemplTOFW[m][i]->SetBinContent(j,w3*HeTemplTOF[m]->GetBinContent(j)/i3*itot);
			}
			else{
                                float i1 = PTemplTOF2[m]->Integral(0,cut1) + PTemplTOF2[m]->Integral(cut2,50);
                                float i2 = DTemplTOF2[m]->Integral(0,cut1) + DTemplTOF2[m]->Integral(cut2,50);
                                float i3=0;
                                float itot= Result->Integral();
				if(He) i3 = HeTemplTOF2[m]->Integral(0,cut1) + HeTemplTOF2[m]->Integral(cut2,50);
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
			DCountsgeoTOF_Dist->SetBinContent(m+1,i+1,0,DTemplTOFW[m][i]->Integral(0,50));

			DCountsgeoTOF_Dist->SetBinContent(m+1,i+1,1,Err[m][i]);	
			if(i==0){
				PCountsTOF_Dist->SetBinContent(m+1,0,PTemplTOFW[m][0]->Integral(0,50));
                        	PCountsTOF_Dist->SetBinContent(m+1,1,Err[m][i]);
			}
		}
		else{
			DataDistTOF[i][m]->SetMarkerStyle(8);
                        DataDistTOF[i][m]->Draw("ep");

			PTemplTOF[m]->SetFillStyle(3001);
	                DTemplTOF[m]->SetFillStyle(3001);
        	        HeTemplTOF[m]->SetFillStyle(3001);
                	PTemplTOF[m]->SetFillColor(2);
                	DTemplTOF[m]->SetFillColor(4);
                	HeTemplTOF[m]->SetFillColor(3);
 
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

	cout<<"******************** Dist NaF TEMPLATE FITS *******************"<<endl;
        TH1F *DataDistNaF[12][18];
        TH1F *DTemplNaF[18];
        TH1F *PTemplNaF[18];
        TH1F *HeTemplNaF[18];
	TH1F *DTemplNaF2[18];
        TH1F *PTemplNaF2[18];
        TH1F *HeTemplNaF2[18];
        for(int m=0;m<18;m++){
                DTemplNaF[m]=new TH1F("","",50,-1,1);
                PTemplNaF[m]=new TH1F("","",50,-1,1);
                HeTemplNaF[m]=new TH1F("","",50,-1,1);
                DTemplNaF2[m]=new TH1F("","",50,-1,1);
                PTemplNaF2[m]=new TH1F("","",50,-1,1);
                HeTemplNaF2[m]=new TH1F("","",50,-1,1);
		for(int l=0;l<12;l++){
                        DataDistNaF[l][m]=new TH1F("","",50,-1,1);
                }
        }
        for(int m=0;m<18;m++){
                for(int i=0;i<50;i++){
                DTemplNaF[m]->SetBinContent(i+1,DTemplatesNaF_Dist->GetBinContent(i+1,m+1));
                PTemplNaF[m]->SetBinContent(i+1,PTemplatesNaF_Dist->GetBinContent(i+1,m+1));
                HeTemplNaF[m]->SetBinContent(i+1,HeTemplatesNaF_Dist->GetBinContent(i+1,m+1));
		DTemplNaF2[m]->SetBinContent(i+1,DTemplatesNaF_Dist2->GetBinContent(i+1,m+1));
                PTemplNaF2[m]->SetBinContent(i+1,PTemplatesNaF_Dist2->GetBinContent(i+1,m+1));
                HeTemplNaF2[m]->SetBinContent(i+1,HeTemplatesNaF_Dist2->GetBinContent(i+1,m+1));
                for(int l=1;l<11;l++)  DataDistNaF[l][m]->SetBinContent(i+1,DhistosgeoNaF_Dist->GetBinContent(l+1,i+1,m+1));
                DataDistNaF[11][m]->SetBinContent(i+1,DhistosNaF_Dist->GetBinContent(i+1,m+1));
		DataDistNaF[0][m]->SetBinContent(i+1,DhistosNaF_Dist2->GetBinContent(i+1,m+1));
		}
        }
	TFractionFitter * fitTNaF[18][12]={{NULL}};
	TObjArray *TplNaF[18][12]={{NULL}};
	int cutNaF=0;
	int s1NaF[18][12]={{0}};
	float ErrNaF[18][12]={{0}};
	TH1F *PTemplNaFW[18][12];
	TH1F *DTemplNaFW[18][12];
	TH1F *HeTemplNaFW[18][12];
	for(int i=0;i<12;i++) for(int m=0;m<18;m++) c40_bis[i][m]=new TCanvas(); 	
		He=true;
	for(int i=0;i<12;i++) for(int m=0;m<18;m++) {
		c40_bis[i][m]->cd();
		cutNaF=0;
		He=false;
		 if(HeTemplNaF[m]->Integral(0,50)<=50) He=false;
		PTemplNaFW[m][i]=new TH1F("","",50,-1,1);
		DTemplNaFW[m][i]=new TH1F("","",50,-1,1);
		HeTemplNaFW[m][i]=new TH1F("","",50,-1,1);
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
		fitTNaF[m][i] = new TFractionFitter( DataDistNaF[i][m], TplNaF[m][i],"q");
		if(He){
			fitTNaF[m][i]->Constrain(0,0.0,1.0);
			fitTNaF[m][i]->Constrain(1,0.0,0.1);
			fitTNaF[m][i]->Constrain(2,0.0,0.01);
		}
		else {
			fitTNaF[m][i]->Constrain(0,0,1);
			fitTNaF[m][i]->Constrain(1,0,1);
			fitTNaF[m][i]->SetRangeX(cutNaF,50);
		}
		if((DataDistNaF[i][m]->Integral(0,50)>50)) s1NaF[m][i]=fitTNaF[m][i]->Fit();
		else s1NaF[m][i]=1;
		cut1=21;
                cut2=29;
                for(int z=cut1;z<=cut2;z++) fitTNaF[m][i]->ExcludeBin(z);
                if(DataDistNaF[i][m]->Integral(0,50)>50) s1NaF[m][i]=fitTNaF[m][i]->Fit();
                else s1NaF[m][i]=1;
                if(s1NaF[m][i]==-1){
                        for(int z=0;z<5;z++){
                                if(s1NaF[m][i]==-1){
                                        cut2+=1;
                                        cut1-=1;
                                        fitTNaF[m][i]->ExcludeBin(cut1);
                                        fitTNaF[m][i]->ExcludeBin(cut2);
                                        s1NaF[m][i]=fitTNaF[m][i]->Fit();
                                }
                                if(s1NaF[m][i]==0) {s1NaF[m][i]=2;break;}
                        }
                }

		double w1,w2,w3=0;
		double e1,e2,e3=0;
		if(s1NaF[m][i]==0||s1NaF[m][i]==2){
			fitTNaF[m][i]->GetResult(0,w1,e1);
			fitTNaF[m][i]->GetResult(1,w2,e2);
			if(He) fitTNaF[m][i]->GetResult(2,w3,e3);
			Result = (TH1F*) fitTNaF[m][i]->GetPlot();
			if(!Result) continue;
			if(i!=0){
				
				float i1 = PTemplNaF[m]->Integral(0,cut1) + PTemplNaF[m]->Integral(cut2,50);
                                float i2 = DTemplNaF[m]->Integral(0,cut1) + DTemplNaF[m]->Integral(cut2,50);
                                float i3=0;
                                float itot= Result->Integral();
                                if(He) i3 = HeTemplNaF[m]->Integral(0,cut1) + HeTemplNaF[m]->Integral(cut2,50);

				if(i1>0) for(int j=0; j<PTemplNaF[m]->GetNbinsX();j++) PTemplNaFW[m][i]->SetBinContent(j,w1*PTemplNaF[m]->GetBinContent(j)/i1*itot);
				if(i2>0) for(int j=0; j<DTemplNaF[m]->GetNbinsX();j++) DTemplNaFW[m][i]->SetBinContent(j,w2*DTemplNaF[m]->GetBinContent(j)/i2*itot);
				if(i3>0) for(int j=0; j<HeTemplNaF[m]->GetNbinsX();j++) HeTemplNaFW[m][i]->SetBinContent(j,w3*HeTemplNaF[m]->GetBinContent(j)/i3*itot);
			}
			else{
				float i1 = PTemplNaF2[m]->Integral(0,cut1) + PTemplNaF2[m]->Integral(cut2,50);
                                float i2 = DTemplNaF2[m]->Integral(0,cut1) + DTemplNaF2[m]->Integral(cut2,50);
                                float i3=0;
                                float itot= Result->Integral();
                                if(He) i3 = HeTemplNaF[m]->Integral(0,cut1) + HeTemplNaF[m]->Integral(cut2,50);
				if(i1>0) for(int j=0; j<PTemplNaF2[m]->GetNbinsX();j++) PTemplNaFW[m][i]->SetBinContent(j,w1*PTemplNaF2[m]->GetBinContent(j)/i1*itot);
                                if(i2>0) for(int j=0; j<DTemplNaF2[m]->GetNbinsX();j++) DTemplNaFW[m][i]->SetBinContent(j,w2*DTemplNaF2[m]->GetBinContent(j)/i2*itot);
                                if(i3>0) for(int j=0; j<HeTemplNaF2[m]->GetNbinsX();j++) HeTemplNaFW[m][i]->SetBinContent(j,w3*HeTemplNaF2[m]->GetBinContent(j)/i3*itot);
			}
			Stack->Add(PTemplNaFW[m][i]);
			Stack->Add(DTemplNaFW[m][i]);
			if(He) Stack->Add(HeTemplNaFW[m][i]);
			DataDistNaF[i][m]->SetMarkerStyle(8);
			Stack->Draw();
			DataDistNaF[i][m]->Draw("epsame");
			Result->SetLineColor(5);
			Result->SetLineWidth(2);
			Result->Draw("SAME");
			
			float Cov01=fitTNaF[m][i]->GetFitter()->GetCovarianceMatrixElement(0,1);
			float Cov02=fitTNaF[m][i]->GetFitter()->GetCovarianceMatrixElement(0,2);
			float Cov12=fitTNaF[m][i]->GetFitter()->GetCovarianceMatrixElement(1,2);
			float Sigma=pow((pow(e2/w2,2)+pow(e1/w1,2))/2,0.5);
			ErrNaF[m][i]= Sigma*DTemplNaFW[m][i]->Integral();
			DCountsgeoNaF_Dist->SetBinContent(m+1,i+1,0,DTemplNaFW[m][i]->Integral(0,50));
			DCountsgeoNaF_Dist->SetBinContent(m+1,i+1,1,ErrNaF[m][i]);	
			if(i==0){
                                PCountsNaF_Dist->SetBinContent(m+1,0,PTemplNaFW[m][0]->Integral(0,50));
                                PCountsNaF_Dist->SetBinContent(m+1,1,Err[m][i]);
                        }

		}
		else{
			DataDistNaF[i][m]->SetMarkerStyle(8);
			DataDistNaF[i][m]->Draw("ep");
			 if(i!=0){PTemplNaF[m]->Draw("same");
                                DTemplNaF[m]->Draw("same");
                                if(He) HeTemplNaF[m]->Draw("same");
                                }
                        if(i==0) {PTemplNaF2[m]->Draw("same");
                                 DTemplNaF2[m]->Draw("same");
                                if(He) HeTemplNaF2[m]->Draw("same");
                                PCountsNaF_Dist->SetBinContent(m+1,0,PTemplNaFW[m][0]->Integral(0,50));
				}

		}

	}
	cout<<"******************** Dist Agl TEMPLATE FITS *******************"<<endl;
        TH1F *DataDistAgl[12][18];
        TH1F *DTemplAgl[18];
        TH1F *PTemplAgl[18];
        TH1F *HeTemplAgl[18];
	TH1F *DTemplAgl2[18];
        TH1F *PTemplAgl2[18];
        TH1F *HeTemplAgl2[18];
        for(int m=0;m<18;m++){
                DTemplAgl[m]=new TH1F("","",50,-1,1);
                PTemplAgl[m]=new TH1F("","",50,-1,1);
                HeTemplAgl[m]=new TH1F("","",50,-1,1);
		DTemplAgl2[m]=new TH1F("","",50,-1,1);
                PTemplAgl2[m]=new TH1F("","",50,-1,1);
                HeTemplAgl2[m]=new TH1F("","",50,-1,1);
                for(int l=0;l<12;l++){
                        DataDistAgl[l][m]=new TH1F("","",50,-1,1);
                }
        }
        for(int m=0;m<18;m++){
                for(int i=0;i<50;i++){
                DTemplAgl[m]->SetBinContent(i+1,DTemplatesAgl_Dist->GetBinContent(i+1,m+1));
                PTemplAgl[m]->SetBinContent(i+1,PTemplatesAgl_Dist->GetBinContent(i+1,m+1));
                HeTemplAgl[m]->SetBinContent(i+1,HeTemplatesAgl_Dist->GetBinContent(i+1,m+1));
                DTemplAgl2[m]->SetBinContent(i+1,DTemplatesAgl_Dist2->GetBinContent(i+1,m+1));
                PTemplAgl2[m]->SetBinContent(i+1,PTemplatesAgl_Dist2->GetBinContent(i+1,m+1));
                HeTemplAgl2[m]->SetBinContent(i+1,HeTemplatesAgl_Dist2->GetBinContent(i+1,m+1));
		for(int l=1;l<11;l++)  DataDistAgl[l][m]->SetBinContent(i+1,DhistosgeoAgl_Dist->GetBinContent(l+1,i+1,m+1));
                DataDistAgl[11][m]->SetBinContent(i+1,DhistosAgl_Dist->GetBinContent(i+1,m+1));
		DataDistAgl[0][m]->SetBinContent(i+1,DhistosAgl_Dist2->GetBinContent(i+1,m+1));
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
        for(int i=0;i<12;i++) for(int m=0;m<18;m++) c40_tris[i][m]=new TCanvas();
                He=false;
        for(int i=0;i<12;i++) for(int m=0;m<18;m++) {
                c40_tris[i][m]->cd();
                cutAgl=0;
		PTemplAglW[m][i]=new TH1F("","",50,-1,1);
                DTemplAglW[m][i]=new TH1F("","",50,-1,1);
                HeTemplAglW[m][i]=new TH1F("","",50,-1,1);
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
		fitTAgl[m][i] = new TFractionFitter( DataDistAgl[i][m], TplAgl[m][i],"q");
                if(He){
                        fitTAgl[m][i]->Constrain(0,0.0,1.0);
                        fitTAgl[m][i]->Constrain(1,0.0,0.1);
                        fitTAgl[m][i]->Constrain(2,0.0,0.015);
                }
                else {
                        fitTAgl[m][i]->Constrain(0,0,1);
                        fitTAgl[m][i]->Constrain(1,0,1);
                        fitTAgl[m][i]->SetRangeX(cutAgl,40);
                }
		if((DataDistAgl[i][m]->Integral(25,50)>50)) s1Agl[m][i]=fitTAgl[m][i]->Fit();
                else s1Agl[m][i]=1;
		cut1=21;
                cut2=29;
                for(int z=cut1;z<=cut2;z++) fitTAgl[m][i]->ExcludeBin(z);
                if(DataDistAgl[i][m]->Integral(0,50)>50) s1Agl[m][i]=fitTAgl[m][i]->Fit();
                else s1Agl[m][i]=1;
                if(s1Agl[m][i]==-1){
                        for(int z=0;z<5;z++){
                                if(s1Agl[m][i]==-1){
                                        cut2+=1;
                                        cut1-=1;
                                        fitTAgl[m][i]->ExcludeBin(cut1);
                                        fitTAgl[m][i]->ExcludeBin(cut2);
                                        s1Agl[m][i]=fitTAgl[m][i]->Fit();
                                }
                                if(s1Agl[m][i]==0) {s1Agl[m][i]=2;break;}
                        }
                }
	
		double w1,w2,w3=0;
                double e1,e2,e3=0;
                if(s1Agl[m][i]==0||s1Agl[m][i]==2){
                        fitTAgl[m][i]->GetResult(0,w1,e1);
                        fitTAgl[m][i]->GetResult(1,w2,e2);
                        if(He) fitTAgl[m][i]->GetResult(2,w3,e3);
                        Result = (TH1F*) fitTAgl[m][i]->GetPlot();
                        if(!Result) continue;
			if(i!=0){
                        	float i1 = PTemplAgl[m]->Integral(0,cut1) + PTemplAgl[m]->Integral(cut2,50);
                                float i2 = DTemplAgl[m]->Integral(0,cut1) + DTemplAgl[m]->Integral(cut2,50);
                                float i3=0;
                                float itot= Result->Integral();
                                if(He) i3 = HeTemplAgl[m]->Integral(0,cut1) + HeTemplAgl[m]->Integral(cut2,50);
				if(i1>0) for(int j=0; j<PTemplAgl[m]->GetNbinsX();j++) PTemplAglW[m][i]->SetBinContent(j,w1*PTemplAgl[m]->GetBinContent(j)/i1*itot);
                        	if(i2>0) for(int j=0; j<DTemplAgl[m]->GetNbinsX();j++) DTemplAglW[m][i]->SetBinContent(j,w2*DTemplAgl[m]->GetBinContent(j)/i2*itot);
                        	if(i3>0) for(int j=0; j<HeTemplAgl[m]->GetNbinsX();j++) HeTemplAglW[m][i]->SetBinContent(j,w3*HeTemplAgl[m]->GetBinContent(j)/i3*itot);
                        }
			else{
                                float i1 = PTemplAgl2[m]->Integral(0,cut1) + PTemplAgl2[m]->Integral(cut2,50);
                                float i2 = DTemplAgl2[m]->Integral(0,cut1) + DTemplAgl2[m]->Integral(cut2,50);
                                float i3=0;
                                float itot= Result->Integral();
                                if(He) i3 = HeTemplAgl2[m]->Integral(0,cut1) + HeTemplAgl2[m]->Integral(cut2,50);
				if(i1>0) for(int j=0; j<PTemplAgl[m]->GetNbinsX();j++) PTemplAglW[m][i]->SetBinContent(j,w1*PTemplAgl2[m]->GetBinContent(j)/i1*itot);
                                if(i2>0) for(int j=0; j<DTemplAgl[m]->GetNbinsX();j++) DTemplAglW[m][i]->SetBinContent(j,w2*DTemplAgl2[m]->GetBinContent(j)/i2*itot);
                                if(i3>0) for(int j=0; j<HeTemplAgl[m]->GetNbinsX();j++) HeTemplAglW[m][i]->SetBinContent(j,w3*HeTemplAgl2[m]->GetBinContent(j)/i3*itot);
			}
			Stack->Add(PTemplAglW[m][i]);
                        Stack->Add(DTemplAglW[m][i]);
                        if(He) Stack->Add(HeTemplAglW[m][i]);
                        DataDistAgl[i][m]->SetMarkerStyle(8);
                        Stack->Draw();
                        DataDistAgl[i][m]->Draw("epsame");
                        Result->SetLineColor(5);
                        Result->SetLineWidth(2);
                        Result->Draw("SAME");
                        float Cov01=fitTAgl[m][i]->GetFitter()->GetCovarianceMatrixElement(0,1);
                        float Cov02=fitTAgl[m][i]->GetFitter()->GetCovarianceMatrixElement(0,2);
                        float Cov12=fitTAgl[m][i]->GetFitter()->GetCovarianceMatrixElement(1,2);
                        float Sigma=pow((pow(e2/w2,2)+pow(e1/w1,2))/2,0.5);
                        ErrAgl[m][i]= Sigma*DTemplAglW[m][i]->Integral();
                        DCountsgeoAgl_Dist->SetBinContent(m+1,i+1,0,DTemplAglW[m][i]->Integral(0,50));
			DCountsgeoAgl_Dist->SetBinContent(m+1,i+1,1,ErrAgl[m][i]);
                	if(i==0){
                                PCountsAgl_Dist->SetBinContent(m+1,0,PTemplAglW[m][0]->Integral(0,50));
                                PCountsAgl_Dist->SetBinContent(m+1,1,Err[m][i]);
                        }

		}
                else{
                        DataDistAgl[i][m]->SetMarkerStyle(8);
                        DataDistAgl[i][m]->Draw("ep");
                        if(i!=0){PTemplAgl[m]->Draw("same");
                        	DTemplAgl[m]->Draw("same");
                        	if(He) HeTemplAgl[m]->Draw("same");
				}
			if(i==0) {PTemplAgl2[m]->Draw("same");
                                 DTemplAgl2[m]->Draw("same");
                                if(He) HeTemplAgl2[m]->Draw("same");
				PCountsAgl_Dist->SetBinContent(m+1,0,PTemplAglW[m][0]->Integral(0,50));
				}
                }

        }

	cout<<"********** TEMPL DIST: FIT RESULTS ********************"<<endl;
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
