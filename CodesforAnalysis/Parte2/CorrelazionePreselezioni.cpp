TH2F * CorrelazionePreselezioni = new TH2F("CorrelazionePreselezioni","CorrelazionePreselezioni",11,0,11,11,0,11);
int Norm[11]= {0};


void Correlazione_Preselezioni(TNtuple *ntupla, int l){
	 ntupla->GetEvent(l);
	for(int S=0;S<10;S++){
		if(((cmask.getMask()>>S)&1)==1){
			for(int F=0;F<10;F++) 
				if(((cmask.getMask()>>F)&1)==1) 
					CorrelazionePreselezioni->Fill(S,F);
			if(Tup.Unbias==0)  
				CorrelazionePreselezioni->Fill(S,10);
		}		
		if(((cmask.getMask()>>S)&1)==1) Norm[S]++;
	}
	if(Tup.Unbias==0){
		for(int F=0;F<10;F++) 
			if(((cmask.getMask()>>F)&1)==1) 
				CorrelazionePreselezioni->Fill(10,F);
		if(Tup.Unbias==0) {
			CorrelazionePreselezioni->Fill(10,10);
			Norm[10]++;
		}
	}
}

void Correlazione_Preselezioni_Write(){
        CorrelazionePreselezioni->Write();
        return;
}


void Correlazione_Preselezioni(TFile * inputHistoFile){
	TH2F * CorrelazionePreselezioni = (TH2F*) inputHistoFile->Get("CorrelazionePreselezioni");
	
	cout<<"***************** CORR. PRESELEZIONI ********************"<<endl;	
	for(int S=0;S<11;S++)
		for(int F=0;F<11;F++) 
			if(Norm[S]>0) CorrelazionePreselezioni->SetBinContent(S+1,F+1,CorrelazionePreselezioni->GetBinContent(S+1,F+1)/(float)Norm[S]);
	
	cout<<"*** Updating P1 file ****"<<endl;
	inputHistoFile->ReOpen("UPDATE");
	if (! inputHistoFile->GetDirectory("Results")) {
		inputHistoFile->mkdir("Results");
	}
	inputHistoFile->cd("Results");
   CorrelazionePreselezioni->Write();
	inputHistoFile-> Write();
	

	TCanvas *c13=new TCanvas("Correlazione Selezioni");
	c13->cd();
	CorrelazionePreselezioni->Draw("col");

	cout<<"*** Updating Results file ***"<<endl;
		cout << (fileFinalPlots->IsOpen()?"Opened  ":"Closed  ") << (const char *)fileFinalPlots->GetOption() << endl;
        fileFinalPlots->cd("MC Results");
        c13->Write();
        fileFinalPlots->Write();


}
