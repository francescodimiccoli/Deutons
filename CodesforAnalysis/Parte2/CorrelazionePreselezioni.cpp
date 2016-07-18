#include "PlottingFunctions/Correlazione_Preselezioni_Plot.h"


TH2F * CorrelazionePreselezioni = new TH2F("CorrelazionePreselezioni","CorrelazionePreselezioni",11,0,11,11,0,11);
int Norm[11]= {0};


void Correlazione_Preselezioni(){
	for(int S=0;S<10;S++){
		if(((cmask.getMask()>>S)&1)==1){
			for(int F=0;F<10;F++) 
				if(((cmask.getMask()>>F)&1)==1) 
					CorrelazionePreselezioni->Fill(S,F);
			if(trgpatt.IsPhysical())  
				CorrelazionePreselezioni->Fill(S,10);
		}		
		if(((cmask.getMask()>>S)&1)==1) Norm[S]++;
	}
	if(trgpatt.IsPhysical()){
		for(int F=0;F<10;F++) 
			if(((cmask.getMask()>>F)&1)==1) 
				CorrelazionePreselezioni->Fill(10,F);
		if(trgpatt.IsPhysical()) {
			CorrelazionePreselezioni->Fill(10,10);
			Norm[10]++;
		}
	}
}

void Correlazione_Preselezioni_Write(){
        CorrelazionePreselezioni->Write();
        return;
}


void Correlazione_Preselezioni(string filename){
	
	cout<<"***************** CORR. PRESELEZIONI ********************"<<endl;

	cout<<"*** Reading  P1 file ****"<<endl;
        TFile * inputHistoFile =TFile::Open(filename.c_str(),"READ");

	TH2F * CorrelazionePreselezioni = (TH2F*) inputHistoFile->Get("CorrelazionePreselezioni");
	
	cout<<"***************** CORR. PRESELEZIONI ********************"<<endl;	
	for(int S=0;S<11;S++)
		for(int F=0;F<11;F++) 
			if(Norm[S]>0) CorrelazionePreselezioni->SetBinContent(S+1,F+1,CorrelazionePreselezioni->GetBinContent(S+1,F+1)/(float)Norm[S]);

	
        finalHistos.Add(CorrelazionePreselezioni);
	finalHistos.writeObjsInFolder("Results");

        cout<<"*** Plotting ...  ****"<<endl;

	Correlazione_Preselezioni_Plot(CorrelazionePreselezioni);
	return;
}

