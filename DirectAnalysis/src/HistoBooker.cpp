#include "HistoBooker.h"

//single histo
SingleHisto::SingleHisto(std::string CollectionName,int nbinsx, float xmin, float xmax,std::string Cut){
		
	TH1F * tmp = new TH1F((CollectionName).c_str(),(CollectionName).c_str(),nbinsx,xmin,xmax);
	Histos.push_back(tmp);
	collectionname = CollectionName;
	cut  = Cut;

};

void SingleHisto::Save(FileSaver finalhisto){

	for(int bin=0;bin<Histos.size();bin++)
		finalhisto.Add(Histos[bin]);
	finalhisto.writeObjsInFolder(collectionname.c_str());
}

void SingleHisto::FillEventByEventData(Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars)){
	if(ApplyCuts(cut,vars)){
     		float weight=1;
               if(ApplyCuts("IsMC",vars)) weight = vars->mcweight;
               else weight = vars->PrescaleFactor;
	    Histos[0]->Fill(var(vars),weight);	
        }
}

//single scatterplot
SingleScatter::SingleScatter(std::string CollectionName, int nbinsx, float xmin, float xmax,int nbinsy, float ymin, float ymax,std::string Cut):SingleHisto(CollectionName,nbinsx,xmin,xmax,Cut){
	Histos.clear();
	TH2F * tmp = new TH2F((CollectionName).c_str(),(CollectionName).c_str(),nbinsx,xmin,xmax,nbinsy,ymin,ymax);
	Histos.push_back(tmp);
};

void SingleScatter::FillEventByEventScatter(Variables * vars, float (*var) (Variables * vars), float (*secondvar) (Variables * vars),float (*discr_var) (Variables * vars)){
	
	 if(ApplyCuts(cut,vars)){
        	 float weight=1;
        	cout<<"ecco"<<endl;
	       if(ApplyCuts("IsMC",vars)) weight = vars->mcweight;
               else weight = vars->PrescaleFactor;
	       ((TH2F*)Histos[0])->Fill(var(vars),secondvar(vars),weight);
        }
}

//binned histogram
BinnedHisto::BinnedHisto(std::string CollectionName, Binning Bins,std::string Cut):SingleHisto(CollectionName,0,0,0,Cut){
	Histos.clear();
	TH1F * tmp = new TH1F((CollectionName).c_str(),(CollectionName).c_str(),Bins.size(),0,Bins.size());
	Histos.push_back(tmp);
	collectionname = CollectionName;
	cut  = Cut;
	bins = Bins;
}
void BinnedHisto::FillEventByEventData(Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars)){

	int kbin;
        kbin =  bins.GetBin(var(vars));
        if(ApplyCuts(cut,vars)&&kbin>0){
		float weight=1;
	       if(ApplyCuts("IsMC",vars)) weight = vars->mcweight;
	       else weight = vars->PrescaleFactor;		 
               Histos[0]->Fill(kbin,weight);
        }
}

//binned scatter
BinnedScatter::BinnedScatter(std::string CollectionName, Binning Bins,std::string Cut):SingleHisto(CollectionName,0,0,0,Cut){
	Histos.clear();
	TH2F * tmp = new TH2F((CollectionName).c_str(),(CollectionName).c_str(),Bins.size(),0,Bins.size(),Bins.size(),0,Bins.size());
	Histos.push_back(tmp);
	collectionname = CollectionName;
	cut  = Cut;
	bins = Bins;	
}
void BinnedScatter::FillEventByEventData(Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars)){

	int kbin1;
        kbin1 =  bins.GetBin(var(vars));
       	int kbin2;
        kbin2 =  bins.GetBin(discr_var(vars));
	 if(ApplyCuts(cut,vars)&&kbin1>0&&kbin2>0){
			float weight=1;
	               	if(ApplyCuts("IsMC",vars)) weight = vars->mcweight;
               		else weight = vars->PrescaleFactor;
			((TH2F*)Histos[0])->Fill(kbin1,kbin2,weight);
        }
}


//collection of histos
HistoBinCollection::HistoBinCollection(std::string CollectionName, Binning Bins, int nbinsx, float xmin, float xmax,std::string Cut):SingleHisto(CollectionName,0,0,0,Cut){
	Histos.clear();
	for(int bin=0;bin<Bins.size();bin++){
		TH1F * tmp = new TH1F((CollectionName+"_bin"+to_string(bin)).c_str(),(CollectionName+"_bin"+to_string(bin)).c_str(),nbinsx,xmin,xmax);
		Histos.push_back(tmp);	
	}
	bins = Bins;
};

void HistoBinCollection::FillEventByEventData(Variables * vars, float (*var) (Variables * vars),float (*discr_var) (Variables * vars)){
	int kbin;
	kbin =  bins.GetBin(discr_var(vars));
	if(ApplyCuts(cut,vars)&&kbin>0){
		float weight=1;
	       if(ApplyCuts("IsMC",vars)) weight = vars->mcweight;
	       else weight = vars->PrescaleFactor;		 
        
		Histos[kbin]->Fill(var(vars),weight);
	}		
}

//collection of scatters
ScatterBinCollection::ScatterBinCollection(std::string CollectionName, Binning Bins, int nbinsx, float xmin, float xmax, int nbinsy, float ymin, float ymax,std::string Cut):SingleHisto(CollectionName,nbinsx,xmin,xmax,Cut){
	Histos.clear();
	for(int bin=0;bin<Bins.size();bin++){
		TH2F * tmp = new TH2F((CollectionName+"_bin"+to_string(bin)).c_str(),(CollectionName+"_bin"+to_string(bin)).c_str(),nbinsx,xmin,xmax,nbinsy,ymin,ymax);
		Histos.push_back(tmp);	
	}
	bins = Bins;
};

void ScatterBinCollection::FillEventByEventScatter(Variables * vars, float (*var) (Variables * vars),float (*secondvar) (Variables * vars),float (*discr_var) (Variables * vars)){
	int kbin;
	kbin =  bins.GetBin(discr_var(vars));
	if(ApplyCuts(cut,vars)&&kbin>0){
		 float weight=1;
               if(ApplyCuts("IsMC",vars)) weight = vars->mcweight;
               else weight = vars->PrescaleFactor;
		((TH2F*) Histos[kbin])->Fill(var(vars),secondvar(vars),weight);
	}		
}

HistoBooker::HistoBooker(){ return;}

void HistoBooker::BookBinCollection(std::string CollectionName, Binning Bins, int nbinsx, float xmin, float xmax,std::string Cut, GetFillinVariable var, GetDiscrimVariable discr_var){
	HistoBinCollection * tmp = new HistoBinCollection(CollectionName,Bins,nbinsx,xmin,xmax,Cut);
	Histos.push_back(tmp);
	Filler.AddObject2beFilled(tmp,var,discr_var);	
};

void HistoBooker::BookSingleHisto(std::string CollectionName, int nbinsx, float xmin, float xmax,std::string Cut, GetFillinVariable var, GetDiscrimVariable discr_var){
	SingleHisto * tmp = new SingleHisto(CollectionName,nbinsx,xmin,xmax,Cut);
	Histos.push_back(tmp);
        Filler.AddObject2beFilled(tmp,var,discr_var);
};

void HistoBooker::BookSingleScatter(std::string CollectionName, int nbinsx, float xmin, float xmax, int nbinsy, float ymin, float ymax, std::string Cut, GetFillinVariable var, GetFillinVariable secondvar,GetDiscrimVariable discr_var){
	SingleScatter * tmp = new SingleScatter(CollectionName,nbinsx,xmin,xmax,nbinsy,ymin,ymax,Cut);
	Histos.push_back(tmp);
        Filler.AddObject2beFilled(tmp,var,discr_var,"",secondvar);
};

void HistoBooker::BookBinnedHisto(std::string CollectionName, Binning Bins,std::string Cut, GetFillinVariable var, GetDiscrimVariable discr_var){
	BinnedHisto * tmp = new BinnedHisto(CollectionName,Bins,Cut);
	Histos.push_back(tmp);
        Filler.AddObject2beFilled(tmp,var,discr_var);
};

void HistoBooker::BookBinnedScatter(std::string CollectionName, Binning Bins, std::string Cut, GetFillinVariable var, GetDiscrimVariable discr_var){
	BinnedScatter * tmp = new BinnedScatter(CollectionName,Bins,Cut);
	Histos.push_back(tmp);
        Filler.AddObject2beFilled(tmp,var,discr_var);
};


void HistoBooker::BookScatterBinCollection(std::string CollectionName, Binning Bins, int nbinsx, float xmin, float xmax,int nbinsy, float ymin, float ymax,std::string Cut,GetFillinVariable var,GetFillinVariable secondvar, GetDiscrimVariable discr_var){
	ScatterBinCollection * tmp = new ScatterBinCollection(CollectionName,Bins,nbinsx,xmin,xmax,nbinsy,ymin,ymax,Cut);
	Histos.push_back(tmp);
        Filler.AddObject2beFilled(tmp,var,discr_var,"",secondvar);
};
void HistoBooker::FillEverything(DBarReader reader){
	Variables * vars = new Variables;	
	Filler.LoopOnGeneric(reader,vars);	
};

void HistoBooker::SaveEverything(FileSaver finalhisto){
	for(int n=0;n<Histos.size();n++)
		Histos[n]->Save(finalhisto);	
};
