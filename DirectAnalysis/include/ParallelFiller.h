#ifndef PARALLELFILER_H
#define PARALLELFILER_H

#include "BadEventSimulator.h"
#include "Livetime.h"
#include "filesaver.h"
#include "DBarReader.h"

template <class T> 
class ParallelFiller{

	typedef float (*GetFillinVariable)  (Variables * vars);
	typedef float (*GetDiscrimVariable) (Variables * vars);

	private:
	std::vector<T> Objects2beFilled;
	std::vector<GetFillinVariable>  FillinVariables;
	std::vector<GetDiscrimVariable> DiscrimVariables;
	std::vector<std::string> 	Cuts;
	bool Refill = true;	

	public:
	void AddObject2beFilled(T Object,GetFillinVariable var=0x0, GetDiscrimVariable discr_var=0x0,std::string cut=""){
		Objects2beFilled.push_back(Object);
		FillinVariables.push_back(var);
		DiscrimVariables.push_back(discr_var);
		Cuts.push_back(cut);
		return;
	}
	void ReinitializeAll(bool refill){
		bool checkifsomeismissing=false;
			for(int nobj=0;nobj<Objects2beFilled.size();nobj++) {
				if(!(Objects2beFilled[nobj]->ReinitializeHistos(refill))) {checkifsomeismissing=true;
					cout<<"missing: "<<nobj<<endl;}
				}
		if(checkifsomeismissing||refill) Refill=true;
		else Refill=false;
		if(Refill) cout<<"Refillo"<<endl;
	}
	
	void LoopOnMC(TTree * treeMC, Variables * vars){
		if(!Refill) return;
		cout<<" MC Filling ..."<< endl;
	        vars->ReadBranches(treeMC);
		for(int i=0;i<treeMC->GetEntries();i++){
			if(i%(int)FRAC!=0) continue;
			UpdateProgressBar(i, treeMC->GetEntries());
			treeMC->GetEvent(i);
			vars->Update();
			for(int nobj=0;nobj<Objects2beFilled.size();nobj++){
				Objects2beFilled[nobj]->LoadEventIntoBadEvSim(vars);
				Objects2beFilled[nobj]->FillEventByEventMC(vars,FillinVariables[nobj],DiscrimVariables[nobj]);
			}
		}
	}

	void LoopOnData(TTree * treeDT, Variables * vars){
		if(!Refill) return;
		cout<<" DATA Filling ..."<< endl;
		vars->ReadBranches(treeDT);
		for(int i=0;i<treeDT->GetEntries()/FRAC;i++){
			UpdateProgressBar(i, treeDT->GetEntries()/FRAC);
			treeDT->GetEvent(i);
			vars->Update();
			for(int nobj=0;nobj<Objects2beFilled.size();nobj++) Objects2beFilled[nobj]->FillEventByEventData(vars,FillinVariables[nobj],DiscrimVariables[nobj]);
		}
	}
	void LoopOnMC(DBarReader treeMC, Variables * vars){
		if(!Refill) return;
		cout<<" MC Filling ..."<< endl;
		for(int i=0;i<treeMC.GetTreeEntries();i++){
			if(i%(int)FRAC!=0) continue;
			UpdateProgressBar(i, treeMC.GetTreeEntries());
			treeMC.FillVariables(i,vars);
			vars->Update();
			for(int nobj=0;nobj<Objects2beFilled.size();nobj++){
				Objects2beFilled[nobj]->LoadEventIntoBadEvSim(vars);
				Objects2beFilled[nobj]->FillEventByEventMC(vars,FillinVariables[nobj],DiscrimVariables[nobj]);
			}
		}
	}

	void LoopOnData(DBarReader treeDT, Variables * vars){
		if(!Refill) return;
		cout<<" DATA Filling ..."<< endl;
		for(int i=0;i<treeDT.GetTreeEntries()/FRAC;i++){
			UpdateProgressBar(i, treeDT.GetTreeEntries()/FRAC);
			treeDT.FillVariables(i,vars);
			vars->Update();
			for(int nobj=0;nobj<Objects2beFilled.size();nobj++) Objects2beFilled[nobj]->FillEventByEventData(vars,FillinVariables[nobj],DiscrimVariables[nobj]);
		}
	}

	void ExposureTimeFilling(TTree * treeDT, Variables * vars,FileSaver finalhistos){
		if(!Refill) return;
                cout<<" Exposure Time Filling ..."<< endl;
                int ActualTime=0;
		vars->ReadBranches(treeDT);
		for(int i=0;i<treeDT->GetEntries()/FRAC;i++){
			UpdateProgressBar(i, treeDT->GetEntries()/FRAC);
			treeDT->GetEvent(i);
			vars->Update();
			if((int)vars->U_time!=ActualTime) {
				for(int nobj=0;nobj<Objects2beFilled.size();nobj++) 
					UpdateZoneLivetime(vars->Livetime,vars->Rcutoff,Objects2beFilled[nobj]->GetExposureTime(),Objects2beFilled[nobj]->GetBins());
				ActualTime=vars->U_time;
			}
		}
		for(int nobj=0;nobj<Objects2beFilled.size();nobj++){ 
			finalhistos.Add(Objects2beFilled[nobj]->GetExposureTime());
			finalhistos.writeObjsInFolder(("Fluxes/"+Objects2beFilled[nobj]->GetName()).c_str());
		}
	}


	void LoopOnMCForGenAcceptance(TTree * treeMC, Variables * vars,FileSaver finalhistos){
		if(!Refill) return;
                cout<<" Geom. Acceptance Filling ..."<< endl;	
		vars->ReadBranches(treeMC);
                int kbin;
		for(int i=0;i<treeMC->GetEntries()/FRAC;i++){
                                if(i%(int)FRAC!=0) continue;
                                UpdateProgressBar(i, treeMC->GetEntries());
                                treeMC->GetEvent(i);
                                vars->Update();
				for(int nobj=0;nobj<Objects2beFilled.size();nobj++) {
					kbin = Objects2beFilled[nobj]->GetBins().GetBin(DiscrimVariables[nobj](vars));
					if(ApplyCuts(Cuts[nobj],vars)&&kbin>0) 
						Objects2beFilled[nobj]->GetGenAcceptance()->Fill(kbin);
					kbin = PRB.GetBin(GetGenMomentum(vars));
					if(ApplyCuts(Cuts[nobj],vars)&&kbin>0)
						 Objects2beFilled[nobj]->GetTriggerCounts()->Fill(kbin);
				}
		}
		for(int nobj=0;nobj<Objects2beFilled.size();nobj++){
			finalhistos.Add(Objects2beFilled[nobj]->GetGenAcceptance());
                        finalhistos.Add(Objects2beFilled[nobj]->GetTriggerCounts());
                        finalhistos.writeObjsInFolder(("Fluxes/"+Objects2beFilled[nobj]->GetName()).c_str());
		}		

	}

};

#endif
