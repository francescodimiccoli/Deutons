#ifndef PARALLELFILER_H
#define PARALLELFILER_H

#include <TF1.h>
#include "Globals.h"
#include "BadEventSimulator.h"
#include "binning.h"
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
    std::vector<GetFillinVariable>  SecondFillinVariables;
    std::vector<GetDiscrimVariable> DiscrimVariables;
    std::vector<std::string> 	Cuts;
	
    bool Refill = true;	

    public:
    bool GetRefillFlag() {return Refill;}
    T GetObject(int i) { return Objects2beFilled[i]; }
    int GetObjectNr() {return Objects2beFilled.size();}
    void AddObject2beFilled(const T Object,GetFillinVariable var=0x0, GetDiscrimVariable discr_var=0x0,std::string cut="",GetFillinVariable second_var=0x0){
        Objects2beFilled.push_back(Object);
        FillinVariables.push_back(var);
        DiscrimVariables.push_back(discr_var);
	SecondFillinVariables.push_back(second_var);
        Cuts.push_back(cut);
        return;
    }
    void ReinitializeAll(bool refill){
        bool checkifsomeismissing=false;
        for(int nobj=0;nobj<Objects2beFilled.size();nobj++) {
            if(!(Objects2beFilled[nobj]->ReinitializeHistos(refill))) {checkifsomeismissing=true;
                cout<<"missing: "<<nobj<<endl;}
        }
        if(checkifsomeismissing||refill){ Refill = true;
        	 cout<<"Refillo"<<endl;
		}
	else Refill = false;
    }


	//// Looping on simple ntuples

    void LoopOnMC(TTree * treeMC, Variables * vars){
        if(!Refill) return;
        cout<<" MC Filling from Ntuples ..."<< endl;
        vars->ReadBranches(treeMC);
	for(int i=0;i<treeMC->GetEntries();i++){
	    if(i%(int)FRAC!=0) continue;
            UpdateProgressBar(i, treeMC->GetEntries());
     	    vars->ResetVariables();
            treeMC->GetEvent(i);
	    //vars->Update();
            for(int nobj=0;nobj<Objects2beFilled.size();nobj++){
                Objects2beFilled[nobj]->LoadEventIntoBadEvSim(vars);
                Objects2beFilled[nobj]->FillEventByEventMC(vars,FillinVariables[nobj],DiscrimVariables[nobj]);
            }
        }
    }

    void LoopOnData(TTree * treeDT, Variables * vars){
        if(!Refill) return;
        cout<<" DATA Filling from NTuples ..."<< endl;
        vars->ReadBranches(treeDT);
        for(int i=0;i<treeDT->GetEntries();i++){
	    if(i%(int)FRAC!=0) continue; // WTF ?!
            UpdateProgressBar(i, treeDT->GetEntries());		
	    vars->ResetVariables();	
            treeDT->GetEvent(i);
	    //vars->PrintCurrentState();	
	    //vars->Update();
            for(int nobj=0;nobj<Objects2beFilled.size();nobj++) Objects2beFilled[nobj]->FillEventByEventData(vars,FillinVariables[nobj],DiscrimVariables[nobj]);
        }
    }


  void LoopOnGeneric(TTree * treeDT, Variables * vars){
        if(!Refill) return;
        cout<<" DATA Filling from NTuples ..."<< endl;
        vars->ReadBranches(treeDT);
        for(int i=0;i<treeDT->GetEntries();i++){
            if(i%(int)FRAC!=0) continue; // WTF ?!
            UpdateProgressBar(i, treeDT->GetEntries());
            vars->ResetVariables();
    		treeDT->GetEvent(i);
            //vars->Update();
           for(int nobj=0;nobj<Objects2beFilled.size();nobj++) {
			if(SecondFillinVariables[nobj]) Objects2beFilled[nobj]->FillEventByEventScatter(vars,FillinVariables[nobj],SecondFillinVariables[nobj],DiscrimVariables[nobj]);
			else Objects2beFilled[nobj]->FillEventByEventData(vars,FillinVariables[nobj],DiscrimVariables[nobj]);
       		}
	 }
    }





   //looping On DBar Ntuples


 
   void LoopOnMCVariables(DBarReader readerMC_P, DBarReader readerMC_D,DBarReader readerMC_He, Variables * vars){
        if(!Refill) return;
        cout<<" MC Filling ..."<< endl;
	if(!readerMC_P.GetTree()) return;
        if(readerMC_P.GetTree()->GetNbranches()>11) {LoopOnMC(readerMC_P.GetTree(),vars); return;}
	else{
	for(int i=0;i<readerMC_P.GetTreeEntries();i++){
              UpdateProgressBar(i, readerMC_P.GetTreeEntries());
	    if(i%(int)FRAC!=0) continue; // WTF ?!
            vars->ResetVariables();
	    readerMC_P.FillVariables(i,vars);	
            vars->Update();
	    //vars->PrintCurrentState();
            for(int nobj=0;nobj<Objects2beFilled.size();nobj++){
                Objects2beFilled[nobj]->LoadEventIntoBadEvSim(vars);
                Objects2beFilled[nobj]->FillEventByEventMC(vars,FillinVariables[nobj],DiscrimVariables[nobj]);
            }
        }
	

	for(int i=0;i<readerMC_D.GetTreeEntries();i++){
            UpdateProgressBar(i, readerMC_D.GetTreeEntries());
           if(i%(int)FRAC!=0) continue; // WTF ?!
             vars->ResetVariables();
	    readerMC_D.FillVariables(i,vars);	
            vars->Update();
	    //vars->PrintCurrentState();
            for(int nobj=0;nobj<Objects2beFilled.size();nobj++){
                Objects2beFilled[nobj]->LoadEventIntoBadEvSim(vars);
                Objects2beFilled[nobj]->FillEventByEventMC(vars,FillinVariables[nobj],DiscrimVariables[nobj]);
            }
        }
	
 
	for(int i=0;i<readerMC_He.GetTreeEntries();i++){
           UpdateProgressBar(i, readerMC_He.GetTreeEntries());
            if(i%(int)FRAC!=0) continue; // WTF ?!
             vars->ResetVariables();
	    readerMC_He.FillVariables(i,vars);	
            vars->Update();
	    //vars->PrintCurrentState();
            for(int nobj=0;nobj<Objects2beFilled.size();nobj++){
                Objects2beFilled[nobj]->LoadEventIntoBadEvSim(vars);
                Objects2beFilled[nobj]->FillEventByEventMC(vars,FillinVariables[nobj],DiscrimVariables[nobj]);
            }
        }
	}
 
   }
   	 
    void LoopOnMC(DBarReader readerMC, Variables * vars){
        if(!Refill) return;
        cout<<" MC Filling ..."<< endl;
	if(!readerMC.GetTree()) return;
        if(readerMC.GetTree()->GetNbranches()>11) {LoopOnMC(readerMC.GetTree(),vars); return;}
	else{
	for(int i=0;i<readerMC.GetCompactEntries();i++){
            
	    if(i%(int)FRAC!=0) continue; // WTF ?!
            //UpdateProgressBar(i, readerMC.GetCompactEntries());
            vars->ResetVariables();
	    readerMC.FillCompact(i,vars);	
            vars->Update();
	   //vars->PrintCurrentState();
				

 	for(int nobj=0;nobj<Objects2beFilled.size();nobj++){
                Objects2beFilled[nobj]->LoadEventIntoBadEvSim(vars);
                Objects2beFilled[nobj]->FillEventByEventMC(vars,FillinVariables[nobj],DiscrimVariables[nobj]);
            }
        }
	}
    }
 
    void LoopOnMC(DBarReader readerMC_P, DBarReader readerMC_D,DBarReader readerMC_He, Variables * vars){
	    if(!Refill) return;
	    cout<<" MC Filling ..."<< endl;
	    if(!readerMC_P.GetTree()) return;
	    if(readerMC_P.GetTree()->GetNbranches()>11) {LoopOnMC(readerMC_P.GetTree(),vars); return;}
	    else{
		    for(int i=0;i<readerMC_P.GetCompactEntries();i++){
			    if(i%(int)FRAC!=0) continue; // WTF ?!
			   // UpdateProgressBar(i, readerMC_P.GetCompactEntries());
			    vars->ResetVariables();
			    readerMC_P.FillCompact(i,vars,0.938);	
			    vars->Update();
			    //vars->PrintCurrentState();
			    for(int nobj=0;nobj<Objects2beFilled.size();nobj++){
				    Objects2beFilled[nobj]->LoadEventIntoBadEvSim(vars);
				    Objects2beFilled[nobj]->FillEventByEventMC(vars,FillinVariables[nobj],DiscrimVariables[nobj]);
			    }
		    }
		
		    for(int i=0;i<readerMC_D.GetCompactEntries();i++){
			    if(i%(int)FRAC!=0) continue; // WTF ?!
			    UpdateProgressBar(i, readerMC_D.GetCompactEntries());
			    vars->ResetVariables();
			    readerMC_D.FillCompact(i,vars,1.875);	
			    vars->Update();
			    //vars->PrintCurrentState();
			    for(int nobj=0;nobj<Objects2beFilled.size();nobj++){
				    Objects2beFilled[nobj]->LoadEventIntoBadEvSim(vars);
				    Objects2beFilled[nobj]->FillEventByEventMC(vars,FillinVariables[nobj],DiscrimVariables[nobj]);
			    }
		    }
		
		    for(int i=0;i<readerMC_He.GetCompactEntries();i++){
			    if(i%(int)FRAC!=0) continue; // WTF ?!
		//	    UpdateProgressBar(i, readerMC_He.GetCompactEntries());
			    vars->ResetVariables();
			    readerMC_He.FillCompact(i,vars,4*0.938);	
			    vars->Update();
			    //vars->PrintCurrentState();
			    for(int nobj=0;nobj<Objects2beFilled.size();nobj++){
				    Objects2beFilled[nobj]->LoadEventIntoBadEvSim(vars);
				    Objects2beFilled[nobj]->FillEventByEventMC(vars,FillinVariables[nobj],DiscrimVariables[nobj]);
			    }
		    }


	    }
    }

    void LoopOnData(DBarReader readerDT, Variables * vars){
        if(!Refill) return;
	if(!readerDT.GetTree() && !readerDT.GetCompactTree()) return;
	if(readerDT.GetTree()) if(readerDT.GetTree()->GetNbranches()>11) {LoopOnData(readerDT.GetTree(),vars); return;}
	cout<<" DATA Filling ..."<< endl;
        for(int i=0;i<readerDT.GetCompactEntries();i++){
	    UpdateProgressBar(i, readerDT.GetCompactEntries());
            if(i%(int)FRAC!=0) continue; // WTF ?!
             vars->ResetVariables();
            readerDT.FillCompact(i,vars); 
            vars->Update();
	    //vars->PrintCurrentState();
            for(int nobj=0;nobj<Objects2beFilled.size();nobj++){ 
	        if(vars->isinsaa==0&&vars->good_RTI==0) {
			Objects2beFilled[nobj]->FillEventByEventData(vars,FillinVariables[nobj],DiscrimVariables[nobj]);
			}
		}
	}
    }

    void LoopOnDataVariables(DBarReader readerDT, Variables * vars){
        if(!Refill) return;
	if(!readerDT.GetTree() && !readerDT.GetCompactTree()) return;
	if(readerDT.GetTree()) if(readerDT.GetTree()->GetNbranches()>11) {LoopOnData(readerDT.GetTree(),vars); return;}
	cout<<" DATA Filling ..."<< endl;
        for(int i=0;i<readerDT.GetTreeEntries();i++){
	    if(i%(int)FRAC!=0) continue; // WTF ?!
            UpdateProgressBar(i, readerDT.GetTreeEntries());
            vars->ResetVariables();
            readerDT.FillVariables(i,vars); 
            vars->Update();
	    //vars->PrintCurrentState();
            for(int nobj=0;nobj<Objects2beFilled.size();nobj++){ 
	        if(vars->isinsaa==0&&vars->good_RTI==0) {
			Objects2beFilled[nobj]->FillEventByEventData(vars,FillinVariables[nobj],DiscrimVariables[nobj]);
			}
		}
	}
    }


    void LoopOnGeneric(DBarReader reader, Variables * vars){
         if(!Refill) return;
        if(!reader.GetTree() && !reader.GetCompactTree()) return;
        if(reader.GetTree()) if(reader.GetTree()->GetNbranches()>11) {LoopOnGeneric(reader.GetTree(),vars); return;}
        else
        cout<<" Generic Filling ..."<< endl;
        for(int i=0;i<reader.GetCompactEntries();i++){
	    if(i%(int)FRAC!=0) continue; // WTF ?!
            UpdateProgressBar(i, reader.GetCompactEntries());
     	    vars->ResetVariables();	
	    //reader.FillVariables(i,vars);
            reader.FillCompact(i,vars);
	    vars->Update();
            //vars->PrintCurrentState();
	    for(int nobj=0;nobj<Objects2beFilled.size();nobj++) 
                if(SecondFillinVariables[nobj]) Objects2beFilled[nobj]->FillEventByEventScatter(vars,FillinVariables[nobj],SecondFillinVariables[nobj],DiscrimVariables[nobj]);
		else Objects2beFilled[nobj]->FillEventByEventData(vars,FillinVariables[nobj],DiscrimVariables[nobj]);
        }
    }

    void LoopOnGeneric(DBarReader reader1,DBarReader reader2,DBarReader reader3, Variables * vars){
	    if(!Refill) return;
	    if(!reader1.GetTree() && !reader1.GetCompactTree()) return;
	    if(reader1.GetTree()) if(reader1.GetTree()->GetNbranches()>11) {LoopOnGeneric(reader1.GetTree(),vars); return;}
	    else
		    cout<<" Generic Filling MC ..."<< endl;
	    for(int i=0;i<reader1.GetCompactEntries();i++){
		    if(i%(int)FRAC!=0) continue; // WTF ?!
		    UpdateProgressBar(i, reader1.GetCompactEntries());
		    vars->ResetVariables();	
		    //reader.FillVariables(i,vars);
		    reader1.FillCompact(i,vars,0.938);
		    vars->Update();
		    //vars->PrintCurrentState();
		    for(int nobj=0;nobj<Objects2beFilled.size();nobj++) 
			    if(SecondFillinVariables[nobj]) Objects2beFilled[nobj]->FillEventByEventScatter(vars,FillinVariables[nobj],SecondFillinVariables[nobj],DiscrimVariables[nobj]);
			    else Objects2beFilled[nobj]->FillEventByEventData(vars,FillinVariables[nobj],DiscrimVariables[nobj]);
	    }

	    for(int i=0;i<reader2.GetCompactEntries();i++){
		    if(i%(int)FRAC!=0) continue; // WTF ?!
		    UpdateProgressBar(i, reader2.GetCompactEntries());
		    vars->ResetVariables();	
		    //reader.FillVariables(i,vars);
		    reader2.FillCompact(i,vars,1.875);
		    vars->Update();
		    //vars->PrintCurrentState();
		    for(int nobj=0;nobj<Objects2beFilled.size();nobj++) 
			    if(SecondFillinVariables[nobj]) Objects2beFilled[nobj]->FillEventByEventScatter(vars,FillinVariables[nobj],SecondFillinVariables[nobj],DiscrimVariables[nobj]);
			    else Objects2beFilled[nobj]->FillEventByEventData(vars,FillinVariables[nobj],DiscrimVariables[nobj]);
	    }
		
	    for(int i=0;i<reader3.GetCompactEntries();i++){
		    if(i%(int)FRAC!=0) continue; // WTF ?!
		    UpdateProgressBar(i, reader3.GetCompactEntries());
		    vars->ResetVariables();	
		    //reader.FillVariables(i,vars);
		    reader3.FillCompact(i,vars,4*0.938);
		    vars->Update();
		    //vars->PrintCurrentState();
		    for(int nobj=0;nobj<Objects2beFilled.size();nobj++) 
			    if(SecondFillinVariables[nobj]) Objects2beFilled[nobj]->FillEventByEventScatter(vars,FillinVariables[nobj],SecondFillinVariables[nobj],DiscrimVariables[nobj]);
			    else Objects2beFilled[nobj]->FillEventByEventData(vars,FillinVariables[nobj],DiscrimVariables[nobj]);
	    }

    }


    void ExposureTimeFilling_RTI(DBarReader reader, Variables * vars){
	cout<<reader.GetTree()<<endl;
	    if(!Refill) return;
	    if(reader.GetTree()->GetNbranches()>11) {ExposureTimeFilling(reader.GetTree(),vars); return;}
	    else
		    cout<<" Exposure Time Filling ..."<< endl;
	    int n_of_seconds =0;	
	    for(int i=0;i<reader.GetTreeEntries();i++){
		    UpdateProgressBar(i, reader.GetTreeEntries());
		    reader.FillVariables(i,vars);
		    vars->Update();
		    for(int nobj=0;nobj<Objects2beFilled.size();nobj++) 
			    if(vars->good_RTI==0&&vars->isinsaa==0) UpdateZoneLivetime(vars->Livetime_RTI,vars->Rcutoff_IGRFRTI,Objects2beFilled[nobj]->GetExposureTime(),Objects2beFilled[nobj]->GetBins());
		    if(vars->good_RTI==0&&vars->isinsaa==0) n_of_seconds++;
		   	
	    }
	    cout<<" Seconds in te RUN: "<<n_of_seconds<<endl;
	    for(int nobj=0;nobj<Objects2beFilled.size();nobj++){ 
		    Objects2beFilled[nobj]->GetOutFileSaver().Add(Objects2beFilled[nobj]->GetExposureTime());
		    Objects2beFilled[nobj]->GetOutFileSaver().writeObjsInFolder(("Fluxes/"+Objects2beFilled[nobj]->GetName()).c_str());
	    }
    }

   void ExposureTimeFilling(TTree * treeDT, Variables * vars){
        if(!Refill) return;
        cout<<" Exposure Time Filling from Ntuples ..."<< endl;
        int ActualTime=0;
        vars->ReadBranches(treeDT);
        for(int i=0;i<treeDT->GetEntries();i++){
            UpdateProgressBar(i, treeDT->GetEntries());
            treeDT->GetEvent(i);
            if((int)vars->U_time!=ActualTime) {
                for(int nobj=0;nobj<Objects2beFilled.size();nobj++)
                   if(vars->R>0) UpdateZoneLivetime(vars->Livetime,vars->Rcutoff_IGRFRTI,Objects2beFilled[nobj]->GetExposureTime(),Objects2beFilled[nobj]->GetBins());
                ActualTime=vars->U_time;
                    }
        }
        for(int nobj=0;nobj<Objects2beFilled.size();nobj++){
            Objects2beFilled[nobj]->GetOutFileSaver().Add(Objects2beFilled[nobj]->GetExposureTime());
            Objects2beFilled[nobj]->GetOutFileSaver().writeObjsInFolder(("Fluxes/"+Objects2beFilled[nobj]->GetName()).c_str());
        }
    }



    void LoopOnMCForGenAcceptance(DBarReader reader, Variables * vars){
	    if(!Refill) return;
	    for(int i=0;i<reader.GetTreeEntries();i++){
		    if(i%(int)20!=0) continue; // WTF ?!
		    UpdateProgressBar(i, reader.GetTreeEntries());
		    vars->ResetVariables();
		    reader.FillVariables(i,vars); 
		    vars->Update();
		    for(int nobj=0;nobj<Objects2beFilled.size();nobj++){ 
			    Objects2beFilled[nobj]->FillTotalMCEvents(vars,FillinVariables[nobj],DiscrimVariables[nobj]);
		    }
	    }
    }




};

#endif
