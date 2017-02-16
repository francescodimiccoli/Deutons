
void CalibrateEdep(Variables *vars){

        float Velocity = vars->Beta;
        if(vars->BetaRICH_new>0) Velocity=vars->BetaRICH_new;
        Velocity=fabs(Velocity);

        if(Velocity>0&&Velocity<1){
                vars->EdepTOFU =((vars->EdepTOFU)*Corr_TOFU->Eval(Velocity));
                vars->EdepTOFD =((vars->EdepTOFD)*Corr_TOFD->Eval(Velocity));
                vars->EdepTrack=((vars->EdepTrack)*Corr_Track->Eval(Velocity));
        }

        return;
}



void ProcessEvent(Variables *vars,bool isMC,Reweighter reweighter){
	vars->joinCutmask=vars->CUTMASK;
	vars->joinCutmask=vars->CUTMASK|(1<<10);
	vars->joinCutmask = vars->joinCutmask|(vars->RICHmask_new<<11);	
	if(isMC){
		vars->mcweight=reweighter.getWeight(fabs(vars->Momento_gen));
                        if(vars->Momento_gen<1) vars->mcweight=1;
			CalibrateEdep(vars);	
		
	}
	if(IsPreselected(vars)){
		Likelihood(vars);
		Eval_Distance(vars);
		Eval_Distance(vars,false);	
	}	
//	vars->PrintCurrentState();
	return;		

}



bool IsMC(TTree * tree){

        Variables * Vars = new Variables;
        Vars->ReadBranches(tree);
        tree->GetEvent(10);

        if (Vars->Massa_gen>0) return true;
        else return false;
}

bool ReadCalibration(string month){

        cout<<"****************** CALIB. READING **************************"<<endl;
        string nomecal=("/storage/gpfs_ams/ams/users/fdimicco/Deutons/CodesforAnalysis/CALIBRAZIONI/"+month+".root");
        TFile *_file2 = TFile::Open(nomecal.c_str());
        if(_file2) cout<<"calibration file found: "<<(month + ".root").c_str()<<endl;
        else { cout<<"calibration file not found"<<endl; return false;}

        Rig = (TSpline3 *) _file2->Get("Fit Results/Splines/Rig");
        beta = (TSpline3 *) _file2->Get("Fit Results/Splines/beta");
        betaNaF = (TF1 *) _file2->Get("Fit Results/Splines/SigmaInvBetaNaF_spl");
        betaAgl = (TF1 *) _file2->Get("Fit Results/Splines/SigmaInvBetaAgl_spl");
        eL1 = (TSpline3 *) _file2->Get("Fit Results/Splines/eL1");
        etofu =  (TSpline3 *) _file2->Get("Fit Results/Splines/etofu");
        etrack =  (TSpline3 *) _file2->Get("Fit Results/Splines/etrack");
        etofd =  (TSpline3 *) _file2->Get("Fit Results/Splines/etofd");
        EdepL1beta =  (TSpline3 *) _file2->Get("Fit Results/Splines/EdepL1beta");
        EdepTOFbeta =  (TSpline3 *) _file2->Get("Fit Results/Splines/EdepTOFbeta");
        EdepTrackbeta =  (TSpline3 *) _file2->Get("Fit Results/Splines/EdepTrackbeta");
        EdepTOFDbeta =  (TSpline3 *) _file2->Get("Fit Results/Splines/EdepTOFDbeta");
        Corr_L1 =  (TSpline3 *) _file2->Get("Fit Results/Splines/Corr_L1");
        Corr_TOFU =  (TSpline3 *) _file2->Get("Fit Results/Splines/Corr_TOFU");
        Corr_Track =  (TSpline3 *) _file2->Get("Fit Results/Splines/Corr_Track");
        Corr_TOFD =  (TSpline3 *) _file2->Get("Fit Results/Splines/Corr_TOFD");
        cout<<Rig<<" "<<beta<<" "<<" "<<betaNaF<<" "<<betaAgl<<" "<<eL1<<" "<<etofu<<" "<<etrack<<" "<<etofd<<" "<<EdepL1beta<<" "<<EdepTOFbeta<<" "<<EdepTrackbeta<<" "<<EdepTOFDbeta<<" "<<Corr_L1<<" "<<Corr_TOFU<<" "<<Corr_Track<<" "<<Corr_TOFD<<endl;

        return true;

}


bool ReadPdfForLikelihood(){
        cout<<"*************** PDF for LIKELIHOOD READING ********************"<<endl;
        string nome;
        TFile *_file1 = TFile::Open("/storage/gpfs_ams/ams/users/fdimicco/Deutons/Ntuple-making/QualityVariables.root");
        TFile *_file3 = TFile::Open("/storage/gpfs_ams/ams/users/fdimicco/Deutons/Ntuple-making/QualityVariables_NaF.root");
        TFile *_file3b = TFile::Open("/storage/gpfs_ams/ams/users/fdimicco/Deutons/Ntuple-making/QualityVariables_Agl.root");
        if(!_file1||!_file3||!_file3b) {cout<<"PDF for likelihood files not found"<<endl; return false;}
        else cout<<"PDF for likelihood files found"<<endl;

        string Variables[9]= {"N. Anti-clusters","Unused TOF Clusters","|Rup-Rdown|:R","Unused Tracker layers","Tracker: Y Hits without X","Track Chi^2","|E.dep(lower TOF) - E.dep(upper TOF)|","|E.dep(layer 2)-E.dep(layer 1)|","|E. dep.(tot)-E.dep.(track)|"};
        for(int u2=0; u2<9; u2++) {
                nome="Splines/Spline: "+Variables[u2]+"_SGNL";
                Signal[u2]=(TSpline3 *) _file1->Get(nome.c_str());
                nome="Splines/Spline: "+Variables[u2]+"_BKGND";
                Bkgnd[u2]=(TSpline3 *) _file1->Get(nome.c_str());
        }
        cout<<_file1<<endl;
        string VariablesRICH[9]= {"N. Anti-clusters","Unused TOF Clusters","|Rup-Rdown|:R","Unused Tracker layers","Tracker: Y Hits without X","Track Chi^2","RICH Hits: tot - used","Rich Photoelectrons","|E. dep.(tot)-E.dep.(track)|"};
        for(int u2=0; u2<9; u2++) {
                nome="Splines/Spline NaF: "+VariablesRICH[u2]+"_SGNL";
                SignalNaF[u2]=(TSpline3 *) _file3->Get(nome.c_str());
                nome="Splines/Spline NaF: "+VariablesRICH[u2]+"_BKGND";
                BkgndNaF[u2]=(TSpline3 *) _file3->Get(nome.c_str());
        }
        cout<<_file3<<endl;
        for(int u2=0; u2<9; u2++) {
                nome="Splines/Spline Agl: "+VariablesRICH[u2]+"_SGNL";
                SignalAgl[u2]=(TSpline3 *) _file3b->Get(nome.c_str());
                nome="Splines/Spline Agl: "+VariablesRICH[u2]+"_BKGND";
                BkgndAgl[u2]=(TSpline3 *) _file3b->Get(nome.c_str());
        }
        cout<<_file3b<<endl;

        for(int qs=0; qs<9; qs++) cout<<Signal[qs]<<" ";
        cout<<endl;
        for(int qs=0; qs<9; qs++) cout<<SignalNaF[qs]<<" ";
        cout<<endl;
        for(int qs=0; qs<9; qs++) cout<<SignalAgl[qs]<<" ";
        cout<<endl;

        return true;
}


Reweighter ReweightInitializer(){
        Histogram   mcFlux = makeLogUniform(500, 0.5, 100);
        Histogram dataFlux = loadGalpropFile("/storage/gpfs_ams/ams/users/fdimicco/Deutons/include/CRDB_ProtonsAMS_R.galprop");
        dataFlux.multiply( mcFlux.at(1.05) / dataFlux.getContent()[0] );
        Reweighter reweighter(mcFlux, dataFlux);
        return reweighter;
}

int AssignMC_type(float Massa_gen)
{
        int MC_type=0;
        int cursor=0;
        if(Massa_gen<1) {
                if(Massa_gen<0.93805) MC_type = MC_type|(1<<(cursor+0));
                if(Massa_gen<0.939  ) MC_type = MC_type|(1<<(cursor+1));
                MC_type = MC_type|(1<<(cursor+2));
        }
        cursor=8;
        if(Massa_gen>1&&Massa_gen<2) {
                MC_type = MC_type|(1<<(cursor+0));
                int moffset=18570;
                int c_s_type=(int)(10000*Massa_gen-moffset + 2);
                MC_type = MC_type|(1<<(cursor+c_s_type));
        }
        cursor=16;
        if(Massa_gen>3) {
                MC_type = MC_type|(1<<(cursor+0));
                MC_type = MC_type|(1<<(cursor+2));
        }
        return  MC_type;
}

