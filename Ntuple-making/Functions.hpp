
void CalibrateEdep(Variables *vars){

        float Velocity = vars->Beta;
        if(vars->BetaRICH_new>0) Velocity=vars->BetaRICH_new;
        Velocity=fabs(Velocity);

        if(Velocity>0&&Velocity<1){
                vars->EdepTOFU =((vars->EdepTOFU)*(1/Corr_TOFU->Eval(Velocity)));
                vars->EdepTOFD =((vars->EdepTOFD)*(1/Corr_TOFD->Eval(Velocity)));
                vars->EdepTrack=((vars->EdepTrack)*(1/Corr_Track->Eval(Velocity)));
        }

        return;
}



void ProcessEvent(Variables *vars,bool isMC,Reweighter reweighter){
	
	if(isMC){
		vars->mcweight=reweighter.getWeight(fabs(vars->Momento_gen));
                        if(vars->Momento_gen<1) vars->mcweight=1;
			//CalibrateEdep(vars);	
		
	}
	if(ApplyCuts("IsPreselected",vars)){
		Likelihood(vars);
		Eval_DistanceFromP(vars);
		Eval_DistanceFromD(vars);	
	}	
	//vars->PrintCurrentState();
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
        string nomecal=("/storage/gpfs_ams/ams/users/fdimicco/Deutons/Calibrations/Calibfiles/"+month+"/Calib.root");
        FileSaver Calibration;	
	Calibration.setName(nomecal.c_str());
	bool checkfile = Calibration.CheckFile();

	if(checkfile) cout<<"calibration file found: "<<(month).c_str()<<endl;
        else { cout<<"calibration file not found"<<endl; return false;}

	Resolution * RigidityResolution_P = new Resolution(Calibration,"RvsR Resolution (P)",PResB);
        Resolution * RigidityResolution_D = new Resolution(Calibration,"RvsR Resolution (D)",PResB);

	Resolution * BetaTOFResolution_P = new Resolution(Calibration,"BetaTOFvsBeta Resolution (P)",ToFResB);
        Resolution * BetaTOFResolution_D = new Resolution(Calibration,"BetaTOFvsBeta Resolution (D)",ToFResB);
	
	Resolution * BetaNaFResolution_P = new Resolution(Calibration,"BetaNaFvsBeta Resolution (P)",NaFResB);
        Resolution * BetaNaFResolution_D = new Resolution(Calibration,"BetaNaFvsBeta Resolution (D)",NaFResB);
	
	Resolution * BetaAglResolution_P = new Resolution(Calibration,"BetaAglvsBeta Resolution (P)",AglResB);
        Resolution * BetaAglResolution_D = new Resolution(Calibration,"BetaAglvsBeta Resolution (D)",AglResB);

	Resolution * EdepUTOFResolution_P= new Resolution(Calibration,"EdepUTOFvsBeta Resolution (P)",ToFResB);
        Resolution * EdepUTOFResolution_D= new Resolution(Calibration,"EdepUTOFvsBeta Resolution (D)",ToFResB);

	Resolution * EdepLTOFResolution_P= new Resolution(Calibration,"EdepLTOFvsBeta Resolution (P)",ToFResB);
	Resolution * EdepLTOFResolution_D= new Resolution(Calibration,"EdepLTOFvsBeta Resolution (D)",ToFResB);

	Resolution * EdepTrackResolution_P= new Resolution(Calibration,"EdepTrackvsBeta Resolution (P)",ToFResB);
	Resolution * EdepTrackResolution_D= new Resolution(Calibration,"EdepTrackvsBeta Resolution (D)",ToFResB);

	Resolution * EdepUTOFMC_P = new Resolution(Calibration,"EdepUTOFvsBeta Measured MC",ToFResB);
	Resolution * EdepUTOFDT_P = new Resolution(Calibration,"EdepUTOFvsBeta Measured DT",ToFResB);

	Resolution * EdepLTOFMC_P = new Resolution(Calibration,"EdepLTOFvsBeta Measured MC",ToFResB);
	Resolution * EdepLTOFDT_P = new Resolution(Calibration,"EdepLTOFvsBeta Measured DT",ToFResB);

	Resolution * EdepTrackMC_P = new Resolution(Calibration,"EdepTrackvsBeta Measured MC",ToFResB);
	Resolution * EdepTrackDT_P = new Resolution(Calibration,"EdepTrackvsBeta Measured DT",ToFResB);


	Rig_p		=(TSpline3 *) RigidityResolution_P  ->Get_SigmasModel(); 
	beta_p		=(TF1 *   )   BetaTOFResolution_P   ->Get_SigmasModel(); 
	betaNaF_p	=(TF1 *	  )   BetaNaFResolution_P   ->Get_SigmasModel(); 
	betaAgl_p	=(TF1 *	  )   BetaAglResolution_P   ->Get_SigmasModel(); 
	Rigmean_p	=(TSpline3 *) RigidityResolution_P  ->Get_MeansModel(); 
	betamean_p	=(TF1 *   )   BetaTOFResolution_P   ->Get_MeansModel(); 
	betaNaFmean_p	=(TF1 *	  )   BetaNaFResolution_P   ->Get_MeansModel(); 
	betaAglmean_p	=(TF1 *	  )   BetaAglResolution_P   ->Get_MeansModel(); 
	etofu_p		=(TSpline3 *) EdepUTOFResolution_P  ->Get_SigmasModel(); 
	etrack_p	=(TSpline3 *) EdepTrackResolution_P ->Get_SigmasModel(); 
	etofd_p		=(TSpline3 *) EdepLTOFResolution_P  ->Get_SigmasModel(); 
	EdepTOFbeta_p	=(TSpline3 *) EdepUTOFResolution_P  ->Get_MeansModel(); 
	EdepTrackbeta_p	=(TSpline3 *) EdepTrackResolution_P ->Get_MeansModel();  
	EdepTOFDbeta_p	=(TSpline3 *) EdepLTOFResolution_P  ->Get_MeansModel(); 

	Rig_d		=(TSpline3 *) RigidityResolution_D  ->Get_SigmasModel(); 
        beta_d		=(TF1 *   )   BetaTOFResolution_D   ->Get_SigmasModel(); 
        betaNaF_d	=(TF1 *	  )   BetaNaFResolution_D   ->Get_SigmasModel(); 
        betaAgl_d	=(TF1 *	  )   BetaAglResolution_D   ->Get_SigmasModel(); 
        Rigmean_d	=(TSpline3 *) RigidityResolution_D  ->Get_MeansModel(); 
        betamean_d	=(TF1 *   )   BetaTOFResolution_D   ->Get_MeansModel(); 
        betaNaFmean_d	=(TF1 *	  )   BetaNaFResolution_D   ->Get_MeansModel(); 
        betaAglmean_d	=(TF1 *	  )   BetaAglResolution_D   ->Get_MeansModel(); 
        etofu_d		=(TSpline3 *) EdepUTOFResolution_D  ->Get_SigmasModel(); 
        etrack_d	=(TSpline3 *) EdepTrackResolution_D ->Get_SigmasModel(); 
        etofd_d		=(TSpline3 *) EdepLTOFResolution_D  ->Get_SigmasModel(); 
        EdepTOFbeta_d	=(TSpline3 *) EdepUTOFResolution_D  ->Get_MeansModel(); 
        EdepTrackbeta_d	=(TSpline3 *) EdepTrackResolution_D ->Get_MeansModel();  
        EdepTOFDbeta_d	=(TSpline3 *) EdepLTOFResolution_D  ->Get_MeansModel(); 

        return checkfile;

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

