class EffCorr{

	private:

	Efficiency * EffMC;
	Efficiency * EffData;

	public:

	EffCorr(FileSaver  File, std::string Basename,std::string Directory, Binning Bins, std::string Cut_before,std::string Cut_after,std::string Cut_Data,std::string Cut_MC){
		EffMC   = new Efficiency(File, (Basename+"_MC" ).c_str(),Directory,Bins, (Cut_before+"&"+Cut_MC  ).c_str(),(Cut_after+"&"+Cut_MC  ).c_str());
		EffData = new Efficiency(File, (Basename+"_lat").c_str(),Directory,Bins, (Cut_before+"&"+Cut_Data).c_str(),(Cut_after+"&"+Cut_Data).c_str(),LatEdges);
		
		cout<<(Cut_before+"&"+Cut_MC  ).c_str()<<endl;
	}	

	void Fill(TNtuple * treeMC,TNtuple * treeDT, Variables * vars, float (*discr_var) (Variables * vars),bool refill=false);
	void Save(FileSaver finalhistos);
	void Eval_Efficiencies();
	void SaveResults(FileSaver finalhistos);
};




void EffCorr::Fill(TNtuple * treeMC,TNtuple * treeDT, Variables * vars, float (*discr_var) (Variables * vars),bool refill){

	EffMC   -> Fill(treeMC,vars,discr_var,refill);
	EffData -> Fill(treeDT,vars,discr_var,refill);
}

void EffCorr::Save(FileSaver finalhistos){
	EffMC  -> Save(finalhistos);
	EffData-> Save(finalhistos);
}

void EffCorr::Eval_Efficiencies(){
	EffMC  -> Eval_Efficiency();
	EffData-> Eval_Efficiency();
}

void EffCorr::SaveResults(FileSaver finalhistos){
	EffMC  -> SaveResults(finalhistos);
	EffData-> SaveResults(finalhistos);
}
