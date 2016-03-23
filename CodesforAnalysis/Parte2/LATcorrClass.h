using namespace std;

class LATcorr
{
public:
    //counts
    TH1 * beforeR,   * afterR;
    TH1 * beforeTOF, * afterTOF;
    TH1 * beforeNaF, * afterNaF;
    TH1 * beforeAgl, * afterAgl;

    //eff	 	
    TH1 * effR;   
    TH1 * effTOF;
    TH1 * effNaF;
    TH1 * effAgl;
    
    //name
    std::string name;				 
    
    //  Creation constructors:

    LATcorr(std::string basename, int n){
        beforeTOF = new TH2F((basename + "1"   ).c_str(),(basename + "1"   ).c_str(),43,0,43, n, 0 ,n);
        afterTOF  = new TH2F((basename + "2"   ).c_str(),(basename + "2"   ).c_str(),43,0,43, n, 0 ,n);
        beforeNaF = new TH2F((basename + "1NaF").c_str(),(basename + "1NaF").c_str(),43,0,43, n, 0 ,n);
        afterNaF  = new TH2F((basename + "2NaF").c_str(),(basename + "2NaF").c_str(),43,0,43, n, 0 ,n);
        beforeAgl = new TH2F((basename + "1Agl").c_str(),(basename + "1Agl").c_str(),43,0,43, n, 0 ,n);
        afterAgl  = new TH2F((basename + "2Agl").c_str(),(basename + "2Agl").c_str(),43,0,43, n, 0 ,n);
        beforeR   = new TH2F((basename + "1_R" ).c_str(),(basename + "1_R" ).c_str(),43,0,43, n, 0 ,n);
        afterR    = new TH2F((basename + "2_R" ).c_str(),(basename + "2_R" ).c_str(),43,0,43, n, 0 ,n);
   	 name = basename; 
   }

  //   Reading constructors

    LATcorr(TFile * file, std::string basename){
        beforeTOF = (TH1 *)file->Get((basename + "1"   ).c_str());
        afterTOF  = (TH1 *)file->Get((basename + "2"   ).c_str());
        beforeNaF = (TH1 *)file->Get((basename + "1NaF").c_str());
        afterNaF  = (TH1 *)file->Get((basename + "2NaF").c_str());
        beforeAgl = (TH1 *)file->Get((basename + "1Agl").c_str());
        afterAgl  = (TH1 *)file->Get((basename + "2Agl").c_str());
        beforeR   = (TH1 *)file->Get((basename + "1_R" ).c_str());
        afterR    = (TH1 *)file->Get((basename + "2_R" ).c_str());
    	name = basename;   
    }

    void Write();
    void UpdateErrorbars();	
    void Eval_Efficiency();		
};


void LATcorr::Write()
{
	if(afterR)	afterR->Write(); 
	if(beforeR)	beforeR->Write();	   
	if(beforeTOF)	beforeTOF->Write(); 
	if(beforeNaF)	beforeNaF->Write();
	if(beforeAgl)	beforeAgl->Write();
	if(afterTOF)	afterTOF->Write();
	if(afterNaF)	afterNaF->Write();
	if(afterAgl)	afterAgl->Write();
}

void LATcorr::UpdateErrorbars()
{

   if(afterR)	 afterR->Sumw2();
   if(beforeR)	 beforeR->Sumw2(); 
   if(beforeTOF) beforeTOF->Sumw2(); 
   if(beforeNaF) beforeNaF->Sumw2();
   if(beforeAgl) beforeAgl->Sumw2();
   if(afterTOF)	 afterTOF->Sumw2();
   if(afterNaF)	 afterNaF->Sumw2();
   if(afterAgl)	 afterAgl->Sumw2();

}
 
void LATcorr::Eval_Efficiency(){
	
	if(afterR)	 effR     = (TH1 *)afterR	->Clone();
	if(afterTOF)	 effTOF   = (TH1 *)afterTOF	->Clone();	
	if(afterNaF)  	 effNaF   = (TH1 *)afterNaF	->Clone();
	if(afterAgl) 	 effAgl   = (TH1 *)afterAgl	->Clone();
	
	LATcorr::UpdateErrorbars();
	if(effR)   	effR   ->Divide( beforeR   );	
        if(effTOF) 	effTOF ->Divide( beforeTOF );
        if(effNaF) 	effNaF ->Divide( beforeNaF );
        if(effAgl) 	effAgl ->Divide( beforeAgl );

	if(effR)  	effR   ->SetName((name	+ "_EffR"  ).c_str()); 
        if(effTOF)	effTOF ->SetName((name	+ "_EffTOF").c_str());
        if(effNaF)	effNaF ->SetName((name	+ "_EffNaF").c_str());
        if(effAgl)	effAgl ->SetName((name	+ "_EffAgl").c_str());
}

