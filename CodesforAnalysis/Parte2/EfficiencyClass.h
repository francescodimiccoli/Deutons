using namespace std;

class Efficiency
{
public:
    TH1 * beforeR,   * afterR;
    TH1 * beforeTOF, * afterTOF;
    TH1 * beforeNaF, * afterNaF;
    TH1 * beforeAgl, * afterAgl;

    //  Creation constructors:
   
     Efficiency(std::string basename){
        beforeTOF = new TH1F((basename + "1"   ).c_str(),(basename + "1"   ).c_str(),18,0,18);
        afterTOF  = new TH1F((basename + "2"   ).c_str(),(basename + "2"   ).c_str(),18,0,18);
        beforeNaF = new TH1F((basename + "1NaF").c_str(),(basename + "1NaF").c_str(),18,0,18);
        afterNaF  = new TH1F((basename + "2NaF").c_str(),(basename + "2NaF").c_str(),18,0,18);
        beforeAgl = new TH1F((basename + "1Agl").c_str(),(basename + "1Agl").c_str(),18,0,18);
        afterAgl  = new TH1F((basename + "2Agl").c_str(),(basename + "2Agl").c_str(),18,0,18);
        beforeR   = new TH1F((basename + "1_R" ).c_str(),(basename + "1_R" ).c_str(),43,0,43);
        afterR    = new TH1F((basename + "2_R" ).c_str(),(basename + "2_R" ).c_str(),43,0,43);
    }

    Efficiency(std::string basename, int n){
        beforeTOF = new TH2F((basename + "1"   ).c_str(),(basename + "1"   ).c_str(),18,0,18, n, 0 ,n);
        afterTOF  = new TH2F((basename + "2"   ).c_str(),(basename + "2"   ).c_str(),18,0,18, n, 0 ,n);
        beforeNaF = new TH2F((basename + "1NaF").c_str(),(basename + "1NaF").c_str(),18,0,18, n, 0 ,n);
        afterNaF  = new TH2F((basename + "2NaF").c_str(),(basename + "2NaF").c_str(),18,0,18, n, 0 ,n);
        beforeAgl = new TH2F((basename + "1Agl").c_str(),(basename + "1Agl").c_str(),18,0,18, n, 0 ,n);
        afterAgl  = new TH2F((basename + "2Agl").c_str(),(basename + "2Agl").c_str(),18,0,18, n, 0 ,n);
        beforeR   = new TH2F((basename + "1_R" ).c_str(),(basename + "1_R" ).c_str(),43,0,43, n, 0 ,n);
        afterR    = new TH2F((basename + "2_R" ).c_str(),(basename + "2_R" ).c_str(),43,0,43, n, 0 ,n);
    }

    Efficiency(std::string basename, int n, int m){
        beforeTOF = new TH3F((basename + "1"   ).c_str(),(basename + "1"   ).c_str(),18,0,18, n, 0 ,n, m, 0, m);
        afterTOF  = new TH3F((basename + "2"   ).c_str(),(basename + "2"   ).c_str(),18,0,18, n, 0 ,n, m, 0, m);
        beforeNaF = new TH3F((basename + "1NaF").c_str(),(basename + "1NaF").c_str(),18,0,18, n, 0 ,n, m, 0, m);
        afterNaF  = new TH3F((basename + "2NaF").c_str(),(basename + "2NaF").c_str(),18,0,18, n, 0 ,n, m, 0, m);
        beforeAgl = new TH3F((basename + "1Agl").c_str(),(basename + "1Agl").c_str(),18,0,18, n, 0 ,n, m, 0, m);
        afterAgl  = new TH3F((basename + "2Agl").c_str(),(basename + "2Agl").c_str(),18,0,18, n, 0 ,n, m, 0, m);
        beforeR   = new TH3F((basename + "1_R" ).c_str(),(basename + "1_R" ).c_str(),43,0,43, n, 0 ,n, m, 0, m);
        afterR    = new TH3F((basename + "2_R" ).c_str(),(basename + "2_R" ).c_str(),43,0,43, n, 0 ,n, m, 0, m);
    }

  //   Reading constructors

    Efficiency(TFile * file, std::string basename){
        beforeTOF = (TH1 *)file->Get((basename + "1"   ).c_str());
        afterTOF  = (TH1 *)file->Get((basename + "2"   ).c_str());
        beforeNaF = (TH1 *)file->Get((basename + "1NaF").c_str());
        afterNaF  = (TH1 *)file->Get((basename + "2NaF").c_str());
        beforeAgl = (TH1 *)file->Get((basename + "1Agl").c_str());
        afterAgl  = (TH1 *)file->Get((basename + "2Agl").c_str());
        beforeR   = (TH1 *)file->Get((basename + "1_R" ).c_str());
        afterR    = (TH1 *)file->Get((basename + "2_R" ).c_str());
    }

    void Write();
    void UpdateErrorbars();	
};


void Efficiency::Write()
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

void Efficiency::UpdateErrorbars()
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
  
