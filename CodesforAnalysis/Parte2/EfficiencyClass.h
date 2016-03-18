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
        beforeTOF = (TH1F *)file->Get((basename + "1"   ).c_str());
        afterTOF  = (TH1F *)file->Get((basename + "2"   ).c_str());
        beforeNaF = (TH1F *)file->Get((basename + "1NaF").c_str());
        afterNaF  = (TH1F *)file->Get((basename + "2NaF").c_str());
        beforeAgl = (TH1F *)file->Get((basename + "1Agl").c_str());
        afterAgl  = (TH1F *)file->Get((basename + "2Agl").c_str());
        beforeR   = (TH1F *)file->Get((basename + "1_R" ).c_str());
        afterR    = (TH1F *)file->Get((basename + "2_R" ).c_str());
    }

    Efficiency(TFile * file, std::string basename, int i){
        beforeTOF = (TH2F *)file->Get((basename + "1"   ).c_str());
        afterTOF  = (TH2F *)file->Get((basename + "2"   ).c_str());
        beforeNaF = (TH2F *)file->Get((basename + "1NaF").c_str());
        afterNaF  = (TH2F *)file->Get((basename + "2NaF").c_str());
        beforeAgl = (TH2F *)file->Get((basename + "1Agl").c_str());
        afterAgl  = (TH2F *)file->Get((basename + "2Agl").c_str());
        beforeR   = (TH2F *)file->Get((basename + "1_R" ).c_str());
        afterR    = (TH2F *)file->Get((basename + "2_R" ).c_str());
    }

    Efficiency(TFile * file, std::string basename, int i, int j){
        beforeTOF = (TH3F *)file->Get((basename + "1"   ).c_str());
        afterTOF  = (TH3F *)file->Get((basename + "2"   ).c_str());
        beforeNaF = (TH3F *)file->Get((basename + "1NaF").c_str());
        afterNaF  = (TH3F *)file->Get((basename + "2NaF").c_str());
        beforeAgl = (TH3F *)file->Get((basename + "1Agl").c_str());
        afterAgl  = (TH3F *)file->Get((basename + "2Agl").c_str());
        beforeR   = (TH3F *)file->Get((basename + "1_R" ).c_str());
        afterR    = (TH3F *)file->Get((basename + "2_R" ).c_str());
    }

    void Write();

};


void Efficiency::Write()
{
    beforeTOF->Write(); 
    beforeNaF->Write();
    beforeAgl->Write();
    afterTOF->Write();
    afterNaF->Write();
    afterAgl->Write();
}  
