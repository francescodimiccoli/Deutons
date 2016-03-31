using namespace std;


void Acceptance(TFile * file1){
	ACCEPTANCE * AcceptanceP = new ACCEPTANCE (file1, bin, BetabinsR_D, BetabinsNaFR_D,BetabinsAglR_D,0.0308232619,"Results","EffpreselMCP","TOTLATCorr","CorrezioneLATp",1);
	ACCEPTANCE * AcceptanceD = new ACCEPTANCE (file1, bin, BetabinsR_D, BetabinsNaFR_D,BetabinsAglR_D,0.0308232619,"Results","EffpreselMCD","TOTLATCorr","CorrezioneLATd",6);
	
        cout<<"****************** ACCEPTANCE CALCULATION ******************"<<endl;

	cout<<  AcceptanceP -> after_TOF <<" \n";
	cout<<  AcceptanceP -> after_NaF <<" \n";
	cout<<  AcceptanceP -> after_Agl <<" \n";
	cout<<  AcceptanceP -> after_R <<" \n";

	AcceptanceP-> Eval_Gen_Acceptance(1);

	cout<<"*** Updating P1 file ****"<<endl;
        string nomefile=percorso + "/Risultati/risultati/"+mese+"_"+frac+"_P1.root";
        file1 =TFile::Open(nomefile.c_str(),"UPDATE");
        if(!file1){
                nomefile=percorso + "/Risultati/"+mese+"/"+mese+"_"+frac+"_P1.root";
                file1 =TFile::Open(nomefile.c_str(),"UPDATE");
        }

        file1->cd("Results");

	AcceptanceP ->Gen_Acceptance_R  ->Write("Gen_Acceptance_R  "); 
	AcceptanceP ->Gen_Acceptance_TOF->Write("Gen_Acceptance_TOF");
	AcceptanceP ->Gen_Acceptance_NaF->Write("Gen_Acceptance_NaF"); 
	AcceptanceP ->Gen_Acceptance_Agl->Write("Gen_Acceptance_Agl");
	
        file1->Write();
        file1->Close();


}
