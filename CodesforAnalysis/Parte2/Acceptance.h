using namespace std;


void Acceptance(){

	string nomefile=percorso + "/Risultati/risultati/"+mese+"_"+frac+"_P1.root";
        TFile * file1 = TFile::Open(nomefile.c_str(),"READ");
        if(!file1){
                nomefile=percorso + "/Risultati/"+mese+"/"+mese+"_"+frac+"_P1.root";
                file1 =TFile::Open(nomefile.c_str(),"READ");
        }
	
	ACCEPTANCE * AcceptanceP = new ACCEPTANCE (file1, bin, BetabinsR_D, BetabinsNaFR_D,BetabinsAglR_D,0.0308232619,"Results","EffpreselMCP","EffpreselMCP","TOTLATCorr","CorrezioneLATp",1);
	ACCEPTANCE * AcceptanceD = new ACCEPTANCE (file1, bin, BetabinsR_D, BetabinsNaFR_D,BetabinsAglR_D,0.0308232619,"Results","EffpreselMCD","EffpreselMCD","TOTLATCorr","CorrezioneLATd",6);
	
        cout<<"****************** ACCEPTANCE CALCULATION ******************"<<endl;

	AcceptanceP-> Eval_Gen_Acceptance(1);
	AcceptanceP-> Eval_MC_Acceptance();
	AcceptanceP-> Eval_Geomag_Acceptance(1);
	
	cout<<"*** Updating P1 file ****"<<endl;
    	nomefile=percorso + "/Risultati/risultati/"+mese+"_"+frac+"_P1.root";
        file1 = TFile::Open(nomefile.c_str(),"UPDATE");
        if(!file1){
                nomefile=percorso + "/Risultati/"+mese+"/"+mese+"_"+frac+"_P1.root";
                file1 =TFile::Open(nomefile.c_str(),"UPDATE");
        }
	file1->cd("Results");
	
	AcceptanceP ->Gen_Acceptance_R  ->Write("Gen_Acceptance_R_P"  ); 
	AcceptanceP ->Gen_Acceptance_TOF->Write("Gen_Acceptance_TOF_P");
	AcceptanceP ->Gen_Acceptance_NaF->Write("Gen_Acceptance_NaF_P"); 
	AcceptanceP ->Gen_Acceptance_Agl->Write("Gen_Acceptance_Agl_P");
	
 	AcceptanceP ->MCAcceptance_R  ->Write("MC_Acceptance_R_P"  );       
 	AcceptanceP ->MCAcceptance_TOF->Write("MC_Acceptance_TOF_P");
 	AcceptanceP ->MCAcceptance_NaF->Write("MC_Acceptance_NaF_P");
 	AcceptanceP ->MCAcceptance_Agl->Write("MC_Acceptance_Agl_P");

	AcceptanceP ->Geomag_Acceptance_R  ->Write("Geomag_Acceptance_R_P"  );
	AcceptanceP ->Geomag_Acceptance_TOF->Write("Geomag_Acceptance_TOF_P");
	AcceptanceP ->Geomag_Acceptance_NaF->Write("Geomag_Acceptance_NaF_P");
	AcceptanceP ->Geomag_Acceptance_Agl->Write("Geomag_Acceptance_Agl_P");

	file1->Write();
        file1->Close();
	
	return;
}
