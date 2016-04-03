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

	AcceptanceD-> Eval_Gen_Acceptance(6);
	AcceptanceD-> Eval_MC_Acceptance();
	AcceptanceD-> Eval_Geomag_Acceptance(6);
	
	cout<<"*** Updating P1 file ****"<<endl;
    	nomefile=percorso + "/Risultati/risultati/"+mese+"_"+frac+"_P1.root";
        file1 = TFile::Open(nomefile.c_str(),"UPDATE");
        if(!file1){
                nomefile=percorso + "/Risultati/"+mese+"/"+mese+"_"+frac+"_P1.root";
                file1 =TFile::Open(nomefile.c_str(),"UPDATE");
        }
	file1->cd("Results");
	
	AcceptanceP ->Gen_Acceptance_R  ->Write("Gen_AcceptanceP_R"  ); 
	AcceptanceP ->Gen_Acceptance_TOF->Write("Gen_AcceptanceP_TOF");
	AcceptanceP ->Gen_Acceptance_NaF->Write("Gen_AcceptanceP_NaF"); 
	AcceptanceP ->Gen_Acceptance_Agl->Write("Gen_AcceptanceP_Agl");
	
	AcceptanceP ->MCAcceptance_R  ->Write("MC_AcceptanceP_R"  );       
	AcceptanceP ->MCAcceptance_TOF->Write("MC_AcceptanceP_TOF");
	AcceptanceP ->MCAcceptance_NaF->Write("MC_AcceptanceP_NaF");
	AcceptanceP ->MCAcceptance_Agl->Write("MC_AcceptanceP_Agl");

	AcceptanceP ->Geomag_Acceptance_R  ->Write("Geomag_AcceptanceP_R"  );
	AcceptanceP ->Geomag_Acceptance_TOF->Write("Geomag_AcceptanceP_TOF");
	AcceptanceP ->Geomag_Acceptance_NaF->Write("Geomag_AcceptanceP_NaF");
	AcceptanceP ->Geomag_Acceptance_Agl->Write("Geomag_AcceptanceP_Agl");


	AcceptanceD ->Gen_Acceptance_R  ->Write("Gen_AcceptanceD_R"  ); 
	AcceptanceD ->Gen_Acceptance_TOF->Write("Gen_AcceptanceD_TOF");
	AcceptanceD ->Gen_Acceptance_NaF->Write("Gen_AcceptanceD_NaF");      	
	AcceptanceD ->Gen_Acceptance_Agl->Write("Gen_AcceptanceD_Agl");

	AcceptanceD ->MCAcceptance_R  ->Write("MC_AcceptanceD_R"  );       
	AcceptanceD ->MCAcceptance_TOF->Write("MC_AcceptanceD_TOF");
	AcceptanceD ->MCAcceptance_NaF->Write("MC_AcceptanceD_NaF");
	AcceptanceD ->MCAcceptance_Agl->Write("MC_AcceptanceD_Agl");

	AcceptanceD ->Geomag_Acceptance_R  ->Write("Geomag_AcceptanceD_R"  );
	AcceptanceD ->Geomag_Acceptance_TOF->Write("Geomag_AcceptanceD_TOF");
	AcceptanceD ->Geomag_Acceptance_NaF->Write("Geomag_AcceptanceD_NaF");
	AcceptanceD ->Geomag_Acceptance_Agl->Write("Geomag_AcceptanceD_Agl");	

	file1->Write();
	file1->Close();
	
	return;
}
