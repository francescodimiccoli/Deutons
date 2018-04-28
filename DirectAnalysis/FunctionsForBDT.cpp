#include <TMVA/Reader.h>
#include <TMVA/Tools.h>

//new members of Analysis class:
TMVA::Reader readerNaF;
TMVA::Reader readerAgl;


void Analysis::InitBDTNaF(){
	TMVA::Tools::Instance();
	readerNaF = new TMVA::Reader( "V:Color:!Silent" );
	readerNaF->AddVariable("Chisquare",&ntpTracker.chisqn[1][0]);
    		//readerNaF->AddVariable("Layernonusati",&);
    		//readerNaF->AddVariable("NTofUsed",&);
    		//readerNaF->AddVariable("diffR",&);
    		//readerNaF->AddVariable("TOF_Up_Down",&);
    		//readerNaF->AddVariable("Richtotused",&);	
    		//readerNaF->AddVariable("RichPhEl",&RichPhEl);
    	readerNaF->AddVariable("RICHprob",&ntpRich.prob);
    		//readerNaF->AddVariable("RICHcollovertotal",&);
 	    	//readerNaF->AddVariable("RICHLipBetaConsistency",&);  
    		//readerNaF->AddVariable("RICHTOFBetaConsistency",&);  
    	readerNaF->AddVariable("RICHChargeConsistency",&ntpRich.q_consistency);
	readerNaF->AddVariable("RICHPmts := RICHPmts/1.40129846432481707e-45",&ntpRich.npmt);
   	readerNaF->AddVariable("RICHgetExpected := RICHgetExpected/1.40129846432481707e-45",&ntpRich.np_exp_uncorr);
    	readerNaF->AddVariable("tot_hyp_p_uncorr",&ntpRich.tot_hyp_p_uncorr[1]);
    		//readerNaF->AddVariable("Bad_ClusteringRICH",&);
    		//readerNaF->AddVariable("NSecondariesRICHrich",&);
	readerNaF->BookMVA("BDTmethod", "QualityNaF_BDT.weights.xml");

}

void Analysis::InitBDTAgl(){
	TMVA::Tools::Instance();
	readerAgl = new TMVA::Reader( "V:Color:!Silent" );
	readerAgl = new TMVA::Reader( "V:Color:!Silent" );
		//readerAgl->AddVariable("Chisquare",&ntpTracker.chisqn[1][0]);
    		//readerAgl->AddVariable("Layernonusati",&);
    		//readerAgl->AddVariable("NTofUsed",&);
    		//readerAgl->AddVariable("diffR",&);
    		//readerAgl->AddVariable("TOF_Up_Down",&);
    		//readerAgl->AddVariable("Richtotused",&);	
    		//readerAgl->AddVariable("RichPhEl",&RichPhEl);
    	readerAgl->AddVariable("RICHprob",&ntpRich.prob);
    		//readerAgl->AddVariable("RICHcollovertotal",&);
 	    	//readerAgl->AddVariable("RICHLipBetaConsistency",&);  
    		//readerAgl->AddVariable("RICHTOFBetaConsistency",&);  
    	readerAgl->AddVariable("RICHChargeConsistency",&ntpRich.q_consistency);
	readerAgl->AddVariable("RICHPmts := RICHPmts/1.40129846432481707e-45",&ntpRich.npmt);
   	readerAgl->AddVariable("RICHgetExpected := RICHgetExpected/1.40129846432481707e-45",&ntpRich.np_exp_uncorr);
    	readerAgl->AddVariable("tot_hyp_p_uncorr",&ntpRich.tot_hyp_p_uncorr[1]);
    		//readerAgl->AddVariable("Bad_ClusteringRICH",&);
    		//readerAgl->AddVariable("NSecondariesRICHrich",&);
	readerAgl->BookMVA("BDTmethod", "QualityAgl_BDT.weights.xml");
}



//to be put in Analysis::FillNtpHeader()
NtpHeader.BDTNaF = readerNaF->EvaluateMVA("BDTmethod");
NtpHeader.BDTAgl = readerAgl->EvaluateMVA("BDTmethod");

