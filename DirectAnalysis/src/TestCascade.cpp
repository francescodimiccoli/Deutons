#include "Analyzer.h"
#include "AllRangesEfficiency.h"
#include "TGraphErrors.h"
#include "filesaver.h"


void Analyzer::BookTestCascade(FileSaver finalhistos, FileSaver finalresults, bool refill)
{

	cout<<"****************************** BINS ***************************************"<<endl;
    	SetUpTOIBinning();

	bool checkfile = finalhistos.CheckFile();
	check_file = checkfile;
	cout<<"****************************** Efficiecny ANALYIS ******************************************"<<endl;



	Efficiency * Cascade1 = new Efficiency(finalhistos,"Cascade1","Cascade1",PRB,"IsProtonMC","IsProtonMC&IsDownGoing&IsLUT2");
	Efficiency * Cascade2 = new Efficiency(finalhistos,"Cascade2","Cascade2",PRB,"IsProtonMC","IsProtonMC&IsDownGoing&IsLUT2&IsGoodTrack");
	Efficiency * Cascade3 = new Efficiency(finalhistos,"Cascade3","Cascade3",PRB,"IsProtonMC","IsProtonMC&IsDownGoing&IsLUT2&IsGoodTrack&IsPhysTrig");
	Efficiency * Cascade4 = new Efficiency(finalhistos,"Cascade4","Cascade4",PRB,"IsProtonMC","IsProtonMC&IsDownGoing&IsLUT2&IsGoodTrack&IsPhysTrig&IsGoodChi2");
	Efficiency * Cascade5 = new Efficiency(finalhistos,"Cascade5","Cascade5",PRB,"IsProtonMC","IsProtonMC&IsDownGoing&IsLUT2&IsGoodTrack&IsPhysTrig&IsGoodChi2&IsCharge1Track");
	Efficiency * Cascade6 = new Efficiency(finalhistos,"Cascade6","Cascade6",PRB,"IsProtonMC","IsProtonMC&IsDownGoing&IsLUT2&IsGoodTrack&IsPhysTrig&IsGoodChi2&IsCharge1Track&IsCharge1UTOF&IsCharge1LTOF");
	Efficiency * Cascade7 = new Efficiency(finalhistos,"Cascade7","Cascade7",PRB,"IsProtonMC","IsProtonMC&IsDownGoing&IsLUT2&IsGoodTrack&IsPhysTrig&IsGoodChi2&IsCharge1Track&IsCharge1UTOF&IsCharge1LTOF&Is1TrTrack");
	Efficiency * Cascade8 = new Efficiency(finalhistos,"Cascade8","Cascade8",PRB,"IsProtonMC","IsProtonMC&IsDownGoing&IsLUT2&IsGoodTrack&IsPhysTrig&IsGoodChi2&IsCharge1Track&IsCharge1UTOF&IsCharge1LTOF&Is1TrTrack&IsMinTOF");

	Efficiency * Cascade1_D = new Efficiency(finalhistos,"Cascade1_D","Cascade1_D",PRB,"IsData","IsData&IsDownGoing&IsLUT2");
	Efficiency * Cascade2_D = new Efficiency(finalhistos,"Cascade2_D","Cascade2_D",PRB,"IsData","IsData&IsDownGoing&IsLUT2&IsGoodTrack");
	Efficiency * Cascade3_D = new Efficiency(finalhistos,"Cascade3_D","Cascade3_D",PRB,"IsData","IsData&IsDownGoing&IsLUT2&IsGoodTrack&IsPhysTrig");
	Efficiency * Cascade4_D = new Efficiency(finalhistos,"Cascade4_D","Cascade4_D",PRB,"IsData","IsData&IsDownGoing&IsLUT2&IsGoodTrack&IsPhysTrig&IsGoodChi2");
	Efficiency * Cascade5_D = new Efficiency(finalhistos,"Cascade5_D","Cascade5_D",PRB,"IsData","IsData&IsDownGoing&IsLUT2&IsGoodTrack&IsPhysTrig&IsGoodChi2&IsCharge1Track");
	Efficiency * Cascade6_D = new Efficiency(finalhistos,"Cascade6_D","Cascade6_D",PRB,"IsData","IsData&IsDownGoing&IsLUT2&IsGoodTrack&IsPhysTrig&IsGoodChi2&IsCharge1Track&IsCharge1UTOF&IsCharge1LTOF");
	Efficiency * Cascade7_D = new Efficiency(finalhistos,"Cascade7_D","Cascade7_D",PRB,"IsData","IsData&IsDownGoing&IsLUT2&IsGoodTrack&IsPhysTrig&IsGoodChi2&IsCharge1Track&IsCharge1UTOF&IsCharge1LTOF&Is1TrTrack");
	Efficiency * Cascade8_D = new Efficiency(finalhistos,"Cascade8_D","Cascade8_D",PRB,"IsData","IsData&IsDownGoing&IsLUT2&IsGoodTrack&IsPhysTrig&IsGoodChi2&IsCharge1Track&IsCharge1UTOF&IsCharge1LTOF&Is1TrTrack&IsMinTOF");

	Cascade1  	 ->SetDefaultOutFile(finalhistos);
	Cascade2  	 ->SetDefaultOutFile(finalhistos);
	Cascade3 	 ->SetDefaultOutFile(finalhistos);
	Cascade4 	 ->SetDefaultOutFile(finalhistos);
	Cascade5 	 ->SetDefaultOutFile(finalhistos);
	Cascade6 	 ->SetDefaultOutFile(finalhistos);
	Cascade7 	 ->SetDefaultOutFile(finalhistos);
	Cascade8 	 ->SetDefaultOutFile(finalhistos);

	Cascade1_D  	 ->SetDefaultOutFile(finalhistos);
	Cascade2_D  	 ->SetDefaultOutFile(finalhistos);
	Cascade3_D 	 ->SetDefaultOutFile(finalhistos);
	Cascade4_D 	 ->SetDefaultOutFile(finalhistos);
	Cascade5_D 	 ->SetDefaultOutFile(finalhistos);
	Cascade6_D 	 ->SetDefaultOutFile(finalhistos);
	Cascade7_D 	 ->SetDefaultOutFile(finalhistos);
	Cascade8_D 	 ->SetDefaultOutFile(finalhistos);


	Filler.AddObject2beFilled(Cascade1,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(Cascade2,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(Cascade3,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(Cascade4,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(Cascade5,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(Cascade6,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(Cascade7,GetGenMomentum,GetGenMomentum);
	Filler.AddObject2beFilled(Cascade8,GetGenMomentum,GetGenMomentum);

	Filler.AddObject2beFilled(Cascade1_D,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(Cascade2_D,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(Cascade3_D,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(Cascade4_D,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(Cascade5_D,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(Cascade6_D,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(Cascade7_D,GetRigidity,GetRigidity);
	Filler.AddObject2beFilled(Cascade8_D,GetRigidity,GetRigidity);


	Filler.ReinitializeAll(refill);

	if(!refill&&checkfile) {	
		Cascade1  	->Eval_Efficiency();
		Cascade2  	->Eval_Efficiency();
		Cascade3 	->Eval_Efficiency();
		Cascade4 	->Eval_Efficiency();
		Cascade5 	->Eval_Efficiency();
		Cascade6 	->Eval_Efficiency();
		Cascade7 	->Eval_Efficiency();
		Cascade8 	->Eval_Efficiency();
	
		Cascade1_D  	->Eval_Efficiency();
		Cascade2_D  	->Eval_Efficiency();
		Cascade3_D 	->Eval_Efficiency();
		Cascade4_D 	->Eval_Efficiency();
		Cascade5_D 	->Eval_Efficiency();
		Cascade6_D 	->Eval_Efficiency();
		Cascade7_D 	->Eval_Efficiency();
		Cascade8_D 	->Eval_Efficiency();
		
		Cascade1  	->SaveResults(finalresults);
		Cascade2  	->SaveResults(finalresults);
		Cascade3 	->SaveResults(finalresults);
		Cascade4 	->SaveResults(finalresults);
		Cascade5 	->SaveResults(finalresults);
		Cascade6 	->SaveResults(finalresults);
		Cascade7 	->SaveResults(finalresults);
		Cascade8 	->SaveResults(finalresults);

		Cascade1_D  	->SaveResults(finalresults);
		Cascade2_D  	->SaveResults(finalresults);
		Cascade3_D 	->SaveResults(finalresults);
		Cascade4_D 	->SaveResults(finalresults);
		Cascade5_D 	->SaveResults(finalresults);
		Cascade6_D 	->SaveResults(finalresults);
		Cascade7_D 	->SaveResults(finalresults);
		Cascade8_D 	->SaveResults(finalresults);


	}

	return ;
}




