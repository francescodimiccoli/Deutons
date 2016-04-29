using namespace std;

class DatavsMC
{

public:
	Efficiency * MCEff;
	Efficiency * DataEff;

	std::string basename;
	
	//creation constructors
	DatavsMC(std::string basename, int n){
		
		MCEff    = new Efficiency((basename + "_MC"  ));
		DataEff  = new Efficiency((basename + "_Data"),n);
	} 
		
	void Write();

};

void DatavsMC::Write(){
	MCEff 	-> Write();
	DataEff -> Write();
	return;
}





