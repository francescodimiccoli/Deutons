using namespace std;

class TriggPatt{

	private:
	int physbpatt = 0;

	public:
	void SetTriggPatt(int   PhysBPatt) {physbpatt = PhysBPatt;};
	void SetTriggPatt(float PhysBPatt) {physbpatt= (int) (PhysBPatt+0.5);};
	
	bool IsPhysical() {return ((physbpatt&0b111110)>0);}
	bool IsUnbias()   {return (((physbpatt>>1)&0b11111)== 0) ; } 
		
};
