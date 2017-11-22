#include "TFile.h"
#include "Ntp.h"

// Forward-declare Variable struct
class Variable;

class DBarReader {

private:
    bool isMC;
    TTree * Tree;	

    NtpHeader     *ntpHeader; 
    NtpMCHeader   *ntpMCHeader;          
    NtpTrd        *ntpTrd;
    NtpTof        *ntpTof;
    NtpTracker    *ntpTracker;
    NtpRich       *ntpRich;
    NtpEcal       *ntpEcal;
    NtpAnti       *ntpAnti;
    NtpStandAlone *ntpStandAlone;


    void Init();
    UInt_t getPackedLayers_1to4();
    bool minTOF();
    bool goldenTOF();
    int RICHmaskConverter();

public:

    DBarReader(TFile * f, bool isMC);
    void FillVariables(int NEvent, Variables * vars);
    int GetTreeEntries(){return Tree->GetEntries();}; 	
};




