#include "TFile.h"
#include "Ntp.h"

// Forward-declare Variable struct
class Variable;

class DBarReader : public Analysis {

private:

    NtpHeader     ntpHeader; 
    NtpMCHeader   ntpMCHeader;          
    NtpTrd        ntpTrd;
    NtpTof        ntpTof;
    NtpTracker    ntpTracker;
    NtpRich       ntpRich;
    NtpEcal       ntpEcal;
    NtpAnti       ntpAnti;
    NtpStandAlone ntpStandAlone;

    UInt_t getPackedLayers_1to4();
    bool minTOF();
    bool goldenTOF();
    int RICHmaskConverter();

public:

    DBarReader(TFile * f);
     
    FillVariables(Variables * vars);
}
