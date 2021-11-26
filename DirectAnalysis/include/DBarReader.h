#ifndef  DBARREAD_H
#define  DBARREAD_H
 
#include "TTree.h"
#include "Ntp.h"
#include "Variables.hpp"

// Forward-declare Variable struct
class Variable;

class DBarReader {

private:
    bool isMC;
    TTree * Tree =0x0;	
    TTree * Tree_RTI=0x0;
    TTree * Tree_Cpct=0x0;

    RTIInfo 	  *rtiInfo; 	
    NtpSHeader    *ntpSHeader; 
    NtpHeader     *ntpHeader; 
    NtpMCHeader   *ntpMCHeader;          
    NtpTrd        *ntpTrd;
    NtpTof        *ntpTof;
    NtpTracker    *ntpTracker;
    NtpRich       *ntpRich;
    NtpEcal       *ntpEcal;
    NtpAnti       *ntpAnti;
    NtpStandAlone *ntpStandAlone;

    NtpCompact    *ntpCompact;

    void Init();
    UInt_t getPackedLayers_1to4();
    bool minTOF();
    bool minTOF_Cpct();
    bool goldenTOF();
    int RICHmaskConverter();
    bool goldenTOF_Cpct();
    int RICHmaskConverter_Cpt();

public:

    DBarReader(TTree * f, bool isMC,TTree * f_RTI, TTree * tree_cpct);
    DBarReader(TTree * f, bool isMC);
    
    void FillVariables(int NEvent, Variables * vars);
    void FillCompact(int NEvent, Variables * vars, float massgen=1);
    
    int GetTreeEntries(){return Tree->GetEntries();}; 	
    int GetCompactEntries(){return Tree_Cpct->GetEntries();}; 	
    
    TTree * GetTree() {return Tree;}
    TTree * GetCompactTree() {return Tree_Cpct;}
    Long64_t ProtonCandidateSelection();

};

#endif


