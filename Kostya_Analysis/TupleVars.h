#ifndef TUPLEVARS__H__
#define TUPLEVARS__H__

#include "TNtuple.h"

struct TupleVars {

    enum MCTYPE { Proton, Deuteron, Helium };

    float Ev_Num, Trig_Num;
    float R, R_pre;
    float Beta, Beta_pre, BetaRICH;
    float EdepECAL, EdepL1, EdepTOFD, EdepTOFU, EdepTrack;
    float BDT_response, LDiscriminant;
    float Dist5D, Dist5D_P, Rmin;
    float Cutmask;
    float Rcutoff, Latitude;
    float Momento_gen;
    float MC_type;
    float Unbias;

    // This checks that the branchname presents int the ntuple
    // if it is not the struct field is set to NAN
    void SetIfExists(TNtuple * ntuple ,const char * brname, float * ptr) {
        if(!ntuple->GetBranch(brname)) {
            *ptr = nanf("");
            return;
        }
        ntuple->SetBranchAddress(brname, ptr);
    }

    TupleVars(TNtuple * ntuple) {
        SetIfExists(ntuple, "Ev_Num",      &Ev_Num     );
        SetIfExists(ntuple, "Trig_Num",    &Trig_Num   );
        SetIfExists(ntuple, "Rcutoff",     &Rcutoff    );
        SetIfExists(ntuple, "Latitude",    &Latitude   );
        SetIfExists(ntuple, "Momento_gen", &Momento_gen);
        SetIfExists(ntuple, "R_pre",       &R_pre      );
        SetIfExists(ntuple, "Beta_pre",    &Beta_pre   );
        SetIfExists(ntuple, "Cutmask",     &Cutmask    );
        SetIfExists(ntuple, "MC_type",     &MC_type    );
        SetIfExists(ntuple, "EdepL1",      &EdepL1     );
        SetIfExists(ntuple, "EdepTOFU",    &EdepTOFU   );
        SetIfExists(ntuple, "EdepTOFD",    &EdepTOFD   );
        SetIfExists(ntuple, "EdepTrack",   &EdepTrack  );
        SetIfExists(ntuple, "EdepECAL",    &EdepECAL   );
        SetIfExists(ntuple, "BetaRICH",    &BetaRICH   );
        SetIfExists(ntuple, "Unbias",      &Unbias     );
        SetIfExists(ntuple, "R",             &R            );
        SetIfExists(ntuple, "Rmin",          &Rmin         );
        SetIfExists(ntuple, "BDT_response",  &BDT_response );
        SetIfExists(ntuple, "LDiscriminant", &LDiscriminant);
        SetIfExists(ntuple, "Cutmask",       &Cutmask      );
        SetIfExists(ntuple, "Dist5D",        &Dist5D       );
        SetIfExists(ntuple, "Dist5D_P",      &Dist5D_P     );
    }


    MCTYPE GetMCType(){
        if ( ( ( (int) MC_type) & 0xFF    ) >0) return Proton;
        if ( ( ( (int) MC_type) & 0xFF00  ) >0) return Deuteron;
        if ( ( ( (int) MC_type) & 0xFF0000) >0) return Helium;
    }

    float GeneratedMass(){
        switch(GetMCType()){
            case Proton:    return 0.938;
            case Deuteron:  return 1.875;
            case Helium:    return 3.725;    
            default: return 0.0;
        }   
    }

    //retrieve MC cross section type
    int MCsubtype()
    {
       int mc_type=-1;
       int cursor=0;
       if (GeneratedMass() <1&&GeneratedMass() >0) cursor=0 ;
       if (GeneratedMass() <2&&GeneratedMass() >1) cursor=8 ;
       if (GeneratedMass() <4&&GeneratedMass() >3) cursor=16;
       for (int i=2; i<8; i++) {
          if ( ( ( ( (int) MC_type) >> (cursor+i) ) & 1 ) ==1) mc_type=i-2;
       }
       return mc_type;
    }
};

#endif
