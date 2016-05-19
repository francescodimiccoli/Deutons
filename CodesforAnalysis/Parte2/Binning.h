

class Particle {
   public:
      Particle (float m) :          mass (m) {}
      Particle (float m, float z) : mass (m), Z (z) {}
      float mass=0;
      int Z=1;
      float ekin  =0;
      float mom =0;
      float rig =0;
      float beta=0;

      void FillFromEk (float);
      void FillFromRig(float);

   protected:
      float BetaFromEk (float ek)  { return sqrt (ek*ek + 2 * ek * mass) / (ek + mass); }
      float GammaFromEk (float ek) { return 1 + ek/mass;                                }
      float RigFromEk (float ek)   { return RigFromMom (MomFromEk (ek) ) ;              }

      float MomFromEk  (float ek)  { return mass * BetaFromEk (ek) * GammaFromEk (ek);  }
      float EkFromMom (float p)  { return sqrt ( mass*mass + p*p) - mass  ;       }

      float RigFromMom (float p) { return p/Z ;     }
      float MomFromRig (float rig) { return rig*Z ;     }

};

void Particle::FillFromEk(float ek)
{
   ekin=ek;
   beta = BetaFromEk(ek);
   mom=MomFromEk(ek);
   rig=RigFromMom(mom);
}

void Particle::FillFromRig( float r)
{
   rig=r;
   mom=MomFromRig(rig);
   ekin=EkFromMom(mom);
   beta=BetaFromEk(ekin);
}


class Binning {
   public:
      Binning () : mass(0) {}   ;                       
      Binning (float m) :           mass(m)       {}
      Binning (float m, float z) :  mass(m), Z(1) {}
      void Setbins (int, float, float, int type=1); // type -- binning in 0: not done, 1 energy, 2 rigidity
      int size() {return ekbin.size(); };

      
      std::vector<float> EkBins  ()   {  return   ekbin;   }
      std::vector<float> MomBins ()   {  return  mombin;   }
      std::vector<float> RigBins ()   {  return  rigbin;   }
      std::vector<float> BetaBins()   {  return betabin;   }

      std::vector<float> EkBinsCent()   { return   ekbincent;  }  // bin centers in log
      std::vector<float> MomBinsCent () { return  mombincent;  }
      std::vector<float> RigBinsCent () { return  rigbincent;  }
      std::vector<float> BetaBinsCent() { return  betabincent; }

      std::vector<float> EkPerMassBins  ();  // returns Ek per mass

      int Type() {return type;}


   protected:
      int Z=1;
      float mass;

      int type=0;       // binning in 0: not done, 1 energy, 2 rigidity

      std::vector<float>   ekbin ;
      std::vector<float>  mombin ;
      std::vector<float>  rigbin ;
      std::vector<float> betabin ;
      std::vector<float>  ekmbin ;

      std::vector<float>   ekbincent ;
      std::vector<float>  mombincent ;
      std::vector<float>  rigbincent ;
      std::vector<float> betabincent ;

};




void Binning::Setbins (int nbins, float min, float max, int typ)
{
   type=typ;
   float logmin=log(min), logmax=log(max);
   float binbeg=logmin;
   float binstep= (logmax-logmin)  / nbins;
   int ibin=0;
   std::vector<float> vbin; // bins
   std::vector<float> vcen; // centers
   Particle Pedge(mass, Z), Pcent(mass, Z);

   // Filling the vectors
   while (binbeg<max) {
      binbeg = exp ( logmin + ibin * binstep);
      float bincent= exp ( logmin + (ibin+0.5) * binstep);
      switch(type) {
      case 1: // Energy
         Pedge.FillFromEk(binbeg);
         Pcent.FillFromEk(bincent);
         break;
      case 2: // Rigidity
         Pedge.FillFromRig(binbeg);
         Pcent.FillFromRig(bincent);
         break;
      default:
         // type not implemented;
         return;
      }
      ekbin.  push_back(Pedge.ekin);
      mombin. push_back(Pedge.mom);
      rigbin. push_back(Pedge.rig);
      betabin.push_back(Pedge.beta);

      ekbincent.  push_back(Pcent.ekin);
      mombincent. push_back(Pcent.mom);
      rigbincent. push_back(Pcent.rig);
      betabincent.push_back(Pcent.beta);

      ibin++;
   }

   switch(type) { // Don't forget the final edge
   case 1:
      Pedge.FillFromEk(binbeg);
      break;
   case 2:
      Pedge.FillFromRig(binbeg);
   }

   ekbin.  push_back(Pedge.ekin);
   mombin. push_back(Pedge.mom);
   rigbin. push_back(Pedge.rig);
   betabin.push_back(Pedge.beta);

}


std::vector<float> Binning::EkPerMassBins() {
   if (ekmbin.size()==0 && mass !=0 )
      for (float e : ekbin)   ekmbin.push_back(e/mass);
   return ekmbin;
}



class PBinning: public Binning {
   public:
      PBinning() : Binning (0.9382720813) {}  // proton mass 938 MeV
};
class DBinning: public Binning {
   public:
      DBinning() : Binning (0.18756129  ) {}  // deuterium mass 1876 MeV
};




