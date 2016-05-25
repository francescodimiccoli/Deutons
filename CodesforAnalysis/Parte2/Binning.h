/** @brief Classed used to convert between kinetic energy, momentum, rigidity and beta
 * @author $Author: L. Basara$
 * @date   $Date: 2016/05/19$
 */


#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <iomanip>

using namespace std;

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











/** @brief Class used for all the binning manipulation
 * @author $Author: L. Basara$
 * @date $Date: 2016/05/19$
 */

class Binning {
   public:
      Binning () : mass(0) {}   ;
      Binning (float m) :           mass(m)       {}
      Binning (float m, float z) :  mass(m), Z(1) {}
      void Setbins (int, float, float, int type=1); ///< type -- binning in 0: not done, 1 energy, 2 rigidity
      int size() {return ekbin.size(); };
      int Type() {return type;}

      struct histo {
         std::vector<float> edges ;
         std::vector<float> content ;
      };


      /** @brief Returns the rigidity bin containing the variable
       *  @param float var : the variable whose location to search
       *  @return int      : the bin number in rigidity of the variable
       */
      int GetRBin(float var);


      std::vector<float> EkBins  ()   {  return   ekbin;   }
      std::vector<float> MomBins ()   {  return  mombin;   }
      std::vector<float> RigBins ()   {  return  rigbin;   }
      std::vector<float> BetaBins()   {  return betabin;   }

      std::vector<float> EkBinsCent()   { return    ekbincent; }  ///< bin centers in log
      std::vector<float> MomBinsCent () { return   mombincent; }
      std::vector<float> RigBinsCent () { return   rigbincent; }
      std::vector<float> BetaBinsCent() { return  betabincent; }

      float EkBin  (int bin)   {  return   ekbin[bin];   }
      float MomBin (int bin)   {  return  mombin[bin];   }
      float RigBin (int bin)   {  return  rigbin[bin];   }
      float BetaBin(int bin)   {  return betabin[bin];   }

      float EkBinCent  (int bin) { return    ekbincent[bin]; }  ///< bin centers in log
      float MomBinCent (int bin) { return   mombincent[bin]; }
      float RigBinCent (int bin) { return   rigbincent[bin]; }
      float BetaBinCent(int bin) { return  betabincent[bin]; }

      std::vector<float> EkPerMassBins  ();  ///< returns Ek per mass

      void SetMatrix(std::vector<float>);  ///< Set the matrix transition from the vector
      void ResetMatrix() { matrix.clear(); } ///< Reset the matrix transition
      std::vector<float> Rebin(histo ); ///< Rebin the histogram passed as parameter to rigbin's bins

      histo LoadCRDB(string filename); ///< Loads a file in the Cosmic-Ray Database data format
      histo LoadData(string filename); ///< Loads a file whose first column is the bin beginning, the second the value
      histo HistoFromVec(std::vector<float> edges,  std::vector<float> content); ///< Make histo from 2 vectors

      void Print(); ///< Print the content of the bins

      void  RFill(TH1* h, float Var); ///< Fill the histogram with var indicating the rigidity bin
      void wRFill(TH1* h, float Var); ///< Fill the histogram with var weighted with the default MC values



   protected:
      float mass;
      int Z=1;
      int type=0;       ///< binning in 0: not done, 1 energy, 2 rigidity

      std::vector<float>   ekbin ;
      std::vector<float>  mombin ;
      std::vector<float>  rigbin ;
      std::vector<float> betabin ;
      std::vector<float>  ekmbin ;

      std::vector<float>   ekbincent ;
      std::vector<float>  mombincent ;
      std::vector<float>  rigbincent ;
      std::vector<float> betabincent ;

      vector< vector<float> > matrix; ///< Migration matrix for fluxes recorded and our binning.
      std::vector<float> mcweight;  ///< MC weights in our binning



};

void Binning::wRFill(TH1* h, float var) {
   if (mcweight.size()==0)
      mcweight=Rebin(LoadData("MCweight.data"));
   int bin=GetRBin(var);
   h->Fill(bin, mcweight[bin]);
}


void Binning::RFill(TH1* h, float var) {
   h->Fill(GetRBin(var));
}



void Binning::Print()
{

   cout << right << setw(3) << setfill(' ') << "Bin";
   cout << right << setw(11) << setfill(' ') << "Ekin";
   cout << right << setw(11) << setfill(' ') << "Momentum";
   cout << right << setw(11) << setfill(' ') << "Rigidity";
   cout << right << setw(11) << setfill(' ') << "Beta";
   cout << endl;

   for (uint i=0; i<ekbin.size(); i++) {
      cout << right << setw(3) << setfill(' ')  << i;
      cout << right << setw(11) << setfill(' ') << setprecision(4) <<   ekbin[i] ;
      cout << right << setw(11) << setfill(' ') << setprecision(4) <<  mombin[i] ;
      cout << right << setw(11) << setfill(' ') << setprecision(4) <<  rigbin[i] ;
      cout << right << setw(11) << setfill(' ') << setprecision(5) << betabin[i] ;
      cout << endl;
   }

   return;

}




Binning::histo Binning::LoadData(string filename)
{

   // Init reading file
   std::fstream Datafile(filename, std::ios_base::in);

   histo histogram;
   float edge, content;
   while(!Datafile.eof()) {
      Datafile >> edge >> content;
      histogram.edges.  push_back(edge);
      histogram.content.push_back(content);
   }
   histogram.edges.  pop_back();
   histogram.content.pop_back();
   histogram.edges.push_back(edge);// last bin

   return histogram;
}


Binning::histo Binning::LoadCRDB(string filename)
{

   // Init reading file
   std::fstream DBfile(filename, std::ios_base::in);
   for (int i=0; i<2; i++) // First 2 lines are descriptors
      DBfile.ignore ( std::numeric_limits<std::streamsize>::max(), '\n' );


   histo histogram;
   float Emoy, Elo, Eup, y, ystat_lo, ystat_up, ysyst_lo, ysyst_up, yerrtot_lo, yerrtot_up; // As described in file
   while(!DBfile.eof()) {
      DBfile >> Emoy >> Elo >> Eup >> y >> ystat_lo >> ystat_up >> ysyst_lo >> ysyst_up >> yerrtot_lo >> yerrtot_up;
      histogram.edges.push_back(Elo);
      histogram.content.push_back(y);
   }
   histogram.edges.pop_back();
   histogram.content.pop_back();
   histogram.edges.push_back(Eup);// last bin

   return histogram;
}


Binning::histo Binning::HistoFromVec(std::vector<float> edges,  std::vector<float> content)
{
   histo h;
   if (edges.size()-content.size() == 1) {
      h.edges=edges;
      h.content=content;
   }
   return h;
}



std::vector<float> Binning::Rebin(histo htorebin)
{
   ResetMatrix();
   std::vector<float> rebinned(rigbin.size()-1); // -1 because rebin is a content, rigbin the edges
   if (htorebin.edges.size()-htorebin.content.size() != 1 ) return rebinned;
   
   SetMatrix(htorebin.edges);

   // What we really do is multiply a matrix by a vector
   for (uint ibin=0; ibin<htorebin.content.size(); ibin++)
      for (uint obin=0; obin<rigbin.size()-1; obin++)
         rebinned[obin] += matrix[ibin][obin] * htorebin.content[ibin];

   return rebinned;
}






void Binning::SetMatrix(std::vector<float> vinput)
{

   //if (matrix.size()!=0) return;     // Already done

   for (uint ib=0; ib<vinput.size()-1; ib++) { // prefix / suffix b for incoming binning
      float bmin=vinput[ib], bmax=vinput[ib+1];
      float bwidth=bmax-bmin;
      vector<float> column;

      for (uint it=0; it<rigbin.size()-1; it++) { // prefix / suffix b for This binning
         float tmin=rigbin[it], tmax=rigbin[it+1];
         float weight=0;
         // Which fraction of the bbin is in tbin? 4 possibilities:
         if      (tmin > bmax || tmax < bmin ) weight=0;                      // either tbin is outside bbin
         else if (tmin > bmin && tmax < bmax ) weight = (tmax-tmin) / bwidth; // or fully included into bbin
         else if (tmin <=bmin && tmax < bmax ) weight = (tmax-bmin) / bwidth; // or partly included to the right
         else if (tmin > bmin && tmax >=bmax ) weight = (bmax-tmin) / bwidth; // ...or to the left
         column.push_back(weight);
      }

      matrix.push_back(column);
   }
   return;
}







void Binning::Setbins (int nbins, float min, float max, int typ)
{
   type=typ;
   float logmin=log(min), logmax=log(max);
   float binbeg=min;
   float binstep= (logmax-logmin)  / nbins;
   int ibin=0;
   std::vector<float> vbin(nbins+1); // bins
   std::vector<float> vcen(nbins); // centers
   Particle Pedge(mass, Z), Pcent(mass, Z);

   // Filling the vectors
   while (binbeg<max) {
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
      binbeg = exp ( logmin + ibin * binstep);
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






std::vector<float> Binning::EkPerMassBins()
{
   if (ekmbin.size()==0 && mass !=0 )
      for (float e : ekbin)   ekmbin.push_back(e/mass);
   return ekmbin;
}






int Binning::GetRBin (float var)
{
   if (var<rigbin[0]) return -1;
   for (uint ib=0; ib<rigbin.size(); ib++)  {
      if (var>rigbin[ib] && var<=rigbin[ib+1])
         return ib;
   }
   return -1;
}










class PBinning: public Binning {
   public:
      PBinning() : Binning (0.9382720813) {}  // proton mass 938 MeV
};
class DBinning: public Binning {
   public:
      DBinning() : Binning (0.18756129  ) {}  // deuterium mass 1876 MeV
};





