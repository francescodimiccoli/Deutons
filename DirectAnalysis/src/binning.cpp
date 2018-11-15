#include "binning.h" 

void Binning::Reset(){
	
	ekbin.clear() ;
	etotbin .clear() ;
	mombin.clear()  ;
	rigbin.clear()  ;
	betabin.clear()  ;
	ekpermassbin.clear()  ;

	ekbincent.clear()  ;
	etotbincent.clear()  ;
	mombincent.clear()  ;
	rigbincent.clear()  ;           	
	betabincent.clear()  ;
	ekpermassbincent.clear()  ;

	ekbin_TOI.clear()  ;
	etotbin_TOI.clear()  ;
	mombin_TOI.clear()  ;
	rigbin_TOI.clear()  ;
	betabin_TOI.clear()  ;
	ekpermassbin_TOI.clear()  ;

	ekbincent_TOI.clear()  ;
	etotbincent_TOI.clear()  ;
	mombincent_TOI.clear()  ;
	rigbincent_TOI.clear()  ;
	betabincent_TOI.clear()  ;
	ekpermassbincent_TOI.clear()  ;

}

void Binning::setBinsFromEk (int nbins, float min, float max, TF1 * SlowDownModel, float alpha_slowdown, float gamma_slowdown)
{
	particle.SetSlowDownModel(SlowDownModel,alpha_slowdown,gamma_slowdown);
	std::vector<float> vedg=computeLogBinEdges(nbins, min, max);
	for (float binedge:vedg) {
		particle.FillFromEk (binedge);
		pushBackVelocities();
	}

	std::vector<float> vcen=computeLogBinCenters(nbins, min, max);
	for (float bincenter:vcen) {
      particle.FillFromEk (bincenter);
      pushBackCentralVelocities();
   }

   return;
}

void Binning::setBinsFromRigidity (int nbins, float min, float max, TF1 * SlowDownModel, float alpha_slowdown, float gamma_slowdown)
{
   particle.SetSlowDownModel(SlowDownModel,alpha_slowdown,gamma_slowdown);
   std::vector<float> vedg=computeLogBinEdges(nbins, min, max);
   for (float binedge:vedg) {
      particle.FillFromRig (binedge);
      pushBackVelocities();
   }
   
   std::vector<float> vcen=computeLogBinCenters(nbins, min, max);
   for (float bincenter:vcen) {
      particle.FillFromRig (bincenter);
      pushBackCentralVelocities();
   }

   return;
}

void Binning::setBinsFromBeta (int nbins, float min, float max, TF1 * SlowDownModel, float alpha_slowdown, float gamma_slowdown)
{
   particle.SetSlowDownModel(SlowDownModel,alpha_slowdown,gamma_slowdown);
   std::vector<float> vedg=computeConstResoBinEdges(nbins, min, max);
   for (float binedge:vedg) {
      particle.FillFromBeta (binedge);
      pushBackVelocities();
   }
   
   std::vector<float> vcen=computeConstResoBinCenters(nbins, min, max);
   for (float bincenter:vcen) {
      particle.FillFromBeta (bincenter);
      pushBackCentralVelocities();
   }

   return;
}



void Binning::setBinsFromEkPerMass (int nbins, float min, float max, TF1 * SlowDownModel,float alpha_slowdown, float gamma_slowdown)
{
   particle.SetSlowDownModel(SlowDownModel,alpha_slowdown,gamma_slowdown);
   std::vector<float> vedg=computeLogBinEdges(nbins, min, max);
   for (float binedge:vedg) {
      particle.FillFromEkPerMass (binedge);
      pushBackVelocities();
   }
   
   std::vector<float> vcen=computeLogBinCenters(nbins, min, max);
   for (float bincenter:vcen) {
      particle.FillFromEkPerMass (bincenter);
      pushBackCentralVelocities();
   }

   return;
}


std::vector<float> Binning::computeLogBinEdges(int nbins, float min, float max) {
   std::vector<float> binEdges(nbins+1);
   float logmin=log (min), logmax=log (max);
   float binstep= (logmax-logmin)  / nbins;
   for (int ibin=0; ibin<nbins; ibin++) {
      float binbeg=min*exp (ibin * binstep);
      binEdges[ibin]=binbeg;
   }
    binEdges[nbins]=max; // So it has the "right" value for sure
   return binEdges;
}

std::vector<float> Binning::computeLogBinCenters(int nbins, float min, float max) {
   std::vector<float> binCent(nbins);
   float logmin=log (min), logmax=log (max);
   float binstep= (logmax-logmin)  / nbins;
   for (int ibin=0; ibin<nbins; ibin++) {
      float bincent=min* exp ( (ibin+0.5) * binstep);
      binCent[ibin]=bincent;
   }
   return binCent;
   
}

std::vector<float> Binning::computeConstResoBinEdges(int nbins, float min, float max) {
   std::vector<float> binEdges(nbins+1);
   binEdges[0]=min;
   for (int ibin=1; ibin<nbins; ibin++) {
       binEdges[ibin]=binEdges[ibin-1]*1.035;
   }
   binEdges[nbins]=max; // So it has the "right" value for sure
   return binEdges;
}

std::vector<float> Binning::computeConstResoBinCenters(int nbins, float min, float max) {
   std::vector<float> binCent(nbins);
   binCent[0]=min*1.0175;
  for (int ibin=1; ibin<nbins; ibin++) {
	binCent[ibin]=binCent[ibin-1]*1.035;
	}
   return binCent;
   
}




void Binning::pushBackVelocities ()
{
   ekbin.  push_back      (particle.getEkin()        );
   etotbin.push_back      (particle.getEtot()        );
   mombin. push_back      (particle.getMom()         );
   rigbin. push_back      (particle.getRig()         );
   betabin.push_back      (particle.getBeta()        );
   ekpermassbin.push_back (particle.getEkinPerMass() );

   ekbin_TOI.  push_back      (particle.getEkin_TOI()        );
   etotbin_TOI.push_back      (particle.getEtot_TOI()        );
   mombin_TOI. push_back      (particle.getMom_TOI()         );
   rigbin_TOI. push_back      (particle.getRig_TOI()         );
   betabin_TOI.push_back      (particle.getBeta_TOI()        );
   ekpermassbin_TOI.push_back (particle.getEkinPerMass_TOI() );
   return;
}

void Binning::pushBackCentralVelocities ()
{
   ekbincent.  push_back      (particle.getEkin()        );
   etotbincent.push_back      (particle.getEtot()        );
   mombincent. push_back      (particle.getMom()         );
   rigbincent. push_back      (particle.getRig()         );
   betabincent.push_back      (particle.getBeta()        );
   ekpermassbincent.push_back (particle.getEkinPerMass() );
  
   ekbincent_TOI.  push_back      (particle.getEkin_TOI()        );
   etotbincent_TOI.push_back      (particle.getEtot_TOI()        );
   mombincent_TOI. push_back      (particle.getMom_TOI()         );
   rigbincent_TOI. push_back      (particle.getRig_TOI()         );
   betabincent_TOI.push_back      (particle.getBeta_TOI()        );
   ekpermassbincent_TOI.push_back (particle.getEkinPerMass_TOI() );
   return;
}


void Binning::Print()
{
   printMatrix::print(
      { ekbin, mombin, rigbin, betabin ,ekbin_TOI, mombin_TOI, rigbin_TOI, betabin_TOI, ekpermassbincent_TOI },

      {"Ekin", "Momentum", "Rigidity", "Beta", "EkinTOI", "MomTOI", "RigTOI", "BetaTOI", "EkpermassCTOI"}
    
      //{betabin},
     // {"Beta"}
   

    );
   return;
}



std::vector<float> Binning::EtotPerMassBins()
{
   std::vector<float> etotpermass (0);
   if (particle.getMass()!=0)
      for (float e : ekbin)   etotpermass.push_back (e/particle.getMass());
   return etotpermass;
}

int Binning::GetBin(float var){
	if(Redges)    return GetRBin    (var);
	if(Betaedges) return GetBetaBin (var);
	if(RTOIedges)    return GetRTOIBin    (var);
	if(BetaTOIedges) return GetBetaTOIBin (var);
}

float Binning::GetBinLowEdge(int bin){
        if(Redges)    return RigBin (bin);
        if(Betaedges) return BetaBin (bin);
}

float Binning::GetBinCenter(int bin){
        if(Redges)    return RigBinCent (bin);
        if(Betaedges) return BetaBinCent (bin);
}

int Binning::GetRBin (float var)
{
   if (var<rigbin[0]||var>rigbin[rigbin.size()-1]) return -1;
   for (uint ib=0; ib<rigbin.size(); ib++)  {
      if (var>rigbin[ib] && var<=rigbin[ib+1])
         return ib;
   }
   return -1;
}

int Binning::GetBetaBin (float var)
{
   if (var<betabin[0]||var>betabin[betabin.size()-1]) return -1;
   for (uint ib=0; ib<betabin.size(); ib++)  {
      if (var>betabin[ib] && var<=betabin[ib+1])
         return ib;
   }
   return -1;
}

int Binning::GetRTOIBin (float var)
{
   if (var<rigbin_TOI[0]||var>rigbin_TOI[rigbin.size()-1]) return -1;
   for (uint ib=0; ib<rigbin.size(); ib++)  {
      if (var>rigbin_TOI[ib] && var<=rigbin_TOI[ib+1])
         return ib;
   }
   return -1;
}

int Binning::GetBetaTOIBin (float var)
{
   if (var<betabin_TOI[0]||var>betabin_TOI[betabin.size()-1]) return -1;
   for (uint ib=0; ib<betabin.size(); ib++)  {
      if (var>betabin_TOI[ib] && var<=betabin_TOI[ib+1])
         return ib;
   }
   return -1;
}


void Binning::RFill (TH1* h, float var)
{
   h->Fill (GetRBin (var) );
}




