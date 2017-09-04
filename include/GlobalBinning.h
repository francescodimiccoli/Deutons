
int nbinsr=43;
int nbinsToF=18;
int nbinsNaF=18;
int nbinsAgl=18;

float ToFsmearSigma= 79.8;
float ToFsmearShift= -1.5;

float FRAC = 1;

Particle proton(0.9382720813, 1, 1);  // proton mass 938 MeV
Particle deuton(1.8756129   , 1, 2);  // deuterium mass 1876 MeV, Z=1, A=2

TRandom3 * Rand= new TRandom3(time(0));


Binning ToFDB(deuton);
Binning ToFPB(proton);
Binning NaFDB(deuton);
Binning NaFPB(proton);
Binning AglDB(deuton);
Binning AglPB(proton);

Binning DRB(deuton);
Binning PRB(proton);


//resolution binning
Binning ToFResB(deuton);
Binning NaFResB(deuton);
Binning AglResB(deuton);

Binning PResB(proton);




void SetBins(){	

	DRB.setBinsFromRigidity(nbinsr, 0.5, 100); 
	PRB.setBinsFromRigidity(nbinsr, 0.5, 100);

	float ekmin=0.15, ekmax=1;
	ToFDB.setBinsFromEkPerMass (nbinsToF, ekmin, ekmax);
	ToFPB.setBinsFromEkPerMass(nbinsToF, ekmin, ekmax);

	ekmin=0.666, ekmax=4.025;
	NaFDB.setBinsFromEkPerMass(nbinsNaF, ekmin, ekmax);
	NaFPB.setBinsFromEkPerMass(nbinsNaF, ekmin, ekmax);

	ekmin=2.57, ekmax=9.01;
	AglDB.setBinsFromEkPerMass(nbinsAgl, ekmin, ekmax);
	AglPB.setBinsFromEkPerMass(nbinsAgl, ekmin, ekmax);


	PResB.setBinsFromRigidity(60, 0.5, 100);	
	ToFResB.setBinsFromRigidity(25, 1,8);
	NaFResB.setBinsFromRigidity(25, 3,13);
	AglResB.setBinsFromRigidity(25, 6,25);



	return;
}


void UpdateProgressBar(int currentevent, int totalentries)
{
        int newratio =(int)100*(currentevent/    (float)totalentries);
        int oldratio =(int)100*((currentevent-1)/(float)totalentries);
        if(newratio>oldratio)
                std::cout<<'\r' << "Progress : "<< newratio+1 << " %"<< std::flush; //+1 pour finir a 100%
}

