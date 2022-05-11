#ifndef COORD_H
#define COORD_H


// Compute Geomagnetic quantities (coords & cutoffs)
// from AMS position and particle direction in GTOD
void coord(int time, float radS, float latS, float phiS, float thetaP, float phiP,
	   float &radSGM, float &latSGM, float &phiSGM, float &rgCutP, float &rgCutN);

// Compute Geomagnetic quantities (coords & cutoffs)
// from AMS position and particle local direction
void coord(int time, float radS, float latS, float phiS, float thetaV, float phiV,
	   float pitch, float yaw, float roll, float thetaL, float phiL,
	   float &radSGM, float &latSGM, float &phiSGM, float &rgCutP, float &rgCutN);

// Geomagnetic dipole quantities from IGRF2010 (1995-2010)                                                                                                                     
void EarthDipole(int time,
		 float &dpValue, float &doGeoLt, float &doGeoLn, 
		 float &dpGeoRd, float &dpGeoLt, float &dpGeoLn);

#endif
