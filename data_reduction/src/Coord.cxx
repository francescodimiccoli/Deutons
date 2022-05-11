#include <iostream>
#include <stdio.h>
#include <math.h>

#include "Coord.h"
#include "point.h"

using namespace std;

void EarthDipole(int time,
		 float &dpValue, float &doGeoLt, float &doGeoLn, 
		 float &dpGeoRd, float &dpGeoLt, float &dpGeoLn)
{
  ///////////////////////////////////////////////////////////////
  //
  //  INPUT:                                     
  //         time    >0 unix time
  //                 <0 year, e.g. -2012
  //                 =0 default  (= AMS values for 2012)
  //                                             
  //  OUTPUT:
  //         dpValue Earth's dipole Value (Am^2)
  //         doGeoLt Orientation Latitude  (deg)
  //         doGeoLn             Longitude (deg)
  //         dpGeoRd Position Radius (cm)
  //         dpGeoLt          Latitude (deg)
  //         dpGeoLn          Longitude (deg)
  //
  ///////////////////////////////////////////////////////////////

  int debug = 1;

  // Constants
  const float pi = 3.1415926536, radeg = 180/pi;
  const float rEarth = 6371.2e3;
  const float dateDflt = 2012;

  ///////////////////////////////////////////////////////////////
  // Geomagnetic dipole quantities from IGRF2010 (1995-2010)
  //     DGRF95  1995.00 10  0  0 1995.00 2000.00   -1.0  600.0           DGRF95   0
  // 1  0 -29692.00      0.00      0.00      0.00                         DGRF95   1
  // 1  1  -1784.00   5306.00      0.00      0.00                         DGRF95   2
  // 2  0  -2200.00      0.00      0.00      0.00                         DGRF95   3
  // 2  1   3070.00  -2366.00      0.00      0.00                         DGRF95   4
  // 2  2   1681.00   -413.00      0.00      0.00                         DGRF95   5
  //   DGRF2000  2000.00 13  0  0 2000.00 2005.00   -1.0  600.0         DGRF2000   0
  // 1  0 -29619.40      0.00      0.00      0.00                       DGRF2000   1
  // 1  1  -1728.20   5186.10      0.00      0.00                       DGRF2000   2
  // 2  0  -2267.70      0.00      0.00      0.00                       DGRF2000   3
  // 2  1   3068.40  -2481.60      0.00      0.00                       DGRF2000   4
  // 2  2   1670.90   -458.00      0.00      0.00                       DGRF2000   5
  //   DGRF2005  2005.00 13  0  0 2005.00 2010.00   -1.0  600.0         DGRF2005   0
  // 1  0 -29554.63      0.00      0.00      0.00                       IGRF2005   1
  // 1  1  -1669.05   5077.99      0.00      0.00                       IGRF2005   2
  // 2  0  -2337.24      0.00      0.00      0.00                       IGRF2005   3
  // 2  1   3047.69  -2594.50      0.00      0.00                       IGRF2005   4
  // 2  2   1657.76   -515.43      0.00      0.00                       IGRF2005   5
  //   IGRF2010  2010.00 13  8  0 2010.00 2015.00   -1.0  600.0         IGRF2010   0
  // 1  0 -29496.50      0.00     11.40      0.00                       IGRF2010   1
  // 1  1  -1585.90   4945.10     16.70    -28.80                       IGRF2010   2
  // 2  0  -2396.60      0.00    -11.30      0.00                       IGRF2010   3
  // 2  1   3026.00  -2707.70     -3.90    -23.00                       IGRF2010   4
  // 2  2   1668.60   -575.40      2.70    -12.90                       IGRF2010   5
  ///////////////////////////////////////////////////////////////

  const int YEARS = 4;
  int year[YEARS] = { 1995, 2000, 2005, 2010 };

  float g10[YEARS] = { -29692.00, -29619.40, -29554.63, -29496.50 };
  float g10_roc = 11.40;

  float g11[YEARS] = {  -1784.00,  -1728.20,  -1669.05,  -1585.90 };
  float g11_roc =  16.70; 

  float h11[YEARS] = {   5306.00,   5186.10,   5077.99,   4945.10 };
  float h11_roc =-28.80;

  float g20[YEARS] = {  -2200.00,  -2267.70,  -2337.24,  -2396.60 };
  float g20_roc =-11.30;

  float g21[YEARS] = {   3070.00,   3068.40,   3047.69,   3026.00 };
  float g21_roc = -3.90;

  float h21[YEARS] = {  -2366.00,  -2481.60,  -2594.50,  -2707.70 };
  float h21_roc =-23.00;

  float g22[YEARS] = {   1681.00,   1670.90,   1657.76,   1668.60 };
  float g22_roc =  2.70;

  float h22[YEARS] = {   -413.00,   -458.00,   -515.43,   -575.40 };
  float h22_roc =-12.90;

  float g10v, g11v, h11v, g20v, g21v, h21v, g22v, h22v;
  float B0, c11, L0, L1, L2, E, xd, yd, zd, rd;
  int i;

  static float dateS=0;
  static float dpValueS, doGeoLtS, doGeoLnS, dpGeoRdS, dpGeoLtS, dpGeoLnS;
#pragma omp threadprivate(dateS,dpValueS,doGeoLtS,doGeoLnS,dpGeoRdS,dpGeoLtS,dpGeoLnS)

  // Set current date 
  struct tm t;
  int leap_year;

  float date = dateDflt;
  if (time>0) {
    time_t ut = time;
    localtime_r(&ut, &t);
    leap_year = (((t.tm_year % 4) == 0) &&
		 (((t.tm_year % 100) != 0) || ((t.tm_year % 400) == 0)));
    date =  1900 + t.tm_year +  (t.tm_yday / (365.0 + leap_year));
  }
  else if (time<0)
    date = -time;

  if (date==dateS) {
    dpValue = dpValueS;
    doGeoLt = doGeoLtS;
    doGeoLn = doGeoLnS;
    dpGeoRd = dpGeoRdS;
    dpGeoLt = dpGeoLtS;
    dpGeoLn = dpGeoLnS;
    return;
  }

  // Compute dipole coefficients for current date
  if (time==0) { // AMS software values for 2012
    dpValue = 7.79585e+22;    // Dipole Value (A m2)  
    doGeoLt = -1*-80.130;     // Dipole Direction Lat (rad) - reversed! 
    doGeoLn = 107.622-180.;   //                  Lon (rad) - reversed!
    dpGeoRd = 567.946e5;      //  Dipole Shift Distance  (cm)
    dpGeoLt = 22.596;         //               Latitude  (rad)
    dpGeoLn = 151.365;        //               Longitude (rad)
  }
  else {
    for (i=0;i<YEARS;i++) if (date < year[i]) break;
    if (i==0) {
      g10v = g10[i];
      g11v = g11[i];
      h11v = h11[i];
      g20v = g20[i];
      g21v = g21[i];
      h21v = h21[i];
      g22v = g22[i];
      h22v = h22[i];
    }
    else if (i<YEARS) {
      g10v = g10[i-1] + (date - year[i-1]) * (g10[i] - g10[i-1]) / (year[i] - year[i-1]);
      g11v = g11[i-1] + (date - year[i-1]) * (g11[i] - g11[i-1]) / (year[i] - year[i-1]);
      h11v = h11[i-1] + (date - year[i-1]) * (h11[i] - h11[i-1]) / (year[i] - year[i-1]);
      g20v = g20[i-1] + (date - year[i-1]) * (g20[i] - g20[i-1]) / (year[i] - year[i-1]);
      g21v = g21[i-1] + (date - year[i-1]) * (g21[i] - g21[i-1]) / (year[i] - year[i-1]);
      h21v = h21[i-1] + (date - year[i-1]) * (h21[i] - h21[i-1]) / (year[i] - year[i-1]);
      g22v = g22[i-1] + (date - year[i-1]) * (g22[i] - g22[i-1]) / (year[i] - year[i-1]);
      h22v = h22[i-1] + (date - year[i-1]) * (h22[i] - h22[i-1]) / (year[i] - year[i-1]);
    }
    else {
      g10v = g10[i-1] + (date - year[i-1]) * g10_roc;
      g11v = g11[i-1] + (date - year[i-1]) * g11_roc;
      h11v = h11[i-1] + (date - year[i-1]) * h11_roc;
      g20v = g20[i-1] + (date - year[i-1]) * g20_roc;
      g21v = g21[i-1] + (date - year[i-1]) * g21_roc;
      h21v = h21[i-1] + (date - year[i-1]) * h21_roc;
      g22v = g22[i-1] + (date - year[i-1]) * g22_roc;
      h22v = h22[i-1] + (date - year[i-1]) * h22_roc;
    } 

    // Compute dipole parameters for current date
    B0 = sqrt(g10v*g10v+g11v*g11v+h11v*h11v);
    c11 = sqrt(g11v*g11v+h11v*h11v);
    L0 = 2*g10v*g20v+sqrt(3.)*(g11v*g21v+h11v*h21v);
    L1 =  -g11v*g20v+sqrt(3.)*(g10v*g21v+g11v*g22v+h11v*h22v);
    L2 =  -h11v*g20v+sqrt(3.)*(g10v*h21v-h11v*g22v+g11v*h22v);
    E = (L0*g10v+L1*g11v+L2*h11v)/(4*B0*B0);
    zd = (L0-g10v*E)/(3*B0*B0);
    xd = (L1-g11v*E)/(3*B0*B0);
    yd = (L2-h11v*E)/(3*B0*B0);
    rd = sqrt(xd*xd+yd*yd+zd*zd);

    // value (A m2)
    dpValue = 1e-2*B0*pow(rEarth,3);

    // position in GTOD (cm,deg,deg)
    dpGeoRd = 1e2*rEarth*rd;
    dpGeoLt = radeg*asin(zd/rd);
    dpGeoLn = radeg*atan2(yd,xd);
    
    // orientation in GTOD (deg,deg)
    doGeoLt = radeg*asin(-g10v/B0);
    doGeoLn = radeg*atan2(-h11v/c11,-g11v/c11);
  }

  if (debug) {
    cout << endl
	 << "********** EARTHs DIPOLE (IGRF2010) **********" << endl
	 << endl
	 << "DATE: " << date               << endl
         << endl
	 << "DIPOLE VALUE:"                << endl
         << dpValue << " (A m2)"           << endl
         << endl 
         << "DIPOLE ORIENTATION:"          << endl
         << "  Lat(deg):" << doGeoLt       << endl
         << "  Lon(deg):" << doGeoLn       << endl
         << endl
         << "DIPOLE POSITION:"             << endl
         << "  R   (cm):" << dpGeoRd       << endl
         << "  Lat(deg):" << dpGeoLt       << endl 
         << "  Lon(deg):" << dpGeoLn       << endl
         << endl
         << "********************************************" << endl
         << endl;
  }

  dateS = date;
  dpValueS = dpValue;
  doGeoLtS = doGeoLt;
  doGeoLnS = doGeoLn;
  dpGeoRdS = dpGeoRd;
  dpGeoLtS = dpGeoLt;
  dpGeoLnS = dpGeoLn;

}


void coord(int time, float radS, float latS, float phiS, float thetaP, float phiP,
	   float &radSGM, float &latSGM, float &phiSGM, float &rgCutP, float &rgCutN)
{
  ///////////////////////////////////////////////////////////////
  //
  //  INPUT:
  //         time    >0  unix time
  //                 <0  year, e.g. -2012
  //                 =0  default  (= AMS values for 2012)
  //         radS     station altitude  GTOD(cm)
  //         latS             latitude  GTOD(rad)
  //         phiS             longitude GTOD(rad)
  //         thetaP  particle polar   angle GTOD(rad)
  //         phiP             azimuth angle GTOD(rad)
  //                                             
  //  OUTPUT:                                    
  //         radSGM  station altitude  GEOM(cm)  
  //         latSGM          latitude  GEOM(rad)
  //         phiSGM          longitude GEOM(rad)
  //         rgCutP  rigidity cutoff for positive charge (GV)
  //         rgCutN                  for negative charge (GV)
  //
  ///////////////////////////////////////////////////////////////

  // Constants
  const float pi=3.1415926536, degrd=pi/180.;
  const float c=2.99792458E8, rEarth=6371.2;

  // Transfer matrices and directions
  float GTOD_LGEOM[3][3], pGTOD[3], pLGEOM[3];

  // Other
  float rssGTOD[3], rdpGTOD[3], rssGEOM[3], uDipDir[3], ussGEOM[3];
  float xLGEOM[3], xLGEOMOD;
  float rg0, gmrd, cl, ct;

  // Eccentric Dipole parameters (IGRF model)
  float dpValue, doGeoLt, doGeoLn, dpGeoRd, dpGeoLt, dpGeoLn;
  EarthDipole(time,
	      dpValue, doGeoLt, doGeoLn, dpGeoRd, dpGeoLt, dpGeoLn);

  //
  ///////////////////////////////////////////////////////////////
  // GEOMAGNETIC ALTITUDE, LATITUDE AND LONGITUDE
  //
  //   Defined from the eccentric dipole coordinate system
  //   based on IGRF model
  //
  //   GEOM axes are centered on the eccentric dipole position.
  //   Directions are defined as
  //   z : along the dipole axis and positive towards 
  //       geographic north
  //   y : z x S(0,0,-1)
  //   x : complets the ortogonal set x = z x y
  //
  ///////////////////////////////////////////////////////////////

  // Station Geomagnetic altitude ( r_geom = r_gtod - r_dip )
  rssGTOD[0] = radS*cos(latS)*cos(phiS);
  rssGTOD[1] = radS*cos(latS)*sin(phiS);
  rssGTOD[2] = radS*sin(latS);

  rdpGTOD[0] = dpGeoRd*cos(dpGeoLt*degrd)*cos(dpGeoLn*degrd);
  rdpGTOD[1] = dpGeoRd*cos(dpGeoLt*degrd)*sin(dpGeoLn*degrd);
  rdpGTOD[2] = dpGeoRd*sin(dpGeoLt*degrd);

  rssGEOM[0] = rssGTOD[0]-rdpGTOD[0]; // cartesian
  rssGEOM[1] = rssGTOD[1]-rdpGTOD[1];
  rssGEOM[2] = rssGTOD[2]-rdpGTOD[2];

  radSGM = sqrt(pow(rssGEOM[0],2) + pow(rssGEOM[1],2) + pow(rssGEOM[2],2) );

  // Station Geomagnetic Latitude & Longitude
  uDipDir[0] = cos(doGeoLt*degrd)*cos(doGeoLn*degrd);
  uDipDir[1] = cos(doGeoLt*degrd)*sin(doGeoLn*degrd);
  uDipDir[2] = sin(doGeoLt*degrd);

  ussGEOM[0] = rssGEOM[0]/radSGM;
  ussGEOM[1] = rssGEOM[1]/radSGM;
  ussGEOM[2] = rssGEOM[2]/radSGM;

  latSGM = asin(uDipDir[0]*ussGEOM[0]+
		uDipDir[1]*ussGEOM[1]+
		uDipDir[2]*ussGEOM[2]);

  phiSGM = atan2( -uDipDir[1]*ussGEOM[0]+uDipDir[0]*ussGEOM[1],
		   uDipDir[0]*uDipDir[2]*ussGEOM[0]+uDipDir[1]*uDipDir[2]*ussGEOM[1]
		  -(pow(uDipDir[0],2)+pow(uDipDir[1],2))*ussGEOM[2]);

  ///////////////////////////////////////////////////////////////
  //   transfer Greenwich true-of-date coordinates to local
  //   geomagnetic coordinate system
  //
  //   LGEOM axes are defined as
  //
  //   z : along the geomagnetic radius vector to the station and
  //       positive toward the zenith
  //   x : pointing to the local geomagnetic east
  //   y : completes the right handed orthogonal system
  ///////////////////////////////////////////////////////////////

  // LGEOMz
  GTOD_LGEOM[2][0] = ussGEOM[0];
  GTOD_LGEOM[2][1] = ussGEOM[1];
  GTOD_LGEOM[2][2] = ussGEOM[2];

  // LGEOMx  ( D x R )
  xLGEOM[0] = uDipDir[1]*ussGEOM[2]-uDipDir[2]*ussGEOM[1];
  xLGEOM[1] = uDipDir[2]*ussGEOM[0]-uDipDir[0]*ussGEOM[2];
  xLGEOM[2] = uDipDir[0]*ussGEOM[1]-uDipDir[1]*ussGEOM[0];
  xLGEOMOD = sqrt( pow(xLGEOM[0],2) + pow(xLGEOM[1],2) + pow(xLGEOM[2],2) );

  GTOD_LGEOM[0][0] = xLGEOM[0]/xLGEOMOD;
  GTOD_LGEOM[0][1] = xLGEOM[1]/xLGEOMOD;
  GTOD_LGEOM[0][2] = xLGEOM[2]/xLGEOMOD;

  // LGEOMy
  GTOD_LGEOM[1][0] = 
     GTOD_LGEOM[2][1]*GTOD_LGEOM[0][2]
    -GTOD_LGEOM[2][2]*GTOD_LGEOM[0][1];
  GTOD_LGEOM[1][1] = 
     GTOD_LGEOM[2][2]*GTOD_LGEOM[0][0]
    -GTOD_LGEOM[2][0]*GTOD_LGEOM[0][2];
  GTOD_LGEOM[1][2] = 
     GTOD_LGEOM[2][0]*GTOD_LGEOM[0][1]
    -GTOD_LGEOM[2][1]*GTOD_LGEOM[0][0];

  // Particle direction in GTOD
  pGTOD[0] = sin(thetaP)*cos(phiP);
  pGTOD[1] = sin(thetaP)*sin(phiP);
  pGTOD[2] = cos(thetaP);

  // Particle coordinates in LGEOM (pLGEOM = GTOD_LGEOM * pGTOD )
  for (int i=0; i<3; i++) {
    pLGEOM[i] = 0;
    for (int j=0; j<3; j++)
      pLGEOM[i]=pLGEOM[i]+GTOD_LGEOM[i][j]*pGTOD[j];
  }

  // Rigidity cutoff from Stormer's formula
  rg0=dpValue*1.E-22*c/pow(rEarth,2);
  gmrd=radSGM*1.E-5/rEarth;
  cl=cos(latSGM);
  ct=pLGEOM[0];
  rgCutP= rg0/pow(gmrd,2)*pow(cl,4)/pow(1+sqrt(1.+ct*pow(cl,3)),2);
  rgCutN=-rg0/pow(gmrd,2)*pow(cl,4)/pow(1+sqrt(1.-ct*pow(cl,3)),2);
}


void coord(int time, float radS, float latS, float phiS, float thetaV, float phiV,
	   float pitch, float yaw, float roll, float thetaL, float phiL,
	   float &radSGM, float &latSGM, float &phiSGM, float &rgCutP, float &rgCutN)
{
  ///////////////////////////////////////////////////////////////
  //
  //  INPUT:
  //         time    >0  unix time
  //                 <0  year, e.g. -2012
  //                 =0  default  (= AMS values for 2012)
  //         radS     station altitude  GTOD(cm)
  //         latS             latitude  GTOD(rad)
  //         phiS             longitude GTOD(rad)
  //         thetaV           velocity vector latitude GTOD(rad)
  //         phiV             velocity vector longitude GTOD(rad)
  //         pitch            pitch (rad)
  //         yaw              yaw   (rad)
  //         roll             roll  (rad)
  //         thetaL   particle polar   angle AMS(rad)
  //         phiL              azimuth angle AMS(rad)
  //                                             
  //  OUTPUT:                                    
  //         radSGM  station altitude  GEOM(cm)  
  //         latSGM          latitude  GEOM(rad)
  //         phiSGM          longitude GEOM(rad)
  //         rgCutP  rigidity cutoff for positive charge (GV)
  //         rgCutN                  for negative charge (GV)
  //
  ///////////////////////////////////////////////////////////////

  // Constants
  const float pi=3.1415926536, degrd=pi/180.;
  const float aAMS=12.0001*degrd; // inclination of AMS-02 

  // Transfer matrices and directions
  float E_BODY_LVLH[3][3], E_LVLH_GTOD[3][3];
  float A_LAMS[3], A_BODY[3], A_LVLH[3], A_GTOD[3];
  float P_LAMS[3], P_BODY[3], P_LVLH[3], P_GTOD[3];

  // Other
  float X_LVLH[3], Y_LVLH[3], Z_LVLH[3];
  float cp, sp, cr, sr, cy, sy;
  float x, y, z, xp, yp, zp, xmod;
  float thetaP, phiP;

  // AMS/Particle coordinates in AMS local system
  A_LAMS[0] = 0;
  A_LAMS[1] = 0;
  A_LAMS[2] = 1;

  P_LAMS[0] = sin(thetaL)*cos(phiL);
  P_LAMS[1] = sin(thetaL)*sin(phiL);
  P_LAMS[2] = cos(thetaL);

  //printf(">> A_LAMS : %10.5E %10.5E %10.5E\n", A_LAMS[0], A_LAMS[1], A_LAMS[2]);
  //printf(">> P_LAMS : %10.5E %10.5E %10.5E\n", P_LAMS[0], P_LAMS[1], P_LAMS[2]);

  ///////////////////////////////////////////////////////////////
  // 
  //  transfer AMS local coordinates to station body coordinates
  //
  //  First rotate -12deg about y azis and then (-y,-x',-z') 
  //
  ///////////////////////////////////////////////////////////////
  A_BODY[0] = - ( A_LAMS[1] );
  A_BODY[1] = - ( A_LAMS[0]*cos(aAMS) + A_LAMS[2]*sin(aAMS) ); 
  A_BODY[2] = - (-A_LAMS[0]*sin(aAMS) + A_LAMS[2]*cos(aAMS) );

  P_BODY[0] = - ( P_LAMS[1] );
  P_BODY[1] = - ( P_LAMS[0]*cos(aAMS) + P_LAMS[2]*sin(aAMS) ); 
  P_BODY[2] = - (-P_LAMS[0]*sin(aAMS) + P_LAMS[2]*cos(aAMS) );

  //printf(">> A_BODY : %10.5E %10.5E %10.5E\n", A_BODY[0], A_BODY[1], A_BODY[2]);
  //printf(">> P_BODY : %10.5E %10.5E %10.5E\n", P_BODY[0], P_BODY[1], P_BODY[2]);

  ///////////////////////////////////////////////////////////////
  //
  //  transfer station body coordinates to the local orbital
  //  coordinate system. 
  //
  //  The Euler angles pitch, yaw and roll define the LVLH axes 
  //  in the station body system
  //
  //  The rotation matrix is Rz(yaw)Ry(pitch)Rx(roll)
  //
  ///////////////////////////////////////////////////////////////
  cr = cos(roll);
  sr = sin(roll);
  cp = cos(pitch);
  sp = sin(pitch);
  cy = cos(yaw);
  sy = sin(yaw); 

  // LVLHx
  E_BODY_LVLH[0][0] = cy*cp;
  E_BODY_LVLH[0][1] = cy*sp*sr - sy*cr;
  E_BODY_LVLH[0][2] = cy*sp*cr + sy*sr;

  // LVLHy
  E_BODY_LVLH[1][0] = sy*cp;
  E_BODY_LVLH[1][1] = sy*sp*sr + cy*cr;
  E_BODY_LVLH[1][2] = sy*sp*cr - cy*sr;

  // LVLHz
  E_BODY_LVLH[2][0] =-sp;
  E_BODY_LVLH[2][1] = cp*sr;
  E_BODY_LVLH[2][2] = cp*cr;

  // AMS/Particle coordinates in LVLH ( P_LVLH = E_BODY_LVLH * P_BODY )
  for (int i=0; i<3; i++) {
    A_LVLH[i] = 0;
    P_LVLH[i] = 0;
    for (int j=0; j<3; j++) {
      A_LVLH[i] += E_BODY_LVLH[i][j]*A_BODY[j];
      P_LVLH[i] += E_BODY_LVLH[i][j]*P_BODY[j];
    }
  }
  //printf(">> A_LVLH : %10.5E %10.5E %10.5E\n", A_LVLH[0], A_LVLH[1], A_LVLH[2]);
  //printf(">> P_LVLH : %10.5E %10.5E %10.5E\n", P_LVLH[0], P_LVLH[1], P_LVLH[2]);

  ///////////////////////////////////////////////////////////////
  //
  // transfer station local orbital coordinates to Greenwich
  // true-of-date coordinate system
  //
  // LVLH axes are defined as
  //
  // z : along the geocentric radius vector to the station and
  //     positive toward the center of the earth
  // x : lies in the station orbit plane, perpendicular to the
  //     z axis and positive in the direction of motion
  // y : completes the right handed orthogonal system
  // 
  ///////////////////////////////////////////////////////////////

  // Station position direction in GTOD
  x = cos(latS)*cos(phiS);
  y = cos(latS)*sin(phiS);
  z = sin(latS);

  // Station speed direction in GTOD
  xp = cos(thetaV)*cos(phiV);
  yp = cos(thetaV)*sin(phiV);
  zp = sin(thetaV);

  // Z in GTOD ( -R )
  Z_LVLH[0] = -x;
  Z_LVLH[1] = -y;
  Z_LVLH[2] = -z;

  // X_LVLH in GTOD ( (RxV)xR ) (and correct for non orthogonality)
  X_LVLH[0] = z*z*xp - zp*x*z - x*y*yp + xp*y*y;
  X_LVLH[1] = x*x*yp - xp*x*y - y*z*zp + yp*z*z;
  X_LVLH[2] = y*y*zp - yp*y*z - x*z*xp + zp*x*x;

  xmod = sqrt(pow(X_LVLH[0],2)+pow(X_LVLH[1],2)+pow(X_LVLH[2],2));

  X_LVLH[0] /= xmod;
  X_LVLH[1] /= xmod;
  X_LVLH[2] /= xmod;

  // Y_LVLH in GTOD ( ZxX )
  Y_LVLH[0] = Z_LVLH[1]*X_LVLH[2] - Z_LVLH[2]*X_LVLH[1];
  Y_LVLH[1] = Z_LVLH[2]*X_LVLH[0] - Z_LVLH[0]*X_LVLH[2];
  Y_LVLH[2] = Z_LVLH[0]*X_LVLH[1] - Z_LVLH[1]*X_LVLH[0];

  E_LVLH_GTOD[0][0]=X_LVLH[0];
  E_LVLH_GTOD[1][0]=X_LVLH[1];
  E_LVLH_GTOD[2][0]=X_LVLH[2];
  E_LVLH_GTOD[0][1]=Y_LVLH[0];
  E_LVLH_GTOD[1][1]=Y_LVLH[1];
  E_LVLH_GTOD[2][1]=Y_LVLH[2];
  E_LVLH_GTOD[0][2]=Z_LVLH[0];
  E_LVLH_GTOD[1][2]=Z_LVLH[1];
  E_LVLH_GTOD[2][2]=Z_LVLH[2];

  // AMS/Particle coordinates in GTOD ( P_GTOD = E_LVLH_GTOD * P_LVLH )
  for (int i=0; i<3; i++) {
    A_GTOD[i] = 0;
    P_GTOD[i] = 0;
    for (int j=0; j<3; j++) {
      A_GTOD[i] += E_LVLH_GTOD[i][j]*A_LVLH[j];
      P_GTOD[i] += E_LVLH_GTOD[i][j]*P_LVLH[j];
    }
  }
  //printf(">> A_GTOD : %10.5E %10.5E %10.5E\n", A_GTOD[0], A_GTOD[1], A_GTOD[2]);
  //printf(">> P_GTOD : %10.5E %10.5E %10.5E\n", P_GTOD[0], P_GTOD[1], P_GTOD[2]);

  // Particle polar & azimuth angles in GTOD
  thetaP = acos(P_GTOD[2]);
  phiP = atan2(P_GTOD[1],P_GTOD[0]);
  //printf(">> thetaP : %10.5E phiP : %10.5E\n", thetaP, phiP);

  ///////////////////////////////////////////////////////////////
  // 
  // Compute Geomagnetic quantities (coords & cutoffs)
  // from AMS position and particle direction in GTOD
  //
  ///////////////////////////////////////////////////////////////
  coord(time, radS, latS, phiS, thetaP, phiP,
	radSGM, latSGM, phiSGM, rgCutP, rgCutN);

}
