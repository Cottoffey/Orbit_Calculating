#include <sofa.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>

#include "Ephemeris.h"

#define EARTH_RADIUS 6378.140 // km
#define SEC_DAY 86400.0 
#define MINSEC_DAY 86400000.0



void ProcessObservingData ()
{
    struct Time
    {
        int y;
        int m;
        int d;
        double fd;
        int h;
        int mn;
        double s;
    };

    std::ifstream input ("Data/observations.txt");
    std::ofstream output ("Data/ObservData.txt");
    output.setf (std::ios::scientific);
    std::cout.setf (std::ios::scientific);

    // For Earth coors
    PlanetEphemeris Earth;
    TimeEphemeris TT_TDB;
    PlanetEphemeris Hubble;
    
    Hubble.init ("Data/Hubble.txt", 0.041666666667);
    Earth.init ("Data/Earth.txt", 0.041666666667);
    Earth.GM = 398600.43552;
    TT_TDB.init ("Data/TDB.txt");

    for (int i = 0; i < 222; i++)
    {
        Time current;
        input >> current.y >> current.m >> current.fd;
        current.d = (int)current.fd;
        current.fd = current.fd - current.d;

        double xp, yp, dut1;
        double l, r, z;
        double coors[3];
        double observCoors[3];
        double radCoors[2];

        input >> observCoors[0] >> observCoors[1] >> observCoors[2];
        radCoors[0] = (observCoors[0] * 15 + observCoors[1] * 0.25 + observCoors[2] * 0.0041666666667) * M_PI / 180;
        input >> observCoors[0] >> observCoors[1] >> observCoors[2];
        if (observCoors[0] < 0.0)
            radCoors[1] = (observCoors[0] - observCoors[1] / 60 - observCoors[2] / 3600) * M_PI / 180;
        else 
            radCoors[1] = (observCoors[0] + observCoors[1] / 60 + observCoors[2] / 3600) * M_PI / 180;

        double tmp;

        // UTC -> TT
        iauDat (current.y, current.m, current.d, current.fd, &tmp);
        Time currentTT = {current.y, current.m, current.d, current.fd + (tmp + 32.184) / SEC_DAY};
        
        if (currentTT.fd >= 1.0)
        {
            currentTT.d = currentTT.d + (int)currentTT.fd;
            currentTT.fd = currentTT.fd - 1.0;
        }

        double TT1, TT2, UT1, UT2;

        // adding h, mn, s, ms
        current.h = (int)(current.fd * 24);
        current.mn = (current.fd * 24.0 - (int)(current.fd * 24)) * 60;
        current.s = ((current.fd * 24.0 - (int)(current.fd * 24)) * 60 - (int)((current.fd * 24.0 - (int)(current.fd * 24)) * 60)) * 60;

        currentTT.h = (int)(currentTT.fd * 24);
        currentTT.mn = (currentTT.fd * 24.0 - (int)(currentTT.fd * 24)) * 60;
        currentTT.s = ((currentTT.fd * 24.0 - (int)(currentTT.fd * 24)) * 60 - (int)((currentTT.fd * 24.0 - (int)(currentTT.fd * 24)) * 60)) * 60;
        
        // to 2-part Julian date
        iauDtf2d("TT", currentTT.y, currentTT.m, currentTT.d, currentTT.h, currentTT.mn, currentTT.s, &TT1, &TT2);
        iauDtf2d("UTC", current.y, current.m, current.d, current.h, current.mn, current.s, &UT1, &UT2);

        std::cout << std::setprecision (15) << i + 1 << ' ' << currentTT.y << ' ' << currentTT.m << ' ' << currentTT.d << ' ' << currentTT.h << ' ' << currentTT.mn << ' ' << currentTT.s << ' ' << TT1+TT2 << std::endl;
        double UT11, UT12;
       
        if ((i > 181 && i < 192) || (i > 196 && i < 202) || i > 206)
        {
            // TT -> TDB
            TT2 -= (TT_TDB.get_dt (TT1 + TT2) / MINSEC_DAY);
            std::cout << std::setprecision (15) << i + 1 << ' ' << currentTT.y << ' ' << currentTT.m << ' ' << currentTT.d << ' ' << currentTT.h << ' ' << currentTT.mn << ' ' << currentTT.s << ' ' << TT1+TT2 << std::endl;
            Hubble.get_coors (TT1 + TT2, coors[0], coors[1], coors[2]);
        }   
        else
        {
            input >> l >> r >> z >> xp >> yp >> dut1;
            // UTC -> UT1
            iauUtcut1(UT1, UT2, dut1, &UT11, &UT12);

            // CylCoors -> DecCoors
            coors[0] = r * cos(l * M_PI / 180.) * EARTH_RADIUS;
            coors[1] = r * sin(l * M_PI / 180.) * EARTH_RADIUS;
            coors[2] = z * EARTH_RADIUS;

            // std::cout << '\n' << coors[0] << ' '<< coors[1] << ' ' << coors[2] << std::endl;
            // celestial to terrestial matrix
            double rc2t[3][3] = {0};
            iauC2t06a(TT1, TT2, UT11, UT12, xp / 206265 , yp / 206265, rc2t);

            // for (int i = 0; i < 3; i++)
            //     std::cout << rc2t[i][0] << ' ' << rc2t[i][1] << ' '<< rc2t[i][2] << std::endl;

            double tmpx = coors[0], tmpy = coors[1], tmpz = coors[2];

            coors[0] = tmpx * rc2t[0][0] + tmpy * rc2t[1][0] + tmpz * rc2t[2][0];
            coors[1] = tmpx * rc2t[0][1] + tmpy * rc2t[1][1] + tmpz * rc2t[2][1];
            coors[2] = tmpx * rc2t[0][2] + tmpy * rc2t[1][2] + tmpz * rc2t[2][2];

            // std::cout << '\n' << coors[0] << ' '<< coors[1] << ' ' << coors[2] << std::endl;
            // TT -> TDB
            TT2 -= (TT_TDB.get_dt (TT1 + TT2) / MINSEC_DAY);
            std::cout << std::setprecision (15) << i + 1 << ' ' << currentTT.y << ' ' << currentTT.m << ' ' << currentTT.d << ' ' << currentTT.h << ' ' << currentTT.mn << ' ' << currentTT.s << ' ' << TT1+TT2 << std::endl;

            // position with Earth
            Earth.get_coors (TT1 + TT2, tmpx, tmpy, tmpz);

            coors[0] = coors[0] + tmpx;
            coors[1] = coors[1] + tmpy;
            coors[2] = coors[2] + tmpz;
        }

        output << std::setprecision (15) << TT1 + TT2 << ' ' << radCoors[0] << ' ' << radCoors[1]
               << ' ' << coors[0] << ' ' << coors[1] << ' ' << coors[2] << std::endl;
    }

    input.close();
    output.close();

    return;
}


int main ()
{
    ProcessObservingData ();

    return -0;
}
