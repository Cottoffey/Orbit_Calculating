#include "creatingModelValues.h"

int main()
{
    PlanetEphemeris sun;
    PlanetEphemeris jupiter;
    PlanetEphemeris venus;
    PlanetEphemeris uranus;
    PlanetEphemeris neptune;
    PlanetEphemeris saturn;
    PlanetEphemeris mars;
    PlanetEphemeris mercury;
    PlanetEphemeris earth;
    PlanetEphemeris RealOrbit;
    PlanetEphemeris moon;

    jupiter.init ("Data/Jupiter.txt", 0.041666666667);
    earth.init ("Data/Earth.txt", 0.041666666667);
    venus.init ("Data/Venus.txt", 0.041666666667);
    uranus.init ("Data/Uranus.txt", 0.041666666667);
    neptune.init ("Data/Neptune.txt", 0.041666666667);
    saturn.init ("Data/Saturn.txt", 0.041666666667);
    sun.init ("Data/Sun.txt", 0.041666666667);
    mars.init ("Data/Mars.txt", 0.041666666667);
    mercury.init ("Data/Mercury.txt", 0.041666666667);
    moon.init ("Data/Moon.txt", 0.041666666667);

    sun.GM = GMS;
    jupiter.GM = GMJ;
    earth.GM = GME;
    venus.GM = GMV;
    uranus.GM = GMU;
    neptune.GM = GMN;
    saturn.GM = GMST;
    mars.GM = GMMS;
    mercury.GM = GMMC;
    moon.GM = GMMN;

    PlanetEphemeris data[10] = {sun, jupiter, earth, venus, uranus, neptune, saturn, mars, mercury, moon};
    // PlanetEphemeris data[4] = {sun, jupiter, earth, moon};

    // start parameters x,y,z,vx,vy,vz
    std::vector<double> X = {1.469591208242925E+08,  7.299762167917201E+07,  2.056299266163284E+07,  3.859428549646102E+06,  3.244525935598258E+05,  1.492020244998816E+06};

    modeling (X, 6, 0.041666666666666667, data, 10);
    std::cout << "Modeling success\n";

    creatingModelingValues ();
    std::cout << "Creating model values success\n";

    return 0;
}