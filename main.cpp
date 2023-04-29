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
    PlanetEphemeris eartht;
    PlanetEphemeris moon;
    // PlanetEphemeris RealOrbit;


    jupiter.init ("Data/Jupiter.txt", 0.041666666667);
    eartht.init ("Data/Earth.txt", 0.041666666667);
    venus.init ("Data/Venus.txt", 0.041666666667);
    uranus.init ("Data/Uranus.txt", 0.041666666667);
    neptune.init ("Data/Neptune.txt", 0.041666666667);
    saturn.init ("Data/Saturn.txt", 0.041666666667);
    sun.init ("Data/Sun.txt", 0.041666666667);
    mars.init ("Data/Mars.txt", 0.041666666667);
    mercury.init ("Data/Mercury.txt", 0.041666666667);
    moon.init ("Data/Moon.txt", 0.041666666667);

    // RealOrbit.init ("Data/RealOrbit.txt", 0.041666666667);

    std::cout << "Initialization success\n";

    PlanetEphemeris datas[10] = {sun, jupiter, eartht, venus, uranus, neptune, saturn, mars, mercury, moon};
    // start parameters x,y,z,vx,vy,vz
    std::vector<double> X = {1.469591208242925E+08,  7.299762167917201E+07,  2.056299266163284E+07,  3.859428549646102E+06,  3.244525935598258E+05,  1.492020244998816E+06,
                             1,0,0,0,0,0,
                             0,1,0,0,0,0,
                             0,0,1,0,0,0,
                             0,0,0,1,0,0,
                             0,0,0,0,1,0,
                             0,0,0,0,0,1};

    modeling (X, 42, 0.041666666666666667, datas, 10);
    std::cout << "Modeling success\n";

    creatingModelingValues ();
    std::cout << "Creating model values success\n";

    return 0;
}