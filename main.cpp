#include "creatingModelValues.h"
#include "Regression/GaussNewton.h"
#include "processingObservingData.h"
#include "Colors.h"

int main() 
{
    std::cout << "\n\t" << FWHT (">> Program start <<\n") << std::endl;

    std::cout << "> Observatories positions and base measures calculations" << std::endl;

    ProcessObservingData ();
    
    std::cout << FGRN ("Success\n") << std::endl;

    std::cout << "> Planets ephemerises initialization" << std::endl;

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

    jupiter.init("Data/Jupiter.txt", EPHEMERISES_STEP);
    eartht.init("Data/Earth.txt", EPHEMERISES_STEP);
    venus.init("Data/Venus.txt", EPHEMERISES_STEP);
    uranus.init("Data/Uranus.txt", EPHEMERISES_STEP);
    neptune.init("Data/Neptune.txt", EPHEMERISES_STEP);
    saturn.init("Data/Saturn.txt", EPHEMERISES_STEP);
    sun.init("Data/Sun.txt", EPHEMERISES_STEP);
    mars.init("Data/Mars.txt", EPHEMERISES_STEP);
    mercury.init("Data/Mercury.txt", EPHEMERISES_STEP);
    moon.init("Data/Moon.txt", EPHEMERISES_STEP);


    std::cout << FGRN ("Success\n") << std::endl;

    PlanetEphemeris datas[10] = {sun, jupiter, eartht, venus, uranus, neptune, saturn, mars, mercury, moon};
    
    // start parameters x,y,z,vx,vy,vz,a
    std::vector<double> X = {1.469591208242925E+08,  7.299762167917201E+07,  2.056299266163284E+07,  3.859428549646102E+06,  3.244525935598258E+05,  1.492020244998816E+06,
                             1, 0, 0, 0, 0, 0,
                             0, 1, 0, 0, 0, 0,
                             0, 0, 1, 0, 0, 0,
                             0, 0, 0, 1, 0, 0,
                             0, 0, 0, 0, 1, 0,
                             0, 0, 0, 0, 0, 1
                             };

    std::cout << "> Proccesing class initialization" << std::endl;

    GaussNewton gn_solver(6, datas);
    gn_solver.init();
    std::cout << FGRN ("Success\n") << std::endl;

    std::cout << "\t>> Computation start <<" << std::endl;

    vec result = gn_solver.fit(X);

    std::cout << "\n\t" << FWHT (">> Finish << \n") << std::endl;

    return 0;
}