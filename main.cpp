#include "creatingModelValues.h"
#include "Regression/GaussNewton.h"

int main() {
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


    jupiter.init("Data/Jupiter.txt", 0.041666666667);
    eartht.init("Data/Earth.txt", 0.041666666667);
    venus.init("Data/Venus.txt", 0.041666666667);
    uranus.init("Data/Uranus.txt", 0.041666666667);
    neptune.init("Data/Neptune.txt", 0.041666666667);
    saturn.init("Data/Saturn.txt", 0.041666666667);
    sun.init("Data/Sun.txt", 0.041666666667);
    mars.init("Data/Mars.txt", 0.041666666667);
    mercury.init("Data/Mercury.txt", 0.041666666667);
    moon.init("Data/Moon.txt", 0.041666666667);

    // RealOrbit.init ("Data/RealOrbit.txt", 0.041666666667);

    std::cout << "Initialization success\n";

    PlanetEphemeris datas[10] = {sun, jupiter, eartht, venus, uranus, neptune, saturn, mars, mercury, moon};
    // start parameters x,y,z,vx,vy,vz
    std::vector<double> X = {1.468787090096414E+08,  7.299085877471100E+07,  2.053190793311784E+07, 
                             3.860105756034324E+06,  3.247863089146682E+05,  1.492113690745696E+06,
                             1, 0, 0, 0, 0, 0,
                             0, 1, 0, 0, 0, 0,
                             0, 0, 1, 0, 0, 0,
                             0, 0, 0, 1, 0, 0,
                             0, 0, 0, 0, 1, 0,
                             0, 0, 0, 0, 0, 1};


    GaussNewton gn_solver(6, datas);
    gn_solver.init();
    gn_solver.fit(X);
    return 0;
}