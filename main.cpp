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


    jupiter.init("Data/Jupiter.txt", 0.0013888889);
    eartht.init("Data/Earth.txt", 0.0013888889);
    venus.init("Data/Venus.txt", 0.0013888889);
    uranus.init("Data/Uranus.txt", 0.0013888889);
    neptune.init("Data/Neptune.txt", 0.0013888889);
    saturn.init("Data/Saturn.txt", 0.0013888889);
    sun.init("Data/Sun.txt", 0.0013888889);
    mars.init("Data/Mars.txt", 0.0013888889);
    mercury.init("Data/Mercury.txt", 0.0013888889);
    moon.init("Data/Moon.txt", 0.0013888889);

    // RealOrbit.init ("Data/RealOrbit.txt", 0.0013888889);

    std::cout << "Initialization success\n";

    PlanetEphemeris datas[10] = {sun, jupiter, eartht, venus, uranus, neptune, saturn, mars, mercury, moon};
    // start parameters x,y,z,vx,vy,vz
    std::vector<double> X = {1.469591208242925E+08,  7.299762167917201E+07,  2.056299266163284E+07,  3.859428549646102E+06,  3.244525935598258E+05,  1.492020244998816E+06,
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