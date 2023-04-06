#pragma once
#include <sofa.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <iomanip>

#include "Ephemeris.h"

#define GMS  132712440043.85333 // Sun
#define GMJ  126712764.13345    // Jupiter
#define GME  398600.43552       // Earth
#define GMV  324858.59200       // Venus
#define GMST 37940585.20000     // Saturn
#define GMU  5794556.46575      // Uranus
#define GMN  6836527.10058      // Neptune
#define GMMS 42828.37521        // Mars
#define GMMC 22031.78000        // Mercury
#define GMMN 4902.80008         // Moon

#define LSD 25902068371.2000 // light speed in km/day
#define LS 299792.458        // light speed in km/s
#define KM_TO_AU 6.68459e-9  // 1 km = .. au

void RK4 (int n, std::vector<double> & x, double & t, double h, PlanetEphemeris *data, int m, void f(int, std::vector<double> &, const double &, double *, PlanetEphemeris *, int));
 
void DP5(int n, std::vector<double> & x, double &t, double h, PlanetEphemeris *data, int m, void f(int, std::vector<double> &, const double &, double *, PlanetEphemeris *, int));

void function(int n, std::vector<double> & X, const double &t, double *result, PlanetEphemeris *data, int m);

void modeling (std::vector<double> X, int n, double h, PlanetEphemeris * data, int m);

double LTC (double time, double * obcoors, PlanetEphemeris & object);

void creatingModelingValues ();

