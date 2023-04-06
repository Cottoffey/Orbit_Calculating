#pragma once
#include <sofa.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <iomanip>

#include "Ephemeris.h"

#define LSD 25902068371.2000 // light speed in km/day
#define LS 299792.458        // light speed in km/s
#define KM_TO_AU 6.68459e-9  // 1 km = .. au

void RK4 (int n, std::vector<double> & x, double & t, double h, PlanetEphemeris *data, int m, void f(int, std::vector<double> &, const double &, double *, PlanetEphemeris *, int));
 
void DP5(int n, std::vector<double> & x, double &t, double h, PlanetEphemeris *data, int m, void f(int, std::vector<double> &, const double &, double *, PlanetEphemeris *, int));

void function(int n, std::vector<double> & X, const double &t, double *result, PlanetEphemeris *data, int m);

void modeling (std::vector<double> X, int n, double h);

double LTC (double time, double * obcoors, PlanetEphemeris & object);

void creatingModelingValues ();

