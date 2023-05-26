#pragma once
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <iomanip>

#include "Ephemeris.h"
#include "sofa/c/src/sofa.h"

void RK4 (int n, std::vector<double> & x, double & t, double h, PlanetEphemeris *data, int m, void f(int, std::vector<double> &, const double &, double *, PlanetEphemeris *, int));
 
void DP5(int n, std::vector<double> & x, double &t, double h, PlanetEphemeris *data, int m, void f(int, std::vector<double> &, const double &, double *, PlanetEphemeris *, int));

void function(int n, std::vector<double> & X, const double &t, double *result, PlanetEphemeris *data, int m);

void modeling (std::vector<double> X, int n, double h, PlanetEphemeris * data, int m);

double LTC (double time, double * obcoors, DataEphemeris & object);

void creatingModelingValues ();

