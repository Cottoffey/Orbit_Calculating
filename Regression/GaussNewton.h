#pragma once

#include "Cholesky.h"
#include "../creatingModelValues.h"


void dgdx(double x, double y, double z, matrix &derivative);

void drdb(const vec &X, matrix &g_deriv, matrix &residual_deriv);

class GaussNewton {
    int p_num; // number of parameters
    vec y; // observations
    vec times; // times
    vec params; // parameters
    vec residual;
    vec result_column; // jacobianT * residual
    matrix grad; // jacobianT * jacobian
    PlanetEphemeris *planets_data;
public:

    GaussNewton(const vec &y, vec &times, int p_num, PlanetEphemeris *datas);

    vec fit(const vec &init_state, double h);

    void residual_calculation(double x, double y, double z, int index);

    void integration(int n, int dim, vec &X, double h);

    void clear_matrices();
};


