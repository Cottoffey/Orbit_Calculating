#pragma once

#include "Cholesky.h"
#include "../creatingModelValues.h"


class GaussNewton {
    int p_num; // number of parameters
    vec y; // observations
    vec times;
    vec params; // parameters
    vec residual;
    vec res_vec; // At * residual
    double delta_ra_sum{};
    double delta_dec_sum{};
    Matrix grad_f; // At * A
    PlanetEphemeris *datas;
public:

    GaussNewton(int p_num, PlanetEphemeris *datas);

    vec fit(const vec &init_state);

    static void dgdx(double x, double y, double z, Matrix &derivative);

    void init();

    void matrices_calculation();

    void clear_matrices();
};


