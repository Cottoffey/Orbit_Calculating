#include "GaussNewton.h"

GaussNewton::GaussNewton(const vec &y, int p_num, PlanetEphemeris *datas) : y(y), p_num(p_num), params(p_num),
                                                                            residual(y.size()),
                                                                            result_column(p_num),
                                                                            grad(p_num),
                                                                            planets_data(datas) {
    for (auto &row: grad)
        row.resize(p_num);
}


void GaussNewton::residual_calculation(double x_coord, double y_coord, double z_coord, int index) {
    double CartesianCoord[3] = {x_coord, y_coord, z_coord};
    double ra, dec;
    iauC2s(CartesianCoord, &ra, &dec);
    residual[index * 2] = y[index * 2] - ra;
    residual[index * 2 + 1] = y[index * 2 + 1] - dec;
}


void GaussNewton::integration(int n, int dim, double h, vec &X) {
    double t = 2458040.937500000;
    matrix r_deriv;    // dr / db
    matrix g_deriv(2); // dg / dx
    g_deriv[0].resize(3);
    g_deriv[1].resize(3);

    for (int i = 0; i < n; ++i) {
        residual_calculation(X[0], X[1], X[2], i);

        // calculation dr/db
        drdb(X, g_deriv, r_deriv);

        for (int j = 0; j < p_num; ++j) {
            // JacobianT * residual
            result_column[j] += r_deriv[0][j] * residual[i * 2] + r_deriv[1][j] * residual[i * 2 + 1];
            // JacobianT * Jacobian
            for (int k = 0; k < p_num; ++k) {
                grad[j][k] += r_deriv[0][j] * r_deriv[0][k] + r_deriv[1][j] + r_deriv[1][k];
            }
        }

        DP5(dim, X, t, h, planets_data, 10, function);
        t += h;
    }


}

vec GaussNewton::fit(const vec &init_state, const double h) {
    // params initialising
    for (int i = 0; i < p_num; ++i)
        params[i] = init_state[i];

    for (int i = 0; i < 2; ++i) {
        // update X
        vec X = init_state;
        clear_matrices();
        for (int j = 0; j < p_num; ++j)
            X[j] = params[j];

        integration(y.size() / 2, X.size(), h, X);
        vec sol = Cholesky::CholeskySLE(grad, result_column);
        printMatrix(grad);
        for (int j = 0; j < p_num; ++j) {
            params[j] = sol[j] + params[j];
            std::cout << params[j] << ' ';
        }
        std::cout << '\n';
    }
    return params;
}

void GaussNewton::clear_matrices() {
    for (int i = 0; i < p_num; ++i) {
        for (int j = 0; j < p_num; ++j) {
            grad[i][j] = 0;
        }
        result_column[i] = 0;
    }
}

void dgdx(double x, double y, double z, matrix &derivative) {
    derivative[0][0] = -y / (x * x + y * y);
    derivative[0][1] = x / (x * x + y * y);
    derivative[0][2] = 0;
    derivative[1][0] = -x * z / ((x * x + y * y + z * z) * sqrt(x * x + y * y));
    derivative[1][1] = -y * z / ((x * x + y * y + z * z) * sqrt(x * x + y * y));
    derivative[1][2] = (x * x + y * y) / ((x * x + y * y + z * z) * sqrt(x * x + y * y));
}

void drdb(const vec &X, matrix &g_deriv, matrix &residual_deriv) {
    dgdx(X[0], X[1], X[2], g_deriv);
    matrix dxdp(3); // d(x,y,z) / db
    for (auto &row: dxdp)
        row.resize(6);
    for (int j = 0; j < 18; ++j)
        dxdp[j / 6][j % 6] = -X[j + 6];

    // dr / db
    residual_deriv = multMatrix(g_deriv, dxdp);
}