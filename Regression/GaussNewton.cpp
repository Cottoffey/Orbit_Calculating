#include "GaussNewton.h"


GaussNewton::GaussNewton(int p_num, PlanetEphemeris *datas) : p_num(p_num), params(p_num), datas(datas), y(444),
                                                              times(222) {
    jacobian.resize(y.size());
    residual.resize(y.size());
    for (auto &row: jacobian)
        row.resize(p_num);


}

void GaussNewton::init() {
    std::cout.setf(std::ios::scientific);
    std::ifstream input_data("Data/ObservData.txt");
    for (int i = 0; i < 222; ++i) {
        char buffer[140];
        input_data >> times[i] >> y[i * 2] >> y[i * 2 + 1];
        input_data.getline(buffer, 140);
    }
    input_data.close();
}

void GaussNewton::matrices_calculation() {
    std::ifstream model_values("Data/ModelingData.txt");
    double x_coord, y_coord, z_coord;
    double tmp, ra, dec;
    // Derivative matrix initialising
    Matrix x_deriv(3);    // d(x, y, z) / db - 3 x p_num
    Matrix g_deriv(2);    // dg / d(x, y, z) - 2 x 3
    for (auto &row: x_deriv)
        row.resize(p_num);
    for (auto &row: g_deriv)
        row.resize(3);

    DataEphemeris object;
    object.init("Data/ModelOrbit.txt", 0.041666666667);

    for (int i = 0; i < 222; ++i) {
        object.get_coors(times[i], x_coord, y_coord, z_coord);

        // derivatives calculation
        dgdx(x_coord, y_coord, z_coord, g_deriv);
        for (int j = 0; j < 18; ++j)
            object.get_dxdx0(times[i], j, x_deriv[j / 6][j % 6]);

        Matrix drdb = multMatrix(g_deriv, x_deriv); // dr / db = dg/df * dx/db

        // Filling jacobian
        for (int j = 0; j < p_num; ++j) {
            jacobian[i * 2][j] = -drdb[0][j];     // d(ra) / db
            jacobian[i * 2 + 1][j] = -drdb[1][j]; // d(dec) / db
        }

        model_values >> tmp >> ra >> dec >> tmp >> tmp;
        residual[i * 2] = (y[i * 2] - ra);                // delta ra
        residual[i * 2 + 1] = (y[i * 2 + 1] - dec);       // delta dec
    }
    model_values.close();

}


vec GaussNewton::fit(const vec &init_state) {
    // params initialising
    std::cout << "Initial parameteres:" << std::endl;
    for (int i = 0; i < p_num; ++i) {
        params[i] = init_state[i];
        std::cout << std::setprecision(15) << params[i] << ' ';
    }
    for (int i = 0; i < 10; ++i) {
        // update X
        vec X = init_state;
        for (int j = 0; j < p_num; ++j)
            X[j] = params[j];

        clear_matrices();
        modeling(X, 42, 0.041666666667, datas, 10);
        creatingModelingValues();
        matrices_calculation();
        Matrix grad = multMatrix(transpose(jacobian), jacobian);
        vec f_b = multMatrix(transpose(jacobian), residual);

        vec sol = Cholesky::CholeskySLE(grad, f_b);

        std::cout << "Iteration #" << i + 1 << std::endl;
        std::cout << "Matrix At * A:" << std::endl;
        printMatrix(grad);

        std::cout << "At * r:" << std::endl;
        for (auto &elem: f_b)
            std::cout << elem << std::endl;

        std::cout << "New_params:" << std::endl;
        for (int j = 0; j < p_num; ++j) {
            params[j] -= sol[j];
            std::cout << params[j] << ' ';
        }
        std::cout << "\n\n";
    }

    return params;
}

void GaussNewton::clear_matrices() {
    for (int j = 0; j < y.size(); ++j) {
        for (int i = 0; i < p_num; ++i) {
            jacobian[j][i] = 0;
        }
    }
}

void GaussNewton::dgdx(double x, double y, double z, Matrix &derivative) {
    derivative.resize(2);
    derivative[0].resize(3);
    derivative[1].resize(3);
    derivative[0][0] = -y / (x * x + y * y);
    derivative[0][1] = x / (x * x + y * y);
    derivative[0][2] = 0;
    derivative[1][0] = -x * z / ((x * x + y * y + z * z) * sqrt(x * x + y * y));
    derivative[1][1] = -y * z / ((x * x + y * y + z * z) * sqrt(x * x + y * y));
    derivative[1][2] = (x * x + y * y) / ((x * x + y * y + z * z) * sqrt(x * x + y * y));
}
