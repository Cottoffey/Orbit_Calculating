#include "creatingModelValues.h"

void RK4(int n, std::vector<double> &x, double &t, double h, PlanetEphemeris *data, int m, void f(int, std::vector<double> &, const double &, double *, PlanetEphemeris *, int))
{
    std::vector<double> xmid(n);
    double tmid;
    auto *k = new double[4 * n];

    f(n, x, t, &k[0], data, m);

    for (int i = 0; i < n; i++)
        xmid[i] = x[i] + h * k[0 + i] / 2.0;
    tmid = t + h / 2.0;
    f(n, xmid, tmid, &k[1 * n], data, m);

    for (int i = 0; i < n; i++)
        xmid[i] = x[i] + h * k[1 * n + i] / 2.0;
    tmid = t + h / 2.0;
    f(n, xmid, tmid, &k[2 * n], data, m);

    for (int i = 0; i < n; i++)
        xmid[i] = x[i] + h * k[2 * n + i];
    tmid = t + h;
    f(n, xmid, tmid, &k[3 * n], data, m);

    for (int i = 0; i < n; i++)
        x[i] = x[i] + h * (k[0 + i] / 6 + k[n + i] / 3 + k[2 * n + i] / 3 + k[3 * n + i] / 6);

    delete[] k;
}

void DP5(int n, std::vector<double> &x, double &t, double h, PlanetEphemeris *data, int m, void f(int, std::vector<double> &, const double &, double *, PlanetEphemeris *, int))
{
    std::vector<double> xmid(n);
    double tmid;
    auto *k = new double[6 * n];

    f(n, x, t, &k[0], data, m);

    for (int i = 0; i < n; i++)
        xmid[i] = x[i] + h * k[0 + i] / 5.0;
    tmid = t + h / 5.0;
    f(n, xmid, tmid, &k[1 * n], data, m);

    for (int i = 0; i < n; i++)
        xmid[i] = x[i] + h * (k[0 + i] * 3.0 / 40.0 + k[1 * n + i] * 9.0 / 40.0);
    tmid = t + h * 3.0 / 10.0;
    f(n, xmid, tmid, &k[2 * n], data, m);

    for (int i = 0; i < n; i++)
        xmid[i] = x[i] + h * (k[0 + i] * 44.0 / 45.0 - k[1 * n + i] * 56.0 / 15.0 + k[2 * n + i] * 32.0 / 9.0);
    tmid = t + h * 4.0 / 5.0;
    f(n, xmid, tmid, &k[3 * n], data, m);

    for (int i = 0; i < n; i++)
        xmid[i] = x[i] + h * (k[0 + i] * 19372.0 / 6561.0 - k[1 * n + i] * 25360.0 / 2187.0 + k[2 * n + i] * 64448.0 / 6561.0 - k[3 * n + i] * 212.0 / 729.0);
    tmid = t + h * 8.0 / 9.0;
    f(n, xmid, tmid, &k[4 * n], data, m);

    for (int i = 0; i < n; i++)
        xmid[i] = x[i] + h * (k[0 + i] * 9017.0 / 3168.0 - k[1 * n + i] * 355.0 / 33.0 + k[2 * n + i] * 46732.0 / 5247.0 + k[3 * n + i] * 49.0 / 176.0 - k[4 * n + i] * 5103.0 / 18656.0);
    tmid = t + h;
    f(n, xmid, tmid, &k[5 * n], data, m);

    for (int i = 0; i < n; i++)
        x[i] = x[i] + h * (k[0 + i] * 35.0 / 384.0 + k[2 * n + i] * 500.0 / 1113.0 + k[3 * n + i] * 125.0 / 192.0 - k[4 * n + i] * 2187.0 / 6784.0 + k[5 * n + i] * 11.0 / 84.0);

    delete[] k;
}

void function(int n, std::vector<double> &X, const double &t, double *result, PlanetEphemeris *data, int m)
{
    result[0] = X[3];
    result[1] = X[4];
    result[2] = X[5];
    result[3] = 0.0;
    result[4] = 0.0;
    result[5] = 0.0;

    double coors[3], R, RC, cc[3];

    double tmp1, tmp2, v[3], a[3], u[3];

    // dF/dX
    double F[36] = {0, 0, 0, 1, 0, 0,
                    0, 0, 0, 0, 1, 0,
                    0, 0, 0, 0, 0, 1,
                    0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0};

    for (int i = 0; i < m; i++)
    {
        data[i].get_coors(t, coors[0], coors[1], coors[2]);
        R = sqrt((X[0] - coors[0]) * (X[0] - coors[0]) + (X[1] - coors[1]) * (X[1] - coors[1]) + (X[2] - coors[2]) * (X[2] - coors[2]));

        // dF/dX
        for (int j = 3; j < 6; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                if (j - 3 == k)
                    F[j * 6 + k] -= (86400.L * 86400.L * data[i].GM  * (1./ (R * R * R) - 3 * (coors[k] - X[k]) * (coors[k] - X[k]) / (R * R * R * R * R)));
                else
                    F[j * 6 + k] += (86400.L * 86400.L * data[i].GM * 3 * (coors[j - 3] - X[j - 3]) * (coors[k] - X[k]) / (R * R * R * R * R));
            }
        }

        

        data[i].get_speed(t, v[0], v[1], v[2]);
        data[i].get_acceleration(t, a[0], a[1], a[2]);

        // unit vector AB
        for (int j = 0; j < 3; j++)
            u[j] = (X[j] - coors[j]) / R;

        // first summation
        tmp1 = X[3] * X[3] + X[4] * X[4] + X[5] * X[5] +
               +2 * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) - 4 * (X[3] * v[0] + X[4] * v[1] + X[5] * v[2]) - 1.5 * (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) * (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) - 0.5 * (u[0] * a[0] + u[1] * a[1] + u[2] * a[2]) * R;

        for (int j = 0; j < m; j++)
        {
            data[j].get_coors(t, cc[0], cc[1], cc[2]);
            RC = sqrt((X[0] - cc[0]) * (X[0] - cc[0]) + (X[1] - cc[1]) * (X[1] - cc[1]) + (X[2] - cc[2]) * (X[2] - cc[2]));
            tmp1 = tmp1 - 4 * data[j].GM * 86400.L * 86400.L / RC;

            if (j != i)
            {
                RC = sqrt((coors[0] - cc[0]) * (coors[0] - cc[0]) + (coors[1] - cc[1]) * (coors[1] - cc[1]) + (coors[2] - cc[2]) * (coors[2] - cc[2]));
                tmp1 = tmp1 - data[j].GM * 86400.L * 86400.L / RC;
            }
        }

        // second summation
        tmp2 = u[0] * (4 * X[3] - 3 * v[0]) + u[1] * (4 * X[4] - 3 * v[1]) + u[2] * (4 * X[5] - 3 * v[2]);

        for (int j = 3; j < 6; j++)
            result[j] = result[j] - 86400.L * 86400.L * data[i].GM * u[j - 3] / (R * R) - 86400.L * 86400.L * data[i].GM * u[j - 3] * tmp1 / (R * R * LSD * LSD) + 86400.L * 86400.L * data[i].GM * tmp2 * (X[j] - v[j - 3]) / (R * R * LSD * LSD) + 7 * 86400.L * 86400.L * data[i].GM * a[j - 3] / (R * 2 * LSD * LSD);
    }
    
    for (int j = 0; j < 6; j++)
    {
        for (int h = 0; h < 6; h++)
        {
            result[6 + j * 6 + h] = 0;
            for (int s = 0; s < 6; s++)
            {
                result[6 + j * 6 + h] += (F[j * 6 + s] * X[6 + s * 6 + h]);
            }
            // std::cout << F[j*6 + h] << ' ';
        }
        // std::cout << std::endl;
    }
    // std::cout << std::endl;

}

void modeling(std::vector<double> X, int n, double h, PlanetEphemeris *data, int m)
{
    std::cout.setf(std::ios::scientific);

    // start time in Julian format
    double t = 2458040.937500000;

    std::ofstream output("Data/ModelOrbit.txt");
    output.setf(std::ios::scientific);

    double x, y, z;

    while (t < 2458123.916666667)
    {
        // RealOrbit.get_coors (t, x, y, z);
        // output << std::setprecision(15) << t << ' ' << X[0] << ' ' << X[1] << ' ' << X[2] << ' ' << X[0] - x << ' ' << X[1] - y << ' ' <<  X[2] - z <<  std::endl;
        // output << std::setprecision(15) << t << ' ' << X[0] << ' ' << X[1] << ' ' << X[2] << ' ' << X[3] << ' ' << X[4] << ' ' << X[5] << std::endl;
        output << std::setprecision (15) << t << ' ';
        for (int i = 0; i < 42; i++)
            output << std::setprecision (15) << X[i] << ' ';

        output << std::endl;
        DP5(n, X, t, h, data, m, function);

        t += h;
    }

    output.close();
}

// light speed correction
double LTC(double time, double *obcoors, DataEphemeris &object)
{
    double X[3];
    double delta = 0;
    double prevdelta = 0, prevprevdelta = 0;

    double R;
    object.get_coors(time - delta, X[0], X[1], X[2]);
    R = sqrt((X[0] - obcoors[0]) * (X[0] - obcoors[0]) + (X[1] - obcoors[1]) * (X[1] - obcoors[1]) + (X[2] - obcoors[2]) * (X[2] - obcoors[2]));
    delta = R / LSD;

    // Fixed-point method
    int i = 0;
    while (delta != prevdelta && delta != prevprevdelta)
    {
        prevprevdelta = prevdelta;
        prevdelta = delta;
        object.get_coors(time - delta, X[0], X[1], X[2]);
        R = sqrt((X[0] - obcoors[0]) * (X[0] - obcoors[0]) + (X[1] - obcoors[1]) * (X[1] - obcoors[1]) + (X[2] - obcoors[2]) * (X[2] - obcoors[2]));
        delta = R / LSD;
        i++;
    }

    
    // std::cout << i << ' ' << delta << std::endl;

    return delta;
}

void creatingModelingValues()
{
    DataEphemeris object;
    PlanetEphemeris sun;
    PlanetEphemeris earth;

    earth.init("Data/Earth.txt", 0.0013888889);
    object.init("Data/ModelOrbit.txt", 0.0013888889);
    sun.init("Data/Sun.txt", 0.0013888889);

    std::cout << "Creating model values: Initialization success\n";

    std::ifstream input("Data/ObservData.txt");
    std::ofstream output("Data/ModelingData.txt");
    output.setf(std::ios::scientific);

    double time;
    double ocoors[3];  // object coors
    double obcoors[3]; // observer coors
    double scoors[3];  // sun coors
    double p[3], q[3], e[3], p1[3], lp, lq, le;
    double em;
    double tmp1, tmp;

    for (int j = 0; j < 222; j++)
    {
        input >> time >> tmp1 >> tmp >> obcoors[0] >> obcoors[1] >> obcoors[2];

        time -= LTC(time, obcoors, object);
        // earth.get_coors (time, obcoors[0], obcoors[1], obcoors[2]);
        object.get_coors(time, ocoors[0], ocoors[1], ocoors[2]);
        sun.get_coors(time, scoors[0], scoors[1], scoors[2]);

        // initialization for calling iauLd
        lp = 0, lq = 0, le = 0;
        for (int i = 0; i < 3; i++)
        {
            p[i] = ocoors[i] - obcoors[i];
            q[i] = ocoors[i] - scoors[i];
            e[i] = obcoors[i] - scoors[i];
            lp += p[i] * p[i];
            lq += q[i] * q[i];
            le += e[i] * e[i];
        }

        // lp = sqrt(lp);
        // lq = sqrt(lq);
        // le = sqrt(le);

        // for (int i = 0; i < 3; i++)
        // {
        //     p[i] = p[i] / lp;
        //     q[i] = q[i] / lq;
        //     e[i] = e[i] / le;
        // }
        // em = sqrt((scoors[0] - obcoors[0]) * (scoors[0] - obcoors[0]) + (scoors[1] - obcoors[1]) * (scoors[1] - obcoors[1]) + (scoors[2] - obcoors[2]) * (scoors[2] - obcoors[2])) * KM_TO_AU;
        
        // iauLd(1, p, q, e, em, 0, p1);
        // // Aberation
        // double v[3], lv = 0;
        // earth.get_speed(time, v[0], v[1], v[2]);
        // for (int i = 0; i < 3; i++)
        // {
        //     v[i] = v[i] / LSD;
        //     lv += v[i] * v[i];
        // }

        // iauAb(p1, v, em, sqrt(1 - lv), p);
        
        // to spherical
        double ra, dec;
        iauC2s(p, &ra, &dec);
        if (ra < 0.0)
            ra = 2 * M_PI + ra;

        output << std::setprecision(15) << time << ' ' << ra << ' ' << dec << ' ' << ((fabs(tmp1 - ra) > M_PI) ? (fabs(tmp1 - ra) - M_PI) : fabs(tmp1 - ra)) << ' ' << ((fabs(tmp - dec) > M_PI) ? (fabs(tmp - dec) - M_PI) : fabs(tmp - dec)) << std::endl;
        // output << std::setprecision (15) << time << ' ' << ra << ' ' << dec << std::endl;
    }
    std::cout << "Creating model values: success\n";

    input.close();
    output.close();
}
