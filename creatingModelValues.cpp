#include "creatingModelValues.h"

void RK4 (int n, std::vector<double> & x, double & t, double h, PlanetEphemeris *data, int m, void f(int, std::vector<double> &, const double &, double *, PlanetEphemeris *, int))  {
    std::vector<double> xmid (n);
    double tmid;
    double * k = new double [4 * n];
    
    f(n, x, t, &k[0], data, m);
    
    for (int i = 0; i < n; i++)
        xmid[i] = x[i] + h * k[0+i] / 2.0;
    tmid = t + h/2.0;
    f(n, xmid, tmid, &k[1 * n], data, m);
        
    for (int i = 0; i < n; i++)
        xmid[i] = x[i] + h * k[1*n + i] / 2.0;
    tmid = t + h/2.0;
    f(n, xmid, tmid, &k[2*n], data, m);

    for (int i = 0; i < n; i++) 
        xmid[i] = x[i] + h * k[2*n + i];
    tmid = t + h;
    f(n, xmid, tmid, &k[3*n], data, m);

    for (int i = 0; i < n; i++) 
        x[i] = x[i] + h*(k[0 + i]/6 + k[n + i]/3 + k[2*n + i]/3 + k[3*n +i]/6);
    
    delete[] k;
}

void DP5(int n, std::vector<double> & x, double &t, double h, PlanetEphemeris *data, int m, void f(int, std::vector<double> &, const double &, double *, PlanetEphemeris *, int))
{
    std::vector<double> xmid (n);
    double tmid;
    double *k = new double[6 * n];

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

void function(int n, std::vector<double> & X, const double &t, double *result, PlanetEphemeris *data, int m)
{
    result[0] = X[3];
    result[1] = X[4];
    result[2] = X[5];
    result[3] = 0.0;
    result[4] = 0.0;
    result[5] = 0.0;

    double x, y, z, R;

    for (int i = 0; i < m; i++)
    {
        data[i].get_coors(t, x, y, z);
        R = sqrt((X[0] - x) * (X[0] - x) + (X[1] - y) * (X[1] - y) + (X[2] - z) * (X[2] - z));

        result[3] = result[3] + 86400. * 86400. * data[i].GM * (x - X[0]) / (R * R * R);
        result[4] = result[4] + 86400. * 86400. * data[i].GM * (y - X[1]) / (R * R * R);
        result[5] = result[5] + 86400. * 86400. * data[i].GM * (z - X[2]) / (R * R * R);
    }
}

void modeling (std::vector<double> X, int n, double h, PlanetEphemeris * data, int m)
{
    PlanetEphemeris RealOrbit;

    RealOrbit.init ("Data/RealOrbit.txt", 0.041666666667);

    std::cout << "Initialization success\n";

    std::cout.setf(std::ios::scientific);
    
    // start time in Julian format
    double t = 2458040.937500000;

    std::ofstream output("Data/ModelOrbit.txt");
    output.setf(std::ios::scientific);

    double x, y, z;    

    while (t < 2458123.916666667)
    {
        RealOrbit.get_coors (t, x, y, z);
        output << std::setprecision(15) << t << ' ' << X[0] << ' ' << X[1] << ' ' << X[2] << ' ' << X[0] - x << ' ' << X[1] - y << ' ' << X[2] - z <<  std::endl;
        // output << std::setprecision(15) << t << ' ' << X[0] << ' ' << X[1] << ' ' << X[2] << std::endl;

        DP5(n, X, t, h, data, m, function);
        // RK4(6, X, t, h, datas, 10, function);

        t += h;
    }

    output.close();
}

// light speed correction
double LTC (double time, double * obcoors, PlanetEphemeris & object)
{
    double X[3];
    double delta = 0;
    double prevdelta = 0;

    double R;
    object.get_coors (time - delta, X[0], X[1], X[2]);
    R = sqrt((X[0] - obcoors[0]) * (X[0] - obcoors[0]) + (X[1] - obcoors[1]) * (X[1] - obcoors[1]) + (X[2] - obcoors[2]) * (X[2] - obcoors[2]));
    delta = R / LSD;

    // Fixed-point method 
    int i = 0;
    while (fabs (delta - prevdelta) / delta > 1e-15)
    {
        prevdelta = delta;
        object.get_coors (time - delta, X[0], X[1], X[2]);
        R = sqrt((X[0] - obcoors[0]) * (X[0] - obcoors[0]) + (X[1] - obcoors[1]) * (X[1] - obcoors[1]) + (X[2] - obcoors[2]) * (X[2] - obcoors[2]));
        delta = R / LSD;
        i++;
    }

    // std::cout << i << ' ' << delta << std::endl;

    return delta;
}

void creatingModelingValues ()
{
    PlanetEphemeris object;
    PlanetEphemeris sun;
    PlanetEphemeris earth;

    earth.init ("Data/Earth.txt", 0.041666666667);
    object.init ("Data/RealOrbit.txt", 0.041666666667);
    sun.init ("Data/Sun.txt", 0.041666666667);

    std::cout << "Iniatilization success\n";

    std::ifstream input ("Data/ObservData.txt");
    std::ofstream output ("Data/ModelingData.txt");
    output.setf (std::ios::scientific);

    double time;
    double ocoors[3]; // object coors
    double obcoors[3]; // observer coors
    double scoors[3];  // sun coors
    double p[3], q[3], e[3], p1[3], lp = 0, lq = 0, le = 0;
    double em;
    double tmp1, tmp;

    for (int j = 0; j < 222; j++)
    {
        input >> time >> tmp1 >> tmp >> obcoors[0] >> obcoors[1] >> obcoors[2];

        time -= LTC (time, obcoors, object);
        // earth.get_coors (time, obcoors[0], obcoors[1], obcoors[2]);
        object.get_coors (time, ocoors[0], ocoors[1], ocoors[2]);
        sun.get_coors (time, scoors[0], scoors[1], scoors[2]);

        // initialization for calling iauLd
        for (int i = 0; i < 3; i++)
        {
            p[i] = ocoors[i] - obcoors[i];
            q[i] = ocoors[i] - scoors[i];
            e[i] = obcoors[i] - scoors[i];
            lp += p[i] * p[i];
            lq += q[i] * q[i];
            le += e[i] * e[i];
        }

        lp = sqrt (lp);
        lq = sqrt (lq);
        le = sqrt (le);

        for (int i = 0; i < 3; i++)
        {
            p[i] = p[i] / lp;
            q[i] = q[i] / lq;
            e[i] = e[i] / le;
        }
        em = sqrt ((scoors[0] - obcoors[0]) * (scoors[0] - obcoors[0]) + (scoors[1] - obcoors[1]) * (scoors[1] - obcoors[1]) + (scoors[2] - obcoors[2]) * (scoors[2] - obcoors[2])) * KM_TO_AU;

        iauLd (1, p, q, e, em, 0, p1);

        // Aberation 
        double v[3], lv = 0;
        earth.get_speed (time, v[0], v[1], v[2]);
        for (int i = 0; i < 3; i++)
        {
            v[i] = v[i] / LSD;
            lv += v[i] * v[i];
        }

        iauAb (p1, v, em, sqrt (1 - lv), p);

        // to spherical
        double ra, dec;
        iauC2s (p, &ra, &dec);
        if (ra < 0.0)
            ra = 2 * M_PI + ra;

        // output << std::setprecision (15) << time << ' ' << ra << ' ' << dec << ' ' << ((fabs(tmp1 - ra) > M_PI) ? (fabs (tmp1 - ra) - M_PI) : fabs (tmp1 - ra)) << ' ' << ((fabs(tmp- dec) > M_PI) ? (fabs (tmp- dec) - M_PI) : fabs (tmp- dec)) <<  std::endl;
        output << std::setprecision (15) << time << ' ' << ra << ' ' << dec << std::endl;
    }

    input.close();
    output.close();
}
