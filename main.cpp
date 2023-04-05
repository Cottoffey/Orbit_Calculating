#include "creatingModelValues.h"

int main()
{
    // start parameters x,y,z,vx,vy,vz
    std::vector<double> X = {1.469591208242925E+08,  7.299762167917201E+07,  2.056299266163284E+07,  3.859428549646102E+06,  3.244525935598258E+05,  1.492020244998816E+06};

    modeling (X, 6, 0.041666666666666667);
    std::cout << "Modeling success\n";

    creatingModelingValues ();
    std::cout << "Creating model values success\n";

    return 0;
}