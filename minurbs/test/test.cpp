#include<iostream>
#include <vector>
#include <cmath>
using namespace std;

#include "../include/Eigen/Dense"
#include "../include/core/curve.hpp"
#include "../include/core/evaluate.hpp"
using namespace Eigen;

int main(){
    cout.precision(16);

    MatrixXd CP(9,2);
    CP <<   0.3, 0.0,
            0.3, -0.3,
            0.0, -0.3,
            -0.3, -0.3,
            -0.3, 0.0,
            -0.3, 0.3,
            0.0, 0.3,
            0.3, 0.3,
            0.3, 0;
    VectorXd K(12);
    K << 0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0;
    VectorXd W(9);
    W << 1.0, sqrt(2.0) / 2.0, 1.0, sqrt(2.0) / 2.0, 1.0, sqrt(2.0) / 2.0, 1.0, sqrt(2.0) / 2.0, 1.0;

    minurbs::RationalCurve crv(2, K, CP, W);

    VectorXd p = curvePoint(crv, 0.75);
    cout << p.transpose() << endl;
    cout << pow(p(0), 2) + pow(p(1), 2) << endl;

	return 0;
}