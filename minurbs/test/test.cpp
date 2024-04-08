#include<iostream>
#include <iomanip>
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

    // cout << "para  x  y  x^2+y^2" << endl;
    // for(int i = 0; i <= 10; i++){
    //     VectorXd p = curvePoint(crv, 0.1 * i);
    //     cout << 0.1*i << " " << p(0) << " " << p(1) << " " << pow(p(0), 2) + pow(p(1), 2) << endl;
    // }
    cout << left << setw(10) << "para" << setw(25) << "x" << setw(25) << "y" << setw(25) << "x^2+y^2" << endl;
    for(int i = 0; i <= 10; i++){
        VectorXd p = curvePoint(crv, 0.1 * i);
        cout << left << setw(10) << setprecision(2) << fixed << 0.1*i
             << setw(25) << setprecision(16) << p(0)
             << setw(25) << setprecision(16) << p(1)
             << setw(25) << setprecision(16) << pow(p(0), 2) + pow(p(1), 2) << endl;
    }

	return 0;
}