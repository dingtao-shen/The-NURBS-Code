#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include<cassert>

#include "tinynurbs/tinynurbs.h"

using namespace std;
using namespace tinynurbs;

int main(){
    cout.precision(16);

    tinynurbs::Curve<float> crv; // Planar curve using float32
    crv.control_points = {glm::vec3(-1.00200000, 0, 0), // std::vector of 3D points
                        glm::vec3(0, 1, 0),
                        glm::vec3(1, 0, 0)
                        };
    crv.knots = {0, 0, 0, 1, 1, 1}; // std::vector of floats
    crv.degree = 2;

    // tinynurbs::RationalCurve3d crv;
    // crv.control_points = {glm::vec<3,double>(0.3, 0.0, 0.0),
    //                     glm::vec<3,double>(0.3, -0.3, 0.0),
    //                     glm::vec<3,double>(0.0, -0.3, 0.0),
    //                     glm::vec<3,double>(-0.3, -0.3, 0.0),
    //                     glm::vec<3,double>(-0.3, 0.0, 0.0),
    //                     glm::vec<3,double>(-0.3, 0.3, 0.0),
    //                     glm::vec<3,double>(0.0, 0.3, 0.0),
    //                     glm::vec<3,double>(0.3, 0.3, 0.0),
    //                     glm::vec<3,double>(0.3, 0.0, 0.0)
    //                     };
    // crv.weights = {1.0, sqrt(2.0) / 2.0, 1.0, sqrt(2.0) / 2.0, 1.0, sqrt(2.0) / 2.0, 1.0, sqrt(2.0) / 2.0, 1.0};
    // // crv.weights = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    // crv.knots = {0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0}; 
    // crv.degree = 2;

    if (tinynurbs::curveIsValid(crv)) {
        // check if degree, knots and control points are configured as per
        // #knots == #control points + degree + 1
        glm::vec3 pt = tinynurbs::curvePoint(crv, 0.f);
        cout << pt[0] << endl;
    }

    glm::vec<3,double> u = {0.3, 0.0, 0.0};
    cout << u[0] << endl;
    
    return 0;
}