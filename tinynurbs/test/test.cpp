#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include<cassert>

#include "tinynurbs/tinynurbs.h"

using namespace std;
using namespace tinynurbs;

int main(){
    cout << "Pass" << endl;
    tinynurbs::Curve<float> crv; // Planar curve using float32
    crv.control_points = {glm::vec3(-1, 0, 0), // std::vector of 3D points
                        glm::vec3(0, 1, 0),
                        glm::vec3(1, 0, 0)
                        };
    crv.knots = {0, 0, 0, 1, 1, 1}; // std::vector of floats
    crv.degree = 2;

    if (tinynurbs::curveIsValid(crv)) {
        // check if degree, knots and control points are configured as per
        // #knots == #control points + degree + 1
        glm::vec3 pt = tinynurbs::curvePoint(crv, 0.f);
        cout << pt[0] << endl;
        // Outputs a point [-1, 0]
        glm::vec3 tgt = tinynurbs::curveTangent(crv, 0.5f);
        // Outputs a vector [1, 0]
    }

    return 0;
}