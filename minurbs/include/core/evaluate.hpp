/**
 * Core functionality for evaluating points, derivatives and related
 * quantities on NURBS curves and surfaces.
 *
 */

#ifndef TINYNURBS_EVALUATE
#define TINYNURBS_EVALUATE

#include <iostream>
using namespace std;
#include "../Eigen/Dense"
#include "basis.hpp"
#include "curve.hpp"

namespace minurbs{

    namespace internal{

        /**
        * Evaluate point on a nonrational NURBS curve
        * @param[in] deg Degree of the given curve.
        * @param[in] knots Knot vector of the curve.
        * @param[in] control_points Control points of the curve.
        * @param[in] u Parameter to evaluate the curve at.
        * @return point Resulting point on the curve at parameter u.
        */
        Eigen::VectorXd curvePoint(int deg, Eigen::VectorXd knots, Eigen::MatrixXd control_points, double u)
        {
            // Initialize result to 0s
            Eigen::VectorXd point = Eigen::VectorXd::Zero(control_points.cols());

            // Find span and corresponding non-zero basis functions
            int span = findSpan(deg, knots, u);
            Eigen::VectorXd N = bsplineBasis(deg, span, knots, u);

            // Compute point
            for (int j = 0; j <= deg; j++)
            {
                for(int k = 0; k < point.size(); k++){
                    point(k) += double(N[j]) * control_points(span - deg + j, k);
                }
            }
            return point;
        }

    } //namespace internal

    /**
    Evaluate point on a nonrational NURBS curve
    @param[in] crv Curve object
    @param[in] u Parameter to evaluate the curve at.
    @return point Resulting point on the curve at parameter u.
    */
    Eigen::VectorXd curvePoint(Curve crv, double u)
    {
        return internal::curvePoint(crv.degree, crv.knots, crv.control_points, u);
    }

    /**
     * Evaluate point on a rational NURBS curve
     * @param[in] crv RationalCurve object
     * @param[in] u Parameter to evaluate the curve at.
     * @return point Resulting point on the curve.
     */
    Eigen::VectorXd curvePoint(RationalCurve crv, double u)
    {
        // Convert cartesian coordinates to homogenous coordinates
        Eigen::MatrixXd Cw = Eigen::MatrixXd::Zero(crv.control_points.rows(), 3);
        for(int i = 0; i < crv.control_points.rows(); i++){
            Cw(i,0) = crv.control_points(i,0) * crv.weights(i);
            Cw(i,1) = crv.control_points(i,1) * crv.weights(i);
            Cw(i,2) = crv.weights(i);
        }
        // Compute point using homogenous coordinates
        Eigen::VectorXd pointw = internal::curvePoint(crv.degree, crv.knots, Cw, u);
        // Convert back to cartesian coordinates
        Eigen::VectorXd point(pointw.size()-1);
        double w = pointw(pointw.size()-1);
        for(int i = 0; i < point.size(); i++){
            point(i) = pointw(i) / w;
        }

        return point;
    }

} // namespace minurbs

#endif