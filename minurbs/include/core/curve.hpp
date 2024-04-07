/**
 * The Curve class represents a non-uniform polynomial B-spline curve, while the
 * RationalCurve class represents a non-uniform rational B-spline (NURBS) curve.
 *
 */

#ifndef MINURBS_CURVE
#define MINURBS_CURVE

#include "../Eigen/Dense"
#include <assert.h>
#include <cmath>

namespace minurbs{

    // Forward declaration
    struct RationalCurve;

    /**
    Struct for holding a polynomial B-spline curve
    @tparam T Data type of control points and knots (float or double)
    */
    struct Curve
    {
        int degree;
        Eigen::VectorXd knots;
        Eigen::MatrixXd control_points; // Np rows by 2 col

        Curve() = default;
        // Curve(RationalCurve crv) : Curve(crv.degree, crv.knots, crv.control_points) {}
        Curve(int degree, Eigen::VectorXd knots, Eigen::MatrixXd control_points)
            : degree(degree), knots(knots), control_points(control_points)
        {
        }
    };

    /**
    Struct for holding a rational B-spline curve
    @tparam T Data type of control points and knots (float or double)
    */
    struct RationalCurve
    {
        int degree;
        Eigen::VectorXd knots;
        Eigen::MatrixXd control_points;
        Eigen::VectorXd weights;

        RationalCurve() = default;
        RationalCurve(Curve crv): RationalCurve(crv, Eigen::VectorXd::Ones(crv.control_points.size()))
        {
        }
        RationalCurve(Curve crv, Eigen::VectorXd weights)
            : RationalCurve(crv.degree, crv.knots, crv.control_points, weights)
        {
        }
        RationalCurve(int degree, Eigen::VectorXd knots, Eigen::MatrixXd control_points, Eigen::VectorXd weights)
            : degree(degree), knots(knots), control_points(control_points), weights(weights)
        {
        }
    };


} // namespace minurbs


#endif