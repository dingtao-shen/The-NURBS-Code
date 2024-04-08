# The NURBS Code

This Project is still under development. 

## Overview

- This tool kits are developed for utilizing 2D NURBS curves with double precision, and mainly involves the following functions:
  - Create 2D NURBS curves with input degree, knot vector, control points and associated weights;
  - Evaluate the points/derivatives on the curve at any parameter.
  - ...
- This  mainly refers to the algorithm in the book and refers to [tinyburns](https://github.com/pradeep-pyro/tinynurbs).
- This tool kits are developed utilizing the [Eigen] library.


## Example

Create a circle origined at $(0,0)$, with the radius of $0.3$, with following input:
```
MatrixXd control_points(9,2);
VectorXd knots(12);
VectorXd weights(9);
control_points << 0.3, 0.0,
                  0.3, -0.3,
                  0.0, -0.3,
                  -0.3, -0.3,
                  -0.3, 0.0,
                  -0.3, 0.3,
                  0.0, 0.3,
                  0.3, 0.3,
                  0.3, 0;
knots << 0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0;
weights << 1.0, sqrt(2.0) / 2.0, 1.0, sqrt(2.0) / 2.0, 1.0, sqrt(2.0) / 2.0, 1.0, sqrt(2.0) / 2.0, 1.0;
```