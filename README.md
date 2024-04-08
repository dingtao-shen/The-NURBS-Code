# The NURBS Code

This Project is still under development. 

## Overview

- This tool kits are developed for utilizing 2D NURBS curves with double precision, and mainly involves the following functions:
  - Create 2D NURBS curves with input degree, knot vector, control points and associated weights;
  - Evaluate the points/derivatives on the curve at any parameter.
  - ...
- This  mainly refers to the algorithm in the book and refers to [tinyburns](https://github.com/pradeep-pyro/tinynurbs).
- This tool kits are developed utilizing the [Eigen](https://libeigen.gitlab.io/docs/) library.


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
Check the evaluation:
```
para      x                        y                        x^2+y^2                  
0.00      0.3000000000000000       0.0000000000000000       0.0900000000000000       
0.10      0.2441478108153225       -0.1743325743344757      0.0900000000000000       
0.20      0.0881435813134764       -0.2867589738320923      0.0900000000000000       
0.30      -0.0881435813134765      -0.2867589738320922      0.0900000000000000       
0.40      -0.2441478108153225      -0.1743325743344756      0.0900000000000000       
0.50      -0.3000000000000000      0.0000000000000000       0.0900000000000000       
0.60      -0.2441478108153224      0.1743325743344758       0.0900000000000000       
0.70      -0.0881435813134762      0.2867589738320923       0.0900000000000000       
0.80      0.0881435813134765       0.2867589738320922       0.0900000000000000       
0.90      0.2441478108153225       0.1743325743344756       0.0900000000000000       
1.00      0.3000000000000000       0.0000000000000000       0.0900000000000000     
```
