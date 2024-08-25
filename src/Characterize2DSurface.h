#ifndef CHARACTERIZE2DSURFACE_H
#define CHARACTERIZE2DSURFACE_H

#include <iostream>
#include <cmath>
#include <limits>

class Characterize2DSurface {
private:
    double c1, c2, c3, c4; // Coefficients for the surface equation

public:
    // Constructor to initialize the coefficients
    Characterize2DSurface(double c1, double c2, double c3, double c4);

	// Method to compute the first derivatives at a given point
    void firstDerivatives(double x, double y, double& fx, double& fy);

    // Method to compute the second derivatives at the center (origin)
    void secondDerivatives(double& fxx, double& fyy, double& fxy);

    // Method to compute the radii and principal curvature lengths
    void computeRadiiAndCurvatureLengths(double x, double y, double& R1, double& R2, double& L1, double& L2);

    // Function to calculate the central angle
    static double calculateCentralAngle(double curveLength, double radius);
};

#endif // CHARACTERIZE2DSURFACE_H
