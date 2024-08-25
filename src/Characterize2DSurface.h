#ifndef CHARACTERIZE2DSURFACE_H
#define CHARACTERIZE2DSURFACE_H

#include <iostream>
#include <cmath>
#include <limits>

class Characterize2DSurface {
private:
    double c1, c2, c3, c4; // Coefficients for the surface equation

    double z(double x, double y) const;
    double dz_dx(double x, double y) const;
    double dz_dy(double x, double y) const;
    double integrand_x(double x, double y) const;
    double integrand_y(double y, double x) const;
    double simpson_integration(double a, double b, double fixed, bool along_x, int n) const;

public:
    // Constructor
    Characterize2DSurface(double c1, double c2, double c3, double c4);

    // Public methods
    void calculate_arc_length(double start, double end, double fixed, bool along_x, double& arc_length, int num_intervals = 1000) const;

    // Method to compute the first derivatives at a given point
    void firstDerivatives(double x, double y, double& fx, double& fy);

    // Method to compute the second derivatives at the center (origin)
    void secondDerivatives(double& fxx, double& fyy, double& fxy);

    // Method to compute the radii and principal curvature lengths
    void computeRadii(double x, double y, double& R1, double& R2);

    // Static method to calculate the central angle
    static double calculateCentralAngle(double curveLength, double radius);
};

#endif // CHARACTERIZE2D2DSURFACE_H