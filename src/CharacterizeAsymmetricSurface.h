#ifndef CHARACTERIZE_ASYMMETRIC_SURFACE_H
#define CHARACTERIZE_ASYMMETRIC_SURFACE_H

#include <iostream>
#include <cmath>

class CharacterizeAsymmetricSurface {
private:
    double c1, c2; // Coefficients shared for both X-axis and Y-axis equations

    double z_x(double x, double y) const; // z for X-axis equation
    double z_y(double x, double y) const; // z for Y-axis equation
    double dz_dx(double x) const;
    double dz_dy(double y) const;
    double integrand(double var, double fixed, bool along_x) const;
    double simpson_integration(double a, double b, double fixed, bool along_x, int n) const;

public:
    // Constructor
    CharacterizeAsymmetricSurface(double c1, double c2);

    // Public methods
    void calculate_arc_length(double start, double end, double fixed, bool along_x, double& arc_length, int num_intervals = 1000) const;
    void print_curve_points(double start, double end, double fixed, bool along_x, int num_points = 10) const;

    // Static method to calculate the central angle
    static double calculateCentralAngle(double curveLength, double radius);
};

#endif 
