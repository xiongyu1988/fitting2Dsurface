#include "CharactDoublySurface.h"
#define _USE_MATH_DEFINES
#include <math.h> 
#include <cmath>
#include <iostream>
#include <iomanip>

// Constructor implementation
CharactDoublySurface::CharactDoublySurface(double c1, double c2, double c3, double c4)
    : c1(c1), c2(c2), c3(c3), c4(c4) {}

// Function to calculate z given x and y
double CharactDoublySurface::z(double x, double y) const {
    return c1 * x * x + c2 * x + c3 * y + c4 * y * y;
}

// Function to calculate dz/dx
double CharactDoublySurface::dz_dx(double x, double y) const {
    return 2 * c1 * x + c2;
}

// Function to calculate dz/dy
double CharactDoublySurface::dz_dy(double x, double y) const {
    return c3 + 2 * c4 * y;
}

// Integrand function for arc length along x
double CharactDoublySurface::integrand_x(double x, double y) const {
    return std::sqrt(1 + std::pow(dz_dx(x, y), 2));
}

// Integrand function for arc length along y
double CharactDoublySurface::integrand_y(double y, double x) const {
    return std::sqrt(1 + std::pow(dz_dy(x, y), 2));
}

// Simpson's rule for numerical integration
double CharactDoublySurface::simpson_integration(double a, double b, double fixed, bool along_x, int n) const {
    double h = (b - a) / n;
    double sum = (along_x ? integrand_x(a, fixed) + integrand_x(b, fixed)
        : integrand_y(a, fixed) + integrand_y(b, fixed));

    for (int i = 1; i < n; i += 2) {
        double var = a + i * h;
        sum += 4 * (along_x ? integrand_x(var, fixed) : integrand_y(var, fixed));
    }
    for (int i = 2; i < n - 1; i += 2) {
        double var = a + i * h;
        sum += 2 * (along_x ? integrand_x(var, fixed) : integrand_y(var, fixed));
    }
    return sum * h / 3;
}

// Function to calculate arc length
void CharactDoublySurface::calculate_arc_length(double start, double end, double fixed, bool along_x, double& arc_length, int num_intervals) const {
    arc_length = simpson_integration(start, end, fixed, along_x, num_intervals);
}

//
//
//
// Method to compute the second derivatives at the center (origin)
void CharactDoublySurface::firstDerivatives(double x, double y, double& fx, double& fy) {
    fx = 2 * c1 * x + c2;
    fy = 2 * c4 * y + c3;
}

// Method to compute second derivatives (constant for this surface)
void CharactDoublySurface::secondDerivatives(double& fxx, double& fyy, double& fxy) {
    fxx = 2 * c1;
    fyy = 2 * c4;
    fxy = 0; // The cross derivative is zero for this surface
}

// Method to compute the radii and principal curvature lengths
void CharactDoublySurface::computeRadii(double x, double y, double& R1, double& R2) {
    double fx, fy, fxx, fyy, fxy;
    firstDerivatives(x, y, fx, fy);
    secondDerivatives(fxx, fyy, fxy);

    // Radii of curvature (principal curvatures are inverse of these)
    R1 = std::abs((fxx != 0) ? 1 / fxx : std::numeric_limits<double>::infinity());
    R2 = std::abs((fyy != 0) ? 1 / fyy : std::numeric_limits<double>::infinity());
}


// Definition of the static method calculateCentralAngle
double CharactDoublySurface::calculateCentralAngle(double curveLength, double radius) {
    double angleInRadians = curveLength / radius;
    double angleInDegrees = angleInRadians * (180.0 / M_PI);
    return angleInDegrees;
}