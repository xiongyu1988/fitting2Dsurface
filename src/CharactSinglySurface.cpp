#include "CharactSinglySurface.h"
#define _USE_MATH_DEFINES
#include <math.h> 
#include <cmath>
#include <iostream>
#include <iomanip>

// Constructor implementation
CharactSinglySurface::CharactSinglySurface(double c1, double c2) : c1(c1), c2(c2) {}

// Function to calculate z for X-axis equation
double CharactSinglySurface::z_x(double x, double y) const {
    return c1 * x * x + c2 * y;
}

// Function to calculate z for Y-axis equation
double CharactSinglySurface::z_y(double x, double y) const {
    return c1 * y * y + c2 * x;
}

// Function to calculate dz/dx
double CharactSinglySurface::dz_dx(double x) const {
    return 2 * c1 * x;
}

// Function to calculate dz/dy
double CharactSinglySurface::dz_dy(double y) const {
    return 2 * c1 * y;
}

// Integrand function for arc length
double CharactSinglySurface::integrand(double var, double fixed, bool along_x) const {
    return std::sqrt(1 + std::pow(along_x ? dz_dx(var) : dz_dy(var), 2));
}

// Simpson's rule for numerical integration
double CharactSinglySurface::simpson_integration(double a, double b, double fixed, bool along_x, int n) const {
    double h = (b - a) / n;
    double sum = integrand(a, fixed, along_x) + integrand(b, fixed, along_x);

    for (int i = 1; i < n; i += 2) {
        double var = a + i * h;
        sum += 4 * integrand(var, fixed, along_x);
    }
    for (int i = 2; i < n - 1; i += 2) {
        double var = a + i * h;
        sum += 2 * integrand(var, fixed, along_x);
    }
    return sum * h / 3;
}

// Function to calculate arc length
void CharactSinglySurface::calculate_arc_length(double start, double end, double fixed, bool along_x, double& arc_length, int num_intervals) const {
    arc_length = simpson_integration(start, end, fixed, along_x, num_intervals);
}

// Function to print curve points
void CharactSinglySurface::print_curve_points(double start, double end, double fixed, bool along_x, int num_points) const {
    std::cout << "\nSome points along the curve:" << std::endl;
    for (int i = 0; i < num_points; ++i) {
        double var = start + (end - start) * i / (num_points - 1);
        double x = along_x ? var : fixed;
        double y = along_x ? fixed : var;
        double z = along_x ? z_x(x, y) : z_y(x, y);
        std::cout << (along_x ? "x = " : "y = ") << std::fixed << std::setprecision(6) << var
            << ", z = " << z << std::endl;
    }
}

// Definition of the static method calculateCentralAngle
double CharactSinglySurface::calculateCentralAngle(double curveLength, double radius) {
    double angleInRadians = curveLength / radius;
    double angleInDegrees = angleInRadians * (180.0 / M_PI);
    return angleInDegrees;
}

// Method to compute the first derivative
void CharactSinglySurface::firstDerivative(double var, double& f_var, bool isXAxis) const {
    if (isXAxis) {
        f_var = 2 * c1 * var; // df/dx for X-axis equation
    }
    else {
        f_var = 2 * c1 * var; // df/dy for Y-axis equation
    }
}

// Method to compute the second derivative
void CharactSinglySurface::secondDerivative(double& f_var_var, bool isXAxis) const {
    if (isXAxis) {
        f_var_var = 2 * c1; // d^2f/dx^2 for X-axis equation
    }
    else {
        f_var_var = 2 * c1; // d^2f/dy^2 for Y-axis equation
    }
}

// Method to compute the radius of curvature
void CharactSinglySurface::computeRadius(double var, double& R, bool isXAxis) const {
    double f_var, f_var_var;
    firstDerivative(var, f_var, isXAxis);
    secondDerivative(f_var_var, isXAxis);

    // Calculate the radius of curvature
    double numerator = std::pow(1 + std::pow(f_var, 2), 1.5);
    double denominator = std::abs(f_var_var);

    R = (denominator != 0) ? numerator / denominator : std::numeric_limits<double>::infinity();
}