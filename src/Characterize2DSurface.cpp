#include "Characterize2DSurface.h"
#define _USE_MATH_DEFINES
#include <math.h> 

// Constructor implementation
Characterize2DSurface::Characterize2DSurface(double c1, double c2, double c3, double c4)
    : c1(c1), c2(c2), c3(c3), c4(c4) {}

// Method to compute the second derivatives at the center (origin)
void Characterize2DSurface::firstDerivatives(double x, double y, double& fx, double& fy) {
    fx = 2 * c1 * x + c2;
    fy = 2 * c4 * y + c3;
}

// Method to compute second derivatives (constant for this surface)
void Characterize2DSurface::secondDerivatives(double& fxx, double& fyy, double& fxy) {
    fxx = 2 * c1;
    fyy = 2 * c4;
    fxy = 0; // The cross derivative is zero for this surface
}

// Method to compute the radii and principal curvature lengths
void Characterize2DSurface::computeRadiiAndCurvatureLengths(double x, double y, double& R1, double& R2, double& L1, double& L2) {
    double fx, fy, fxx, fyy, fxy;
    firstDerivatives(x, y, fx, fy);
    secondDerivatives(fxx, fyy, fxy);

    // Radii of curvature (principal curvatures are inverse of these)
    R1 = (fxx != 0) ? 1 / fxx : std::numeric_limits<double>::infinity();
    R2 = (fyy != 0) ? 1 / fyy : std::numeric_limits<double>::infinity();

    // Principal curvature lengths, assuming significant curvature occurs at z = 1
    double z0 = 1.0;
    L1 = std::sqrt(2 * std::abs(R1) * z0);
    L2 = std::sqrt(2 * std::abs(R2) * z0);
}

// Function to calculate the central angle
// Definition of the static method calculateCentralAngle
double Characterize2DSurface::calculateCentralAngle(double curveLength, double radius) {
    // Calculate the angle in radians
    double angleInRadians = curveLength / radius;

    // Convert radians to degrees
    double angleInDegrees = angleInRadians * (180.0 / M_PI);

    return angleInDegrees;
}
