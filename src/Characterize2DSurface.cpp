#include "Characterize2DSurface.h"
#include <iostream>
#include <cmath>

// QuadraticSurface methods
Characterize2DSurface::QuadraticSurface::QuadraticSurface(double c, double a1, double a2, double a3, double a4, double a5)
    : c(c), a1(a1), a2(a2), a3(a3), a4(a4), a5(a5) {}

double Characterize2DSurface::QuadraticSurface::Z(double X, double Y) const {
    return c + a1 * X + a2 * Y + a3 * X * X + a4 * X * Y + a5 * Y * Y;
}

double Characterize2DSurface::QuadraticSurface::Zx(double X, double Y) const {
    return a1 + 2 * a3 * X + a4 * Y;
}

double Characterize2DSurface::QuadraticSurface::Zy(double X, double Y) const {
    return a2 + a4 * X + 2 * a5 * Y;
}

double Characterize2DSurface::QuadraticSurface::Zxx(double, double) const {
    return 2 * a3;
}

double Characterize2DSurface::QuadraticSurface::Zyy(double, double) const {
    return 2 * a5;
}

double Characterize2DSurface::QuadraticSurface::Zxy(double, double) const {
    return a4;
}

// SinglyCurvedShell methods
Characterize2DSurface::SinglyCurvedShell::SinglyCurvedShell(double c1, double c2, double c3, double c4, bool rotate_along_x)
    : c1(c1), c2(c2), c3(c3), c4(c4), rotate_along_x(rotate_along_x) {}

double Characterize2DSurface::SinglyCurvedShell::Z(double X, double Y) const {
    if (rotate_along_x) {
        return c1 * pow(Y - c2, 2) + c3 * X + c4;
    }
    else {
        return c1 * pow(X - c2, 2) + c3 * Y + c4;
    }
}

double Characterize2DSurface::SinglyCurvedShell::Zx(double X, double Y) const {
    return rotate_along_x ? c3 : 2 * c1 * (X - c2);
}

double Characterize2DSurface::SinglyCurvedShell::Zy(double X, double Y) const {
    return rotate_along_x ? 2 * c1 * (Y - c2) : c3;
}

double Characterize2DSurface::SinglyCurvedShell::Zxx(double, double) const {
    return rotate_along_x ? 0 : 2 * c1;
}

double Characterize2DSurface::SinglyCurvedShell::Zyy(double, double) const {
    return rotate_along_x ? 2 * c1 : 0;
}

double Characterize2DSurface::SinglyCurvedShell::Zxy(double, double) const {
    return 0;
}

// Characterize2DSurface methods
Characterize2DSurface::Characterize2DSurface(std::pair<double, double> x_range, std::pair<double, double> y_range, int num_x, int num_y)
    : x_range(x_range), y_range(y_range), num_x(num_x), num_y(num_y) {}

Characterize2DSurface::SurfaceProperties Characterize2DSurface::characterize(const Surface& surface) {
    return discretize_surface(surface, x_range, y_range, num_x, num_y);
}

void Characterize2DSurface::print_properties(const SurfaceProperties& props) {
    std::cout << "Total Surface Area: " << props.area << std::endl;
    std::cout << "Average Radius of Curvature R1: " << props.avg_R1 << std::endl;
    std::cout << "Average Radius of Curvature R2: " << props.avg_R2 << std::endl;
    std::cout << "Theta span 1: " << props.span_Theta1 << std::endl;
    std::cout << "Theta span 2: " << props.span_Theta2 << std::endl;
}

Characterize2DSurface::SurfaceProperties Characterize2DSurface::calculate_patch_properties(const Surface& surface, double X, double Y, double dX, double dY) {
    double Zx = surface.Zx(X, Y);
    double Zy = surface.Zy(X, Y);
    double Zxx = surface.Zxx(X, Y);
    double Zyy = surface.Zyy(X, Y);
    double Zxy = surface.Zxy(X, Y);

    double E = 1 + Zx * Zx;
    double F = Zx * Zy;
    double G = 1 + Zy * Zy;
    double L = Zxx / sqrt(1 + Zx * Zx + Zy * Zy);
    double M = Zxy / sqrt(1 + Zx * Zx + Zy * Zy);
    double N = Zyy / sqrt(1 + Zx * Zx + Zy * Zy);

    double k1 = (L * G - 2 * F * M + E * N) / (E * G - F * F);
    double k2 = (L * G - E * N) / (E * G - F * F);

    double R1 = 1 / std::abs(k1);
    double R2 = 1 / std::abs(k2);

    double Theta1 = acos(sqrt(E / (E + G)));
    double Theta2 = acos(sqrt(G / (E + G)));

    double patch_area = sqrt(E * G - F * F) * dX * dY;

    return { patch_area, R1, R2, Theta1, Theta2 };
}

Characterize2DSurface::SurfaceProperties Characterize2DSurface::discretize_surface(const Surface& surface, std::pair<double, double> x_range,
    std::pair<double, double> y_range, int num_x, int num_y) {
    double x_min = x_range.first, x_max = x_range.second;
    double y_min = y_range.first, y_max = y_range.second;

    double dX = (x_max - x_min) / (num_x - 1);
    double dY = (y_max - y_min) / (num_y - 1);

    double total_area = 0, total_R1 = 0, total_R2 = 0, total_Theta1 = 0, total_Theta2 = 0;

    for (int i = 0; i < num_x; ++i) {
        double x = x_min + i * dX;
        for (int j = 0; j < num_y; ++j) {
            double y = y_min + j * dY;
            auto props = calculate_patch_properties(surface, x, y, dX, dY);
            total_area += props.area;
            total_R1 += props.avg_R1;
            total_R2 += props.avg_R2;
            total_Theta1 += props.span_Theta1;
            total_Theta2 += props.span_Theta2;
        }
    }

    int num_patches = num_x * num_y;
    return {
        total_area,
        total_R1 / num_patches,
        total_R2 / num_patches,
        total_Theta1 / num_x,  
        total_Theta2 / num_y   
    };
}
