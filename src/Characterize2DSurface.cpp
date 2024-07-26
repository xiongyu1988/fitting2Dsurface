#include "Characterize2DSurface.h"
#include <iostream>
#include <cmath>

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

Characterize2DSurface::SinglyCurvedShell::SinglyCurvedShell(double c1, double c2, double c3, double c4)
    : c1(c1), c2(c2), c3(c3), c4(c4) {}

double Characterize2DSurface::SinglyCurvedShell::Z(double X, double Y) const {
    return c1 * pow(X - c2, 2) + c3 * Y + c4;
}

double Characterize2DSurface::SinglyCurvedShell::Zx(double X, double) const {
    return 2 * c1 * (X - c2);
}

double Characterize2DSurface::SinglyCurvedShell::Zy(double, double) const {
    return c3;
}

double Characterize2DSurface::SinglyCurvedShell::Zxx(double, double) const {
    return 2 * c1;
}

double Characterize2DSurface::SinglyCurvedShell::Zyy(double, double) const {
    return 0;
}

double Characterize2DSurface::SinglyCurvedShell::Zxy(double, double) const {
    return 0;
}

Characterize2DSurface::Characterize2DSurface(std::pair<double, double> x_range, std::pair<double, double> y_range, int num_x, int num_y)
    : x_range(x_range), y_range(y_range), num_x(num_x), num_y(num_y) {}

Characterize2DSurface::SurfaceProperties Characterize2DSurface::characterize(const Surface& surface) {
    return discretize_surface(surface, x_range, y_range, num_x, num_y);
}

void Characterize2DSurface::print_properties(const SurfaceProperties& props) {
    std::cout << "Total Surface Area: " << props.area << std::endl;
    std::cout << "Average R: " << props.avg_R << std::endl;
    std::cout << "Theta span: " << props.span_Theta << std::endl;
}

Characterize2DSurface::SurfaceProperties Characterize2DSurface::calculate_patch_properties(const Surface& surface, double X, double Y, double dX, double dY) {
    double Zx = surface.Zx(X, Y);
    double Zy = surface.Zy(X, Y);
    double Zxx = surface.Zxx(X, Y);

    double E = 1 + Zx * Zx;
    double F = Zx * Zy;
    double G = 1 + Zy * Zy;
    double L = Zxx / sqrt(1 + Zx * Zx + Zy * Zy);

    double k = L / (E * sqrt(1 + Zx * Zx + Zy * Zy));

    double R = 1 / std::abs(k);
    double Theta = acos(sqrt(E / (E + G)));

    double patch_area = sqrt(E * G - F * F) * dX * dY;

    return { patch_area, R, Theta };
}

Characterize2DSurface::SurfaceProperties Characterize2DSurface::discretize_surface(const Surface& surface, std::pair<double, double> x_range,
    std::pair<double, double> y_range, int num_x, int num_y) {
    double x_min = x_range.first, x_max = x_range.second;
    double y_min = y_range.first, y_max = y_range.second;

    double dX = (x_max - x_min) / (num_x - 1);
    double dY = (y_max - y_min) / (num_y - 1);

    double total_area = 0, total_R = 0, total_Theta = 0;

    for (int i = 0; i < num_x; ++i) {
        double x = x_min + i * dX;
        for (int j = 0; j < num_y; ++j) {
            double y = y_min + j * dY;
            auto props = calculate_patch_properties(surface, x, y, dX, dY);
            total_area += props.area;
            total_R += props.avg_R;
            total_Theta += props.span_Theta;
        }
    }

    int num_patches = num_x * num_y;
    return {
        total_area,
        total_R / num_patches,
        total_Theta / num_x  // Assuming Theta varies along X
    };
}