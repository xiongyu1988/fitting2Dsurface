#include <iostream>
#include <cmath>
#include <vector>
#include <memory>

class Characterize2DSurface {
public:
    struct SurfaceProperties {
        double area;
        double avg_R;
        double span_Theta;
    };

    class Surface {
    public:
        virtual double Z(double X, double Y) const = 0;
        virtual double Zx(double X, double Y) const = 0;
        virtual double Zy(double X, double Y) const = 0;
        virtual double Zxx(double X, double Y) const = 0;
        virtual double Zyy(double X, double Y) const = 0;
        virtual double Zxy(double X, double Y) const = 0;
        virtual ~Surface() = default;
    };

    class QuadraticSurface : public Surface {
    private:
        double c, a1, a2, a3, a4, a5;

    public:
        QuadraticSurface(double c, double a1, double a2, double a3, double a4, double a5)
            : c(c), a1(a1), a2(a2), a3(a3), a4(a4), a5(a5) {}

        double Z(double X, double Y) const override {
            return c + a1 * X + a2 * Y + a3 * X * X + a4 * X * Y + a5 * Y * Y;
        }
        double Zx(double X, double Y) const override {
            return a1 + 2 * a3 * X + a4 * Y;
        }
        double Zy(double X, double Y) const override {
            return a2 + a4 * X + 2 * a5 * Y;
        }
        double Zxx(double, double) const override { return 2 * a3; }
        double Zyy(double, double) const override { return 2 * a5; }
        double Zxy(double, double) const override { return a4; }
    };

    class SinglyCurvedShell : public Surface {
    private:
        double a, b, c;

    public:
        SinglyCurvedShell(double a, double b, double c)
            : a(a), b(b), c(c) {}

        double Z(double X, double Y) const override {
            return a * pow(X - b, 2) + c;
        }
        double Zx(double X, double) const override {
            return 2 * a * (X - b);
        }
        double Zy(double, double) const override {
            return 0;
        }
        double Zxx(double, double) const override { return 2 * a; }
        double Zyy(double, double) const override { return 0; }
        double Zxy(double, double) const override { return 0; }
    };

    Characterize2DSurface(std::pair<double, double> x_range, std::pair<double, double> y_range, int num_x, int num_y)
        : x_range(x_range), y_range(y_range), num_x(num_x), num_y(num_y) {}

    SurfaceProperties characterize(const Surface& surface) {
        return discretize_surface(surface, x_range, y_range, num_x, num_y);
    }

    static void print_properties(const SurfaceProperties& props) {
        std::cout << "Total Surface Area: " << props.area << std::endl;
        std::cout << "Average R: " << props.avg_R << std::endl;
        std::cout << "Theta span: " << props.span_Theta << std::endl;
    }

private:
    std::pair<double, double> x_range;
    std::pair<double, double> y_range;
    int num_x;
    int num_y;

    SurfaceProperties calculate_patch_properties(const Surface& surface, double X, double Y, double dX, double dY) {
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

    SurfaceProperties discretize_surface(const Surface& surface, std::pair<double, double> x_range,
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
};
