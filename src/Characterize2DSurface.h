#ifndef CHARACTERIZE_2D_SURFACE_H
#define CHARACTERIZE_2D_SURFACE_H

#include <utility>
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
        QuadraticSurface(double c, double a1, double a2, double a3, double a4, double a5);

        double Z(double X, double Y) const override;
        double Zx(double X, double Y) const override;
        double Zy(double X, double Y) const override;
        double Zxx(double X, double Y) const override;
        double Zyy(double X, double Y) const override;
        double Zxy(double X, double Y) const override;
    };

    class SinglyCurvedShell : public Surface {
    private:
        double c1, c2, c3, c4;

    public:
        SinglyCurvedShell(double c1, double c2, double c3, double c4);

        double Z(double X, double Y) const override;
        double Zx(double X, double Y) const override;
        double Zy(double X, double Y) const override;
        double Zxx(double X, double Y) const override;
        double Zyy(double X, double Y) const override;
        double Zxy(double X, double Y) const override;
    };

    Characterize2DSurface(std::pair<double, double> x_range, std::pair<double, double> y_range, int num_x, int num_y);

    SurfaceProperties characterize(const Surface& surface);

    static void print_properties(const SurfaceProperties& props);

private:
    std::pair<double, double> x_range;
    std::pair<double, double> y_range;
    int num_x;
    int num_y;

    SurfaceProperties calculate_patch_properties(const Surface& surface, double X, double Y, double dX, double dY);
    SurfaceProperties discretize_surface(const Surface& surface, std::pair<double, double> x_range,
        std::pair<double, double> y_range, int num_x, int num_y);
};

#endif // CHARACTERIZE_2D_SURFACE_H