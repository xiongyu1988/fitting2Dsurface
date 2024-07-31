#include "FittingAlgorithms.h"
#include "ReadParseData.h"
#include <vector>
#include <map>
#include <Eigen/Dense>

Eigen::VectorXd FittingAlgorithms::fitDoublyCurvedShell() const {
    std::vector<Vector3D> points;
    std::map<int, Vector3D> gridData = readParseData.getGridData();
    for (const auto& point : gridData) {
        points.push_back(point.second);
    }

    size_t n = points.size();
    Eigen::MatrixXd A(n, 6);
    Eigen::VectorXd z(n);

    for (int i = 0; i < n; ++i) {
        double x = points[i].x;
        double y = points[i].y;
        A.row(i) << 1, x, y, x* x, x* y, y* y;
        z(i) = points[i].z;
    }

    return A.colPivHouseholderQr().solve(z);
}

Eigen::VectorXd FittingAlgorithms::fitCylindricalParaboloid() const {
    std::vector<Vector3D> points;
    std::map<int, Vector3D> gridData = readParseData.getGridData();
    for (const auto& point : gridData) {
        points.push_back(point.second);
    }

    size_t n = points.size();
    Eigen::MatrixXd A(n, 4);
    Eigen::VectorXd z(n);

    for (int i = 0; i < n; ++i) {
        double x = points[i].x;
        double y = points[i].y;
        A.row(i) << x * x, x, y, 1;
        z(i) = points[i].z;
    }

    Eigen::VectorXd params = A.colPivHouseholderQr().solve(z);

    double p = params(0);
    double q = params(1);
    double c = params(2);
    double r = params(3);

    double a = p;
    double b = -q / (2 * p);
    double d = r - (q * q) / (4 * p);

    return Eigen::Vector4d(a, b, c, d);
}

Eigen::VectorXd FittingAlgorithms::fitFlatPanel() const {
    std::vector<Vector3D> points;
    std::map<int, Vector3D> gridData = readParseData.getGridData();
    for (const auto& point : gridData) {
        points.push_back(point.second);
    }

    size_t n = points.size();
    Eigen::MatrixXd A(n, 3);
    Eigen::VectorXd z(n);

    for (int i = 0; i < n; ++i) {
        double x = points[i].x;
        double y = points[i].y;
        A.row(i) << 1, x, y;
        z(i) = points[i].z;
    }

    return A.colPivHouseholderQr().solve(z);
}