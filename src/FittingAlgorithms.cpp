#include "FittingAlgorithms.h"
#include "ReadParseData.h"
#include <vector>
#include <map>
#include <Eigen/Dense>
#include <stdexcept>
#include <cmath>
#include <limits>

// Constructor definition
FittingAlgorithms::FittingAlgorithms(const ReadParseData& data)
    : readParseData(data) {
    // Optionally, you can initialize or perform other setup here.
}

Eigen::VectorXd FittingAlgorithms::fitDoublyCurvedShell() const {
    std::vector<Vector3D> points;
    std::map<int, Vector3D> gridData = readParseData.getGridData();

    if (gridData.empty()) {
        throw std::runtime_error("No points available for fitting.");
    }

    for (const auto& point : gridData) {
        points.push_back(point.second);
    }

    size_t n = points.size();
    Eigen::MatrixXd A(n, 6);
    Eigen::VectorXd z(n);

    for (size_t i = 0; i < n; ++i) {
        double x = points[i].x;
        double y = points[i].y;
        A.row(i) << 1, x, y, x* x, x* y, y* y;
        z(i) = points[i].z;
    }

    // Perform the decomposition once and check for rank deficiency
    Eigen::FullPivHouseholderQR<Eigen::MatrixXd> qr(A);
    if (qr.rank() == A.cols()) { // Ensure the matrix is full rank
        return qr.solve(z);
    }
    else {
        throw std::runtime_error("Matrix is singular or near-singular, fitting failed.");
    }
}


Eigen::VectorXd FittingAlgorithms::fitSinglyCurvedShell(Axis axis) const {
    std::vector<Vector3D> points;
    std::map<int, Vector3D> gridData = readParseData.getGridData();

    if (gridData.empty()) {
        throw std::runtime_error("No points available for fitting.");
    }

    for (const auto& point : gridData) {
        points.push_back(point.second);
    }

    size_t n = points.size();
    Eigen::MatrixXd A(n, 4);
    Eigen::VectorXd z(n);

    for (size_t i = 0; i < n; ++i) {
        double x = points[i].x;
        double y = points[i].y;
        if (axis == Axis::X) {
            A.row(i) << x * x, x, y, 1;
        }
        else {
            A.row(i) << y * y, y, x, 1;
        }
        z(i) = points[i].z;
    }

    Eigen::VectorXd params = A.colPivHouseholderQr().solve(z);

    double p = params(0);
    double q = params(1);
    double c = params(2);
    double r = params(3);

    if (std::abs(p) < std::numeric_limits<double>::epsilon()) {
        throw std::runtime_error("Fitting failed due to division by near-zero value.");
    }

    double a = p;
    double b = -q / (2 * p);
    double d = r - (q * q) / (4 * p);

    return Eigen::Vector4d(a, b, c, d);
}

Eigen::VectorXd FittingAlgorithms::fitFlatPanel() const {
    std::vector<Vector3D> points;
    std::map<int, Vector3D> gridData = readParseData.getGridData();

    if (gridData.empty()) {
        throw std::runtime_error("No points available for fitting.");
    }

    for (const auto& point : gridData) {
        points.push_back(point.second);
    }

    size_t n = points.size();
    Eigen::MatrixXd A(n, 3);
    Eigen::VectorXd z(n);

    for (size_t i = 0; i < n; ++i) {
        double x = points[i].x;
        double y = points[i].y;
        A.row(i) << 1, x, y;
        z(i) = points[i].z;
    }

    return A.colPivHouseholderQr().solve(z);
}

Eigen::VectorXd FittingAlgorithms::fitCloseSinglyCurvedShell(Axis axis) const {
    std::vector<Vector3D> points;
    std::map<int, Vector3D> gridData = readParseData.getGridData();

    if (gridData.empty()) {
        throw std::runtime_error("No points available for fitting.");
    }

    for (const auto& point : gridData) {
        points.push_back(point.second);
    }

    size_t n = points.size();
    Eigen::MatrixXd A(n, 3);
    Eigen::VectorXd result(n);

    for (size_t i = 0; i < n; ++i) {
        double x = points[i].x;
        double y = points[i].y;
        if (axis == Axis::X) {
            A(i, 0) = 1;
            A(i, 1) = x;
            A(i, 2) = y;
            result(i) = -x * x - y * y;
        }
        else {
            A(i, 0) = 1;
            A(i, 1) = y;
            A(i, 2) = x;
            result(i) = -y * y - x * x;
        }
    }

    Eigen::VectorXd params = A.colPivHouseholderQr().solve(result);

    double a = -params(1) / 2;
    double b_center = -params(2) / 2;
    double R = std::sqrt(std::max(0.0, a * a + b_center * b_center - params(0)));

    return Eigen::Vector3d(a, b_center, R);
}
