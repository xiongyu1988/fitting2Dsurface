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

Eigen::VectorXd FittingAlgorithms::fitDoublyCurvedShell2() const {
    std::vector<Vector3D> points;
    std::map<int, Vector3D> gridData = readParseData.getGridData();
    for (const auto& point : gridData) {
        points.push_back(point.second);
    }

    size_t n = points.size();
    Eigen::MatrixXd A(n, 4); // 4 coefficients: a, b, c, d
    Eigen::VectorXd z(n);

    for (size_t i = 0; i < n; ++i) {
        double x = points[i].x;
        double y = points[i].y;
        A.row(i) << x * x, x, y, y* y;
        z(i) = points[i].z;
    }

    Eigen::VectorXd params = A.colPivHouseholderQr().solve(z);

    // The output vector contains [a, b, c, d] corresponding to z = a*x^2 + b*x + c*y + d*y^2
    return params;
}

Eigen::VectorXd FittingAlgorithms::fitSinglyCurvedShell(Axis axis) const {
    std::vector<Vector3D> points;
    std::map<int, Vector3D> gridData = readParseData.getGridData();
    for (const auto& point : gridData) {
        points.push_back(point.second);
    }

    size_t n = points.size();
    Eigen::MatrixXd A(n, 2); // 2 coefficients: a for quadratic term and b for linear term
    Eigen::VectorXd z(n);

    for (size_t i = 0; i < n; ++i) {
        double x = points[i].x;
        double y = points[i].y;
        if (axis == Axis::X) {
            // Fit z = a * x^2 + b * y
            A.row(i) << x * x, y;
        }
        else if (axis == Axis::Y) {
            // Fit z = a * y^2 + b * x
            A.row(i) << y * y, x;
        }
        z(i) = points[i].z;
    }

    Eigen::VectorXd params = A.colPivHouseholderQr().solve(z);

    // The output vector contains [a, b] corresponding to either:
    // - z = a * x^2 + b * y (if axis == Axis::X)
    // - z = a * y^2 + b * x (if axis == Axis::Y)
    return params;
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
