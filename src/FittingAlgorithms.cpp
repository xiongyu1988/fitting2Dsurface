#include "FittingAlgorithms.h"
#include "ReadParseData.h"
#include <vector>
#include <map>
#include <Eigen/Dense>
#include <stdexcept>

#include <cmath>
#include <limits>

//Eigen::VectorXd FittingAlgorithms::fitDoublyCurvedShell() const {
//    std::vector<Vector3D> points;
//    std::map<int, Vector3D> gridData = readParseData.getGridData();
//    for (const auto& point : gridData) {
//        points.push_back(point.second);
//    }
//
//    size_t n = points.size();
//    Eigen::MatrixXd A(n, 6);
//    Eigen::VectorXd z(n);
//
//    for (int i = 0; i < n; ++i) {
//        double x = points[i].x;
//        double y = points[i].y;
//        A.row(i) << 1, x, y, x* x, x* y, y* y;
//        z(i) = points[i].z;
//    }
//
//    return A.colPivHouseholderQr().solve(z);
//}
//
//Eigen::VectorXd FittingAlgorithms::fitSinglyCurvedShell(Axis axis) const {
//    std::vector<Vector3D> points;
//    std::map<int, Vector3D> gridData = readParseData.getGridData();
//    for (const auto& point : gridData) {
//        points.push_back(point.second);
//    }
//
//    size_t n = points.size();
//    Eigen::MatrixXd A(n, 4);
//    Eigen::VectorXd z(n);
//
//    for (size_t i = 0; i < n; ++i) {
//        double x = points[i].x;
//        double y = points[i].y;
//        if (axis == Axis::X) {
//            A.row(i) << x * x, x, y, 1;
//        }
//        else {
//            A.row(i) << y * y, y, x, 1;
//        }
//        z(i) = points[i].z;
//    }
//
//    Eigen::VectorXd params = A.colPivHouseholderQr().solve(z);
//
//    double p = params(0);
//    double q = params(1);
//    double c = params(2);
//    double r = params(3);
//
//    double a = p;
//    double b = -q / (2 * p);
//    double d = r - (q * q) / (4 * p);
//
//    return Eigen::Vector4d(a, b, c, d);
//}
//
//Eigen::VectorXd FittingAlgorithms::fitFlatPanel() const {
//    std::vector<Vector3D> points;
//    std::map<int, Vector3D> gridData = readParseData.getGridData();
//    for (const auto& point : gridData) {
//        points.push_back(point.second);
//    }
//
//    size_t n = points.size();
//    Eigen::MatrixXd A(n, 3);
//    Eigen::VectorXd z(n);
//
//    for (int i = 0; i < n; ++i) {
//        double x = points[i].x;
//        double y = points[i].y;
//        A.row(i) << 1, x, y;
//        z(i) = points[i].z;
//    }
//
//    return A.colPivHouseholderQr().solve(z);
//}
//
//Eigen::VectorXd FittingAlgorithms::fitCloseSinglyCurvedShell(Axis axis) const {
//    std::vector<Vector3D> points;
//    std::map<int, Vector3D> gridData = readParseData.getGridData();
//    for (const auto& point : gridData) {
//        points.push_back(point.second);
//    }
//
//    size_t n = points.size();
//    Eigen::MatrixXd A(n, 3);
//    Eigen::VectorXd b(n);
//
//    for (size_t i = 0; i < n; ++i) {
//        double x = points[i].x;
//        double y = points[i].y;
//        if (axis == Axis::X) {
//            double z = points[i].z;
//            A(i, 0) = 1;
//            A(i, 1) = x;
//            A(i, 2) = y;
//            b(i) = -x * x - y * y;
//        }
//        else {
//            double z = points[i].z;
//            A(i, 0) = 1;
//            A(i, 1) = y;
//            A(i, 2) = x;
//            b(i) = -y * y - x * x;
//        }
//    }
//
//    Eigen::VectorXd params = A.colPivHouseholderQr().solve(b);
//
//    double a = -params(1) / 2;
//    double b_center = -params(2) / 2;
//    double R = std::sqrt(a * a + b_center * b_center - params(0));
//
//    return Eigen::Vector3d(a, b_center, R);
//}

// ------------------------------------------------------------
// New aproach for fitting SynclasticShell
// ------------------------------------------------------------

FittingAlgorithms::FittingAlgorithms(const ReadParseData& data)
    : m_data(data), m_fittedParams(Eigen::VectorXd::Zero(5)), m_center(Eigen::Vector3d::Zero()), m_scale(1.0) {}

void FittingAlgorithms::fitSynclasticShell() {
    centerAndScaleData();
    Eigen::MatrixXd X = buildDesignMatrix();
    Eigen::VectorXd y = buildTargetVector();

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(X, Eigen::ComputeThinU | Eigen::ComputeThinV);
    m_fittedParams = svd.solve(y);
}

Eigen::VectorXd FittingAlgorithms::getFittedParameters() const {
    return m_fittedParams;
}

std::pair<double, double> FittingAlgorithms::getRadiiOfCurvature() const {
    double R1 = 0.0, R2 = 0.0;
    if (std::abs(m_fittedParams[0]) > std::numeric_limits<double>::epsilon() &&
        std::abs(m_fittedParams[1]) > std::numeric_limits<double>::epsilon()) {
        R1 = 1.0 / (2 * std::abs(m_fittedParams[0]));
        R2 = 1.0 / (2 * std::abs(m_fittedParams[1]));
    }
    return std::make_pair(R1, R2);
}

double FittingAlgorithms::calculateFitError() const {
    Eigen::MatrixXd X = buildDesignMatrix();
    Eigen::VectorXd y = buildTargetVector();
    Eigen::VectorXd residuals = X * m_fittedParams - y;
    return residuals.squaredNorm() / residuals.size();
}

Eigen::Vector3d FittingAlgorithms::getDataCenter() const {
    return m_center;
}

double FittingAlgorithms::getDataScale() const {
    return m_scale;
}

void FittingAlgorithms::centerAndScaleData() {
    const std::map<int, Vector3D>& gridData = m_data.getGridData();
    Eigen::Vector3d sum = Eigen::Vector3d::Zero();
    double maxDist = 0.0;

    for (const auto& point : gridData) {
        Eigen::Vector3d p(point.second.x, point.second.y, point.second.z);
        sum += p;
        maxDist = std::max(maxDist, p.norm());
    }

    m_center = sum / gridData.size();
    m_scale = maxDist > 0 ? 1.0 / maxDist : 1.0;
}

Eigen::MatrixXd FittingAlgorithms::buildDesignMatrix() const {
    const std::map<int, Vector3D>& gridData = m_data.getGridData();
    Eigen::MatrixXd X(gridData.size(), 5);

    int row = 0;
    for (const auto& point : gridData) {
        Eigen::Vector3d p(point.second.x, point.second.y, point.second.z);
        p = (p - m_center) * m_scale;
        X(row, 0) = p.x() * p.x();
        X(row, 1) = p.y() * p.y();
        X(row, 2) = p.x();
        X(row, 3) = p.y();
        X(row, 4) = 1;
        ++row;
    }

    return X;
}

Eigen::VectorXd FittingAlgorithms::buildTargetVector() const {
    const std::map<int, Vector3D>& gridData = m_data.getGridData();
    Eigen::VectorXd y(gridData.size());
    int row = 0;
    for (const auto& point : gridData) {
        Eigen::Vector3d p(point.second.x, point.second.y, point.second.z);
        p = (p - m_center) * m_scale;
        y(row) = p.z();
        ++row;
    }
    return y;
}