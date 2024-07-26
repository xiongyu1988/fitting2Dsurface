#include "ReadParseData.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <limits>

Vector3D::Vector3D() : x(0), y(0), z(0) {}
Vector3D::Vector3D(double x, double y, double z) : x(x), y(y), z(z) {}

Vector3D Vector3D::operator-(const Vector3D& other) const {
    return Vector3D(x - other.x, y - other.y, z - other.z);
}

void ReadParseData::parseGridLine(const std::string& line) {
    if (line.substr(0, 4) == "GRID") {
        int nodeId = std::stoi(line.substr(8, 8));
        double x = std::stod(line.substr(24, 8));
        double y = std::stod(line.substr(32, 8));
        double z = std::stod(line.substr(40, 8));
        gridData[nodeId] = Vector3D(x, y, z);
    }
}

void ReadParseData::parseMeshLine(const std::string& line) {
    if (line.substr(0, 6) == "CTRIA3") {
        std::istringstream lineStream(line);
        std::string dummy;
        int elementId, node1, node2, node3;

        lineStream >> dummy >> elementId >> dummy >> node1 >> node2 >> node3;
        meshElements[elementId] = std::vector<int>{ node1, node2, node3 };
    }
}

void ReadParseData::readFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open file: " + filename);
    }

    std::string line;
    while (std::getline(file, line)) {
        parseGridLine(line);
        parseMeshLine(line);
    }

    file.close();
}

double ReadParseData::calculateSurfaceArea() const {
    double totalArea = 0;
    for (const auto& element : meshElements) {
        const auto& nodes = element.second;
        if (nodes.size() == 3) {
            const Vector3D& p1 = gridData.at(nodes[0]);
            const Vector3D& p2 = gridData.at(nodes[1]);
            const Vector3D& p3 = gridData.at(nodes[2]);

            Vector3D v1 = p2 - p1;
            Vector3D v2 = p3 - p1;
            Vector3D cross(
                v1.y * v2.z - v1.z * v2.y,
                v1.z * v2.x - v1.x * v2.z,
                v1.x * v2.y - v1.y * v2.x
            );
            double elementArea = 0.5 * std::sqrt(cross.x * cross.x + cross.y * cross.y + cross.z * cross.z);
            totalArea += elementArea;
        }
    }
    return totalArea;
}

Eigen::VectorXd ReadParseData::fitDoublyCurvedShell() const {
    std::vector<Vector3D> points;
    for (const auto& point : gridData) {
        points.push_back(point.second);
    }

    int n = points.size();
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

Eigen::VectorXd ReadParseData::fitCylindricalParaboloid() const {
    std::vector<Vector3D> points;
    for (const auto& point : gridData) {
        points.push_back(point.second);
    }

    int n = points.size();
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

void ReadParseData::calculateCoordinateRanges() {
    if (gridData.empty()) {
        std::cout << "No data available." << std::endl;
        return;
    }

    ranges.x = { std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest() };
    ranges.y = { std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest() };
    ranges.z = { std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest() };

    for (const auto& point : gridData) {
        const Vector3D& coord = point.second;
        ranges.x.min = std::min(ranges.x.min, coord.x);
        ranges.x.max = std::max(ranges.x.max, coord.x);
        ranges.y.min = std::min(ranges.y.min, coord.y);
        ranges.y.max = std::max(ranges.y.max, coord.y);
        ranges.z.min = std::min(ranges.z.min, coord.z);
        ranges.z.max = std::max(ranges.z.max, coord.z);
    }
}

void ReadParseData::printCoordinateRanges() const {
    std::cout << "\nCoordinate Ranges:" << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "x range: [" << ranges.x.min << ", " << ranges.x.max << "]" << std::endl;
    std::cout << "y range: [" << ranges.y.min << ", " << ranges.y.max << "]" << std::endl;
    std::cout << "z range: [" << ranges.z.min << ", " << ranges.z.max << "]" << std::endl;
}

const ReadParseData::CoordinateRanges& ReadParseData::getCoordinateRanges() const {
    return ranges;
}

void ReadParseData::printAllGridData() const {
    std::cout << "All Grid Data:" << std::endl;
    for (const auto& point : gridData) {
        std::cout << "Node " << point.first << ": ("
            << std::setprecision(6) << std::fixed << point.second.x << ", "
            << std::setprecision(6) << std::fixed << point.second.y << ", "
            << std::setprecision(6) << std::fixed << point.second.z << ")" << std::endl;
    }
}

void ReadParseData::printAllMeshElements() const {
    std::cout << "\nAll Mesh Elements:" << std::endl;
    for (const auto& element : meshElements) {
        std::cout << "Element " << element.first << ": ";
        for (int node : element.second) {
            std::cout << node << " ";
        }
        std::cout << std::endl;
    }
}