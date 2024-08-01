#define _USE_MATH_DEFINES

#include "ReadParseData.h"
#include <cmath> 
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <regex>


Vector3D::Vector3D() : x(0), y(0), z(0) {}
Vector3D::Vector3D(double x, double y, double z) : x(x), y(y), z(z) {}

Vector3D Vector3D::operator-(const Vector3D& other) const {
    return Vector3D(x - other.x, y - other.y, z - other.z);
}

std::string trim(const std::string& str) {
    size_t start = str.find_first_not_of(" \t");
    size_t end = str.find_last_not_of(" \t");
    return (start == std::string::npos || end == std::string::npos) ? "" : str.substr(start, end - start + 1);
}

std::string fixScientificNotation(const std::string& str) {
    static const std::regex sci_notation_pattern(R"(([-+]?\d*\.?\d+)([-+]\d+))");
    std::smatch match;
    std::string result = str;
    if (std::regex_search(str, match, sci_notation_pattern)) {
        result = match[1].str() + "e" + match[2].str();
    }
    return result;
}

void ReadParseData::parseGridLine(const std::string& line) {
    if (trim(line.substr(0, 4)) == "GRID") {
        int nodeId = std::stoi(trim(line.substr(8, 8)));
        std::string xStr = trim(line.substr(24, 8));
        std::string yStr = trim(line.substr(32, 8));
        std::string zStr = trim(line.substr(40, 8));

        xStr = fixScientificNotation(xStr);
        yStr = fixScientificNotation(yStr);
        zStr = fixScientificNotation(zStr);

        double x = std::stod(xStr);
        double y = std::stod(yStr);
        double z = std::stod(zStr);

        gridData[nodeId] = Vector3D(x, y, z);
    }
}

void ReadParseData::parseMeshLine(const std::string& line) {
    if (trim(line.substr(0, 8)) == "CTRIA3") {
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

void ReadParseData::estimateAndAlignAxis() {
    if (gridData.empty()) {
        std::cout << "No data available for axis estimation." << std::endl;
        return;
    }

    // Step 1: Compute the centroid
    Eigen::Vector3d centroid(0, 0, 0);
    for (const auto& point : gridData) {
        centroid += Eigen::Vector3d(point.second.x, point.second.y, point.second.z);
    }
    centroid /= gridData.size();

    // Step 2: Build the covariance matrix
    Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
    for (const auto& point : gridData) {
        Eigen::Vector3d p(point.second.x, point.second.y, point.second.z);
        Eigen::Vector3d centered = p - centroid;
        covariance += centered * centered.transpose();
    }
    covariance /= gridData.size();

    // Step 3: Compute eigenvectors and eigenvalues
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(covariance);
    Eigen::Vector3d main_axis = eigensolver.eigenvectors().col(2); // The eigenvector with the largest eigenvalue

    // Step 4: Check if alignment is needed
    Eigen::Vector3d z_axis(0, 0, 1);
    double angle = std::acos(main_axis.dot(z_axis));

    std::cout << "Estimated main axis: (" << main_axis.x() << ", " << main_axis.y() << ", " << main_axis.z() << ")" << std::endl;
    std::cout << "Angle with z-axis: " << angle * 180 / M_PI << " degrees" << std::endl;

    // Step 5: Align with z-axis if needed
    if (std::abs(angle) > 1e-6) { // If angle is not close to 0 or 180 degrees
        Eigen::Vector3d rotation_axis = main_axis.cross(z_axis).normalized();
        Eigen::AngleAxisd rotation(angle, rotation_axis);
        Eigen::Matrix3d rotation_matrix = rotation.toRotationMatrix();

        // Apply rotation to all points
        for (auto& point : gridData) {
            Eigen::Vector3d p(point.second.x, point.second.y, point.second.z);
            p = rotation_matrix * (p - centroid) + centroid;
            point.second = Vector3D(p.x(), p.y(), p.z());
        }

        std::cout << "Geometry aligned with z-axis." << std::endl;
    }
    else {
        std::cout << "Geometry is already aligned with z-axis." << std::endl;
    }

    // Recalculate coordinate ranges after alignment
    calculateCoordinateRanges();
}

std::map<int, Vector3D> ReadParseData::getGridData() const {return gridData;}

std::map<int, std::vector<int>> ReadParseData::getMeshElements() const {return meshElements;}