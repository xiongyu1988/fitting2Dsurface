#ifndef READ_PARSE_DATA_H
#define READ_PARSE_DATA_H

#include <map>
#include <vector>
#include <string>
#include <Eigen/Dense>

class Vector3D {
public:
    double x, y, z;

    Vector3D();
    Vector3D(double x, double y, double z);

    Vector3D operator-(const Vector3D& other) const;
};

class ReadParseData {
private:
    std::map<int, Vector3D> gridData;
    std::map<int, std::vector<int>> meshElements;

    struct Range {
        double min;
        double max;
    };

    struct CoordinateRanges {
        Range x;
        Range y;
        Range z;
    };

    CoordinateRanges ranges;

    void parseGridLine(const std::string& line);
    void parseMeshLine(const std::string& line);

public:
    void readFromFile(const std::string& filename);
    double calculateSurfaceArea() const;
    /// <summary>
    /// The equation of an elliptic paraboloid is given by:
    /// z = c1 + c2 * x + c3 * y + c4 * x^2 + c5 * x * y + c6 * y^2
    /// </summary>
    /// <returns></returns>
    Eigen::VectorXd fitDoublyCurvedShell() const;
    /// <summary>
    /// The equation of a cylindrical paraboloid is given by:
    /// z = c1 * (x - c2)^2 + c3 * y + c4
    /// a simpler approach to solve this non - linear least squares fitting problem.
    /// One such approach is to transform the problem into a linear least squares problem 
    /// by rearranging the equation.Here's how we can do it:
    /// Start with the equation : z = a * (x - b)² + c * y + d
    /// Expand it : z = ax² - 2abx + ab² + c * y + d
    /// Let p = a, q = -2ab, r = ab² + d
    /// Now we have : z = px² + qx + c * y + r
    /// This is now linear in terms of the parameters p, q, c, and r.We can solve this using simple linear least squares.
    /// </summary>
    /// <returns></returns>
    Eigen::VectorXd fitCylindricalParaboloid() const;
    void calculateCoordinateRanges();
    void printCoordinateRanges() const;
    const CoordinateRanges& getCoordinateRanges() const;
    void printAllGridData() const;
    void printAllMeshElements() const;
};

#endif // READ_PARSE_DATA_H