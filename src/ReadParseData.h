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

	struct CenterCoordinates {
		double centerX;
		double centerY;
		double centerZ;
	};

    struct CoordinateRanges {
        Range x;
        Range y;
        Range z;
    };

    CoordinateRanges ranges;
    CenterCoordinates center;

    void parseGridLine(const std::string& line);
    void parseMeshLine(const std::string& line);

public:
    void readFromFile(const std::string& filename);
    double calculateSurfaceArea() const;
    void calculateCoordinateRanges();
    void printCoordinateRanges() const;
    void calculateCenterCoordantes();
    const CoordinateRanges& getCoordinateRanges() const;
	const CenterCoordinates& getCenterCoordinates() const;
    void printAllGridData() const;
    void printAllMeshElements() const;
    std::map<int, Vector3D> getGridData() const;
    std::map<int, std::vector<int>> getMeshElements() const;
};

#endif // READ_PARSE_DATA_H