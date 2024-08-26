#include "ReadParseData.h"
#include "CharactDoublySurface.h"
#include "CharactSinglySurface.h"
#include "FittingAlgorithms.h"
#include <iostream>
#include <iomanip>

int main() {

    // ReadParseData
    ReadParseData readParseData;

    try {
        readParseData.readFromFile("data/doublyCurvedShell_2.fem");
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    // FittingAlgorithms
    FittingAlgorithms fittingAlgorithms(readParseData);

    //readParseData.printAllGridData();
    //readParseData.printAllMeshElements();

    std::cout << "Read and Analyze FEM Data" << std::endl;
    std::cout << "--------------------------------------------" << std::endl;
    double meshArea = readParseData.calculateSurfaceArea();
    std::cout << "Mesh surface area: " << std::setprecision(6) << std::fixed << meshArea << std::endl;

    std::cout << "\n\n";
    std::cout << "Coordinate Range and Center" << std::endl;
    std::cout << "--------------------------------------------" << std::endl;
    readParseData.calculateCoordinateRanges();
    readParseData.printCoordinateRanges();
    // Access ranges for later use
    const auto& ranges = readParseData.getCoordinateRanges();
    readParseData.calculateCenterCoordantes();

  //  std::cout << "\nFitting doubly curved shell (elliptic paraboloid):" << std::endl;
  //  Eigen::VectorXd doublyCurvedCoeffs = fittingAlgorithms.fitDoublyCurvedShell();
  //  std::cout << "Shell equation: z = " << doublyCurvedCoeffs(0) << " + " << doublyCurvedCoeffs(1) << " * x + "
		//<< doublyCurvedCoeffs(2) << " * y + " << doublyCurvedCoeffs(3) << " * x^2 + " << doublyCurvedCoeffs(4) << " * x * y + "
		//<< doublyCurvedCoeffs(5) << " * y^2 " << std::endl;
    std::cout << "\n\n";
	std::cout << "1) Fitting doubly curved shell (elliptic paraboloid) using a NEW approach:" << std::endl;
    std::cout << "--------------------------------------------" << std::endl;
    std::cout << "Fitting Algorithms: " << std::endl;
    Eigen::VectorXd coefficients = fittingAlgorithms.fitDoublyCurvedShell2();
    std::cout << "Fitted equation: z = " << coefficients(0) << " * x^2 + " << coefficients(1) 
              << " * x + " << coefficients(2) << " * y + " << coefficients(3) << " * y^2" << std::endl;
    std::cout << "Characterize Fitted 2D Surface: " << std::endl;
    CharactDoublyCurved2DSurface surface(coefficients(0), coefficients(1), coefficients(2), coefficients(3));
    const auto& center = readParseData.getCenterCoordinates();
    double R1X, R2Y, L1X, L2Y;
    // For fixed y, varying x
    surface.calculate_arc_length(ranges.x.min, ranges.x.max, center.centerY, true, L1X);
    // For fixed x, varying y
    surface.calculate_arc_length(ranges.y.min, ranges.y.max, center.centerX, false, L2Y);
    surface.computeRadii(center.centerX, center.centerY, R1X, R2Y);
    std::cout << "Principal Radius of Curvature R1: " << R1X << std::endl;
    std::cout << "Principal Radius of Curvature R2: " << R2Y << std::endl;
    std::cout << "Principal Curvature Length L1: " << L1X << std::endl;
    std::cout << "Principal Curvature Length L2: " << L2Y << std::endl;
    double centerAngle1, centerAngle2;
    centerAngle1 = CharactDoublySurface::calculateCentralAngle(L1X, R1X);
    centerAngle2 = CharactDoublySurface::calculateCentralAngle(L2Y, R2Y);
    std::cout << "Central Angle 1: " << centerAngle1 << std::endl;
    std::cout << "Central Angle 2: " << centerAngle2 << std::endl;

    std::cout << "\n\n";
    std::cout << "2) Fitting singly curved shell (cylindrical paraboloid) using a NEW approach:" << std::endl;
    std::cout << "--------------------------------------------" << std::endl;
    Eigen::VectorXd coefitSinglyCurvedShell2X = fittingAlgorithms.fitSinglyCurvedShell(Axis::X);
    std::cout << "Fitted equation along X-axis: z = " << coefitSinglyCurvedShell2X(0) << " * x^2 + " << coefitSinglyCurvedShell2X(1) << " * y" << std::endl;
    Eigen::VectorXd coefitSinglyCurvedShell2Y = fittingAlgorithms.fitSinglyCurvedShell(Axis::Y);
    std::cout << "Fitted equation along Y-axis: z = " << coefitSinglyCurvedShell2Y(0) << " * y^2 + " << coefitSinglyCurvedShell2Y(1) << " * x" << std::endl;
    CharactSinglySurface surface2(coefitSinglyCurvedShell2Y(0), coefitSinglyCurvedShell2Y(1));
    double arc_length_y;
    surface2.calculate_arc_length(ranges.y.min, ranges.y.max, center.centerX, false, arc_length_y); // Along Y-axis
    double R_x, R_y;
    surface2.computeRadius(center.centerY, R_y, false);
    std::cout << "Radius of Curvature along Y: " << R_y << std::endl;
    std::cout << "Principal Curvature Length L2: " << arc_length_y << std::endl;
	double centerAngleY = surface2.calculateCentralAngle(arc_length_y, R_y);
	std::cout << "Central Angle: " << centerAngleY << std::endl;

    std::cout << "\n\n";
    std::cout << "3) Fitting flat panel surface:" << std::endl;
    std::cout << "--------------------------------------------" << std::endl;
    Eigen::VectorXd flatPanelCoeffs = fittingAlgorithms.fitFlatPanel();
    std::cout << "Panel equation: z = " << flatPanelCoeffs(0) << " + " << flatPanelCoeffs(1) << " * x + "
        << flatPanelCoeffs(2) << " * y " << std::endl;



    return 0;
}