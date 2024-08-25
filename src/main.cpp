#include "ReadParseData.h"
#include "characterize2DSurface.h"
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
    std::cout << "Fitting Algorithms: " << std::endl;
    std::cout << "--------------------------------------------" << std::endl;
    Eigen::VectorXd coefficients = fittingAlgorithms.fitDoublyCurvedShell2();
    std::cout << "Fitted equation: z = " << coefficients(0) << " * x^2 + " << coefficients(1) 
              << " * x + " << coefficients(2) << " * y + " << coefficients(3) << " * y^2" << std::endl;
    std::cout << "Characterize Fitted 2D Surface: " << std::endl;
    std::cout << "--------------------------------------------" << std::endl;
    Characterize2DSurface surface(coefficients(0), coefficients(1), coefficients(2), coefficients(3));
    const auto& center = readParseData.getCenterCoordinates();
    double R1X, R2Y, L1X, L2Y;
    // For fixed y, varying x
    double y_fixed = center.centerY;
    double x_start = ranges.x.min;
    double x_end = ranges.x.max;
    surface.calculate_arc_length(x_start, x_end, y_fixed, true, L1X);
    // For fixed x, varying y
    double x_fixed = center.centerX;
    double y_start = ranges.y.min;
    double y_end = ranges.y.max;
    surface.calculate_arc_length(y_start, y_end, x_fixed, false, L2Y);
    surface.computeRadii(center.centerX, center.centerY, R1X, R2Y);
    std::cout << "Principal Radius of Curvature R1: " << R1X << std::endl;
    std::cout << "Principal Radius of Curvature R2: " << R2Y << std::endl;
    std::cout << "Principal Curvature Length L1: " << L1X << std::endl;
    std::cout << "Principal Curvature Length L2: " << L2Y << std::endl;
    double centerAngle1, centerAngle2;
    centerAngle1 = Characterize2DSurface::calculateCentralAngle(L1X, R1X);
    centerAngle2 = Characterize2DSurface::calculateCentralAngle(L2Y, R2Y);
    std::cout << "Central Angle 1: " << centerAngle1 << std::endl;
    std::cout << "Central Angle 2: " << centerAngle2 << std::endl;

    std::cout << "\n\n";
    std::cout << "2) Fitting singly curved shell (cylindrical paraboloid) using a NEW approach:" << std::endl;
    Eigen::VectorXd coefitSinglyCurvedShell2X = fittingAlgorithms.fitSinglyCurvedShell(Axis::X);
    std::cout << "Fitted equation along X-axis: z = " << coefitSinglyCurvedShell2X(0) << " * x^2 + " << coefitSinglyCurvedShell2X(1) << " * y" << std::endl;
    Eigen::VectorXd coefitSinglyCurvedShell2Y = fittingAlgorithms.fitSinglyCurvedShell(Axis::Y);
    std::cout << "Fitted equation along Y-axis: z = " << coefitSinglyCurvedShell2Y(0) << " * y^2 + " << coefitSinglyCurvedShell2Y(1) << " * x" << std::endl;

    std::cout << "\n\n";
    std::cout << "3) Fitting flat panel surface:" << std::endl;
    Eigen::VectorXd flatPanelCoeffs = fittingAlgorithms.fitFlatPanel();
    std::cout << "Panel equation: z = " << flatPanelCoeffs(0) << " + " << flatPanelCoeffs(1) << " * x + "
        << flatPanelCoeffs(2) << " * y " << std::endl;



    return 0;
}