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
    std::cout << "\nMesh surface area: " << std::setprecision(6) << std::fixed << meshArea << std::endl;

    std::cout << "\n\n";
    std::cout << "Fitting Algorithms" << std::endl;
    std::cout << "--------------------------------------------" << std::endl;
    std::cout << "\nFitting doubly curved shell (elliptic paraboloid):" << std::endl;
    Eigen::VectorXd doublyCurvedCoeffs = fittingAlgorithms.fitDoublyCurvedShell();
    std::cout << "Shell equation: z = " << doublyCurvedCoeffs(0) << " + " << doublyCurvedCoeffs(1) << " * x + "
		<< doublyCurvedCoeffs(2) << " * y + " << doublyCurvedCoeffs(3) << " * x^2 + " << doublyCurvedCoeffs(4) << " * x * y + "
		<< doublyCurvedCoeffs(5) << " * y^2 " << std::endl;

	std::cout << "\nFitting doubly curved shell (elliptic paraboloid) using a different approach:" << std::endl;
    Eigen::VectorXd coefficients = fittingAlgorithms.fitDoublyCurvedShell2();
    std::cout << "Fitted equation: z = " << coefficients(0) << " * x^2 + " << coefficients(1) 
              << " * x + " << coefficients(2) << " * y + " << coefficients(3) << " * y^2" << std::endl;

    std::cout << "\nFitting singly curved shell (cylindrical paraboloid) using a different approach:" << std::endl;
    Eigen::VectorXd coefitSinglyCurvedShell2X = fittingAlgorithms.fitSinglyCurvedShell(Axis::X);
    std::cout << "Fitted equation along X-axis: z = " << coefitSinglyCurvedShell2X(0) << " * x^2 + " << coefitSinglyCurvedShell2X(1) << " * y" << std::endl;
    Eigen::VectorXd coefitSinglyCurvedShell2Y = fittingAlgorithms.fitSinglyCurvedShell(Axis::Y);
    std::cout << "Fitted equation along Y-axis: z = " << coefitSinglyCurvedShell2Y(0) << " * y^2 + " << coefitSinglyCurvedShell2Y(1) << " * x" << std::endl;

    std::cout << "\nFitting flat panel surface:" << std::endl;
    Eigen::VectorXd flatPanelCoeffs = fittingAlgorithms.fitFlatPanel();
    std::cout << "Panel equation: z = " << flatPanelCoeffs(0) << " + " << flatPanelCoeffs(1) << " * x + "
        << flatPanelCoeffs(2) << " * y " << std::endl;

    std::cout << "\n\n";
    std::cout << "Coordinate Range and Center" << std::endl;
    std::cout << "--------------------------------------------" << std::endl;
    readParseData.calculateCoordinateRanges();
    readParseData.printCoordinateRanges();
    // Access ranges for later use
    const auto& ranges = readParseData.getCoordinateRanges();
    readParseData.calculateCenterCoordantes();

    std::cout << "\n\n";
    std::cout << "Characterize Fitted 2D Surface" << std::endl;
    std::cout << "--------------------------------------------" << std::endl;
    Characterize2DSurface surface(coefficients(0), coefficients(1), coefficients(2), coefficients(3));
    const auto& center = readParseData.getCenterCoordinates();
    double x = center.centerX, y = center.centerY; // The point at which to compute the curvature
    double R1, R2, L1, L2;
    surface.computeRadiiAndCurvatureLengths(x, y, R1, R2, L1, L2);
    std::cout << "Principal Radius of Curvature R1: " << R1 << std::endl;
    std::cout << "Principal Radius of Curvature R2: " << R2 << std::endl;
    std::cout << "Principal Curvature Length L1: " << L1 << std::endl;
    std::cout << "Principal Curvature Length L2: " << L2 << std::endl;
	double centerAngle1, centerAngle2;
	centerAngle1 = Characterize2DSurface::calculateCentralAngle(L1, R1);
	centerAngle2 = Characterize2DSurface::calculateCentralAngle(L2, R2);
	std::cout << "Central Angle 1: " << centerAngle1 << std::endl;
	std::cout << "Central Angle 2: " << centerAngle2 << std::endl;

    return 0;
}