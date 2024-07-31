#include "ReadParseData.h"
#include "characterize2DSurface.h"
#include "FittingAlgorithms.h"
#include <iostream>
#include <iomanip>


int main() {

    // ReadParseData
    ReadParseData readParseData;

    try {
        readParseData.readFromFile("data/openCylinder.fem");
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

    std::cout << "\nFitting doubly curved shell (elliptic paraboloid):" << std::endl;
    Eigen::VectorXd doublyCurvedCoeffs = fittingAlgorithms.fitDoublyCurvedShell();
    std::cout << "Shell equation: z = " << doublyCurvedCoeffs(0) << " + " << doublyCurvedCoeffs(1) << " * x + "
		<< doublyCurvedCoeffs(2) << " * y + " << doublyCurvedCoeffs(3) << " * x^2 + " << doublyCurvedCoeffs(4) << " * x * y + "
		<< doublyCurvedCoeffs(5) << " * y^2 " << std::endl;

    std::cout << "\nFitting singly curved shell (cylindrical paraboloid):" << std::endl;
    Eigen::VectorXd singlyCurvedCoeffs = fittingAlgorithms.fitCylindricalParaboloid();
    std::cout << "Shell equation: z = " << singlyCurvedCoeffs(0) << " * (x - " << singlyCurvedCoeffs(1) << ")^2 + "
        << singlyCurvedCoeffs(2) << " * y + " << singlyCurvedCoeffs(3) << std::endl;

    std::cout << "\nFitting flat panel surface:" << std::endl;
    Eigen::VectorXd flatPanelCoeffs = fittingAlgorithms.fitFlatPanel();
    std::cout << "Panel equation: z = " << flatPanelCoeffs(0) << " + " << flatPanelCoeffs(1) << " * x + "
        << flatPanelCoeffs(2) << " * y " << std::endl;

    readParseData.calculateCoordinateRanges();
    readParseData.printCoordinateRanges();

    // Access ranges for later use
    const auto& ranges = readParseData.getCoordinateRanges();
    double x_min = ranges.x.min;
    double x_max = ranges.x.max;
    double y_min = ranges.y.min;
    double y_max = ranges.y.max;
    double z_min = ranges.z.min;
    double z_max = ranges.z.max;
    // Similarly for y and z


    std::cout << "\n\n";
    std::cout << "Characterize Fitted 2D Surface" << std::endl;
    std::cout << "--------------------------------------------" << std::endl;

    // Characterize2DSurface
    std::pair<double, double> x_range = { x_min, x_max };
    std::pair<double, double> y_range = { y_min, y_max };
    int num_x = 50;
    int num_y = 50;

    Characterize2DSurface characterizer(x_range, y_range, num_x, num_y);

    // Quadratic surface: Z = c + a1*X + a2*Y + a3*X^2 + a4*X*Y + a5*Y^2
    Characterize2DSurface::QuadraticSurface surface1(
        doublyCurvedCoeffs(0), doublyCurvedCoeffs(1), doublyCurvedCoeffs(2),
        doublyCurvedCoeffs(3), doublyCurvedCoeffs(4), doublyCurvedCoeffs(5));

    // Singly curved shell: Z = c1*(X - c2)^2 + c3*Y + c4
    Characterize2DSurface::SinglyCurvedShell surface2(
        singlyCurvedCoeffs(0), singlyCurvedCoeffs(1), 
        singlyCurvedCoeffs(2), singlyCurvedCoeffs(3));

    auto props1 = characterizer.characterize(surface1);
    auto props2 = characterizer.characterize(surface2);

    std::cout << "Quadratic Surface properties:" << std::endl;
    Characterize2DSurface::print_properties(props1);

    std::cout << "\nSingly Curved Shell properties:" << std::endl;
    Characterize2DSurface::print_properties(props2);



    return 0;
}