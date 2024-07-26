#include "readParseData.h"
#include "characterize2DSurface.h"



int main() {

    // ReadParseData
    ReadParseData analyzer;

    try {
        analyzer.readFromFile("data/geo1.fem");
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    //analyzer.printAllGridData();
    //analyzer.printAllMeshElements();

    std::cout << "Read and Analyze FEM Data" << std::endl;
    std::cout << "--------------------------------------------" << std::endl;

    double meshArea = analyzer.calculateSurfaceArea();
    std::cout << "\nMesh surface area: " << std::setprecision(6) << std::fixed << meshArea << std::endl;

    std::cout << "\nFitting doubly curved shell (elliptic paraboloid):" << std::endl;
    Eigen::VectorXd doublyCurvedCoeffs = analyzer.fitDoublyCurvedShell();
    std::cout << "Shell coefficients: " << doublyCurvedCoeffs.transpose() << std::endl;

    std::cout << "\nFitting singly curved shell (cylindrical paraboloid):" << std::endl;
    Eigen::VectorXd singlyCurvedCoeffs = analyzer.fitCylindricalParaboloid();
    std::cout << "Shell coefficients: " << singlyCurvedCoeffs.transpose() << std::endl;

    analyzer.calculateCoordinateRanges();
    analyzer.printCoordinateRanges();

    // Access ranges for later use
    const auto& ranges = analyzer.getCoordinateRanges();
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
    Characterize2DSurface::QuadraticSurface surface1(64.49757509292367, -26.730064479096814, 0.02179036522711368,
        2.8019565720091357, 0.0014880729562687031, 1.283908862065956);

    // Singly curved shell: Z = a*(X - b)^2 + c
    Characterize2DSurface::SinglyCurvedShell surface2(1.051454, 4.790226, 0.880074);

    auto props1 = characterizer.characterize(surface1);
    auto props2 = characterizer.characterize(surface2);

    std::cout << "Quadratic Surface properties:" << std::endl;
    Characterize2DSurface::print_properties(props1);

    std::cout << "\nSingly Curved Shell properties:" << std::endl;
    Characterize2DSurface::print_properties(props2);



    return 0;
}