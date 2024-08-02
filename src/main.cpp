#include "ReadParseData.h"
#include "characterize2DSurface.h"
#include "FittingAlgorithms.h"
#include <iostream>
#include <iomanip>


// Helper function to calculate z value from the fitted equation
double calculateZ(const Eigen::VectorXd& params, double x, double y, const Eigen::Vector3d& center, double scale) {
    double x_scaled = (x - center[0]) * scale;
    double y_scaled = (y - center[1]) * scale;
    double z_scaled = params[0] * x_scaled * x_scaled + params[1] * y_scaled * y_scaled +
        params[2] * x_scaled + params[3] * y_scaled + params[4];
    return z_scaled / scale + center[2];
}

int main() {

  //  // ReadParseData
  //  ReadParseData readParseData;

  //  try {
  //      readParseData.readFromFile("data/testGeo.fem");
  //  }
  //  catch (const std::exception& e) {
  //      std::cerr << "Error: " << e.what() << std::endl;
  //      return 1;
  //  }

  //  // FittingAlgorithms
  //  FittingAlgorithms fittingAlgorithms(readParseData);

  //  //readParseData.printAllGridData();
  //  //readParseData.printAllMeshElements();

  //  std::cout << "Read and Analyze FEM Data" << std::endl;
  //  std::cout << "--------------------------------------------" << std::endl;

  //  double meshArea = readParseData.calculateSurfaceArea();
  //  std::cout << "\nMesh surface area: " << std::setprecision(6) << std::fixed << meshArea << std::endl;

  //  std::cout << "\nFitting doubly curved shell (elliptic paraboloid):" << std::endl;
  //  Eigen::VectorXd doublyCurvedCoeffs = fittingAlgorithms.fitDoublyCurvedShell();
  //  std::cout << "Shell equation: z = " << doublyCurvedCoeffs(0) << " + " << doublyCurvedCoeffs(1) << " * x + "
		//<< doublyCurvedCoeffs(2) << " * y + " << doublyCurvedCoeffs(3) << " * x^2 + " << doublyCurvedCoeffs(4) << " * x * y + "
		//<< doublyCurvedCoeffs(5) << " * y^2 " << std::endl;

  //  std::cout << "\nFitting singly curved shell (cylindrical paraboloid):" << std::endl;
  //  Eigen::VectorXd singlyCurvedCoeffs = fittingAlgorithms.fitSinglyCurvedShell(Axis::Y);
  //  std::cout << "Shell equation: z = " << singlyCurvedCoeffs(0) << " * (x - " << singlyCurvedCoeffs(1) << ")^2 + "
  //      << singlyCurvedCoeffs(2) << " * y + " << singlyCurvedCoeffs(3) << std::endl;

  //  std::cout << "\nFitting closed-form singly curved shell (cylinder):" << std::endl;
  //  Eigen::VectorXd closedSinglyCurvedCoeffs = fittingAlgorithms.fitCloseSinglyCurvedShell(Axis::Y);
  //  std::cout << "Cylinder equation: (x - " << closedSinglyCurvedCoeffs(0) << ")^2 + (y - "
  //      << closedSinglyCurvedCoeffs(1) << ")^2 = " << closedSinglyCurvedCoeffs(2) << "^2" << std::endl;

  //  std::cout << "\nFitting flat panel surface:" << std::endl;
  //  Eigen::VectorXd flatPanelCoeffs = fittingAlgorithms.fitFlatPanel();
  //  std::cout << "Panel equation: z = " << flatPanelCoeffs(0) << " + " << flatPanelCoeffs(1) << " * x + "
  //      << flatPanelCoeffs(2) << " * y " << std::endl;

  //  readParseData.calculateCoordinateRanges();
  //  readParseData.printCoordinateRanges();

  //  // Access ranges for later use
  //  const auto& ranges = readParseData.getCoordinateRanges();
  //  double x_min = ranges.x.min;
  //  double x_max = ranges.x.max;
  //  double y_min = ranges.y.min;
  //  double y_max = ranges.y.max;
  //  double z_min = ranges.z.min;
  //  double z_max = ranges.z.max;
  //  // Similarly for y and z


  //  std::cout << "\n\n";
  //  std::cout << "Characterize Fitted 2D Surface" << std::endl;
  //  std::cout << "--------------------------------------------" << std::endl;

  //  // Characterize2DSurface
  //  std::pair<double, double> x_range = { x_min, x_max };
  //  std::pair<double, double> y_range = { y_min, y_max };
  //  int num_x = 50;
  //  int num_y = 50;

  //  Characterize2DSurface characterizer(x_range, y_range, num_x, num_y);

  //  // Create surfaces
  //  Characterize2DSurface::QuadraticSurface surface1(
  //      doublyCurvedCoeffs(0), doublyCurvedCoeffs(1), doublyCurvedCoeffs(2),
  //      doublyCurvedCoeffs(3), doublyCurvedCoeffs(4), doublyCurvedCoeffs(5));

  //  Characterize2DSurface::SinglyCurvedShell surface2(
  //      singlyCurvedCoeffs(0), singlyCurvedCoeffs(1),
  //      singlyCurvedCoeffs(2), singlyCurvedCoeffs(3), true);

  //  // Characterize the surfaces
  //  auto props1 = characterizer.characterize(surface1);
  //  auto props2 = characterizer.characterize(surface2);

  //  // Print results
  //  std::cout << "Quadratic Surface properties:" << std::endl;
  //  Characterize2DSurface::print_properties(props1);

  //  std::cout << "\nSingly Curved Shell properties:" << std::endl;
  //  Characterize2DSurface::print_properties(props2);


    ReadParseData data;
    data.readFromFile("data/testGeo.fem");

    FittingAlgorithms fitter(data);
    fitter.fitSynclasticShell();

    Eigen::VectorXd params = fitter.getFittedParameters();
    Eigen::Vector3d center = fitter.getDataCenter();
    double scale = fitter.getDataScale();
    std::pair<double, double> radii = fitter.getRadiiOfCurvature();
    double error = fitter.calculateFitError();

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Data center: " << center.transpose() << std::endl;
    std::cout << "Data scale: " << scale << std::endl;

    std::cout << "\nFitted parameters [a, b, c, d, e] (after centering and scaling):\n"
        << params.transpose() << std::endl;

    std::cout << "Radii of curvature: R1 = " << radii.first << ", R2 = " << radii.second << std::endl;
    std::cout << "Fit error: " << error << std::endl;

    std::cout << "\nFitted surface equation (after centering and scaling):" << std::endl;
    std::cout << "z = " << params[0] << "x^2 + " << params[1] << "y^2 + "
        << params[2] << "x + " << params[3] << "y + " << params[4] << std::endl;

    std::cout << "\nFitted surface equation (original coordinates):" << std::endl;
    std::cout << "z = " << params[0] * scale * scale << "(x - " << center[0] << ")^2 + "
        << params[1] * scale * scale << "(y - " << center[1] << ")^2 + "
        << params[2] * scale << "(x - " << center[0] << ") + "
        << params[3] * scale << "(y - " << center[1] << ") + "
        << (params[4] / scale + center[2]) << std::endl;

    // Calculate range for the fitted equation
    double x_min = std::numeric_limits<double>::max();
    double x_max = std::numeric_limits<double>::lowest();
    double y_min = std::numeric_limits<double>::max();
    double y_max = std::numeric_limits<double>::lowest();

    const std::map<int, Vector3D>& gridData = data.getGridData();
    for (const auto& point : gridData) {
        x_min = std::min(x_min, point.second.x);
        x_max = std::max(x_max, point.second.x);
        y_min = std::min(y_min, point.second.y);
        y_max = std::max(y_max, point.second.y);
    }

    double z_min = std::min({
        calculateZ(params, x_min, y_min, center, scale),
        calculateZ(params, x_min, y_max, center, scale),
        calculateZ(params, x_max, y_min, center, scale),
        calculateZ(params, x_max, y_max, center, scale)
        });

    double z_max = std::max({
        calculateZ(params, x_min, y_min, center, scale),
        calculateZ(params, x_min, y_max, center, scale),
        calculateZ(params, x_max, y_min, center, scale),
        calculateZ(params, x_max, y_max, center, scale)
        });

    std::cout << "\nRange for the fitted equation:" << std::endl;
    std::cout << "x range: [" << x_min << ", " << x_max << "]" << std::endl;
    std::cout << "y range: [" << y_min << ", " << y_max << "]" << std::endl;
    std::cout << "z range: [" << z_min << ", " << z_max << "]" << std::endl;


    return 0;
}