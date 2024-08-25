#ifndef FITTINGALGORITHMS_H
#define FITTINGALGORITHMS_H

#include <Eigen/Dense>
#include "ReadParseData.h"

enum class Axis {
    X,
    Y
};

class FittingAlgorithms {
public:
    // Constructor
    explicit FittingAlgorithms(const ReadParseData& data);

    /**
     * @brief Fits a doubly curved shell (elliptic paraboloid).
     * The equation of an elliptic paraboloid is given by:
     * z = c1 + c2 * x + c3 * y + c4 * x^2 + c5 * x * y + c6 * y^2
     *
     * @return Eigen::VectorXd - The coefficients [c1, c2, c3, c4, c5, c6] for the fitted surface.
     */
    Eigen::VectorXd fitDoublyCurvedShell() const;

	/**
	 * @brief Fits a doubly curved shell (elliptic paraboloid) using a different approach.
	 * The equation of an elliptic paraboloid is given by:
	 * z = c1 * x^2 + c2 * x + c3 * y + c4 * y^2
	 *
	 * This function uses a simplified approach to transform a non-linear least squares fitting problem
	 * into a linear one by expanding the equation into:
	 * z = p * x^2 + q * x + c * y + r * y^2,
	 * where p and q are the parameters to be solved.
	 *
	 * @return Eigen::VectorXd - The coefficients [p, q, c, r] for the fitted surface.
	 */
	Eigen::VectorXd fitDoublyCurvedShell2() const;

    /**
	* @brief Fits a singly curved shell (cylindrical paraboloid) using a different approach.
	* The equation of a cylindrical paraboloid is given by:
	* z = c1 * x^2 + c2 * y.
    */
    Eigen::VectorXd fitSinglyCurvedShell(Axis axis) const;

    /**
     * @brief Fits a flat panel described by a plane equation.
     * The equation of a flat panel (plane) is given by:
     * z = c1 + c2 * x + c3 * y.
     *
     * @return Eigen::VectorXd - The coefficients [c1, c2, c3] for the fitted plane.
     */
    Eigen::VectorXd fitFlatPanel() const;

private:
    const ReadParseData& readParseData;
};

#endif // FITTINGALGORITHMS_H
