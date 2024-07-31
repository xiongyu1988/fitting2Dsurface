#ifndef FittingAlgorithms_H
#define FittingAlgorithms_H

#include <Eigen/Dense>
#include "ReadParseData.h"

enum class Axis {
	X,
	Y
};

class FittingAlgorithms {
public:

	FittingAlgorithms(const ReadParseData& data) : readParseData(data) {}

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
	/// This is now linear in terms of the parameters p, q, c, and r.
	/// We can solve this using simple linear least squares.
	/// </summary>
	/// <returns></returns>
	Eigen::VectorXd fitSinglyCurvedShell(Axis axis) const;

	/// <summary>
	/// A flat panel can be described by a plane equation of the form:
	/// z = c1 + c2 * x + c3 * y
	/// </summary>
	/// <returns></returns>
	Eigen::VectorXd fitFlatPanel() const;

	/// <summary>
	/// The general cylinder equation with the axis aligned with the z-axis
	/// is given by:
	/// (x - c1)^2 + (y - c2)^2 = r^2
	/// </summary>
	Eigen::VectorXd fitCloseSinglyCurvedShell(Axis axis) const;

private:
	const ReadParseData& readParseData;
};

#endif // FittingAlgorithms_H
