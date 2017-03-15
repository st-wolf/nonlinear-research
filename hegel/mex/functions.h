#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include "Point.h"
using namespace NR;

// Lyuba Hegel task.
// NLS with presence of quadratic linear potential and cosine nonlinear potential
// Equation for nonlinear stationary modes
// u_{xx} + (\omega - x^2) u - \cos{2 \Omega x} u^3 = 0
// Parameters: [\omega \Omega]
Point<2> hegel(double x, Point<2> u, double* param)
{
	Point<2> du;
	double omega = param[0], Omega = param[1];
	
	du[0] = u[1];
	du[1] = -(omega - pow(x, 2)) * u[0] + (cos(2 * Omega * x)) * pow(u[0], 3);
	return du;
}

#endif