#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include "Point.h"
using namespace NR;

// Kolesnikova task.
// NLS with presence of quadratic linear potential and quadratic nonlinear potential
// Equation for nonlinear stationary modes
// u_{xx} + (\omega - x^2) u + (1 + \alpha x^2) u^3 = 0
// Parameters: [\omega \alpha]
Point<2> hegel(double x, Point<2> u, double* param)
{
	Point<2> du;
	double omega = param[0], alpha = param[1];
	
	du[0] = u[1];
	du[1] = -(omega - pow(x, 2)) * u[0] - (1 + alpha * pow(x, 2)) * pow(u[0], 3);
	return du;
}

#endif