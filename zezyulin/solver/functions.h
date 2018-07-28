#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include <cmath>
#include "Point.h"
using namespace NR;

// Kolesnikova task.
// NLS with presence of quadratic linear potential and quadratic nonlinear potential
// Equation for nonlinear stationary modes
// u_{xx} + (\omega - x^2) u + (1 + \alpha x^2) u^3 = 0
// Parameters: [\omega \alpha]
Point<2> kolesnikova(double x, Point<2> u, double* param)
{
	Point<2> du;
	double omega = param[0], alpha = param[1];
	
	du[0] = u[1];
	du[1] = -(omega - pow(x, 2)) * u[0] + (1 + alpha * pow(x, 2)) * pow(u[0], 3);
	return du;
}

// Zezyulin task
Point<2> zezyulin(double x, Point<2> u, double* param)
{
    Point<2> du;
	double mu = param[0], omega = param[1], A = param[2];
	
	du[0] = u[1];
	du[1] = -(mu - 0.5 * pow(omega, 2) * pow(x, 2)) * u[0]
            -(1 + A * pow(std::tanh(x), 2)) * pow(u[0], 3);
	return du;
}

#endif