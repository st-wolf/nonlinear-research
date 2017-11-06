#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include "Point.h"
using namespace NR;

// Lyuba Hegel task.
// NLS with presence of quadratic linear potential and cosine nonlinear potential
// Equation for nonlinear stationary modes
// u_{xx} + (\mu - x^2) u + (1 + P1 \cos{2 \Omega x}) u^3 = 0
// Parameters: [\mu \Omega P1]
Point<2> hegel(double x, Point<2> u, double* param)
{
	Point<2> du;
	double mu = param[0], Omega = param[1], P1 = param[2];
	
	du[0] = u[1];
	du[1] = -(mu - pow(x, 2)) * u[0] - (1 + P1*cos(Omega * x)) * pow(u[0], 3);
	return du;
}

// Lyuba Hegel task.
// We study approximation for the equation:
// u_{xx} + (\mu - x^2) u + \sigma_1 \cos(\Omega x) u^3 = 0
Point<2> hegel_sigma(double x, Point<2> u, double* param)
{
	Point<2> du;
	double mu = param[0], Omega = param[1], sigma_1 = param[2];
	
	du[0] = u[1];
	du[1] = -(mu - pow(x, 2)) * u[0] - sigma_1 * cos(Omega * x) * pow(u[0], 3);
	return du;
}

// Lyuba Hegel approximation task.
// We study approximation for the equation:
// u_{xx} + (\mu - x^2) u + \sigma_1 \cos(\Omega x) u^3 = 0
// Now we need to solve: U_{xx} + (\mu - x^2) U + (3/2) U^5 = 0
Point<2> hegel_approx(double x, Point<2> u, double* param)
{
	Point<2> du;
	double mu = param[0], Omega = param[1], sigma_1 = param[2];
	
	du[0] = u[1];
	du[1] = -(mu - pow(x, 2)) * u[0] - (3.0 / 2.0) * pow(u[0], 5);
	return du;
}

#endif