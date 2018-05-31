#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include "Point.h"
using namespace NR;

Point<2> test_task(double x, Point<2> u, ComplexNumber* param)
{
	Point<2> du;

	// Parse complex parameters
	ComplexNumber factor = param[0];

	// Use them in equation
	du[0] = u[1];
	du[1] = factor * u[0];
	return du;
}

Point<2> kucheryavy(double x, Point<2> u, ComplexNumber* param)
{
	Point<2> du;
	// Parse compel parameters
	// double mu = param[0], Omega = param[1], P1 = param[2];

	du[0] = u[1];

	// Use them in equation
	// du[1] = -(mu - pow(x, 2)) * u[0] - (1 + P1*cos(Omega * x)) * pow(u[0], 3);
	return du;
}

#endif