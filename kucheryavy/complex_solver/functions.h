#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include "Point.h"
using namespace NR;

Point<2> test_task(double x, Point<2> u, double* param)
{
	Point<2> du;

	// Parse parameters
	double factor = param[0];

	// Use them in equation
	du[0] = u[1];
	du[1] = u[0] * factor;
	return du;
}

/*
 * ODE problem for Kucheryavy equation:
 * u_{xx} - \omega u + (\alpha + \cos \phi \cos 2x + i * \sin \phi \sin 2x) |u|^2 u = 0 
 */
Point<2> kucheryavy(double x, Point<2> u, double* param)
{
	const complex imaginary_unit = complex(0, 1);

	Point<2> du;

	// Parse parameters
	double omega = param[0], alpha = param[1], phi = param[2];

	du[0] = u[1];

	// Use them in equation
	du[0] = u[1];
	du[1] = u[0] * omega - (
		alpha + cos(phi) * cos(2 * x) +
		imaginary_unit * sin(phi) * sin(2 * x)) * 
		(u[0].abs() * u[0].abs()) * u[0];
	return du;
}

#endif