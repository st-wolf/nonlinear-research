#ifndef ODEPROBLEM_H_
#define ODEPROBLEM_H_

namespace NR
{

template <int dim>
class Problem
{
	ComplexNumber* p;
	Point<dim>(*f)(double t, Point<dim> y, ComplexNumber* param);
public:
	Problem(Point<dim>(*function)(double, Point<dim>, ComplexNumber*), ComplexNumber* param)
	{
		f = function;
		p = param;
	}

	Point<dim> operator() (double t, Point<dim> y)
	{
		return f(t, y, p);
	}
};

template <int dim>
class OdeProblem
{
public:
	Problem<dim> f;
	Point<dim> Init;
	double tspan[2]; //!!!
	//double prec;
	//int intCount;

	OdeProblem(Point<dim>(*function)(double, Point<dim>, ComplexNumber*),
				ComplexNumber* param, Point<dim> i_Init, double* i_tspan)
		:f(function, param)
	{
		Init = i_Init;
		tspan[0] = i_tspan[0];
		tspan[1] = i_tspan[1];
	}
};

}

#endif