#include "mex.h"
#include "functions.h"
#include "Solver.h"
using namespace NR;

// MATLAB CALL:
// [xgrid, solution] =  f_solve(params, xspan, u0, intCount)

// INPUT:
// :params:   [\omega \alpha \phi]
// :xspan:	  spatial interval of solution
// :u0:		  initial values for Cauchy problem
// :intCount: number of intervals

// OUTPUT:
// :xgrid:	   massive of spatial coordinates
// :solutions: n-dimentioonal massive of solution values on xgrid 

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *i_param = mxGetPr(prhs[0]), *i_xspan = mxGetPr(prhs[1]),
		*i_U0_real = mxGetPr(prhs[2]), *i_U0_imag = new double[2];

	// Validation
	if (!mxIsComplex(prhs[2]))
	{
		// Set zero imaginary part
		i_U0_imag[0] = 0;
		i_U0_imag[1] = 0;
	} else {
		i_U0_imag = mxGetPi(prhs[2]);
	}
	
	const double defIntCount = 1024;
	const double* pIntCount;

	if (nrhs == 4)
		pIntCount = mxGetPr(prhs[3]);
	else
		pIntCount = &defIntCount;

	// Create initial complex point
	complex *i_U0 = new complex[2];
	i_U0[0] = complex(i_U0_real[0], i_U0_imag[0]);
	i_U0[1] = complex(i_U0_real[1], i_U0_imag[1]);

	// Preparing
	Point<2> U0(i_U0);
	int intCount = (int)(*pIntCount);

	// Solving
	OdeProblem<2> OP(kucheryavy, i_param, U0, i_xspan);
	Solver<2> Slv;
	Solution<2> Sol = Slv.Solve(OP, intCount);

	// Packing
	int M = Sol.nodesCount;
	plhs[0] = mxCreateDoubleMatrix(M, 1, mxREAL); // xgrid
	plhs[1] = mxCreateDoubleMatrix(M, 2, mxCOMPLEX); // solution
	plhs[2] = mxCreateDoubleScalar((double)Sol.isInfinit);

	double *X = mxGetPr(plhs[0]),
		*U_real = mxGetPr(plhs[1]), *U_imag = mxGetPi(plhs[1]);

	for (int i = 0; i < M; i++)
	{
		X[i] = Sol.tspan[0] + i*Sol.step;
		
		U_real[i] = Sol[i][0].real;
		U_real[i + M] = Sol[i][1].real;

		U_imag[i] = Sol[i][0].imag;
		U_imag[i + M] = Sol[i][1].imag;
	}
	return;
}
