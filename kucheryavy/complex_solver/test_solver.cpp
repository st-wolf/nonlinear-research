#include "functions.h"
#include "Solver.h"
using namespace NR;

int main()
{
    /* test RK4 complex solver */

    // Parameters
    ComplexNumber* param = new ComplexNumber(0, 1);

    // Preparing
    ComplexNumber* U0 = new ComplexNumber[2];
    U0[0] = ComplexNumber(1, 0);
    U0[1] = ComplexNumber(0, 0);

    int intCount = 1024;
    double xspan[2] = {0.0, 1.0};

    // Solving
    OdeProblem<2> OP(test_task, param, U0, xspan);
    Solver<2> Slv;
    Solution<2> Sol = Slv.Solve(OP, intCount);

    // Extracting
    const Point<2>* Data = Sol.GetData();
    int M = Sol.nodesCount;

    for (int i = 0; i < M; i++)
	{
        std::cout << Sol[i][0] << "\t" << Sol[i][1] << "\n";
	}

    return 0;
}
