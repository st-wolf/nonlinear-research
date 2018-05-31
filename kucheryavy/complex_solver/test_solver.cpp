#include "functions.h"
#include "Solver.h"
using namespace NR;

int main()
{
    /* test RK4 complex solver */

    // Parameters
    // double* param = new double(1);
    double param = 1.0;

    // Preparing
    complex* U0 = new complex[2];
    U0[0] = complex(1, 0);
    U0[1] = complex(0, 0);

    int intCount = 1024;
    double xspan[2] = {0.0, 1.0};

    // Solving
    OdeProblem<2> OP(test_task, &param, U0, xspan);
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
