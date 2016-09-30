#include "CommonSolversTest.h"
#include "StressSolverTest.h"

#include "../Solvers/Stress/StressStrainCppSolver.h"
#include "../Solvers/Stress/StressStrainSolverExports.h"
#include "../Solvers/Stress/FTimer.h"

#include <sstream>
#include <vector>
#include "TestFactory.h"


using std::string;
using std::vector;
using std::stringstream;

namespace SpecialSolversTest
{

	namespace StressStrainStuff
	{

		void Test10x1x1(int solverType, ECode code)
		{
			TestFactory factory;
			SolverHandler _hsolver = factory
				.E(2.1e11f)
				.Density(7900.f)
				.Damping(2.f)
				.ScaleFactor(1e6f)

				.IterationsCount(500)
				.SubIterationsCount(1000)
				.TimeStep(0.00001f)

				.GridStep(0.01f)
				.Dims(10, 1, 1, code)

				.Force(1000)

				.SolverType(solverType)

				.BuildBeam();

			if (_hsolver != nullptr)
			{
				PerformanceCounter pc;
				pc.Start();
				Solve
					(
					_hsolver,
					factory.IntegrationParams()
					);
				pc.Print("Solving time: ", true);

				Stress::ReleaseMemory((void* &)_hsolver);
			}
		}
	}
}