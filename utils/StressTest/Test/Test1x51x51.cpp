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
		void Test1x51x51(int solverType, ECode code)
		{
			TestFactory factory;
			SolverHandler _hsolver = factory
				.E(2.1e12f)
				.Density(7900.f)
				.Damping(1.f)
				.ScaleFactor(1e9f)

				.IterationsCount(200)
				.SubIterationsCount(1000)
				.TimeStep(0.0001f)

				.GridStep(0.01f)
				.Dims(1, 51, 51, code)

				.Force(1000)
				.ForceDof(dof_x) // для пластины!!!
				.SolverType(solverType)

				.BuildPlate();

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