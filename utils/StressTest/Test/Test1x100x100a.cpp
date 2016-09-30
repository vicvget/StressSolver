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
		void Test1x100x100a(int solverType, ECode code)
		{
			TestFactory factory;
			SolverHandler _hsolver = factory
				.E(2.1e12f)
				.Density(7900.f)
				.Damping(1.f)
				.ScaleFactor(1e11f)

				.IterationsCount(1000)
				.SubIterationsCount(1)
				.TimeStep(0.001f)

				.GridStep(0.01f)
				.Dims(100, 100, 1, code)

				.Force(1)
				.ForceDof(dof_x) // для пластины!!!
				.SolverType(solverType)

				.BuildQuarterPlate();

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