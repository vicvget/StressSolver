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
				.E(1e6)
				.Density(1953125)
				.Damping(0.1f)
				.ScaleFactor(1e5)

				.IterationsCount(1000)
				.SubIterationsCount(1)
				.TimeStep(0.00001f)

				.GridStep(0.008f)
				.Dims(100, 100, 1, code)

				.Force(1.f)
				.ForceDof(dof_x) // для пластины!!!
				.SolverType(solverType)
				
				.BuildQuarterPlate();

			
			if (_hsolver != nullptr)
			{
				OverrideStiffness(_hsolver,
					1000000.,
					1000000.,
					900.,
					900.,
					1.);
				OverrideInertia(_hsolver, 1., 1.);

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