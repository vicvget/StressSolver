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
		void Test1x22x22a(int solverType, ECode code)
		{
			TestFactory factory;
			SolverHandler _hsolver = factory
				.E(1e6)
				.Density(7900)
				.Damping(900.f)
				.ScaleFactor(1.f)

				.IterationsCount(20000)
				.SubIterationsCount(1)
				.TimeStep(0.0002f)

				.GridStep(0.036f)
				.Dims(22, 22, 1, code)

				.Force(1000.f)
				.ForceDof(dof_x) // для пластины!!!
				.SolverType(solverType)
				
				.BuildQuarterPlate();

			
			if (_hsolver != nullptr)
			{
				OverrideStiffness(_hsolver,
					1e6,
					1.,
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