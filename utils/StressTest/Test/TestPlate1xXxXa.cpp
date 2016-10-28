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
		void Test1xXxXa(int solverType, float plateSideLength, float plateWidth, int sideElements)
		{
			float gridStep = (float)(plateSideLength / sideElements);
			float stressScalingFactor = (float)(gridStep/plateWidth);

			TestFactory factory;
			SolverHandler _hsolver = factory
				.E(1e6)
				.Density(7900)
				.Damping(0.1f)
				.ScaleFactor(1.f)

				.IterationsCount(500)
				.SubIterationsCount(100)
				.TimeStep(0.0003f)

				.GridStep(gridStep)
				.Dims(sideElements, sideElements, 1, SpecialSolversTest::StressStrainStuff::ECode::xlr)

				.Force(1000.f)
				.ForceDof(dof_x) // для пластины!!!
				.SolverType(solverType)
				
				.BuildQuarterPlate();

			
			if (_hsolver != nullptr)
			{
				OverrideStiffness(_hsolver,
					1e7,
					1.,
					1800.,
					1000.,
					1.);
				OverrideInertia(_hsolver, 1., 1.);
				OverrideScalingFactors(
					_hsolver,
					stressScalingFactor,
					stressScalingFactor,
					1);

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