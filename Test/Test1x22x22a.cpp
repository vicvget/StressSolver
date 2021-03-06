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
			int side = 22;
			float sideLength = 0.8;
			TestFactory factory;
			float gridStep = sideLength / side;
			SolverHandler _hsolver = factory
				.E(1e6)
				.Density(7900)
				.Damping(900.f)
				.ScaleFactor(1.f)

				.IterationsCount(200)
				.SubIterationsCount(100)
				.TimeStep(0.0003f)

				.GridStep(sideLength/side)
				.Dims(side, side, 1, code)

				.Force(1000.f)
				.ForceDof(dof_x) // ��� ��������!!!
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
				double factor = gridStep / 0.01;
				OverrideScalingFactors(_hsolver, factor, factor, factor);
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