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
			int sideElements = 51;
			double plateSideLength = 0.8;
			double plateWidth = 0.01;
			float gridStep = (float)(plateSideLength / sideElements);
			float stressScalingFactor = (float)(gridStep / plateWidth);

			double timeStep = 1e-4;
			double stiff = 1e6;
			double stiffa = 1e3;
			double ldamp = 1e3;
			double adamp = 1e2;

			//double timeStep = 1e-6;
			//double stiff = 1e11;
			//double stiffa = 1e5;
			//double ldamp = 1e5;
			//double adamp = 1e2;

			TestFactory factory;
			SolverHandler _hsolver = factory
				.E(2.1e12f)
				.Density(7900.f)
				.Damping(1.f)
				.ScaleFactor(1.f)

				.IterationsCount(1000)
				.SubIterationsCount(100)
				.TimeStep(timeStep)

				.GridStep(gridStep)
				.Dims(1, sideElements, sideElements, code)

				.Force(1)
				.ForceDof(dof_x) // для пластины!!!
				.SolverType(solverType)

				.BuildPlate();

			if (_hsolver != nullptr)
			{
				OverrideStiffness(_hsolver,
					stiff,
					stiffa,
					ldamp,//1800.,
					adamp,//1000.,
					1.);
				OverrideInertia(_hsolver, 1., 1.);// *(gridStep*gridStep*0.25));
				OverrideScalingFactors(
					_hsolver,
					1,//stressScalingFactor,
					1,//stressScalingFactor,
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