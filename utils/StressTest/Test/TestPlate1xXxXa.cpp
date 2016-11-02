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

			//double timeStep = 1e-4;
			//double stiff = 1e8;
			//double ldamp = 1e3;
			//double adamp = 1e2;

			double timeStep = 1e-6;
			double stiff = 1e11;
			double ldamp = 1e5;
			double adamp = 1e2;

			TestFactory factory;
			SolverHandler _hsolver = factory
				.E(2.1e11)
				.Density(7900)
				.Damping(0.1f)
				.ScaleFactor(1.f)

				.IterationsCount(100)
				.SubIterationsCount(100)
				.TimeStep(timeStep)

				.GridStep(gridStep)
				.Dims(sideElements, sideElements, 1, SpecialSolversTest::StressStrainStuff::ECode::xlr)

				.Force(1000.f)
				.ForceDof(dof_x) // для пластины!!!
				.SolverType(solverType)
				
				.BuildQuarterPlate();

			
			if (_hsolver != nullptr)
			{
				OverrideStiffness(_hsolver,
					stiff,
					0.,
					ldamp,//1800.,
					adamp,//1000.,
					1.);
				OverrideInertia(_hsolver, 1., 1.);// *(gridStep*gridStep*0.25));
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