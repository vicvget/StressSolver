#include "CommonSolversTest.h"
#include "StressSolverTest.h"

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

		/** Stress solver testing for grid
		* @param solverType - type of solver
		* solverType = 0 - default CPP
		* solverType = 1 - AVX
		* solverType = 2 - FMA
		* solverType = 3 - KNC
		* solverType = 4 - Unaligned
		*/
		void Test1x1x3(int solverType, ECode code)
		{
			TestFactory factory;
			SolverHandler _hsolver = factory
				.E(100.f)
				.Density(1000.f)
				.Damping(1.f)
				.ScaleFactor(1.f)
				.IterationsCount(100)
				.SubIterationsCount(1000)
				.TimeStep(0.00001f)
				.GridStep(0.1f)
				.Force(10)
				.SolverType(solverType)
				.Dims(3, 1, 1, code)
				.Build();
	
			if (_hsolver != nullptr)
			{
				// “ест дл€ сравнени€ с эталонной моделью
				OverrideStiffness(_hsolver,
					1000.,
					10.,
					10.,
					10.,
					1.);

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