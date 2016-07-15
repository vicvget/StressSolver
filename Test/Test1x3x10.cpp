#include "CommonSolversTest.h"
#include "StressSolverTest.h"

#include "../Solvers/Stress/StressStrainCppSolver.h"
#include "../Solvers/Stress/StressStrainSolverExports.h"
#include "../Solvers/Stress/FTimer.h"

#include <sstream>
#include <vector>


using std::string;
using std::vector;
using std::stringstream;

namespace SpecialSolversTest
{

	namespace StressStrainStuff
	{
		/** Stress solver testing for grid 1x3x10
		* @param solverType - type of solver
		* solverType = 0 - default, fortran translated solver
		* solverType = 1 - static old
		* solverType = 2 - static 2 solver (mkl)
		* solverType = other - cpp, rewritten solver
		*/
		void Test1x3x10(int solverType)
		{
			stringstream str;
			str << "test_stress_type_" << solverType;

			string fileRlc = str.str() + ".rlc";
			string fileResult = str.str() + ".mpr";

			SpecialParams specialParams;
			IntegrationParams integrationParams;
			GridParams gridParams;

			specialParams._scaleFactor = 1e6;
			//specialParams._dampingRatio = 1e-4;

			//integrationParams._nIterations = 120;// 600;
			//integrationParams._nSubIterations = 60;//300;
			integrationParams._nIterations = 600;
			integrationParams._nSubIterations = 300;
			integrationParams._timeStep = 0.00005f;


			gridParams._nz = 1;
			gridParams._ny = 3;
			gridParams._nx = 10;
			gridParams._gridStep = 0.005;

			SolverHandler _hsolver = MakeSolver
				(
				gridParams,
				specialParams,
				integrationParams,
				fileRlc,
				solverType
				);
			if (_hsolver == NULL)
				return;
			PerformanceCounter pc;
			pc.Start();
			Solve
				(
				_hsolver,
				fileResult,
				gridParams,
				integrationParams
				);
			pc.Print("Solving time: ", true);
			Stress::ReleaseMemory((void* &)_hsolver);
		}
	}
}