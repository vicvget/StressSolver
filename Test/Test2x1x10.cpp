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

		/** Stress solver testing for grid 3x3x10
		* @param solverType - type of solver
		* solverType = 0 - default, fortran translated solver
		* solverType = 1 - static old
		* solverType = 2 - static 2 solver (mkl)
		* solverType = other - cpp, rewritten solver
		*/
		void Test2x1x10(int solverType)
		{
			stringstream str;
			str << "test_stress_type_" << solverType;

			string fileRlc = str.str() + ".rlc";
			string fileResult = str.str() + ".mpr";

			SpecialParams specialParams;
			IntegrationParams integrationParams;
			GridParams gridParams;
			specialParams._E = 100;
			specialParams._density = 1000;
			specialParams._dampingRatio = 1;
			specialParams._scaleFactor = 1;// 1e8;

			//integrationParams._nIterations = 500;

			integrationParams._nIterations = 100;
			integrationParams._nSubIterations = 1000;
			integrationParams._timeStep = 0.00001f;

			if (solverType == 2)
			{
				integrationParams._nIterations = 1;
				integrationParams._nSubIterations = 1;
			}

			enum ECode
			{
				xlr,
				xrl,
				yfb,
				ybf,
				ztb,
				zbt
			};

			SolverHandler _hsolver;
			ECode code = xlr;
			switch (code)
			{
			case xlr:
				gridParams._nx = 10;
				gridParams._ny = 2;
				gridParams._nz = 1;
				gridParams._gridStep = 0.1;

				_hsolver = MakeSolver
					(
					gridParams,
					specialParams,
					integrationParams,
					fileRlc,
					face_left,
					face_right,
					10,
					dof_y,
					solverType
					);
				break;
			case ztb:
				gridParams._nx = 1;
				gridParams._ny = 1;
				gridParams._nz = 3;
				gridParams._gridStep = 0.1;

				_hsolver = MakeSolver
					(
					gridParams,
					specialParams,
					integrationParams,
					fileRlc,
					face_top,
					face_bottom,
					1,
					dof_x,
					solverType
					);
				break;
			}
			PerformanceCounter pc;
			pc.Start();
			if (_hsolver != nullptr)
			{
				OverrideStiffness(_hsolver,
					1000.,
					10.,
					10.,
					10.,
					1.);

				Solve
					(
					_hsolver,
					fileResult,
					gridParams,
					integrationParams
					);
				pc.Print("Solving time: ", true);
			}
			Stress::ReleaseMemory((void* &)_hsolver);
		}
	}
}