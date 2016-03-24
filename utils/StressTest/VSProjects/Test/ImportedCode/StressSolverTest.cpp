#include "StressSolverTest.h"

#include "../Stress/StressStrainSolverExports.h"
#include "CommonSolversTest.h"
#include "PerformanceCounter.h"


#include <sstream>
#include <vector>


using std::string;
using std::vector;
using std::stringstream;

namespace SpecialSolversTest
{

	namespace StressStrainStuff
	{


		StressStrainSolver MakeSolver
			(
				const GridParams& gridParams,
				const SpecialParams& specialParams,
				const IntegrationParams& integrationParams,
				const string& fileRlc,
				const int solverType
			)
		{
			double* nodes = nullptr;
			int* links = nullptr;
			int nLinks;

			CreateTestGrid
				(
					nodes,
					links,
					nLinks,
					gridParams._nx,
					gridParams._ny,
					gridParams._nz,
					gridParams._gridStep,
					fileRlc
				);

			double params[5];
			specialParams.GetParams(params);
			StressStrainSolver _hsolver = Stress::Init
				(
					params,
					links,
					nLinks,
					nodes,
					gridParams.NodesCount(),
					gridParams._gridStep,
					integrationParams._timeStep,
					0,
					0,
					false,
					solverType // -1 = ÑPP // 0 - fortr // 2 = static2
				);

			if (nodes != nullptr)
				delete [] nodes;
			if (links != nullptr)
				delete [] links;
			SetLeftSealRightForceBc(_hsolver, 10, gridParams);

			return _hsolver;
		}

		void SetLeftSealRightForceBc
			(
				StressStrainSolver hStressSolver,
				double force,
				const GridParams& gridParams
			)
		{
			double bcParams1[6] = {-1, -1, -1, -1, -1, -1};
			double bcParams2[6] = { 0, force, 0, 0, 0, 0 };

			vector<int> bcIndices1;

			for (int i = 0; i < gridParams._nz * gridParams._ny; i++)
			{
				bcIndices1.push_back(i + 1);
			}

			vector<int> bcIndices2;

			//bcIndices2.push_back(nNodes);
			for (int i = 0; i < gridParams._nz * gridParams._ny; i++)
			{
				bcIndices2.push_back(gridParams.NodesCount() - i);
			}

			Stress::AddBoundary
				(
					hStressSolver,
					&bcIndices1[0],
					static_cast<int>(bcIndices1.size()),
					4,
					bcParams1
				);
			Stress::AddBoundary
				(
					hStressSolver,
					&bcIndices2[0],
					static_cast<int>(bcIndices2.size()),
					3,
					bcParams2
				);
		}

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

			if (solverType == 2)
			{
				integrationParams._nIterations = 1;
				integrationParams._nSubIterations = 1;
			}

			gridParams._nz = 1;
			gridParams._ny = 3;
			gridParams._nx = 10;
			gridParams._gridStep = 0.005;

			StressStrainSolver _hsolver = MakeSolver
				(
				gridParams,
				specialParams,
				integrationParams,
				fileRlc,
				solverType
				);

			SolverPerformanceCounter pc;

			pc.Reset();
			Solve
				(
				_hsolver,
				fileResult,
				gridParams,
				integrationParams
				);
			pc.Print("Solving time: ");
			Stress::ReleaseMemory((void* &)_hsolver);
		}

		/** Stress solver testing for grid 3x3x10
		* @param solverType - type of solver
		* solverType = 0 - default, fortran translated solver
		* solverType = 1 - static old
		* solverType = 2 - static 2 solver (mkl)
		* solverType = other - cpp, rewritten solver
		*/
		void Test3x3x10(int solverType)
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

			if (solverType == 2)
			{
				integrationParams._nIterations = 1;
				integrationParams._nSubIterations = 1;
			}

			gridParams._nz = 3;
			gridParams._ny = 3;
			gridParams._nx = 10;
			gridParams._gridStep = 0.005;

			StressStrainSolver _hsolver = MakeSolver
				(
					gridParams,
					specialParams,
					integrationParams,
					fileRlc,
					solverType
				);
	
			SolverPerformanceCounter pc;

			pc.Reset();
			Solve
				(
					_hsolver,
					fileResult,
					gridParams,
					integrationParams
				);
			pc.Print("Solving time: ");
			Stress::ReleaseMemory((void* &)_hsolver);
		}


		/** Stress solver testing for grid 3x3x10
		* @param solverType - type of solver
		* solverType = 0 - default, fortran translated solver
		* solverType = 1 - static old
		* solverType = 2 - static 2 solver (mkl)
		* solverType = other - cpp, rewritten solver
		*/
		void Test3x1x1(int solverType)
		{
			stringstream str;
			str << "test_stress_type_" << solverType;

			string fileRlc = str.str() + ".rlc";
			string fileResult = str.str() + ".mpr";

			SpecialParams specialParams;
			IntegrationParams integrationParams;
			GridParams gridParams;
			specialParams._E = 1000;
			specialParams._density = 1000;			
			specialParams._dampingRatio = 10;
			specialParams._scaleFactor = 1;// 1e8;
			//specialParams._dampingRatio = 1e-4;

			//integrationParams._nIterations = 120;// 600;
			//integrationParams._nSubIterations = 60;//300;
			integrationParams._nIterations = 1000;
			integrationParams._nSubIterations = 100;
			integrationParams._timeStep = 0.00005f;

			if (solverType == 2)
			{
				integrationParams._nIterations = 1;
				integrationParams._nSubIterations = 1;
			}

			gridParams._nx = 3;
			gridParams._ny = 1;
			gridParams._nz = 1;
			gridParams._gridStep = 0.1;

			StressStrainSolver _hsolver = MakeSolver
				(
				gridParams,
				specialParams,
				integrationParams,
				fileRlc,
				solverType
				);

			SolverPerformanceCounter pc;

			pc.Reset();
			Solve
				(
				_hsolver,
				fileResult,
				gridParams,
				integrationParams
				);
			pc.Print("Solving time: ");
			Stress::ReleaseMemory((void* &)_hsolver);
		}

		void TestSolveSystemOfLinearEquationsForStiffness()
		{
			Stress::SolveSystemOfLinearEquationsForStiffness();
		}

	}

}