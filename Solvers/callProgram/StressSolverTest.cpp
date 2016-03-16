#include "StressSolverTest.h"

#include "../Stress/StressStrainSolverExports.h"
#include "CommonSolversTest.h"
#include "PerformanceCounter.h"

#include <vector>


using std::string;
using std::vector;


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
			SetLeftSealRightForceBc(_hsolver, -10.01, gridParams);

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
			double bcParams2[6] = { 0, 0, force, 0, 0, 0 };

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

		void Test()
		{
			const string& fileRlc = "test_stress.rlc";
			const string& fileResult = "test_stress.mpr";
			SpecialParams specialParams;
			IntegrationParams integrationParams;
			GridParams gridParams;

			specialParams._scaleFactor = 1e6;
			//specialParams._dampingRatio = 1e-4;


			integrationParams._nIterations = 600;// 600;
			integrationParams._nSubIterations = 300;//300;
			integrationParams._timeStep = 0.00001f;
			gridParams._nz = 3;
			gridParams._ny = 3;
			gridParams._nx = 10;
			gridParams._gridStep = 0.005;
			const int solverType = 0; 

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