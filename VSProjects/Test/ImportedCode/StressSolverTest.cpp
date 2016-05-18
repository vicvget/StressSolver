#include "StressSolverTest.h"

#include "../Stress/StressStrainSolverExports.h"
#include "CommonSolversTest.h"
#include "PerformanceCounter.h"


#include <sstream>
#include <vector>
#include <StressStrainCppSolver.h>


using std::string;
using std::vector;
using std::stringstream;

namespace SpecialSolversTest
{

	namespace StressStrainStuff
	{
		SolverHandler MakeSolver(
			const GridParams& gridParams,
			const SpecialParams& specialParams,
			const IntegrationParams& integrationParams,
			const string& fileRlc,
			EFACE sealedFace,
			EFACE forcedFace,
			double force,
			EDOF dof,
			const int solverType)
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
			SolverHandler _hsolver = Stress::Init
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
				solverType // -1 = СPP // 0 - fortr // 2 = static2
				);

			if (nodes != nullptr)
				delete[] nodes;
			if (links != nullptr)
				delete[] links;
			SetSealedForceBc(_hsolver, sealedFace, forcedFace, force, dof, gridParams);
			//SetSealedForceForceBc(_hsolver, sealedFace, forcedFace, force, dof, (EDOF)((dof+2)%3), gridParams);
			//SetElementForceBc(_hsolver, force*gridParams._gridStep*gridParams._gridStep*10, 119, EDOF::dof_ry);
			

			return _hsolver;
		}

		SolverHandler MakeSolver(const GridParams& gridParams, const SpecialParams& specialParams, const IntegrationParams& integrationParams, const std::string& fileRlc, const int solverType)
		{
			return MakeSolver(gridParams, specialParams, integrationParams, fileRlc, face_left, face_right, 10, dof_y, solverType);
		}

		//void OverrideStiffness(SolverHandler solver_handler, int i, int i1, int i2, int i3, int i4);

		void OverrideStiffness(
			SolverHandler hStressSolver,
			double elasticFactorLinear,
			double elasticFactorAngular,
			double dampingFactorLinear,
			double dampingFactorAngular,
			double stiffnessScale)
		{
			((Stress::StressStrainCppSolver*)hStressSolver)->OverrideStiffness(
				elasticFactorLinear,
				elasticFactorAngular,
				dampingFactorLinear,
				dampingFactorAngular,
				stiffnessScale);
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
			specialParams._E = 100;
			specialParams._density = 1000;
			specialParams._dampingRatio = 1;
			specialParams._scaleFactor = 1;// 1e8;

			//integrationParams._nIterations = 500;
			
			integrationParams._nIterations = 10000;
			integrationParams._nSubIterations = 10;
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
				gridParams._nx = 3;
				gridParams._ny = 1;
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
					dof_z,
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
			SolverPerformanceCounter pc;
			pc.Reset();
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
				pc.Print("Solving time: ");
			}
			Stress::ReleaseMemory((void* &)_hsolver);
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

			SolverHandler _hsolver = MakeSolver
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

			SolverHandler _hsolver = MakeSolver
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


		void FillFaceIndices(vector<int>& bcIndices, EFACE face, const GridParams& gridParams)
		{
			switch (face)
			{
			case face_left:
				for (int i = 0; i < gridParams._nz * gridParams._ny; i++)
				{
					bcIndices.push_back(i + 1);
				}
				break;
			case face_right:
				for (int i = 0; i < gridParams._nz * gridParams._ny; i++)
				{
					bcIndices.push_back(gridParams.NodesCount() - i);
				}
				break;
			case face_bottom:
				for (size_t i = 0; i < gridParams._nx; i++)
					for (size_t j = 0; j < gridParams._ny; j++)
					{
						bcIndices.push_back(i*gridParams._nx*gridParams._ny + j * gridParams._nx + 1);
					}
				break;
			case face_top:
				for (size_t i = 0; i < gridParams._nx; i++)
					for (size_t j = 0; j < gridParams._ny; j++)
					{
						bcIndices.push_back(i*gridParams._nx*gridParams._ny + j * gridParams._nx + gridParams._nz);
					}
				break;
			case face_front:
				for (size_t i = 0; i < gridParams._nx; i++)
					for (size_t j = 0; j < gridParams._nz; j++)
					{
						bcIndices.push_back(i*gridParams._ny*gridParams._nz + j + 1);
					}
				break;
			case face_back:
				for (size_t i = 0; i < gridParams._nx; i++)
					for (size_t j = 0; j < gridParams._nz; j++)
					{
						bcIndices.push_back((i+1)*gridParams._ny*gridParams._nz - j);
					}
				break;
			}
		}

		void AddForceBoundary(
			SolverHandler hStressSolver,
			double force,
			EDOF dof,
			EFACE face,
			const GridParams& gridParams)
		{
			double bcParams[6] = { 0 };
			bcParams[dof] = force;
			vector<int> bcIndices;
			FillFaceIndices(bcIndices, face, gridParams);

			Stress::AddBoundary
				(
				hStressSolver,
				&bcIndices[0],
				static_cast<int>(bcIndices.size()),
				3,
				bcParams
				);
		}

		void AddElementForceBoundary(
			SolverHandler hStressSolver,
			double force,
			EDOF dof,
			size_t nodeId)
		{
			double bcParams[6] = { 0 };
			bcParams[dof] = force;
			vector<int> bcIndices;			
			bcIndices.push_back(nodeId);


			Stress::AddBoundary
				(
				hStressSolver,
				&bcIndices[0],
				static_cast<int>(bcIndices.size()),
				3,
				bcParams
				);
		}


		void AddSealedBoundary(
			SolverHandler hStressSolver,
			EFACE face,
			const GridParams& gridParams)
		{
			double bcParams[6] = { -1, -1, -1, -1, -1, -1 };
			vector<int> bcIndices;
			FillFaceIndices(bcIndices, face, gridParams);

			Stress::AddBoundary
				(
				hStressSolver,
				&bcIndices[0],
				static_cast<int>(bcIndices.size()),
				4,
				bcParams
				);
		}


		void SetSealedForceBc(
			SolverHandler hStressSolver,
			EFACE faceSealed,
			EFACE faceForced,
			double force,
			EDOF dof,
			const GridParams& gridParams)
		{
			AddSealedBoundary(hStressSolver, faceSealed, gridParams);
			AddForceBoundary(hStressSolver, force, dof, faceForced, gridParams);
		}

		// добавляет 2 силы и одну заделку
		void SetSealedForceForceBc(
			SolverHandler hStressSolver,
			EFACE faceSealed,
			EFACE faceForced,
			double force,
			EDOF dof,
			EDOF dof2,
			const GridParams& gridParams)
		{
			AddSealedBoundary(hStressSolver, faceSealed, gridParams);
			AddForceBoundary(hStressSolver, force, dof, faceForced, gridParams);
			AddForceBoundary(hStressSolver, force, dof2, faceForced, gridParams);
		}

		void SetSealedBc(
			SolverHandler hStressSolver,
			EFACE faceSealed,
			const GridParams& gridParams)
		{
			AddSealedBoundary(hStressSolver, faceSealed, gridParams);
		}

		void SetElementForceBc(
			SolverHandler hStressSolver,
			double torque,
			size_t nodeId,
			EDOF dof)
		{
			AddElementForceBoundary(hStressSolver, torque, dof, nodeId);
		}

		void SetLeftSealRightForceBc
			(
			SolverHandler hStressSolver,
			double force,
			EDOF dof,
			const GridParams& gridParams
			)
		{
			SetSealedForceBc(hStressSolver, face_left, face_right, force, dof, gridParams);
		}

		void SetLeftSealRightForceBc1
			(
			SolverHandler hStressSolver,
			double force,
			const GridParams& gridParams
			)
		{
			double bcParams1[6] = { -1, -1, -1, -1, -1, -1 };
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
	}

}