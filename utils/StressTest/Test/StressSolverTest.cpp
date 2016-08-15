#include "CommonSolversTest.h"
#include "StressSolverTest.h"

#include "../Solvers/Stress/StressStrainCppSolver.h"
#include "../Solvers/Stress/StressStrainSolverExports.h"
#include "../Solvers/Stress/FTimer.h"

#include <sstream>
#include <vector>
#include <omp.h>


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
			const string& solverUid,
			EFACE sealedFace,
			EFACE forcedFace,
			double force,
			EDOF dof,
			const int solverType)
		{
			double* nodes = nullptr;
			int* links = nullptr;
			int nLinks;
			string fileRlc = solverUid + ".rlc";
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
				solverUid,
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
				solverType
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

		SolverHandler MakePlateSolver(
			const GridParams& gridParams,
			const SpecialParams& specialParams,
			const IntegrationParams& integrationParams,
			const string& solverUid,
			double force,
			EDOF dof,
			const int solverType)
		{
			double* nodes = nullptr;
			int* links = nullptr;
			int nLinks;
			string fileRlc = solverUid + ".rlc";
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
			int numThreads = omp_get_max_threads();

			SolverHandler _hsolver = Stress::Init
				(
				solverUid,
				params,
				links,
				nLinks,
				nodes,
				gridParams.NodesCount(),
				gridParams._gridStep,
				integrationParams._timeStep,
				0,
				numThreads,
				false,
				solverType
				);

			if (nodes != nullptr)
				delete[] nodes;
			if (links != nullptr)
				delete[] links;
			size_t size = dof == dof_x ? gridParams._ny : gridParams._nx;
			SetSealedFullForcePlateBc(_hsolver, size, force, dof);
			//SetSealedForcePlateBc(_hsolver, size, force, dof);
			return _hsolver;
		}


		SolverHandler MakeSolver(const GridParams& gridParams, const SpecialParams& specialParams, const IntegrationParams& integrationParams, const std::string& solverUid, const int solverType)
		{
			return MakeSolver(gridParams, specialParams, integrationParams, solverUid, face_left, face_right, 10, dof_y, solverType);
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
						bcIndices.push_back((int)(i*gridParams._nx*gridParams._ny + j * gridParams._nx + 1));
					}
				break;
			case face_top:
				for (size_t i = 0; i < gridParams._nx; i++)
					for (size_t j = 0; j < gridParams._ny; j++)
					{
						bcIndices.push_back((int)(i*gridParams._nx*gridParams._ny + j * gridParams._nx + gridParams._nz));
					}
				break;
			case face_front:
				for (size_t i = 0; i < gridParams._nx; i++)
					for (size_t j = 0; j < gridParams._nz; j++)
					{
						bcIndices.push_back((int)(i*gridParams._ny*gridParams._nz + j + 1));
					}
				break;
			case face_back:
				for (size_t i = 0; i < gridParams._nx; i++)
					for (size_t j = 0; j < gridParams._nz; j++)
					{
						bcIndices.push_back((int)((i+1)*gridParams._ny*gridParams._nz - j));
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
			bcIndices.push_back((int)nodeId);

			//std::cout << "Force applied to node " << nodeId << std::endl;
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

		void AddSealedPlateBoundary(
			SolverHandler hStressSolver,
			size_t side)
		{
			double bcParams[6] = { -1, -1, -1, -1, -1, -1 };
			vector<int> bcIndices;
			bcIndices.push_back(1);
			bcIndices.push_back(side);
			bcIndices.push_back(side*side);
			bcIndices.push_back(side*(side - 1)+1);

			Stress::AddBoundary
				(
				hStressSolver,
				&bcIndices[0], 
				static_cast<int>(bcIndices.size()),
				4,
				bcParams
				);
		}

		void AddSealedFullPlateBoundary(
			SolverHandler hStressSolver,
			size_t side)
		{
			double bcParams[6] = { -1, -1, -1, -1, -1, -1 };
			vector<int> bcIndices;
			for (int i = 1; i <= side; i++)
			{
				bcIndices.push_back(i);
				bcIndices.push_back(side*side-i+1);
			}
			for (int i = 1; i < side; i++)
			{
				bcIndices.push_back(side*i+1);
				bcIndices.push_back(side*i+side);
			}

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

		void SetSealedForcePlateBc(
			SolverHandler hStressSolver,
			size_t side,
			double force,
			EDOF dof)
		{
			AddSealedPlateBoundary(hStressSolver, side);
			AddElementForceBoundary(hStressSolver, force, dof, side*(side / 2) + side / 2 + 1);
		}

		void SetSealedFullForcePlateBc(
			SolverHandler hStressSolver,
			size_t side,
			double force,
			EDOF dof)
		{
			AddSealedFullPlateBoundary(hStressSolver, side);
			AddElementForceBoundary(hStressSolver, force, dof, side*(side / 2) + side / 2 + 1);
		}


		// ��������� 2 ���� � ���� �������
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