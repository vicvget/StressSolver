#include "OilMistSolverTest.h"

#include "PerformanceCounter.h"


namespace SpecialSolversTest
{

	namespace OilMistStuff
	{

		void SetLeftHeatingBc
			(
				IOilMistSolver* oilMistSolver,
				const GridParams& gridParams,
				const SpecialParams& specialParams
			)
		{
			std::vector<double> bcParams1;

			bcParams1.push_back(0.2/*101*/);

			int bcIndicesCount = gridParams._nz * gridParams._ny;
			std::vector<int> bcIndices1(bcIndicesCount);

			for (int i = 0; i < bcIndicesCount; i++)
			{
				bcIndices1[i] = i * gridParams._nx;
			}

			oilMistSolver->AddBoundaryCondition(bcIndices1, 0, bcParams1);
		}

		void SetFreeCoolingBc
			(
				IOilMistSolver* oilMistSolver,
				double oilMistFlow,
				const GridParams& gridParams,
				const SpecialParams& specialParams
			)
		{
			std::vector<double> bcParams1;

			bcParams1.push_back(1000);

			std::vector<double> bcParams2;

			bcParams2.push_back(47);
			bcParams2.push_back(specialParams._tBody);

			std::vector<int> bcIndices1;

			for (int i = 0; i < gridParams._nz * gridParams._ny; i++)
			{
				bcIndices1.push_back(i * gridParams._nx);
			}

			//std::vector<int> bcIndices2;

			//bcIndices2.push_back(nNodes);
			//for (int i = 0; i < gridParams._nz * gridParams._ny; i++)
			//{
			//	bcIndices2.push_back(gridParams.NodesCount() - i);
			//}
			oilMistSolver->AddBoundaryCondition(bcIndices1, 1, bcParams1);
			//OilMist::AddBoundaryCondition(hOilMistSolver, bcIndices2, 2, bcParams2);
		}
		
		void Test()
		{

			// SET DATA
			//const string& fileRlc = "test_oilMist2.rlc";
			//const string& fileRlc = "step_4_part_axle_case_inner.result_good.rlc";
			//const string& fileRlc = "65601-2502011-50_asm_inner.rlc";
			const string& fileRlc = "65601-2502011-50_asm_inner_step_6_180k.rlc";
			//const string& fileRlc = "65601-2502011-50_asm_inner_step_4_600k.rlc";
			//const string& fileRlc = "step_4_part_shaft.rlc";
			//const string& fileRlc = "65601-differential_assembly_488k.rlc";
			//const string& fileRlc = "65601-2502011-50_asm_inner_step_1_38.6kk.rlc";
			//const string& fileRlc = "sphere_plain_step_2.4_52k.rlc";
			//const string& fileRlc = "sphere_plain_step_1.1_557k.rlc";
			//const string& fileRlc = "sphere_plain_step_1.1_557k_v1.rlc";

			const string& fileResult = "test_oilMist2.mpr";
			SpecialParams specialParams;
			IntegrationParams integrationParams;
			GridParams gridParams;

			specialParams._tBody = 0;
			gridParams._gridStep = 0.006;
			integrationParams._nIterations = 100*5;// 500 * 4 * 5;// *4 * 5 * 3;// *2 * 9;
			integrationParams._nSubIterations = 1;
			integrationParams._timeStep = 0.001f;// 25;

			/*gridParams._nx = 25;
			gridParams._ny = 25;
			gridParams._nz = 25;
			CreateTestGrid(gridParams, fileRlc);*/

			OilMistSolver oilMistSolver = MakeSolver
				(
					specialParams,
					integrationParams,
					gridParams,
					fileRlc
				);

			if (!oilMistSolver)
			{
				std::cout << "Oil mist solver was not created: error!" << std::endl;

				return;
			}
			/*oilMistSolver->FillGridWithScalarValue(0.0);

			double point[]{0.0, 0.05, 0.0};

			oilMistSolver->FillGridWithScalarValueToPlane
				(
					point,
					CoordinateDirection::YMinus,
					10.0
				);*/

			SetLeftHeatingBc
				(
					oilMistSolver.get(),
					gridParams,
					specialParams
				);

			SolverPerformanceCounter pc;

			pc.Reset();
			Solve
				(
					oilMistSolver.get(),
					fileResult,
					gridParams,
					integrationParams
				);
			pc.Print("test oilMist: ");
		}

	}

}