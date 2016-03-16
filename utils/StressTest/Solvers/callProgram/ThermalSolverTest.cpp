#include "ThermalSolverTest.h"

#include "../Thermal_new/ThermalSolverExports.h"
#include "CommonSolversTest.h"
#include "PerformanceCounter.h"
#include "../../FormatProviders/GridProvider/MeshDataSolverFormatter.h"


using std::string;


namespace SpecialSolversTest
{

	namespace ThermalStuff
	{

		ThermalSolver MakeSolver
			(
				const GridParams& gridParams,
				const SpecialParams& specialParams,
				const IntegrationParams& integrationParams,
				const string& fileRlc
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

			int* nodesLinks;			// массив связей, содержит по 6 индексов узлов, с которыми связан текущий узел
			MeshDataSolverFormatter* gridProvider = new MeshDataSolverFormatter(fileRlc);

			gridProvider->ExportGrid(nodesLinks);

			double params[5];

			specialParams.GetParams(params);
			ThermalSolver hSolver = Thermal::InitSolver
				(
					gridParams.NodesCount(),
					nodesLinks,
					gridParams._gridStep,
					integrationParams._timeStep,
					params,
					//true,	// use float
					true,	// use float
					0,		// processor type
					1		// threads count
				);

			//Thermal::AddThermalLinks(hSolver, nodesLinks);
			//hSolver->AddLinks(nodesLinks);// size = numberOfNodes*6
			delete[] nodes;
			delete[] links;
			delete[] nodesLinks;

			return hSolver;
		}

		void SetLeftHeatingBc
			(
				ThermalSolver hThermalSolver,
				double thermalFlow,
				const GridParams& gridParams,
				const SpecialParams& specialParams
			)
		{

			std::vector<int> bcIndices1;
			std::vector<double> bcParams1;
			bcParams1.push_back(1000);
			std::vector<double> bcParams2;
			bcParams2.push_back(47);
			bcParams2.push_back(specialParams._tBody);

			for (int i = 0; i < gridParams._nz * gridParams._ny; i++)
			{
				bcIndices1.push_back(i * gridParams._nx);
			}

			std::vector<int> bcIndices2;

			//bcIndices2.push_back(nNodes);
			//for (int i = 0; i < gridParams._nz * gridParams._ny; i++)
			//{
			//	bcIndices2.push_back(gridParams.NodesCount() - i);
			//}
			Thermal::AddBoundaryCondition(hThermalSolver,bcIndices1, 1, bcParams1);
			//Thermal::AddBoundaryCondition(hThermalSolver,bcIndices2, 2, bcParams2);
		}

		void Test()
		{
			const string& fileRlc = "test_thermal.rlc";
			const string& fileResult = "test_thermal.mpr";
			SpecialParams specialParams;
			IntegrationParams integrationParams;
			GridParams gridParams;

			gridParams._gridStep = 0.001;
			integrationParams._nIterations = 100;
			integrationParams._nSubIterations = 100;

			gridParams._nx = 256;
			gridParams._ny = 256;
			gridParams._nz = 256;

			ThermalSolver _hsolver = MakeSolver
				(
					gridParams,
					specialParams,
					integrationParams,
					fileRlc
				);
	
			SetLeftHeatingBc
				(
					_hsolver,
					1000,
					gridParams,
					specialParams
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
			pc.Print("test thermal: ");

			Thermal::ReleaseSolver(_hsolver);
		}

	}

}