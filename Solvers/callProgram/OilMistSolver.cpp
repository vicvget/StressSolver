#include "OilMistSolver.h"

#include "../OilMist/OilMistSolverExports.h"
#include "../../FormatProviders/GridProvider/MeshDataSolverFormatter.h"
#include "../../FormatProviders/ResProvider/ProviderMpr.h"


using std::ofstream;


namespace SpecialSolvers
{
	
	namespace OilMistStuff
	{

		// struct SpecialParams

		SpecialParams::SpecialParams()
			:
				_tConduction(46.5),
				_tCapacity(460),
				_density(7800),
				_tBody(20)
		{
		}

		SpecialParams::SpecialParams
			(
				double tConduction,
				double tCapacity,
				double density,
				double tBody
			)
			:
				_tConduction(tConduction),
				_tCapacity(tCapacity),
				_density(density),
				_tBody(tBody)
		{
		}

		void SpecialParams::GetParams
			(
				double (&params)[4]
			)	const
		{
			params[0] = _tConduction;
			params[1] = _tCapacity;
			params[2] = _density;
			params[3] = _tBody;
		}

		void ExportNormals
			(
				const vector<BoundaryNormal*>& boundaryNormals,
				size_t nodesCount,
				double (*&normals)[3]
			)
		{
			normals = new double[nodesCount][3];
			for (size_t nodeIndex = 0; nodeIndex < nodesCount; ++nodeIndex)
			{
				std::fill_n(normals[nodeIndex], 3, 0.0);
			}
			for (const BoundaryNormal* boundaryNormal : boundaryNormals)
			{
				size_t nodeNumber = boundaryNormal->GetPointNumber();
				size_t nodeIndex = nodeNumber - 1;
				const VertexPoint& unitaryNormal = boundaryNormal->GetUnitaryNormal();
				double(&normal)[3] = normals[nodeIndex];

				normal[0] = unitaryNormal.x;
				normal[1] = unitaryNormal.y;
				normal[2] = unitaryNormal.z;
			}
		}

		OilMistSolver MakeSolver
			(
				const SpecialParams& specialParams,
				const IntegrationParams& integrationParams,
				GridParams& gridParams,
				const string& rlcGridFileName
			)
		{
			auto gridProvider = MeshDataSolverFormatter::Create(rlcGridFileName);
			const vector<BoundaryNormal*>& boundaryNormals = gridProvider->GetNormals();
			int gridNodesCount = gridProvider->GetCount();

			gridParams._gridNodesCount = gridNodesCount;
			gridProvider->GetGridSize(gridParams._nx, gridParams._ny, gridParams._nz);

			double ex, ey, ez;

			gridProvider->GetElementSize(ex, ey, ez);
			if ((ex != ey) && (ex != ez))
			{
				return {nullptr, DeleteSolver};
			}
			gridParams._gridStep = ex / 1000.0;

			double* gridNodes{}; // массив связей, содержит по 6 индексов узлов, с которыми связан текущий узел
			int* gridLinks{}; // массив связей, содержит по 6 индексов узлов, с которыми связан текущий узел

			gridProvider->ExportGrid(gridLinks, gridNodes);

			std::unique_ptr<double[]> gridNodesPointer(gridNodes);
			std::unique_ptr<int[]> gridLinksPointer(gridLinks);
			double (*gridNormals)[3]{};

			ExportNormals(boundaryNormals, gridNodesCount, gridNormals);
			gridProvider.reset();

			std::unique_ptr<double[][3]> gridNormalsPointer(gridNormals);

			double params[4];

			specialParams.GetParams(params);

			IOilMistSolver* oilMistSolver = CreateSolver
				(
					gridNodesCount,
					gridNodesPointer.get(),
					gridLinksPointer.get(),
					gridNormalsPointer.get(),
					gridParams._gridStep,
					integrationParams._timeStep,
					params,
					true,	// use float
					0,		// processor type
					1		// threads count
				);

			//OilMist::AddOilMistLinks(hSolver, nodesLinks);
			//hSolver->AddLinks(nodesLinks);// size = numberOfNodes*6

			return {oilMistSolver, DeleteSolver};
		}

		void Solve
			(
				IOilMistSolver* oilMistSolver,
				const string& fileResults,
				const GridParams& gridParams,
				const IntegrationParams& integrationParams
			)
		{
			MultiphysicsResultsHeader mprHeader
				(
					SolverTypes::ST_Thermal,
					gridParams._gridNodesCount,
					integrationParams._initialTime,
					integrationParams.TotalTime()
				);

			oilMistSolver->UpdateReturnedBuffer();

			const float* data = oilMistSolver->GetReturnedBuffer();
			int dataSize = oilMistSolver->GetReturnedBufferSize();
			int subIteration = 0;
			ProviderMpr writer;
			writer.InitWriter(fileResults, &mprHeader);

			writer.WriteFrame(data, dataSize, (subIteration++) * integrationParams._timeStep);

			int percentCounter = 0;

			for (int iteration = 0; iteration < integrationParams._nIterations; ++iteration)
			{
				oilMistSolver->Solve();
				if ((iteration + 1) % integrationParams._nSubIterations == 0)
				{
					oilMistSolver->UpdateReturnedBuffer();
					writer.WriteFrame(data, dataSize, subIteration * integrationParams._timeStep);
				}
				if (iteration * 100.0 / integrationParams._nIterations > percentCounter)
				{
					std::cout << percentCounter << '%' << std::endl;
					++percentCounter;
				}
				++subIteration;
			}
		}

		/**
		* Заполнить скалярные параметры узлов области сетки данным значением
		* @param oilMistSolver - решатель масляного тумана
		* @param selectedArea - выделенная область сетки
		* @param scalarValue - значение скалярного параметра
		*/
		void FillAreaWithScalarValue
			(
				IOilMistSolver* oilMistSolver,
				const GridSelectedArea& selectedArea,
				double scalarValue
			)
		{
			const auto& areaPoints = selectedArea.GetAreaPoints();

			if (areaPoints.empty())
			{
				return;
			}
			oilMistSolver->FillAreaWithScalarValue
				(
					areaPoints.size(),
					areaPoints.data(),
					scalarValue
				);
		}

		/**
		* Заполнить векторные параметры узлов области сетки данным значением
		* @param oilMistSolver - решатель масляного тумана
		* @param selectedArea - выделенная область сетки
		* @param vectorValue - значение векторного параметра
		*/
		void FillAreaWithVectorValue
			(
				IOilMistSolver* oilMistSolver,
				const GridSelectedArea& selectedArea,
				const double(&vectorValue)[3]
			)
		{
			const auto& areaPoints = selectedArea.GetAreaPoints();

			if (areaPoints.empty())
			{
				return;
			}
			oilMistSolver->FillAreaWithVectorValue
				(
					areaPoints.size(),
					areaPoints.data(),
					vectorValue
				);
		}
	}

}