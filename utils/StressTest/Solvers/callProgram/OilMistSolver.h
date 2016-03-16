#pragma once


#include "CommonSolvers.h"
#include "../OilMist/IOilMistSolver.h"
#include "../../FormatProviders/GridProvider/GridSelectedArea.h"

#include <memory>


namespace SpecialSolvers
{

	namespace OilMistStuff
	{

		using OilMistSolver = std::unique_ptr<IOilMistSolver, void(*)(IOilMistSolver*)>;


		struct SpecialParams
		{

			double _tConduction;	// коэф. теплопроводности
			double _tCapacity;		// удельна€ теплоемкость
			double _density;		// плотность материала
			double _tBody;			// температура тела


			SpecialParams();

			SpecialParams
				(
					double tConduction,
					double tCapacity,
					double density,
					double tBody
				);

			void GetParams
				(
					double (&params)[4]
				)	const;

		};


		void ExportNormals
			(
				const std::vector<BoundaryNormal*>& boundaryNormals,
				std::size_t nodesCount,
				double (*&normals)[3]
			);

		OilMistSolver MakeSolver
			(
				const SpecialParams& specialParams,
				const IntegrationParams& integrationParams,
				GridParams& gridParams,
				const std::string& rlcGridFileName
			);

		void Solve
			(
				IOilMistSolver* oilMistSolver,
				const std::string& fileResults,
				const GridParams& gridParams,
				const IntegrationParams& integrationParams
			);


		// –абота с начальными услови€ми

		/**
		* «аполнить скал€рные параметры узлов области сетки данным значением
		* @param oilMistSolver - решатель масл€ного тумана
		* @param selectedArea - выделенна€ область сетки
		* @param scalarValue - значение скал€рного параметра
		*/
		void FillAreaWithScalarValue
			(
				IOilMistSolver* oilMistSolver,
				const GridSelectedArea& selectedArea,
				double scalarValue
			);

		/**
		* «аполнить векторные параметры узлов области сетки данным значением
		* @param oilMistSolver - решатель масл€ного тумана
		* @param selectedArea - выделенна€ область сетки
		* @param vectorValue - значение векторного параметра
		*/
		void FillAreaWithVectorValue
			(
				IOilMistSolver* oilMistSolver,
				const GridSelectedArea& selectedArea,
				const double(&vectorValue)[3]
			);

	}

}