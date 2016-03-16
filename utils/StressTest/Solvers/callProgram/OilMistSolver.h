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

			double _tConduction;	// ����. ����������������
			double _tCapacity;		// �������� ������������
			double _density;		// ��������� ���������
			double _tBody;			// ����������� ����


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


		// ������ � ���������� ���������

		/**
		* ��������� ��������� ��������� ����� ������� ����� ������ ���������
		* @param oilMistSolver - �������� ��������� ������
		* @param selectedArea - ���������� ������� �����
		* @param scalarValue - �������� ���������� ���������
		*/
		void FillAreaWithScalarValue
			(
				IOilMistSolver* oilMistSolver,
				const GridSelectedArea& selectedArea,
				double scalarValue
			);

		/**
		* ��������� ��������� ��������� ����� ������� ����� ������ ���������
		* @param oilMistSolver - �������� ��������� ������
		* @param selectedArea - ���������� ������� �����
		* @param vectorValue - �������� ���������� ���������
		*/
		void FillAreaWithVectorValue
			(
				IOilMistSolver* oilMistSolver,
				const GridSelectedArea& selectedArea,
				const double(&vectorValue)[3]
			);

	}

}