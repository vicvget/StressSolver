#pragma once


#include "CommonSolvers.h"

#include <string>


namespace SpecialSolvers
{

	namespace StressStrainStuff
	{

		typedef void* SolverHandler;

		struct SpecialParams
		{
			double _E;			  // ������ ���������
			double _dampingRatio; // ����������� �������������
			double _density;	  // ��������� ���������
			double _scaleFactor;  // �������� ������ ���������

			SpecialParams();

			SpecialParams
				(
					double E,
					double dampingRatio,
					double density,
					double scaleFactor
				);

			void GetParams
				(
					double* params
				)	const;

		};


		void Solve
			(
				SolverHandler hStressSolver,
				const std::string& fileResults,
				const GridParams& gridParams,
				const IntegrationParams& integrationParams
			);

	}

}