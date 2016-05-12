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
			double _E;			  // модуль упругости
			double _dampingRatio; // коэффициент демпфирования
			double _density;	  // плотность материала
			double _scaleFactor;  // делитель модуля упругости

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