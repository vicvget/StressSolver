#pragma once


#include "CommonSolvers.h"

#include <string>


namespace SpecialSolvers
{

	namespace ThermalStuff
	{

		typedef void* ThermalSolver;


		struct SpecialParams
		{

			double _tConduction; // коэф. теплопроводности
			double _tCapacity;	// удельная теплоемкость
			double _density;	// плотность материала
			double _tBody;		// температура тела

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
					double* params
				)	const;

		};


		void Solve
			(
				ThermalSolver hStressSolver,
				const std::string& fileResults,
				const GridParams& gridParams,
				const IntegrationParams& integrationParams
			);

	}

}