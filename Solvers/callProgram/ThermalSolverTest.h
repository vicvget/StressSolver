#pragma once


#include "ThermalSolver.h"


namespace SpecialSolversTest
{

	namespace ThermalStuff
	{

		using namespace SpecialSolvers;
		using namespace SpecialSolvers::ThermalStuff;


		ThermalSolver MakeSolver
			(
				const GridParams& gridParams,
				const SpecialParams& specialParams,
				const IntegrationParams& integrationParams,
				const std::string& fileResults
			);

		void SetLeftHeatingBc
			(
				ThermalSolver hStressSolver,
				double thermalFlow,
				const GridParams& gridParams,
				const SpecialParams& specialParams
			);

		void Test();

	}

}