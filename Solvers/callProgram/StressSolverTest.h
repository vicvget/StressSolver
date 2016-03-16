#pragma once

#include "StressSolver.h"


namespace SpecialSolversTest
{

	namespace StressStrainStuff
	{

		using namespace SpecialSolvers;
		using namespace SpecialSolvers::StressStrainStuff;


		StressStrainSolver MakeSolver
			(
				const GridParams& gridParams,
				const SpecialParams& specialParams,
				const IntegrationParams& integrationParams,
				const std::string& fileRes,
				const int solverType
			);

		void SetLeftSealRightForceBc
			(
				StressStrainSolver hStressSolver,
				double force,
				const GridParams& gridParams
			);

		void Test();

		void TestSolveSystemOfLinearEquationsForStiffness();

	}

}