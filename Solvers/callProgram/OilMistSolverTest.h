#pragma once


#include "OilMistSolver.h"


namespace SpecialSolversTest
{

	namespace OilMistStuff
	{

		using namespace SpecialSolvers;
		using namespace SpecialSolvers::OilMistStuff;


		void SetLeftHeatingBc
			(
				IOilMistSolver* oilMistSolver,
				const GridParams& gridParams,
				const SpecialParams& specialParams
			);

		void SetFreeCoolingBc
			(
				IOilMistSolver* oilMistSolver,
				double oilMistFlow,
				const GridParams& gridParams,
				const SpecialParams& specialParams
			);

		void Test();

	}

}